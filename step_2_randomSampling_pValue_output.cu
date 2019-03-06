
/* Execution Format : ./<exe> <drug_result_1_dict_compounds.txt> <drug_result_2_dict_compounds.txt> <drug_result_1_dict_proteins.txt> <drug_result_2_dict_proteins.txt> <para.txt> <drug name>
*/

#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/dir.h>
#include <stdbool.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

//Structure for key-value pairs in dictionary
struct kvpair {
    char *key;
    int value;
    struct kvpair *next;
};

//Structure for dictionary
typedef struct dictionary {
    int size;           // size of the pointer table 
    int n;              // number of elements stored 
    struct kvpair **table;
}*Dictionary;

//Structure for nodes in the CUDA hashtable
typedef struct node {
	char key[80];
	int index;
	struct node *next;
} Node;

//Function to compute hash value
unsigned long computeHash(const char *s)
{
    unsigned const char *us;
    unsigned long h;
    h = 0;
    for(us = (unsigned const char *) s; *us; us++) {
        h = h * 401 + *us;
    }
    return h;
}

// Function to create empty dictionary
Dictionary createDictionary()
{
    Dictionary d;
    int i;
    d = (Dictionary)malloc(sizeof(*d));
	if(d==NULL){
		printf("dictionary d malloc failed\n");
		exit(0);
	}
    assert(d != 0);
    d->size = 102397;
    d->n = 0;
    d->table = (kvpair **)malloc(sizeof(struct kvpair *) * d->size);
	if(d->table==NULL){
		printf("d->table malloc failed\n");
		exit(0);
	}
    assert(d->table != 0);
    for(i = 0; i < d->size; i++) 
    	d->table[i] = 0;
    return d;
}

// Function to insert a new key-value pair into the dictionary 
void insertDictionary(Dictionary d, const char *key, int val)
{
    struct kvpair *e;
    unsigned long h;
    assert(key);
    assert(val);
    e = (kvpair*)malloc(sizeof(*e));
	if(e==NULL){
		printf("e kvpair malloc failed\n");
		exit(0);
	}
    assert(e);
    e->key = strdup(key);
    e->value = val;
    h = computeHash(key) % d->size;
    e->next = d->table[h];
    d->table[h] = e;
    d->n++;
	return;
}

//Function to search for a key in the dictionary, returns NULL or the Node of the key if found in the dictionary
struct kvpair * searchDictionary(Dictionary d, const char *key)
{
    struct kvpair *e;
    for(e = d->table[computeHash(key) % d->size]; e != 0; e = e->next) {
        if(!strcmp(e->key, key)) {
            return e;
        }
    }
    return NULL;
}

// Function to delete key-value pair in dictionary --this is currently not used
void deleteDictionary(Dictionary d, const char *key)
{
    struct kvpair **prev;          
    struct kvpair *e;              
    for(prev = &(d->table[computeHash(key) % d->size]); 
        *prev != 0; 
        prev = &((*prev)->next)) {
        if(!strcmp((*prev)->key, key)) {
            e = *prev;
            *prev = e->next;
            free(e->key);
            free(e);
            return;
        }
    }
}

//Function to free dictionary
void destroyDictionary(Dictionary d)
{
    int i;
    struct kvpair *e;
    struct kvpair *next;
    for(i = 0; i < d->size; i++) {
        for(e = d->table[i]; e != 0; e = next) {
            next = e->next;
            free(e->key);
            free(e);
        }
    }
    free(d->table);
    free(d);
    return;
}

//Function to compute normal distribution of a value, equivalent to Python's CDF.norm from NVIDIA CUDA samples
__device__ float CND(float d)
{
    const double       A1 = 0.31938153;
    const double       A2 = -0.356563782;
    const double       A3 = 1.781477937;
    const double       A4 = -1.821255978;
    const double       A5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;
    double
    K = 1.0 / (1.0 + 0.2316419 * fabs(d));
    double
    cnd = RSQRT2PI * exp(- 0.5 * d * d) *
          (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));
    if (d > 0)
        cnd = 1.0 - cnd;
    return cnd;
}

//Function to remove specific characters from input string, used to remove carriage returns
void removeChar(char *str, char garbage) {
    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        *dst = *src;
        if (*dst != garbage) dst++;
    }
    *dst = '\0';
	return;
}

//Function used by qsort to sort the records based on number of tokens
int sort(const void* a, const void* b)
{

	 char *ia = strdup(*(const char **)a);
     char *ib = strdup(*(const char **)b);
	 char *split1, *saveptr, *saveptr1;
	 split1 = strtok_r(ia, ";", &saveptr);
	 split1 = strtok_r(NULL, ";", &saveptr);
	 int x = atoi(split1);
	 split1 = strtok_r(ib, ";", &saveptr1);
	 split1 = strtok_r(NULL, ";", &saveptr1);
	 return (x-atoi(split1));
}

/*Kernel function performs random sampling 
It is designed in this way: One block does one sampling and every thread processes one record. In case number of records to be processed exceed 1024, then some threads will take more than one stride. That is, some threads process more than one record during one sampling.
Stages in kernel function are: 
1) Build the dictionary 's_hashtab' is shared memory, for O(1) time lookup of keyword while sampling, first thread will ensure all the keys 'd_r1_dict_keys' are linked in the hashtable.
2) Generate the random numbers and sort them. Every thread will generate a random number and first thread in every block will sort the random numbers using iterative quick sort function. In case number of random numbers required are higher than 1024, then some threads will take more strides to generate the required number of random numbers. 
3) Shared memory initialization for sampling. This is required, as during sampling if keyword is found then we increment the count.
4) Sampling, every block performs one sampling. And, every thread will process atleast one record. That is, thread will extract the keywords/tokens in the record and then finds for the keyword in the dictionary, if found then shared memory is incremented.
5) Copy data to global memory from shared memory for Z-score and P-value calculation.

Note: Shared memory s_hashtab - is the hash table in the shared memory. 
	  Shared memory s_r1_value_list - is the value list in the shared memory.

Arguments passed to kernel function:
	  * d_r2_str - is the list of records to be used for sampling.
	  * d_r1_dict_value - Global values for vector produced from sampling.
	  * sampleTimes - number of samples.
	  * sampleSize - size of the sample.
	  * randomRange - maximum value of each random number.
	  * r1_dict_cnt - number of keywords in dictionary 1.
	  * d_r1_dict_keys - keywords of dictionary 1, to populate hashtable in kernel function.
	  * d_hashtab - global memory hash table.
	  * sampleStrides - maximum number of strides every thread will take for sampling.
	  * threadCount - number of threads per block.
	  * samplesCompleted - number of samplings completed before this kernel launch.
	  * relaunch - to decide whether kernel is launched for the first time or relaunched.
*/
__global__ void deviceDDI(char * d_r2_str, int * d_r1_dict_value, int sampleTimes, int sampleSize, int randomRange, int r1_dict_cnt, Node *d_r1_dict_keys, Node **d_hashtab, int sampleStrides, int threadCount, int samplesCompleted, bool relaunch){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j, k, x = 0, ind1, random;
	unsigned hashval;
 	char str_split[80];
	extern __shared__ int s_r1_value_list[];
    __shared__ Node *s_hashtab[5003];
	//s_r1_value_list array has array for keywords found during sampling.
	//This is used for incrementing array of keywords for sampling.
	//Build the dictionary
	if(threadIdx.x==0){
		if(!relaunch){
			//Initialize the hash table
			for(j=0;j<5003;j++){
            	s_hashtab[j] = 0;
			}
			j=0;
			k=0;
			//Build the dictionary
			for(j=0;j<r1_dict_cnt;j++){
				x=0;
				hashval = 0;
				while(d_r1_dict_keys[j].key[x]!='\0'){
					hashval = ((int)d_r1_dict_keys[j].key[x++] + 401*hashval)%5003;
				}
				d_r1_dict_keys[j].next = (s_hashtab[hashval]==0?NULL:s_hashtab[hashval]);
				s_hashtab[hashval] = &d_r1_dict_keys[j];
			}
 			for(j=0;j<5003;j++){
            	d_hashtab[j] = s_hashtab[j];
			}
		}
		else{
			for(j=0;j<5003;j++){
            	s_hashtab[j] = d_hashtab[j];
			}
		}
		//Initialize the shared memory
		for(j=0;j<(r1_dict_cnt);j++){
 			s_r1_value_list[j] = 0;
 		}
	}	
	__syncthreads();		
	//Generate the random numbers
	curandState_t state;
 	curand_init(clock64(), i, 0, &state);		
 	//Sampling
 	for(j=0;(j<sampleStrides)&&((threadIdx.x+(j*threadCount))<sampleSize);j++){
 		//char str_split[80];
 		x=0;
 		hashval = 0;
		random = curand(&state)%randomRange;
		for(k=random*1000;k<random*1000+1000;k++){
			if(d_r2_str[k] == '^')
				break;
 			if(d_r2_str[k] != '~'){
 		 		str_split[x++] = d_r2_str[k];
 		 		hashval = ((int)d_r2_str[k] + 401*hashval)%5003;
 			}
 			else{
 				//str_split[x] = '\0';
				Node *np = s_hashtab[hashval];
				//np = s_hashtab[hashval];
 				while((np!= NULL)&&(np!=0)){
 		 			ind1 = 0;
 		 			while((np->key[ind1] != '\0')&&(ind1<x)){
 		 				if(np->key[ind1] == str_split[ind1])
 		 					ind1++;
 		 				else
 		 					break;
 		 			}
 		 			if((np->key[ind1] == '\0')&&(ind1==x)){
						atomicAdd(&s_r1_value_list[(np->index)],1);
 		 				break;
 		 			}
 		 			if(np->next == NULL||np->next==0)
 		 				break;
 		 			np = np->next;
 				}
 				x=0;
 				hashval = 0;
 			}
 		}
 	}
	__syncthreads(); 

	//Copy to global memory from shared memory
	if(threadIdx.x==0){
		for(j=0;j<(r1_dict_cnt);j++){
 			d_r1_dict_value[(j)*(sampleTimes+1)+blockIdx.x+1+samplesCompleted] = s_r1_value_list[j];
 		}
	}
}

/*	deviceZP: Function to calculate P-value
	One thread operates on one array to compute Z-score and P-value

	Arguments passed to the kernel function:
		* d_r1_dict_value -  array holding the 
		* d_z_score - array to hold Z-scores.
	  	* d_p_value - array to hold P-values.
*/
__global__ void deviceZP(int * d_r1_dict_value, int sampleTimes,int r1_dict_cnt, float * d_z_score, float * d_p_value){
    int x, j,  i = blockDim.x * blockIdx.x + threadIdx.x;
    float mean =0,sd;
    if(i<r1_dict_cnt){
            x = 0;
            sd = 0;
            for(j=1;j<=sampleTimes;j++){
                x += d_r1_dict_value[(i*(sampleTimes+1))+j];
            }
            mean = x/(sampleTimes);
            for(j=1;j<=sampleTimes;j++){
                sd += (d_r1_dict_value[(i*(sampleTimes+1))+j]-mean)*(d_r1_dict_value[(i*(sampleTimes+1))+j]-mean);
            }
            sd = sqrt(sd/(sampleTimes));

            if(fabs(sd)>pow(10.0,-7))
                d_z_score[i] = (d_r1_dict_value[i*(sampleTimes+1)] - mean)/sd;
            else{
                if(d_r1_dict_value[i*(sampleTimes+1)] != (int)mean)
                    d_z_score[i] = d_r1_dict_value[i*(sampleTimes+1)]*100;
                else
                    d_z_score[i] = -100;
            }
            d_p_value[i] = 1-CND(d_z_score[i]);
    }
}

//Function to partition records while sorting based on Z-score, called by quickSort
int partition( float a[], int index[], int l, int r) {
   int i, j, t;
   float temp;
   float pivot = a[l];
   i = l; 
   j = r+1;		
   while( 1)
   {
   	do ++i; 
   	while( a[i] >= pivot && i <= r );
   	do --j; 
   	while( a[j] < pivot );
   	if( i >= j ) 
   		break;
   	temp = a[i]; 
   	a[i] = a[j]; 
   	a[j] = temp;
   	t = index[i];
   	index[i] = index[j];
   	index[j] = t;
   }
   temp = a[l]; 
   a[l] = a[j]; 
   a[j] = temp;
   t = index[l];
   index[l] = index[j];
   index[j] = t;
   return j;
}

//Function to quicksort the records based on Z-score
void quickSort(float a[], int index[], int l, int r)
{
   int j;
   if( l < r ) 
   {
       j = partition( a, index, l, r);
       quickSort( a, index, l, j-1);
       quickSort( a, index, j+1, r);
   }	
   return;
}


int main(int argc, char *argv[])
{
	if(argc!=7){
		printf("\nIncorrect arguments passed, Please pass <Compounds with interactions>, <Compounds without interactions>, <Proteins with interactions>, <Proteins without interactions>, <PMID Substances>, <para.txt>, <Drug Name> as arguments\n");
		exit(1);
	}
	FILE *inp_r1, *inp_r2, *inp_para, *op1, *op2;
	char *split0,*split1, *saveptr, *saveptr1, *saveptr2, *inp2_list[100000];
	char filename1[100], filename2[100], cutoffstr[20], pvaluestr[20], str1[10000], rmode[2] = "r";
	size_t len = 0;
	Dictionary d_cinp1;
	int cutoff, sampleTimes, i=0, j=0, k=0, r1_cnt, r2_cnt, r1_dict_cnt, threadCount, sampleStrides;
	float p_value, elapsedTime, totalTime=0;
	cudaEvent_t start, stop;

	printf("Drug name = %s\n",argv[6]);
	printf("Read input files\n");
	d_cinp1 = createDictionary();			

	//Read the parameters from para.txt - 4th argument
	inp_para = fopen(argv[5],rmode);
	if (inp_para == NULL) 
	{
		fprintf(stderr, "Can't open input file %s!\n", argv[5]);
		exit(1);
	}	
	while(1)
	{	
		fscanf(inp_para,"%[^\n]%*c", str1);
		if(feof(inp_para)) break;
		split0 = strtok_r(str1, "\t", &saveptr);
		split1 = strtok_r(NULL, "\t", &saveptr);
	 	removeChar(split0,'\r');
		removeChar(split1,'\r');
		if( strcmp(split0,"sampleTimes") == 0)
		{
			char temp[20];
			strcpy(temp, split1);
			sampleTimes = atoi(temp);	
		}
		
		 else if( strcmp(split0,"cutoff") == 0)
		{
			char temp[20];
			strcpy(temp, split1);
			strcpy(cutoffstr,temp);
			cutoff = atoi(temp);	
		}
		else if( strcmp(split0,"p_value") == 0)
		{
			char temp[20];
			strcpy(temp, split1);
			strcpy(pvaluestr,temp);
			p_value = atof(temp);	
		}
	}
	fclose(inp_para);
	printf("Number of Samples = %d\n",sampleTimes);
	if(sampleTimes <=0){
		printf("Incorrect number of samples specified = %d, value of atleast 1 is expected\n", sampleTimes);
		exit(0);
	}
	// Reading the dictionary of compounds of result 1 - 1st argument
	// Create and populate dictionary 'd_cinp1' while reading the records
	inp_r1 = fopen(argv[1], rmode);
	if (inp_r1 == NULL) 
	{
		fprintf(stderr, "Can't open input file %s!\n", argv[1]);
		exit(1);
	}	
	
	r1_cnt = 0;
	r1_dict_cnt = 0;

	while(1){
		fscanf(inp_r1, "%[^\n]%*c", str1);
		if( feof(inp_r1)) break;
		removeChar(str1,'\r');
		r1_cnt++;
		len = strlen(str1);
		for(i=0;(i<len);i++){
			char *newstr = (char*)malloc(len+1);
			if(newstr==NULL){
				printf("malloc to newstr failed\n");
				exit(0);
			}
			j=0;
			while(str1[i] != '~'){
				newstr[j++] = str1[i++];
			}	
			newstr[j] = '\0';
			struct kvpair * e = searchDictionary(d_cinp1,newstr);
			if(e!=NULL){
				e->value++;
			}
			else{
				insertDictionary(d_cinp1,newstr,1);
				r1_dict_cnt++;
			}
			free(newstr);
		}
	}	
	fclose(inp_r1);

	// Reading the list of result 2- 2nd argument
	inp_r2 = fopen(argv[2], rmode);
	if (inp_r2 == NULL) 
	{
		fprintf(stderr, "Can't open input file %s!\n", argv[2]);
		exit(1);
	}	
	r2_cnt = 0;
	while (1) 
	{
		fscanf(inp_r2, "%[^\n]%*c", str1);
		if( feof(inp_r2)) break;
		removeChar(str1,'\r');
		inp2_list[r2_cnt] = (char*)malloc(strlen(str1)+1);
		if(inp2_list[r2_cnt]==NULL){
			printf("malloc to inp2_list[r2_cnt] failed\n");
			exit(0);
		}
		strcpy(inp2_list[r2_cnt++],str1);		
	}
	fclose(inp_r2);
	printf("Input files read completed\n");
	printf("Sample size = %d\n", r1_cnt);
	//Sort inp2_list based on the number of tokens or length
	qsort(inp2_list,r2_cnt,sizeof(char *), sort);

	printf("Pre-process records for kernel launch\n");
	
	cudaSetDevice(0);
	//populate value list for dictionary 1
	cudaError_t err = cudaSuccess;
	int * r1_dict_value;
	//pinned memory for optimized usage of memory transfer bandwidth
	err = cudaMallocHost((void**)&r1_dict_value, sizeof(int)*r1_dict_cnt*(sampleTimes+1));
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate r1_dict_value host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	//Device value list for dictionary 1
	int *d_r1_dict_value;
	Node * r1_dict_keys = (Node*)malloc(sizeof(Node)*r1_dict_cnt);
	if(r1_dict_keys==NULL){
		printf("malloc to r1_dict_keys failed\n");
		exit(0);
	}
	j=0;
	k=0;
	for(i=0;i<d_cinp1->size;i++)
    {
    	if(d_cinp1->table[i]!=0){
			while(1)
			{
				strcpy(r1_dict_keys[j].key,d_cinp1->table[i]->key);
				r1_dict_keys[j].index = j;
				r1_dict_keys[j].next = NULL;
				r1_dict_value[j*(sampleTimes+1)] = d_cinp1->table[i]->value;
				j++;
				if(d_cinp1->table[i]->next!= NULL)
					d_cinp1->table[i] =  d_cinp1->table[i]->next;
				else
					break;
			}
		}
	}
	destroyDictionary(d_cinp1);

	//Strip off the number of tokens from every record in list 2.
	for(i=0;i<r2_cnt;i++){
		split0 = strtok_r(inp2_list[i], ";", &saveptr1);
	}

	//Process the records for shipping to kernel
	char * temp1 = (char*) malloc(1000*r2_cnt*sizeof(char));
	if(temp1==NULL){
		printf("temp1 malloc failed\n");
		exit(0);
	}
	char * d_r2_str;
	j=0;
	for(i=0;i<r2_cnt;i++){
		for(k=0;k<1000;k++){
			if(k<strlen(inp2_list[i])){
				temp1[j++] = inp2_list[i][k];
			}
			else
				temp1[j++] = '^';
		}
		free(inp2_list[i]);
	}

	//char* d_r1_dict_list;
	//Allocate global memory for dictionary 1 keywords
	/*err = cudaMalloc((void **)&d_r1_dict_list,80*sizeof(char)*r1_dict_cnt);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_r1_dict_list (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r1_dict_list,r1_dict_list,80*sizeof(char)*r1_dict_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r1_dict_list (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}*/
	
	//Allocate global memory for input list 2 records
	err = cudaMalloc((void **)&d_r2_str,1000*sizeof(char)*r2_cnt);
	if(err != cudaSuccess){
			fprintf(stderr,"Failed to allocate device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r2_str,temp1,1000*sizeof(char)*r2_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	//Allocate global memory for dictionary 1 value list 
	err = cudaMalloc((void **)&d_r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r1_dict_value,r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	Node * d_r1_dict_keys = NULL;
	err = cudaMalloc((void **)&d_r1_dict_keys,sizeof(Node)*r1_dict_cnt);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	err = cudaMemcpy(d_r1_dict_keys,r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	Node **hashtab;
	Node **d_hashtab;
	err = cudaMalloc((void **)&d_hashtab,sizeof(Node*)*5003);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_hashtab (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	hashtab = (Node**)malloc(sizeof(Node*)*5003);
	if(hashtab==NULL){
		printf("hashtab malloc failed\n");
		exit(0);
	}
	printf("Sampling for compounds begin\n");
	for(i=0;i<sampleTimes;i=i+256){

		threadCount = (r1_cnt>1024)?1024:r1_cnt;
		//sampleStrides: maximum number of strides every thread needs to take for sampling
		sampleStrides = ceil(r1_cnt/threadCount);
		printf("Kernel deviceDDI launched with %d blocks of %d threads each\n", (sampleTimes-i)>256?256:(sampleTimes-i), threadCount);

		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord( start, 0 );

		deviceDDI<<<(sampleTimes-i)>256?256:(sampleTimes-i), threadCount, (sizeof(int)*(r1_dict_cnt))>>>(d_r2_str, d_r1_dict_value, sampleTimes, r1_cnt, r2_cnt, r1_dict_cnt, d_r1_dict_keys, d_hashtab, sampleStrides, threadCount,i,(i==0)?false:true);
		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );
		cudaEventElapsedTime( &elapsedTime, start, stop );
		cudaEventDestroy( start );
		cudaEventDestroy( stop );
		err = cudaDeviceSynchronize();
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to launch deviceDDI kernel device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}
		err = cudaGetLastError();
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to launch deviceDDI kernel device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}

		err = cudaMemcpy(r1_dict_value,d_r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_value from device to Host(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}

		err = cudaMemcpy(d_r1_dict_value, r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_value from host to device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}
		
		err = cudaMemcpy(r1_dict_keys,d_r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_keys to host (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_r1_dict_keys,r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(hashtab,d_hashtab,sizeof(Node*)*5003,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_hashtab to host (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_hashtab,hashtab,sizeof(Node*)*5003,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_hashtab (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_r2_str,temp1,1000*sizeof(char)*r2_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

	}
	printf("Sampling for compounds completed\n");

    err = cudaFree(d_r2_str);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to free from device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    err = cudaFree(d_r1_dict_keys);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to free from device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }
    free(temp1);

    float *d_z_score;
    float *z_score_arr;
    //Allocate array for Z-score, pinned memory for optimized usage of memory transfer bandwidth
    err = cudaMallocHost((void**)&z_score_arr, sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate z-score host (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    err = cudaMalloc((void **)&d_z_score,sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate device d_z_score (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    float *d_p_value;
    float *p_value_arr;
    //Allocate array for P-value, pinned memory for optimized usage of memory transfer bandwidth
    err = cudaMallocHost((void**)&p_value_arr, sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate p-value host (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    err = cudaMalloc((void **)&d_p_value,sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate device d_p_value (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }
	
	printf("Z-score and P-value calculation for Compounds begin\n");
	printf("Kernel deviceZP launched with %d blocks of %d threads each\n", (int)ceil(r1_dict_cnt/256.0), 256);
    totalTime += elapsedTime;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord( start, 0 );

    deviceZP<<<ceil(r1_dict_cnt/256.0), 256>>>(d_r1_dict_value, sampleTimes, r1_dict_cnt, d_z_score, d_p_value);
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );

	cudaEventElapsedTime( &elapsedTime, start, stop );
	cudaEventDestroy( start );
	cudaEventDestroy( stop );
	err = cudaDeviceSynchronize();
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to launch deviceZP kernel device(error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE); 
	}
	err = cudaGetLastError();
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to launch deviceZP kernel device(error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	printf("Z-score, P-value calculation completed\n");
	err = cudaMemcpy(z_score_arr,d_z_score,sizeof(float)*r1_dict_cnt,cudaMemcpyDeviceToHost);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy from z-score device to host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(p_value_arr,d_p_value,sizeof(float)*r1_dict_cnt,cudaMemcpyDeviceToHost);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy from p-value device to host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	
	err = cudaFree(d_r1_dict_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	err = cudaFree(d_p_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_p_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaFree(d_z_score);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_z_score (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	
	//Initialize the sortedIndex array, as sortedIndex will have the values sorted with quickSort based on descending order of Z-score
	//After sorting sortedIndex contains the new index of Z-score.
	int * sortedIndex = (int*) malloc(sizeof(int)*r1_dict_cnt);
	if(sortedIndex == NULL){
		printf("malloc error for sortedIndex\n");
	}
	for(i=0;i<r1_dict_cnt;i++){
		sortedIndex[i] = i;
	}

	quickSort(z_score_arr, sortedIndex, 0, r1_dict_cnt-1);
	
	printf("Write extracted compounds to output files\n");
	//Write to output files
	strcpy(filename1, argv[6]);
	strcat(filename1, "_temp_result1_Substance_compounds_cutoff_");
	strcat(filename1,cutoffstr);
	strcat(filename1,"_p_");
	strcat(filename1,pvaluestr);
	strcat(filename1,".txt");
	strcpy(filename2, argv[6]);
	strcat(filename2, "_temp_result1_Substance_compounds_cutoff_");
	strcat(filename2,cutoffstr);
	strcat(filename2,".txt");
	
	
	op1 = fopen(filename1, "w");
	fprintf(op1,"Term Pair\tMeSHID\tDistribution\tZ-Score\tP-value\n");
	op2 = fopen(filename2, "w");
	fprintf(op2,"Term Pair\tMeSHID\tDistribution\tZ-Score\tP-value\n");
	k=0;
	for(i=0;i<r1_dict_cnt;i++){
		if((r1_dict_value[(sortedIndex[i]*(sampleTimes+1))]>=cutoff) && (p_value_arr[sortedIndex[i]]<=p_value)){
			fprintf(op1,"%s;%s\t[",argv[6],r1_dict_keys[sortedIndex[i]].key);
			for(j=0;j<sampleTimes;j++)
				fprintf(op1,"%d, ",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j]);
			fprintf(op1,"%d]\t%f\t%f\n",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j], z_score_arr[i],p_value_arr[sortedIndex[i]]);
		}
		if((r1_dict_value[(sortedIndex[i]*(sampleTimes+1))]>=cutoff) && (p_value_arr[sortedIndex[i]]<=1.0)){
			fprintf(op2,"%s;%s\t[",argv[6],r1_dict_keys[sortedIndex[i]].key);
			for(j=0;j<=sampleTimes;j++)
				fprintf(op2,"%d, ",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j]);
			fprintf(op2,"%d]\t%f\t%f\n",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j], z_score_arr[i],p_value_arr[sortedIndex[i]]);
		}
		k++;
	}
	fclose(op1);
	fclose(op2);
	printf("Compounds output files written\n");

	free(r1_dict_keys);

	err = cudaFreeHost(p_value_arr);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host p_value_arr (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaFreeHost(z_score_arr);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host z_score_arr (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	
	err = cudaFreeHost(r1_dict_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	free(sortedIndex);
	d_cinp1 = createDictionary();
	printf("Processing proteins\n");
	printf("Read input files\n");
	// Reading the dictionary of proteins of result 1 - 3rd argument
	inp_r1 = fopen(argv[3], rmode);
	if (inp_r1 == NULL) 
	{
		fprintf(stderr, "Can't open input file %s!\n", argv[3]);
		exit(1);
	}	
	
	r1_dict_cnt = 0;

	while(1){
		fscanf(inp_r1, "%[^\n]%*c", str1);
		if( feof(inp_r1)) break;
		removeChar(str1,'\r');

		len = strlen(str1);
		for(i=0;(i<len);i++){
			char *newstr = (char*)malloc(len+1);
			if(newstr==NULL){
				printf("newstr malloc failed\n");
				exit(0);
			}
			j=0;
			while(str1[i] != '~'){
				newstr[j++] = str1[i++];
			}	
			newstr[j] = '\0';
			struct kvpair * e = searchDictionary(d_cinp1,newstr);
			if(e!=NULL){
				e->value++;
			}
			else{
				insertDictionary(d_cinp1,newstr,1);
				r1_dict_cnt++;
			}
			free(newstr);
		}
	}	
	fclose(inp_r1);

	// Reading the list of result 2- 4th argument
	inp_r2 = fopen(argv[4], rmode);
	if (inp_r2 == NULL) 
	{
		fprintf(stderr, "Can't open input file %s!\n", argv[4]);
		exit(1);
	}	
	r2_cnt = 0;
	while (1) 
	{
		fscanf(inp_r2, "%[^\n]%*c", str1);
		if( feof(inp_r2)) break;
		removeChar(str1,'\r');
		inp2_list[r2_cnt] = (char*)malloc(strlen(str1)+1);
		if(inp2_list[r2_cnt]==NULL){
			printf("inp2_list[r2_cnt] malloc failed\n");
			exit(0);
		}
		strcpy(inp2_list[r2_cnt++],str1);		
	}
	fclose(inp_r2);
	printf("Input files read completed\n");
	//Sort inp2_list based on the number of tokens
	qsort(inp2_list,r2_cnt,sizeof(char *), sort);
	printf("Pre-process records for kernel launch\n");
	//pinned memory for optimized usage of memory transfer bandwidth
	err = cudaMallocHost((void**)&r1_dict_value, sizeof(int)*r1_dict_cnt*(sampleTimes+1));
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate r1_dict_value host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	//Process the input data for shipping
	free(r1_dict_keys);
	r1_dict_keys = (Node*)malloc(sizeof(Node)*r1_dict_cnt);
	if(r1_dict_keys==NULL){
		printf("r1_dict_keys malloc failed\n");
		exit(0);
	}
	j=0;
	k=0;
	for(i=0;i<d_cinp1->size;i++)
    {
    	if(d_cinp1->table[i]!=0){
			while(1)
			{
				strcpy(r1_dict_keys[j].key,d_cinp1->table[i]->key);
				r1_dict_keys[j].index = j;
				r1_dict_keys[j].next = NULL;
				r1_dict_value[j*(sampleTimes+1)] = d_cinp1->table[i]->value;
				j++;
				if(d_cinp1->table[i]->next!= NULL)
					d_cinp1->table[i] =  d_cinp1->table[i]->next;
				else
					break;
			}
		}
	}
	destroyDictionary(d_cinp1);
	//Strip off the number of tokens from list 2 records
	for(i=0;i<r2_cnt;i++){
		split0 = strtok_r(inp2_list[i], ";", &saveptr2);
	}

	temp1 = (char*) malloc(1000*r2_cnt*sizeof(char));
	if(temp1 == NULL){
		printf("temp1 malloc failed\n");
		exit(0);
	}
	j=0;
	for(i=0;i<r2_cnt;i++){
		for(k=0;k<1000;k++){
			if(k<strlen(inp2_list[i])){
				temp1[j++] = inp2_list[i][k];
			}
			else
				temp1[j++] = '^';
		}
		free(inp2_list[i]);
	}

	err = cudaMalloc((void **)&d_r1_dict_keys,sizeof(Node)*r1_dict_cnt);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r1_dict_keys,r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	//Allocate global memory for list 2 records
	err = cudaMalloc((void **)&d_r2_str,1000*sizeof(char)*r2_cnt);
	if(err != cudaSuccess){
			fprintf(stderr,"Failed to allocate device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r2_str,temp1,1000*sizeof(char)*r2_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	//Allocate global memory for index of dictionary 1 value list
	err = cudaMalloc((void **)&d_r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to allocate device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(d_r1_dict_value,r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyHostToDevice);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}


	printf("Sampling for proteins begin\n");
	for(i=0;i<sampleTimes;i=i+256){
		threadCount = (r1_cnt>1024)?1024:r1_cnt;
		//sampleStrides: maximum number of strides every thread need to take for sampling
		sampleStrides = ceil(r1_cnt/threadCount);
		printf("Kernel deviceDDI launched with %d blocks of %d threads each\n", (sampleTimes-i)>256?256:(sampleTimes-i), threadCount);

		totalTime += elapsedTime;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord( start, 0 );
	
		deviceDDI<<<(sampleTimes-i)>256?256:(sampleTimes-i), threadCount, (sizeof(int)*(r1_dict_cnt))>>>( d_r2_str, d_r1_dict_value, sampleTimes, r1_cnt, r2_cnt, r1_dict_cnt, d_r1_dict_keys, d_hashtab, sampleStrides, threadCount,i,(i==0)?false:true);

		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );
		cudaEventElapsedTime( &elapsedTime, start, stop );
		cudaEventDestroy( start );
		cudaEventDestroy( stop );
    	err = cudaDeviceSynchronize();
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to launch deviceDDI kernel device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}
		err = cudaGetLastError();
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to launch deviceDDI kernel device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}
		
		err = cudaMemcpy(r1_dict_keys,d_r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_keys to host (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_r1_dict_keys,r1_dict_keys,sizeof(Node)*r1_dict_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(r1_dict_value,d_r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_value from device to Host(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}

		err = cudaMemcpy(d_r1_dict_value, r1_dict_value,sizeof(int)*(sampleTimes+1)*r1_dict_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r1_dict_value from host to device(error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE); 
		}

		err = cudaMemcpy(hashtab,d_hashtab,sizeof(Node*)*5003,cudaMemcpyDeviceToHost);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_hashtab to host (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_hashtab,hashtab,sizeof(Node*)*5003,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_hashtab (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}

		err = cudaMemcpy(d_r2_str,temp1,1000*sizeof(char)*r2_cnt,cudaMemcpyHostToDevice);
		if(err != cudaSuccess){
			fprintf(stderr,"Failed to copy device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);	
		}
	}
	printf("Sampling completed\n");
    err = cudaFree(d_r1_dict_keys);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to free from device d_r1_dict_keys (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }
    err = cudaFree(d_r2_str);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to free from device d_r2_str (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }
	free(temp1);
    //Allocate Z-score array pinned memory for optimized usage of memory transfer bandwidth
    err = cudaMallocHost((void**)&z_score_arr, sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate z-score host (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

	//float * d_z_score_p;
    err = cudaMalloc((void **)&d_z_score,sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate device d_z_score (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    //Allocate P-value array pinned memory for optimized usage of memory transfer bandwidth
    err = cudaMallocHost((void**)&p_value_arr, sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate p-value host (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

    err = cudaMalloc((void **)&d_p_value,sizeof(float)*r1_dict_cnt);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to allocate device d_p_value (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }

	printf("Z-score, P-value calculation for proteins begin\n");
	printf("Kernel deviceZP launch with %d blocks of %d threads each\n", (int)ceil(r1_dict_cnt/256.0),256);
    totalTime += elapsedTime;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord( start, 0 );

    deviceZP<<<ceil(r1_dict_cnt/256.0), 256>>>(d_r1_dict_value, sampleTimes, r1_dict_cnt, d_z_score, d_p_value);

	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &elapsedTime, start, stop );
	cudaEventDestroy( start );
	cudaEventDestroy( stop );
	err = cudaDeviceSynchronize();
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to launch deviceZP kernel device(error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE); 
	}
	err = cudaGetLastError();
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to launch deviceZP kernel device(error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	printf("Z-score, P-value calculation for proteins completed\n");
	printf( "\n******** Total Running Time of Kernel = %0.5f seconds ******* \n", (elapsedTime+totalTime)/1000);
	printf("Copy output data to host memory\n");
	err = cudaMemcpy(p_value_arr,d_p_value,sizeof(float)*r1_dict_cnt,cudaMemcpyDeviceToHost);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy from p-value device to host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	err = cudaFree(d_r1_dict_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	err = cudaFree(d_p_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_p_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaMemcpy(z_score_arr,d_z_score,sizeof(float)*r1_dict_cnt,cudaMemcpyDeviceToHost);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to copy from z-score device to host (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}
	
	err = cudaFree(d_z_score);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free from device d_z_score (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	//Initialize sortedIndex, this will hold correct index of the dictionary 1 records after sorting based on descending order of Z-score
	
	sortedIndex = (int*)malloc(sizeof(int)*r1_dict_cnt);
	if(sortedIndex == NULL){
		printf("sortedIndex malloc error\n");
	}
	for(i=0;i<r1_dict_cnt;i++){
		sortedIndex[i] = i;
	}
	//Sort the array based on descending order of Z-score
	quickSort(z_score_arr, sortedIndex, 0, r1_dict_cnt-1);
	//Write to output files
	strcpy(filename1, argv[6]);
	strcat(filename1, "_temp_result1_Substance_proteins_cutoff_");
	strcat(filename1,cutoffstr);
	strcat(filename1,"_p_");
	strcat(filename1,pvaluestr);
	strcat(filename1,".txt");
	strcpy(filename2, argv[6]);
	strcat(filename2, "_temp_result1_Substance_proteins_cutoff_");
	strcat(filename2,cutoffstr);
	strcat(filename2,".txt");

	printf("Write extracted proteins to output files\n");
	op1 = fopen(filename1, "w");
	fprintf(op1,"Term Pair\tMeSHID\tDistribution\tZ-Score\tP-value\n");
	op2 = fopen(filename2, "w");
	fprintf(op2,"Term Pair\tMeSHID\tDistribution\tZ-Score\tP-value\n");
	k=0;
	for(i=0;i<r1_dict_cnt;i++){
		
		if((r1_dict_value[(sortedIndex[i]*(sampleTimes+1))]>=cutoff) && (p_value_arr[sortedIndex[i]]<=p_value)){
			fprintf(op1,"%s;%s\t[",argv[6],r1_dict_keys[sortedIndex[i]].key);
			for(j=0;j<sampleTimes;j++)
				fprintf(op1,"%d, ",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j]);
			fprintf(op1,"%d]\t%f\t%f\n",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j], z_score_arr[i],p_value_arr[sortedIndex[i]]);
		}
		if((r1_dict_value[(sortedIndex[i]*(sampleTimes+1))]>=cutoff) && (p_value_arr[sortedIndex[i]]<=1.0)){
			fprintf(op2,"%s;%s\t[",argv[6],r1_dict_keys[sortedIndex[i]].key);
			for(j=0;j<=sampleTimes;j++)
				fprintf(op2,"%d, ",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j]);
			fprintf(op2,"%d]\t%f\t%f\n",r1_dict_value[(sortedIndex[i]*(sampleTimes+1))+j], z_score_arr[i],p_value_arr[sortedIndex[i]]);
		}
		k++;
	}
	fclose(op1);
	fclose(op2);

	printf("Processing completed\n");
	free(r1_dict_keys);
	err = cudaFreeHost(p_value_arr);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host p_value_arr (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaFreeHost(z_score_arr);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host z_score_arr (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	err = cudaFreeHost(r1_dict_value);
	if(err != cudaSuccess){
		fprintf(stderr,"Failed to free pinned host r1_dict_value (error code %s) !\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);	
	}

	free(sortedIndex);
	free(hashtab);

	err = cudaFree(d_hashtab);
    if(err != cudaSuccess){
        fprintf(stderr,"Failed to free from device d_hashtab (error code %s) !\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE); 
    }
}
//End of program