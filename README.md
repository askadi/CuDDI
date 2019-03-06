CuDDI
================================================

Instructions for CuDDI:<br />
Module 1 is performed in step 1, Module 2 is performed in step2 and Module 3 is performed in step3<br />
The detailed commands to compile and execute source code are as follows:<br />

Here, we use Simvastatin as an example. The related files are located in the directory of example<br />
1. Run the python script:<br />
	`python step_1_download_Substance.py drugNameList_original.txt para.txt`<br />
	In this step, please make sure that the Biopython package is installed.<br />
	Pay attention: drugNameList_original.txt and para.txt must be manually configured in advance.
	There are five output files which are:<br />
	simvastatin_PMID_Substances.txt<br />
	simvastatin_result_1_compounds.txt (DDI-related)<br />
	simvastatin_result_1_proteins.txt (DDI-related)<br />
	simvastatin_result_2_compounds.txt (DDI-unrelated)<br />
	simvastatin_result_2_proteins.txt (DDI-unrelated)<br />

2. <br />
	2.1 Random sampling and calculate the p value.<br />
	2.1.1 CUDA-based implementation for Module 2.<br />
		Compile the source code in the computer which has NVIDIA GPU by using the command:<br />
		nvcc step_2_1_randomSampling_pValue_output.cu -o step_2_1_randomSampling_pValue_output -Wno-deprecated-gpu-targets
		Run the executable file:<br />
		`step_2_1_randomSampling_pValue_output simvastatin_result_1_compounds.txt simvastatin_result_2_compounds.txt simvastatin_result_1_proteins.txt simvastatin_result_2_proteins.txt para.txt "Simvastatin"`<br />

	2.1.2 Python-based implementation for Module 2.<br />
		In this step, please make sure that the scipy and numpy packages are installed.<br />
		Run the python script:<br />
		`python step_2_1_randomSampling_pValue_output.py simvastatin_result_1_compounds.txt simvastatin_result_2_compounds.txt simvastatin_result_1_proteins.txt simvastatin_result_2_proteins.txt para.txt "Simvastatin"`<br />
	
	There are four output files containing the compounds and proteins information for Ibuprofen related DDI which are:<br />
	Simvastatin_temp_result1_Substance_compounds_cutoff_4.txt<br />
	Simvastatin_temp_result1_Substance_compounds_cutoff_4_p_0.05.txt<br />
	Simvastatin_temp_result1_Substance_proteins_cutoff_4.txt<br />
	Simvastatin_temp_result1_Substance_proteins_cutoff_4_p_0.05.txt<br />
	
	2.2 Calculate the Benjamini-Hochberg adjusted p value.<br />
		Run the python scripts:<br />
		`python step_2_2_adjusted_pValue.py Simvastatin_temp_result1_Substance_compounds_cutoff_4_p_0.05.txt`<br />
		`python step_2_2_adjusted_pValue.py Simvastatin_temp_result1_Substance_proteins_cutoff_4_p_0.05.txt`<br />
	
		There are two output files containing the compounds and proteins information for Ibuprofen related DDI ranked by Benjamini-Hochberg adjusted p value which are:<br />
		Simvastatin_temp_result1_Substance_compounds_cutoff_4_adjustedPvalue_0.05.txt<br />
		Simvastatin_temp_result1_Substance_proteins_cutoff_4_adjustedPvalue_0.05.txt<br />

3. Run the python scripts:<br />
	`python step_3_network_file.py Simvastatin_temp_result1_Substance_compounds_cutoff_4_adjustedPvalue_0.05.txt Simvastatin_temp_result1_Substance_proteins_cutoff_4_adjustedPvalue_0.05.txt simvastatin_PMID_Substances.txt 4 0.05`<br />
	In this step,<br />
	firstly, generate the vertices and substance-PMID matrix file.<br />
	secondly, Rscript will be called to run Plot_network.R. The social network where the substances are nodes is generated. The color of the target compounds is blue, the color of the interactive compounds is yellow, and the color of the proteins related to DDIs is red.<br />
	Pay attention: termList_to_keep.txt must be manually created so that the terms in termList_to_keep.txt will be eventually used to build the social network.<br />