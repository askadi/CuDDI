import os
import sys
import random
from scipy.stats import norm
import numpy
import math
from time import time

# Calculate P value , the input must be a list of int
# Return a new list with the last two element is z_score and p_value
def calPValue(list):
	outlist = []
	temparray = numpy.array(list)
	mean = numpy.mean(temparray[1:])
	std = numpy.std(temparray[1:])
	n = numpy.shape(temparray[1:])[0]
	if std != 0:
		z_score_temp = (temparray[0]-mean)/(std)
	else:
		if temparray[0] != mean:
			z_score_temp = temparray[0]*100
		else:
			z_score_temp = -100
	p_value_temp = 1- norm.cdf(z_score_temp)
	outlist = list + [z_score_temp] + [p_value_temp]
	return outlist

# Sort the calculated results based on the z_score, from largest to smallest  
# Return the sorted list
def sortResult(filtered_drug_term_pair_dic):
	sorted_drug_tem_pair_list = []
	inputtempList = []
	for temppair in filtered_drug_term_pair_dic.keys():
		inputtempList.append([filtered_drug_term_pair_dic[temppair][-2],temppair])
	inputtempList.sort(reverse = True)
	for subitem in inputtempList:
		sorted_drug_tem_pair_list.append(subitem[1])
	return sorted_drug_tem_pair_list

def randomSample(inputList, sampleSize, sampleRange):
	d_cr2 = {}
	for i in range(sampleSize):
		index = random.randint(0,sampleRange-1)
		for key in inputList[index]:
			if key not in d_cr2:
				d_cr2.update({key:1})
			else:
				d_cr2[key] = d_cr2[key]+1
	return d_cr2

def main(arg = sys.argv):
	if len(arg) != 7:
		print >> sys.stderr,"Please pass all the 7 arguments in the format <Drug result 1 compounds>, <Drug result 2 compounds>, <Drug result 1 proteins>, <Drug result 2 proteins>, <para.txt> <Drug name>"
		sys.exit(1)
	inp_r1 = open(arg[1],"r")
	d_cinp1 = {}
	sampleSize = 0
	for line in inp_r1:
		sampleSize += 1
		keywords = 	line.split("~")
		#print keywords
		for key in keywords:
			if key != '\n':
				if key not in d_cinp1:
					d_cinp1.update({key:[1]})
				else:
					d_cinp1[key][0] = d_cinp1[key][0]+1
	inp_r1.close()

	print '\nsample size = ', sampleSize
	inp_r2 = open(arg[2],"r")
	list_r2 = [[]]
	sampleRange = 0
	for line in inp_r2:
		sampleRange += 1
		keywords = line.split(";")
		tokens = keywords[0].split("~")
		temp = []
		for key in tokens:
			if key != '\n' and key != '':
				temp.append(key)
		list_r2.append(temp)
	inp_r2.close()
	#print list_r2
	print "\nsample range = ", sampleRange
	inp_para = open(arg[5],"r")
	for line in inp_para:
		content = line.strip().split("\t")
		if content[0] == "sampleTimes":
			sampleTimes = int(content[1])
		elif content[0] == "cutoff":
			cutoff = int(content[1])
		elif content[0] == "p_value":
			p_value = float(content[1])
		elif content[0] == "z_score":
			corresponding_z_score = float(content[1])
		elif content[0] == "endDate":
			endDate = content[1]
	inp_para.close()
	print "\nNumber of samples = ", sampleTimes
	print "\nDrug name = ", arg[6]
	print "Sampling and Z-score, P-value calculation for compounds begin"
	totalTime = 0.0
	start_time = time()
	for i in range(sampleTimes):		
		d_sample = randomSample(list_r2, sampleSize, sampleRange)
		for key, value in d_cinp1.iteritems():
			if key in d_sample:
				value.append(d_sample[key])
			else:
				value.append(0)

	for key, value in d_cinp1.iteritems():
		new_value = calPValue(value)
		d_cinp1[key] = new_value

	end_time = time()
	totalTime += end_time - start_time
	#print "total time for compounds = ", totalTime
	print "Compounds processing completed"
	inp_r1 = open(arg[3],"r")
	d_cinp2 = {}
	for line in inp_r1:
		keywords = 	line.split("~")
		#print keywords
		for key in keywords:
			if key != '\n':
				if key not in d_cinp2:
					d_cinp2.update({key:[1]})
				else:
					d_cinp2[key][0] = d_cinp2[key][0]+1
	inp_r1.close()

	inp_r2 = open(arg[4],"r")
	list2_r2 = [[]]
	sampleRange = 0
	for line in inp_r2:
		sampleRange += 1
		keywords = line.split(";")
		tokens = keywords[0].split("~")
		temp = []
		for key in tokens:
			if key != '\n' and key != '':
				temp.append(key)
		list2_r2.append(temp)
	inp_r2.close()
	print "Sampling and Z-score, P-value calculation for proteins begin"
	start_time = time()	
	for i in range(sampleTimes):	
		d_sample = randomSample(list2_r2, sampleSize, sampleRange)
		for key, value in d_cinp2.iteritems():
			if key in d_sample:
				value.append(d_sample[key])
			else:
				value.append(0)

	for key, value in d_cinp2.iteritems():
		new_value = calPValue(value)
		d_cinp2[key] = new_value

	end_time = time()
	totalTime += end_time - start_time
	print "\nTotal time for compounds and proteins= ", totalTime, "seconds"
	print "Processing proteins completed"
	sorted_result_compound_list = sortResult(d_cinp1);
	sorted_result_protein_list = sortResult(d_cinp2);
	p_value2 = 1.0
	'''inp_r5 = open(arg[5],"r")
	d_pmid = {}
	for line in inp_r5:
		line = line.strip()
		pmids = line.split(";")
		if pmids[-1] not in d_pmid:
			d_pmid.update({pmids[-1]:pmids[0]})
	inp_r5.close()
	'''
	outf1 = open(arg[6]+"_temp_result1_Substance_compounds_cutoff_"+str(cutoff)+"_p_"+str(p_value)+".txt","w")
	outf2 = open(arg[6]+"_temp_result1_Substance_proteins_cutoff_"+str(cutoff)+"_p_"+str(p_value)+".txt","w")
	outf3 = open(arg[6]+"_temp_result1_Substance_compounds_cutoff_"+str(cutoff)+".txt","w")
	outf4 = open(arg[6]+"_temp_result1_Substance_proteins_cutoff_"+str(cutoff)+".txt","w")

	outf1.write("Term Pair\tDistribution\tz score\tp value\n")
	outf3.write("Term Pair\tDistribution\tz score\tp value\n")
	for key in sorted_result_compound_list:
		if d_cinp1[key][0]>=cutoff and d_cinp1[key][-1]<=p_value:
			outf1.write(arg[6]+";"+key+"\t"+str(d_cinp1[key][:-2])+"\t"+str(d_cinp1[key][-2])+"\t"+str(d_cinp1[key][-1])+"\n")
		if d_cinp1[key][0]>=cutoff and d_cinp1[key][-1]<=p_value2:
			outf3.write(arg[6]+";"+key+"\t"+str(d_cinp1[key][:-2])+"\t"+str(d_cinp1[key][-2])+"\t"+str(d_cinp1[key][-1])+"\n")
	outf1.close()
	outf3.close()

	outf2.write("Term Pair\tDistribution\tz score\tp value\n")
	outf4.write("Term Pair\tDistribution\tz score\tp value\n")
	for key in sorted_result_protein_list:
		if d_cinp2[key][0]>=cutoff and d_cinp2[key][-1]<=p_value:
			outf2.write(arg[6]+";"+key+"\t"+str(d_cinp2[key][:-2])+"\t"+str(d_cinp2[key][-2])+"\t"+str(d_cinp2[key][-1])+"\n")
		if d_cinp2[key][0]>=cutoff and d_cinp2[key][-1]<=p_value2:
			outf4.write(arg[6]+";"+key+"\t"+str(d_cinp2[key][:-2])+"\t"+str(d_cinp2[key][-2])+"\t"+str(d_cinp2[key][-1])+"\n")
	outf2.close()
	outf4.close()

if __name__ == "__main__":
	main()