import sys
import os
import platform

def main(arg = sys.argv):
	if len(arg) != 6:
		print >> sys.stderr,"Please pass all the 6 arguments in the format <Drug result compounds>, <Drug result proteins>, <Drug PMID substances>, <Frequency Threshold> <adjusted P value Threshold>"
		sys.exit(1)
	dir = os.path.dirname(os.path.abspath(arg[1]))
	threshold = int(arg[4])
	pThreshold = float(arg[5])
	targetName = ''
	compoundList = []
	proteinList = []
	substanceList = []
	with open(arg[1], 'r') as inf:
		line = inf.readline()
		line = inf.readline()
		while line != '':
			pValue = float(line.strip().split('\t')[-1])
			targetName = line.strip().split('\t')[0].split(';')[0]
			compound = line.strip().split('\t')[0].split(';')[1]
			count = int(line.strip().split('\t')[1].split(',')[0][1:].strip())
			if count >= threshold and pValue <= pThreshold:
				if compound not in compoundList:
					compoundList.append(compound)
			line = inf.readline()
	with open(arg[2], 'r') as inf:
		line = inf.readline()
		line = inf.readline()
		while line != '':
			pValue = float(line.strip().split('\t')[-1])
			protein = line.strip().split('\t')[0].split(';')[1]
			if pValue <= pThreshold:
				if protein not in proteinList:
					proteinList.append(protein)
			line = inf.readline()
	substanceList = compoundList + proteinList
	substanceList.append(targetName)
	uniqueSub = []
	uniquePMID = []
	PMIDSubDic = {}
	SubPMIDDic = {}
	with open(arg[3], 'r') as inf:
		line = inf.readline()
		line = inf.readline()
		while line != '':
			PMID = line.strip().split(";")[0]
			Sub = line.strip().split(";")[-1]
			for sub in substanceList:
				if sub not in uniqueSub:
					uniqueSub.append(sub)
				else:
					pass
				if PMID not in uniquePMID:
					uniquePMID.append(PMID)
				else:
					pass
				if PMID not in PMIDSubDic.keys():
					newPair = {PMID:[Sub]}
					PMIDSubDic.update(newPair)
				else:
					PMIDSubDic[PMID].append(Sub)
				if Sub not in SubPMIDDic.keys():
					newPair2 = {Sub:[PMID]}
					SubPMIDDic.update(newPair2)
				else:
					SubPMIDDic[Sub].append(PMID)
			else:
				pass
			line = inf.readline()
	subCount1List = []
	for iterm in uniqueSub:
		itermCount = 0
		for PMID_1 in PMIDSubDic.keys():
			if iterm in PMIDSubDic[PMID_1]:
				itermCount = itermCount + 1
		subCount1List.append([itermCount,iterm])
	subCount1List.sort()
	subCount1List.reverse()
	outf1 = open(dir+"\\compound_protein_Counts.txt","w")
	outf1.write("Substances\tCount\n")
	for iterm in subCount1List:
		outf1.write(iterm[-1]+"\t"+str(iterm[0])+"\n")
	outf1.close()
	outf_new1 = open(dir+"\\matrix.txt","w")
	outf_new1.write("Substances\t")
	for i in range(len(uniquePMID)):
		if i != len(uniquePMID)-1:
			outf_new1.write(uniquePMID[i]+"\t")
		else:
			outf_new1.write(uniquePMID[i]+"\n")
	for item in uniqueSub:
		outf_new1.write(item+"\t")
		for k in range(len(uniquePMID)):
			if k != len(uniquePMID)-1:
				if uniquePMID[k] in SubPMIDDic[item]:
					outf_new1.write("1\t")
				else:
					outf_new1.write("0\t")
			else:
				if uniquePMID[k] in SubPMIDDic[item]:
					outf_new1.write("1\n")
				else:
					outf_new1.write("0\n")
	outf_new1.close()
	outf2 = open(dir+'\\vertex_attribute.txt', 'w')
	outf2.write('Substances\tType\n')
	for item in compoundList:
		outf2.write(item+'\tCompounds\n')
	outf2.write(targetName+'\tTarget\n')
	for item in proteinList:
		outf2.write(item+'\tProteins\n')
	outf2.close()
	if os.path.isfile('Plot_network.R'):
		if platform.system() == "Windows":
			os.system('Rscript.exe Plot_network.R')
		elif platform.system() == "Linux" or platform.system() == "Darwin":
			os.system('Rscript Plot_network.R')
		else:
			print "Unknown Operation System. plot_network.R can't run."
	else:
		print "Plot_network.R is missing. Please check it"

if __name__ == "__main__":
	main()