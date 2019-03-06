import sys
import os

def is_number(in_string):
	try:
		float(in_string)
		return True
	except ValueError:
		return False

def main(arg = sys.argv):
	if len(arg) != 2:
		print >> sys.stderr,"Please pass 1 argument in the format <compounds or proteins output file with p value>"
		sys.exit(1)
	# Check the file name.
	filePath = os.path.split(os.path.abspath(arg[1]))[0]
	bName = os.path.splitext(os.path.basename(arg[1]))[0]
	itemList = bName.split("_")
	if not (itemList[-2] == "p" and is_number(itemList[-1])):
		print >> sys.stderr, "The file name of input argument need to be checked."
		sys.exit(1)
	itemList[-2] = "adjustedPvalue"
	outfName = os.path.join(filePath, "_".join(itemList)+".txt")
	outf = open(outfName, "w")
	inf = open(arg[1],"r")
	line = inf.readline()
	outf.write(line.strip()+'\tadjusted p value\n')
	line = inf.readline()
	content = []
	while line != "":
		content.append(line.strip().split("\t"))
		line = inf.readline()
	len1 = len(content)
	contentUpdated = []
	for i, item in enumerate(content):
		rank = i+1
		adjusted_p_value = float(item[-1])*len1/rank
		item.append(str(adjusted_p_value))
		contentUpdated.append(item)
	for item in contentUpdated:
		outf.write("\t".join(item)+'\n')
	outf.close()

if __name__ == "__main__":
	main()