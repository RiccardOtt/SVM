import sys
import numpy as np


def from_class_to_str(dssp,pred):
	str = ''
	pred_class = []
	start = 0
	pred_str = ''
	for pred_lines in pred:
		pred_lines = pred_lines.rstrip()
		pred_class.append(pred_lines)
	for dssp_lines in dssp:
		dssp_lines = dssp_lines.rstrip()
		if '>' in dssp_lines:
			id = dssp_lines
		else:
			structure = dssp_lines
			str += structure
			len_str = len(structure)
			end = start + len_str
#			print(start,end)
#			print(id+'\n'+structure,len_str)
			pred_str += '\n'+id+'\n'
			for i in range(start, end):
				if pred_class[i] == '1':
					pred_str += 'H'
				if pred_class[i] == '2':
					pred_str += 'E'
				if pred_class[i] == '3':
					pred_str += '-'
			start += len_str
	print(pred_str)
#	print(len(str))
#	print(len(pred_class))


if __name__ == '__main__':
	dssp_file = sys.argv[1]
	predclass = sys.argv[2]
	with open(dssp_file) as dssp, open(predclass) as predclas:
		from_class_to_str(dssp,predclas)
