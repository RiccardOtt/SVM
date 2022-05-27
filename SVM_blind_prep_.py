import sys
import os
import numpy as np
np.set_printoptions(threshold=np.inf)

def matrix_pssm(filename):
	seq_prof = []
	for j in filename:
		L = j.split()
		try: L[0]
		except: pass
		else:
			if L[0].isdigit():
				seq_prof.append(L[22:-2])
	x = np.array(seq_prof, dtype=np.float64)
	x /= 100

	if x.sum() != float(0):
		return(x), True
	else:
		return(None), False

def SVM_input(profile1, dssp1):
	P = len(profile1)
#	print(P)
	pad = np.zeros((8,20))
	SVM_profile = np.zeros(shape=(P,340))
	padding = np.vstack((pad,profile1,pad))
	Lp = len(padding)
	dssp_seq = dssp1.read().splitlines()[1]
	dssp_ok='-'*8 + dssp_seq + '-'*8
	for e,i in zip(range(8,Lp-8), range(0,P)):
		j=e-8
		k=e+8
		window = padding[j:k+1]
		SVM_line = window.flatten() #it is a 60 items element
		SVM_profile[i] += SVM_line
	class_column = np.zeros(shape=(P,1))
	SVM_profile_w_class = np.append(class_column,SVM_profile,axis=1)
	SVM_profile_list_w_class = SVM_profile_w_class.tolist()
	for e in range(0,P):
		if dssp_ok[e] == 'H':
			SVM_profile_list_w_class[e][0] = 1
		if dssp_ok[e] == 'E':
			SVM_profile_list_w_class[e][0] = 2
		if dssp_ok[e] == '-':
			SVM_profile_list_w_class[e][0] = 3
	for line in SVM_profile_list_w_class:
#		print(SVM_profile_list_w_class)
#		print(line)
#		print(line[0])					#the error is the fact that the first residue has value 0.0
		for indice1 in range(0,341):
			line[indice1] = indice1,':', line[indice1]			#indice1, line[indice1]
		filtered_line = [i for i in line if i[2] > 0]				#i[1] > 0
#		print(filtered_line[0][2])
#		print(line)
#		print(filtered_line)
#		indice = 0
#		print(filtered_line[0])
#GOOD FOR NOW
#		if filtered_line[indice][0] == 0:
#			print(filtered_line[0])
		filtered_line[0] = filtered_line[0][2]		#=filtered_line[indice][1]
		print(filtered_line)

if __name__ == '__main__':
	fileid = sys.argv[1]
	with open(fileid) as filein:
		for id in filein:
#			print('>'+id)
#			output_file= open(id+"_SVM_in.txt",'w+')
			id=id.rstrip()
			profile_file = '/home/riccardo/Documents/Documents/LB2/Castrense/project/blindset.pssm/' + id + '.pssm'
			dssp_file = '/home/riccardo/Documents/Documents/LB2/Castrense/project/blindset_dssp/' + id + '.dssp'
			try:
				prof=open(profile_file)
				dssp=open(dssp_file)
			except: continue
			else:
				profile, ret = matrix_pssm(prof)
				if ret:
#					print(id)
					SVM_input(profile, dssp)
#					results =  SVM_input(profile, dssp)
#			output_file.write(results)
#			output_file.close()
