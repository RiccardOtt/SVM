import sys
import math
import re


def conf_matrix(dssp,predstr):
	cm = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
	dssp_dict = {}
	dssp_ids = []
	dssp_strs = []
	pred_dict = {}
	pred_ids = []
	pred_strs = []
	for lines in dssp:
		lines = lines.rstrip()
		if '>' in lines:
			dssp_ids.append(lines)
		else:
			dssp_strs.append(lines)
	for i in range(len(dssp_ids)):     #create a dictionary from 2 parallel list
		dssp_dict[dssp_ids[i]] = dssp_strs[i]
#	print(dssp_dict)
	for lines in predstr:
		lines = lines.rstrip()
		if '>' in lines:
			pred_ids.append(lines)
		else:
			pred_strs.append(lines)
	for i in range(len(pred_ids)):
		pred_dict[pred_ids[i]] = pred_strs[i]
#	print(pred_dict)
	for keys_1,values_1 in dssp_dict.items():
		for keys_2,values_2 in pred_dict.items():
			if keys_1 == keys_2:
				if len(values_1) == len(values_2):
#					print(keys_1,'\n',values_1,'\n',values_2)
					for i in range(len(values_1)):
							if values_1[i] == values_2[i]:
								if values_1[i] == 'H':
									cm[0][0] += 1
								if values_1[i] == 'E':
									cm[1][1] += 1
								if values_1[i] == '-':
									cm[2][2] += 1
							elif values_1[i] != values_2[i]:
								if values_2[i] == 'H':
									if values_1[i] == 'E':
										cm[0][1] += 1
									elif values_1[i] == '-':
										cm[0][2] += 1
								if values_2[i] == 'E':
									if values_1[i] == 'H':
										cm[1][0] += 1
									elif values_1[i] == '-':
										cm[1][2] += 1
								if values_2[i] == '-':
									if values_1[i] == 'H':
										cm[2][0] += 1
									elif values_1[i] == 'E':
										cm[2][1] += 1



	for h in cm:
		print(h)
	return(cm,dssp_dict,pred_dict)

def performance(cm):
	C_H = cm[0][0]
	C_E = cm[1][1]
	C_C = cm[2][2]

	U_H = cm[1][0]+cm[2][0]
	U_E = cm[0][1]+cm[2][1]
	U_C = cm[0][2]+cm[1][2]

	O_H = cm[0][1]+cm[0][2]
	O_E = cm[1][0]+cm[1][2]
	O_C = cm[2][0]+cm[2][1]

	N_H = cm[1][1]+cm[1][2]+cm[2][1]+cm[2][2]
	N_E = cm[0][0]+cm[0][2]+cm[2][0]+cm[2][2]
	N_C = cm[0][0]+cm[0][1]+cm[1][0]+cm[1][1]

	SEN_H = C_H/(C_H+U_H)
	SEN_E = C_E/(C_E+U_E)
	SEN_C = C_C/(C_C+U_C)

	PPV_H = C_H/(C_H+O_H)
	PPV_E = C_E/(C_E+O_E)
	PPV_C = C_C/(C_C+O_C)

	ACC_H = (C_H+N_H)/(C_H+N_H+U_H+O_H)
	ACC_E = (C_E+N_E)/(C_E+N_E+U_E+O_E)
	ACC_C = (C_C+N_C)/(C_C+N_C+U_C+O_C)
	TOT_ACC = (C_H+C_E+C_C)/(C_H+U_H+O_H+N_H)

	MCC_H = ((C_H*N_H)-(O_H*U_H))/math.sqrt((C_H+O_H)*(C_H+U_H)*(N_H+O_H)*(N_H+U_H))
	MCC_E = ((C_E*N_E)-(O_E*U_E))/math.sqrt((C_E+O_E)*(C_E+U_E)*(N_E+O_E)*(N_E+U_E))
	MCC_C = ((C_C*N_C)-(O_C*U_C))/math.sqrt((C_C+O_C)*(C_C+U_C)*(N_C+O_C)*(N_C+U_C))

	print('sens_helix =',SEN_H,'\n','sens_beta =',SEN_E,'\n','sens_coil =',SEN_C,'\n','ppv_helix =',PPV_H,'\n','ppv_beta =',PPV_E,'\n','ppv_coil =',PPV_C,'\n','Acc_H =',ACC_H,'\n','Acc_E =',ACC_E,'\n','Acc_C =',ACC_C,'\n', 'Tot_Acc =',TOT_ACC, '\n', 'MCC_helix =',MCC_H,'\n','MCC_beta =',MCC_E,'\n','MCC_coil =',MCC_C)




def SegmentsOVerlap(dssp_dict,pred_dict):
	SOV_H = []
	SOV_E = []
	SOV_C = []
	SOV_tot = []
#######	FIND SEGMENTS ####################################################
	for keys_1,values_1 in dssp_dict.items():
		for keys_2,values_2 in pred_dict.items():
			if keys_1 == keys_2:
				if len(values_1) == len(values_2):
					seqs = [values_1,values_2]
					frag = {'H':[[],[]], 'E':[[], []], '-':[[], []]}
					rr = ['H', 'E', '-']
					for c in rr:
						my_regex = re.escape(c) + r'+'
						for seq, num in zip(seqs,range(2)):
							idx = 0
							while c in seq[idx:]:
								sr = seq[idx:]
								sing_frag = []
								x = re.search(my_regex,sr)
								for i in range(x.span()[0]+idx, x.span()[1]+idx):
									sing_frag.append(i)
								frag[c][num].append(sorted(sing_frag))
								idx += x.span()[1]
#	print(keys_1)
#	print(values_1)
#	print(values_2)
#	print(values_1)
#	print(values_2)
#					print(frag['H'][0], frag['H'][1], frag['E'][0], frag['E'][1], frag['-'][0], frag['-'][1])
#					print(frag)


####### SEGMS OVERLAP #########################################################
					for strc,seq in frag.items():
#						print(strc)
#						print(seq)
						N_obs = 0
						N_pred = 0
						N_struct = 0
						SOV_segm = []
						minov = []
						maxov = []
						lensegs_obs = []
						lensegs_pred = []
						delta = []
						obs_seqs = seq[0]
						pred_seqs = seq[1]
#						print(strc)
#						print(obs_seqs)
#						print(pred_seqs)
						for obs_seq in obs_seqs:
							for pred_seq in pred_seqs:
								N_obs += len(obs_seq)
								intersect = set(obs_seq) & set(pred_seq)
								intersect = sorted(intersect)
								if len(intersect) > 0:
#									print(intersect)
									minov.append(len(intersect))
									union = set(obs_seq) | set(pred_seq)
									union = sorted(union)
#									print(union)
									maxov.append(len(union))
									lensegs_obs.append(len(obs_seq))
									lensegs_pred.append(len(pred_seq))



						N_struct += sum(lensegs_obs)


#						print(N_obs)
#						print('MAXOV','\n',maxov)
#						print('MINOV','\n',minov)
#						print('Lenght Observe Segments','\n',lensegs_obs)
#						print('Lenght Predicted Segments','\n',lensegs_pred)
#						print(N_struct)


						for a,b,c,d in zip(range(len(maxov)),range(len(minov)),range(len(lensegs_obs)),range(len(lensegs_pred))):
							dlt = [(maxov[a]-minov[b]), (minov[b]), ((lensegs_obs[c])/2), ((lensegs_pred[d])/2)]
#							print(a)
#							print(dlt)
							delta.append(min(dlt))
#						print('DELTA','\n',delta)




						for a,b,c,d in zip(range(len(maxov)),range(len(minov)),range(len(lensegs_obs)),range(len(delta))):
							SOV_segm.append(((minov[b]+delta[d])/maxov[a])*lensegs_obs[c])
						if N_struct == 0:	continue
						if strc == 'H':
							SOV_H.append(sum(SOV_segm)*(100/N_struct))
						if strc == 'E':
							SOV_E.append(sum(SOV_segm)*(100/N_struct))
						if strc == '-':
							SOV_C.append(sum(SOV_segm)*(100/N_struct))
						SOV_tot.append(sum(SOV_segm)*(100/N_struct))
	print('SOV_H =',sum(SOV_H)/len(SOV_H))
	print('SOV_E =',sum(SOV_E)/len(SOV_E))
	print('SOV_C =',sum(SOV_C)/len(SOV_C))
	print('SOV_tot =',sum(SOV_tot)/len(SOV_tot))






if __name__ == '__main__':
	dssp_str = sys.argv[1]
	pred_str = sys.argv[2]
	with open(dssp_str) as dssp, open(pred_str) as pred:
		cm,dsspdict,preddict = conf_matrix(dssp,pred)
		performance(cm)
		SegmentsOVerlap(dsspdict,preddict)





