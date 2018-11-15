#!/usr/bin/env python
'''
usage: 

python database_search.py -sig [signature|pattern] ....

parameters:
    -sig       signature or pattern to search
    -v         v gene
    -j         j gene
    -usage     out put the v gene or j gene usage, defualt is 0: no, 1 for yes.
    -vj        specific recombination of V and J frequency and usage, defualt is 0. 1 need both V and J input
    -wsr       wheather write signature matched reads, defualt is 0. 1 is yes
    -chain     HKL H is heavy chain, L is lambda chain, K is kappa chain, defualt is heavy chain
    -cdr3_len  length of signature cdr3 defualt is no limit, or use x~y means range of cdr3 length
    -subtype   For heavy chain search: IGHG,IGHM,IGHA, defualt is all subtype.
    -mode      cdr3: cdr3 mode (defualt), pos: position model.
								Example signature for cdr3 mode: '[A-Z]{3}[AFILMYWV][EQ][A-Z]{2}'.
								Example signature for position mode: 22,A,25,G,50,K,60,[EQ]
    -scheme    Which numbering scheme should be used. defualt is kabat. imgt and chothia is optional.

A script to search signature and usage etc from the database. Writen by yicheng guo at 27th march 2018 in columbia university.
All rights reseverd to NIH and CUMC.
'''

import os
import sys
import numpy
import re
import pickle
global dic_seq
import time

##command line and paramter setting
q = lambda x: x in sys.argv
if any([q(x) for x in ["h", "-h", "--h", "help", "-help", "--help"]]):
    print __doc__
    sys.exit(0)

p = sys.argv
vgene = ''
jgene = ''
signature = '[A-Z]{3}[AFILMYWV][EQ][A-Z]{2}'
wsr = 0
chain = 'H'
cdr3_len = ''
vj_need = 0
subtype = 'all'
mode = 'cdr3'
scheme = 'kabat'
usage = 1
for i in range(len(p)):
    if p[i] == '-sig':
        signature = p[i+1]
    if p[i] == '-v':
        vgene = p[i+1]
    if p[i] == '-j':
        jgene = p[i+1]
    if p[i] == '-vj':
        vj_need = p[i+1]
    if p[i] == '-wsr':
        wsr = p[i+1]
    if p[i] == '-chain':
        chain = p[i+1]
    if p[i] == '-cdr3_len':
        cdr3_len = p[i+1]
    if p[i] == '-subtype':
        subtype = p[i+1]
    if p[i] == '-mode':
        mode = p[i+1]
    if p[i] == '-scheme':
        scheme = p[i+1]
    if p[i] == '-usage':
        scheme = p[i+1]

path = os.getcwd()

if chain == 'H':
	data_file = 'all_unique_reads_HEAVY.fasta_good.fasta'
	if subtype != 'all':
		print subtype
if chain == 'K':
	data_file = 'all_unique_reads_light.fasta_good.fasta'
if chain == 'L':
	data_file = 'all_unique_reads_light.fasta_good.fasta'
if chain == 'LK' or chain == 'KL':
	data_file = 'all_unique_reads_light.fasta_good.fasta'
if chain == 'SLE':
	data_file = 'all_unique_reads_HEAVY.fasta_good_sle.fasta'


#load pos dic
if chain == 'H': 
	if mode == 'pos':
		dic_pos = pickle.load(open('heavy_position_'+scheme+'_dic.txt','r'))
		print 'load position complete....'
else:
	if mode == 'pos':
		dic_pos = pickle.load(open('light_position_'+scheme+'_dic.txt','r'))


##sub functions
def get_codon_dic(x):
    dic = {'AAA':'K','AAT':'N','AAC':'N','AAG':'K','ATA':'I','ATT':'I','ATC':'I','ATG':'M','ACA':'T','ACT':'T','ACC':'T','ACG':'T','AGA':'R','AGT':'S','AGC':'S','AGG':'R'\
    ,'TAA':'*','TAT':'Y','TAC':'Y','TAG':'*','TTA':'L','TTT':'F','TTC':'F','TTG':'L','TCA':'S','TCT':'S','TCC':'S','TCG':'S','TGA':'*','TGT':'C','TGC':'C','TGG':'W'\
    ,'CAA':'Q','CAT':'H','CAC':'H','CAG':'Q','CTA':'L','CTT':'L','CTC':'L','CTG':'L','CCA':'P','CCT':'P','CCC':'P','CCG':'P','CGA':'R','CGT':'R','CGC':'R','CGG':'R'\
    ,'GAA':'E','GAT':'D','GAC':'D','GAG':'E','GTA':'V','GTT':'V','GTC':'V','GTG':'V','GCA':'A','GCT':'A','GCC':'A','GCG':'A','GGA':'G','GGT':'G','GGC':'G','GGG':'G'}
    return dic
dic_codon = get_codon_dic(1)

def translate_seq(seq):
	seq = seq.strip()
	trans_seq = ''
	for i in range(0,len(seq),3):
		if dic_codon.has_key(seq[i:i+3]):
			trans_seq += dic_codon[seq[i:i+3]]
		else:
			trans_seq += 'X'
	return trans_seq

def subtype_identify(line,subtype):
	x = 0
	constant = line.split( )[4].split('=')[1]
	if subtype != 'all':
		if subtype in constant:
			x = 1
	else:
		x = 1
	return x

def get_sequences(file_name,subtype):
	dic = {}
	f = open(file_name,'r').readlines()
	for i in range(len(f)):
		if f[i][0] == '>':
			if subtype_identify(f[i],subtype) == 1:
				name = f[i][1:]
				dic[name] = f[i+1]
	return dic

def detransform_pos(lst):
	lst = lst.split(',')
	final_lst = []
	for ele in lst:
		if '~' in ele:
			for i in range(int(ele.split('~')[0]),int(ele.split('~')[1])+1):
				final_lst.append(str(i))
		else:
			final_lst.append(ele)
	return final_lst

#NEED determine#################################################################
#dic_pos = pickle.dump(open('xx.txt'),'r')
def search_pos(dic_seq,key,signature,abs_pos):#need complete
	#print signature, abs_pos
	control = 0
	reads = key.split( )[0]
	search_pos = []
	target_seq = ''
	search_seq = ''
	signature = signature.split(',')
	#print signature,abs_pos,translate_seq(dic_seq[key])
	for i in range(0,len(signature),2):
		if str(signature[i]) in abs_pos:
			search_pos = str(signature[i])
			#print search_pos,[abs_pos.index(search_pos)],str(signature[i+1])
			target_seq+=translate_seq(dic_seq[key])[abs_pos.index(search_pos)]
			search_seq+= str(signature[i+1])
		else:
			control = 2
			break
	#print control,target_seq,search_seq
	if control != 2:
		#print target_seq,search_seq
		if re.search(search_seq,target_seq):
			control = 1
	return control

def judge_cdr3_len(cdr3_len,seq_cdr3_len):
	if '~' in cdr3_len:
		if seq_cdr3_len in range(int(cdr3_len.split('~')[0])+2,int(cdr3_len.split('~')[1])+3):
			return 1
		else:
			return 0
	else:
		if int(cdr3_len)+2 == seq_cdr3_len:
			return 1

def get_big_dic(dic_seq,signature,mode,subtype):
	if wsr != 0:
		f_sig = open(signature+'_signature_reads.fasta','w')
	v_count = 0
	j_count = 0
	vj_count = 0
	dic_sample_statistic = {}## like SRR440767:{total_reads:10000, v1~vn gene:10000, j1~jn gene:10000}
	search_count = 0
	total_count = float(len(dic_seq.keys()))
	time_control = 0
	for key in dic_seq.keys():
		search_count += 1
		if int((search_count*100/total_count))//2 >= time_control:
			str_time = '>'*(int((search_count*100/total_count))//2)+' '*(int((100-(search_count*100/total_count)))//2)
			sys.stdout.write('\rsignature searching: '+str_time+'[%s%%]'%(str(round(search_count*100/total_count,2))))
			sys.stdout.flush()
			time_control+=1
		key_sp = key.strip().split( )
		CDR3 = key_sp[-2].split('=')[1]
		#CDR3_start = translate_seq(dic_seq[key]).index(CDR3,50)+1
		#CDR3_end = CDR3_start + len(CDR3)
		sample = key_sp[0].split('_')[0]
		if dic_sample_statistic.has_key(sample):
			dic_sample_statistic[sample]['total_reads'] += 1
		else:
			dic_sample_statistic[sample] = {'total_reads':1,'signature_reads':0}
		V = key_sp[1].split('=')[1].split(',')[0:1]
		#print V
		J = key_sp[2].split('=')[1].split(',')[0:1]
		for vv in V:
			if dic_sample_statistic[sample].has_key(vv):
				dic_sample_statistic[sample][vv] += 1
			else:
				dic_sample_statistic[sample][vv] = 1
		for jj in J:
			if dic_sample_statistic[sample].has_key(jj):
				dic_sample_statistic[sample][jj] += 1
			else:
				dic_sample_statistic[sample][jj] = 1
		if vgene != '':
			if vgene in V[0]:
				v_count += 1
				if jgene != '':
					if jgene in J[0]:
						vj_count +=1
		if jgene != '':
			if jgene in J[0]:
				j_count += 1
		#print v_count,j_count,vj_count
		if signature != '':
			if vgene != '':
				if mode == 'cdr3' and vgene in V[0]:
					if re.search(signature,CDR3):#r'[A-Z]{3}[AFILMYWV][EQ][A-Z]{2}'
						if cdr3_len != '':
							if judge_cdr3_len(cdr3_len,len(CDR3))==1:
								dic_sample_statistic[sample]['signature_reads'] += 1
						else:
							dic_sample_statistic[sample]['signature_reads'] += 1
							#print len(CDR3)
						if wsr != 0:
							f_sig.write('>'+key+dic_seq[key])
				if mode == 'pos' and vgene in V[0]:##need complete
					if dic_pos.has_key(key.split( )[0]):
						if search_pos(dic_seq,key,signature,detransform_pos(dic_pos[key.split( )[0]][0])) == 1:
							dic_sample_statistic[sample]['signature_reads'] += 1 
							if wsr != 0:
								f_sig.write('>'+key+dic_seq[key])
			else:
				if mode =='cdr3':
					if re.search(signature,CDR3):#r'[A-Z]{3}[AFILMYWV][EQ][A-Z]{2}'
						if cdr3_len != '':
							if judge_cdr3_len(cdr3_len,len(CDR3))==1:
								dic_sample_statistic[sample]['signature_reads'] += 1
						else:
							dic_sample_statistic[sample]['signature_reads'] += 1
						if wsr != 0:
							f_sig.write('>'+key+dic_seq[key])
				if mode == 'pos':##need complete
					if dic_pos.has_key(key.split( )[0]):
						#print dic_pos[key.split( )[0]]
						if search_pos(dic_seq,key,signature,detransform_pos(dic_pos[key.split( )[0]][0])) == 1:
							dic_sample_statistic[sample]['signature_reads'] += 1 
							if wsr != 0:
								f_sig.write('>'+key+dic_seq[key])
	for sample in dic_sample_statistic.keys():
		kkk = 1
		#print sample
		#print sample + '\t' + str(float(dic_sample_statistic[sample]['signature_reads']*1000000/dic_sample_statistic[sample]['total_reads']))
		#print dic_sample_statistic[sample][vgene]*1000000/dic_sample_statistic[sample]['total_reads']
	return dic_sample_statistic
##		print re.search(CDR3,translate_seq(dic_seq[key])).group()##if re.search: do something

##save donor_class_dic
f = open('donor_class.txt','r').readlines()
fo = open('donor_class_dic.txt','w')
dic_donor = {}
for line in f[1:]:
	line = line.strip().split('\t')
	dic_donor[line[2]] = line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[4]
pickle.dump(dic_donor,fo)
fo.close()


def main():
	print 'loading database....'
	dic_seq = get_sequences(data_file,subtype)
	print 'load database complete...'
	dic_sample_statistic = get_big_dic(dic_seq,signature,mode,subtype)
#	for sample in dic_sample_statistic.keys():
#		print sample,dic_sample_statistic[sample]['signature_reads']
	if usage == 1:
		f = open('v_useage_for_r_plot.txt','w')
		dic_donor = pickle.load(open('donor_class_dic.txt','r'))
		f.write('sample\tusage\tvgene\tcell_type\tdonor_class\tproject\tdonor_id\n')
		for sample in dic_sample_statistic.keys():
			#print sample
			for vgene in dic_sample_statistic[sample].keys():
				if vgene != 'total_reads' and vgene != 'signature_reads':
					if int(dic_sample_statistic[sample]['total_reads']) >= 1700 and dic_donor.has_key(sample):
						f.write(sample+'\t'+str(dic_sample_statistic[sample][vgene]*100/float(dic_sample_statistic[sample]['total_reads']))+'\t'+vgene+'\t'+dic_donor[sample]+'\n')
		f.close()
	f = open('v_useage_for_r_plot_signature.txt','w')
	f.write('sample\tusage\tvgene\tcell_type\tdonor_class\tproject\tdonor_id\n')
	if signature != '':
		for sample in dic_sample_statistic.keys():
			for vgene in dic_sample_statistic[sample].keys():
				if vgene != 'total_reads' and vgene == 'signature_reads':
					if int(dic_sample_statistic[sample]['total_reads']) >= 1700 and dic_donor.has_key(sample):
						f.write(sample+'\t'+str(dic_sample_statistic[sample][vgene]*100/float(dic_sample_statistic[sample]['total_reads']))+'\t'+vgene+'\t'+dic_donor[sample]+'\n')
		f.close()


if __name__ == '__main__':
	main()



#use for subtype count
if 1 == 0:
	dic_subcount = {}
	fx = open('subtype_count_for_each_sample.txt','w')
	subtypes = ['IGHG','IGHA','IGHE','IGHM','IGHD','not_found']
	fx.write('sample\t'+'\t'.join(subtypes)+'\n')
	for subtype in subtypes:
		dic_sample_statistic = get_big_dic(get_sequences(data_file,subtype),signature,mode,subtype)
		for sample in dic_sample_statistic.keys():
			if dic_donor.has_key(sample):
				if dic_subcount.has_key(sample):
					dic_subcount[sample][subtype] = dic_sample_statistic[sample]['total_reads']
				else:
					dic_subcount[sample] = {}
					dic_subcount[sample][subtype] = dic_sample_statistic[sample]['total_reads']
	for key in dic_subcount.keys():
		aa = key
		for subtype in subtypes:
			if dic_subcount[key].has_key(subtype):
				aa += '\t' + str(dic_subcount[key][subtype])
			else:
				aa += '\t0'
	fx.write(aa+'\n')
	fx.close()
