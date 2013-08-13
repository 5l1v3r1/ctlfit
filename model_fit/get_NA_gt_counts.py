#!/ebio/ag-neher/share/programs/EPD/bin/python
'''
author:     Taylor Kessinger
date:       19/06/2013
content:    Pull out epitope mutations from FASTA files.
'''

from Bio import AlignIO
import re
import numpy as np

patient = 'CH40'

AA_dict = {}
AA_dict['A'] = ['GCU', 'GCC', 'GCA', 'GCG']
AA_dict['C'] = ['UGU', 'UGC']
AA_dict['D'] = ['GAU', 'GAC']
AA_dict['E'] = ['GAA', 'GAG']
AA_dict['F'] = ['UUU', 'UUC']
AA_dict['G'] = ['GGU', 'GGC', 'GGA', 'GGG']
AA_dict['H'] = ['CAU', 'CAC']
AA_dict['I'] = ['AUU', 'AUC', 'AUA']
AA_dict['K'] = ['AAA', 'AAG']
AA_dict['L'] = ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG']
AA_dict['M'] = ['AUG']
AA_dict['N'] = ['AAU', 'AAC']
#AA_dict['O'] = ['UAG'] #pyrrlosine
AA_dict['P'] = ['CCU', 'CCC', 'CCA', 'CCG']
AA_dict['Q'] = ['CAA', 'CAG']
AA_dict['R'] = ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
AA_dict['S'] = ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']
AA_dict['T'] = ['ACU', 'ACC', 'ACA', 'ACG']
AA_dict['V'] = ['GUU', 'GUC', 'GUA', 'GUG']
AA_dict['W'] = ['UGG']
AA_dict['Y'] = ['UAU', 'UAC']
AA_dict['*'] = ['UAA', 'UAG', 'UGA']

if patient == 'CH40':
    epitopes = ['nef191', 'gag125', 'gag405', 'vpr84', 'vif167', 'env145']
    NA_locations = ['nef571', 'gag373', 'gag1216', 'vpr252', 'vif499', 'env433']
    tp = [0, 16, 45, 111, 181]
    foundername = 'CH40_con'

elif patient == 'CH58':
    epitopes = ['env578', 'env829', 'nef112', 'gag246']
    NA_locations = ['env1732', 'env2488', 'nef334', 'gag736']
    tp = [0, 9, 45, 85]
    foundername = 'CH58_transmitted'

elif patient == 'CH77':
    epitopes = ['tat61', 'env354', 'nef26', 'nef92']
    NA_locations = ['tat181', 'env1060', 'nef76', 'nef274']
    tp = [0, 14, 32]
    foundername = 'CH77_transmitted'

sequence_gt = {}
gt_counts = {}
for ei, epitope in enumerate(epitopes):
    if not epitope in gt_counts.keys():
        gt_counts[epitope] = {}
    gene = epitope[:3]
    loc = int(epitope[3:])
    NA_loc = int(NA_locations[ei][3:])
    
    alignment_loc = 'gt_data/alignments/' + patient + '_all_proteins/' + gene.upper() + '.AA.FASTA'
    NA_alignment_loc = 'gt_data/alignments/' + patient + '_all_proteins_NA/' + gene.upper() + '.NA.FASTA'
    alignment = AlignIO.read(alignment_loc, 'fasta')
    NA_alignment = AlignIO.read(NA_alignment_loc, 'fasta')
    for seq in NA_alignment.get_all_seqs():
        seq = seq.upper()
        desc = seq.description
        temp = desc[4:].split('_')[0]
        if patient == 'CH77' and desc != 'CH77_transmitted':
            temp = desc[5:].split('_')[0]
        if temp == '':
            day = 0
        elif temp == 'E':
            day = 16
        else:
            day = int(temp[1:])
        if not day in gt_counts[epitope].keys() and day in tp:
            gt_counts[epitope][day] = {}
        sequence = seq.seq.tostring()[NA_loc-18:NA_loc+20]
        if day in tp:
            #if re.match('[A-Z]', sequence):
            print day, epitope, sequence, desc
            #desctime = 
            if not desc in sequence_gt.keys():
                sequence_gt[desc] = {}
            sequence_gt[desc][epitope] = sequence
            if not sequence in gt_counts[epitope][day].keys():
                gt_counts[epitope][day][sequence] = 1
            else:
                gt_counts[epitope][day][sequence] += 1
            #print alignment[:,loc-5:loc+5]


seqfile = open('gt_data/' + patient + '_NA_sequences.dat', 'wb')
seqfile.write('sequence\t' + '\t'.join(epitopes) + '\n')
seqfile.write(foundername + '\t' + '\t'.join([sequence_gt[foundername][epitope] for epitope in sequence_gt[foundername].keys()]) + '\n')
founder_seqs = {}
for epitope in sequence_gt[foundername].keys():
    founder_seqs[epitope] = sequence_gt[foundername][epitope]
A = sequence_gt.keys()
A.remove(foundername)
A.sort()
for desc in A:
    if patient == 'CH40':
        if 'nef191' not in sequence_gt[desc].keys():
            sequence_gt[desc]['nef191'] = '------------'
    temp_seqs = []
    for epitope in epitopes:
        if (sequence_gt[desc][epitope] == founder_seqs[epitope]):
            temp_seqs.append('wildtype')
        else:
            temp_seqs.append(sequence_gt[desc][epitope])
    seqfile.write(desc + '\t' + '\t'.join([temp_seqs[i] for i in range(len(temp_seqs))]) + '\n')
seqfile.close()


transitions = 0
transversions = 0


origin = {}
names = {}
changes = {}

dumpfile = open('gt_data/' + patient + '_NA_genotypes_raw.dat', 'wb')
B = gt_counts.keys()
for epitope in epitopes:
    observed_seqs = []
    number_observed = 0
    for day in tp:
        mutant_sum = 0
        for sequence in gt_counts[epitope][day].keys():
            sequence = sequence.upper()
            changes[sequence] = []
            if re.match('[A-Z]', sequence):
                if not sequence in observed_seqs:
                    if len(observed_seqs) == 0:
                        origin[sequence] = 'wildtype'
                    elif len(observed_seqs) != 0:
                        hamming_dist = np.zeros(len(observed_seqs))
                        for si, seq in enumerate(observed_seqs):
                            for ch1, ch2 in zip(sequence, seq):
                                if ch1 != ch2:
                                    hamming_dist[si] += 1
                        parent_seq = observed_seqs[np.argmin(hamming_dist)]
                        origin[sequence] = "ancestor:" + names[parent_seq]
                        position = 0
                        for ch1, ch2 in zip(parent_seq, sequence):
                            position += 1
                            if sorted([ch1,ch2]) == ['A', 'G'] or sorted([ch1,ch2]) == ['C', 'T']:
                                transitions += 1
                                changes[sequence].append(ch1 + str(position) + ch2)
                            elif sorted([ch1,ch2]) == ['A', 'C'] or sorted([ch1,ch2]) == ['A', 'T'] or sorted([ch1,ch2]) == ['C', 'G'] or sorted([ch1,ch2]) == ['C', 'T']:
                                transversions += 1
                                changes[sequence].append(ch1 + str(position) + ch2)
                    if not sequence in names.keys():
                        names[sequence] = epitope + '_' + str(number_observed)
                        number_observed += 1
                    observed_seqs.append(sequence)
                temp_sequence = sequence
                if day != 0 and temp_sequence == gt_counts[epitope][0].keys()[0]:
                    temp_sequence = 'wildtype'
                dumpfile.write(str(day) +'\t' + epitope +'\t' + temp_sequence + '\t' + str(gt_counts[epitope][day][sequence])+ '\t' + names[sequence] + '\t' + origin[sequence] + '\t' + ' '.join(changes[sequence]) + '\n')
                if temp_sequence != 'wildtype':
                    mutant_sum += gt_counts[epitope][day][sequence]
        if day != 0:
            dumpfile.write(str(day) +'\t' + epitope +'\t' + 'all mutants' + '\t' + str(mutant_sum)+ '\n')
        dumpfile.write('\n')
    #for seq in observed_seqs:
    #    print seq
    
    print transitions, transversions
dumpfile.write('total transitions: ' + str(transitions) + '\n')
dumpfile.write('total transversions: ' + str(transversions) + '\n')
dumpfile.close()