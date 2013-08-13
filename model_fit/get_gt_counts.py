from Bio import AlignIO
import re

patient = 'CH40'

if patient == 'CH40':
    epitopes = ['nef191', 'gag125', 'gag405', 'vpr84', 'vif167', 'env145']
    tp = [0, 16, 45, 111, 181]
    foundername = 'CH40_con'

elif patient == 'CH58':
    epitopes = ['env578', 'env829', 'nef112', 'gag246']
    tp = [0, 9, 45, 85]
    foundername = 'CH58_transmitted'

elif patient == 'CH77':
    epitopes = ['tat61', 'env354', 'nef26', 'nef92']
    tp = [0, 14, 32]
    foundername = 'CH77_transmitted'


sequence_gt = {}
gt_counts = {}
for epitope in epitopes:
    if not epitope in gt_counts.keys():
        gt_counts[epitope] = {}
    gene = epitope[:3]
    loc = int(epitope[3:])
    alignment_loc = 'gt_data/alignments/' + patient + '_all_proteins/' + gene.upper() + '.AA.FASTA'
    alignment = AlignIO.read(alignment_loc, 'fasta')
    for seq in alignment.get_all_seqs():
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
        sequence = seq.seq.tostring()[loc-6:loc+6]
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

'''
seqfile = open('gt_data/' + patient + '_sequences.dat', 'wb')
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


dumpfile = open('gt_data/' + patient + '_genotypes_raw.dat', 'wb')
B = gt_counts.keys()
for epitope in epitopes:
    for day in tp:
        mutant_sum = 0
        for sequence in gt_counts[epitope][day].keys():
            if re.match('[A-Z]', sequence):
                temp_sequence = sequence
                if day != 0 and temp_sequence == gt_counts[epitope][0].keys()[0]:
                    temp_sequence = 'wildtype'
                dumpfile.write(str(day) +'\t' + epitope +'\t' + temp_sequence + '\t' + str(gt_counts[epitope][day][sequence])+ '\n')
                if temp_sequence != 'wildtype':
                    mutant_sum += gt_counts[epitope][day][sequence]
        if day != 0:
            dumpfile.write(str(day) +'\t' + epitope +'\t' + 'all mutants' + '\t' + str(mutant_sum)+ '\n')
        dumpfile.write('\n')
dumpfile.close()
'''