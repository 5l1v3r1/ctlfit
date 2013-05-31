'''
authors:  Taylor Kessinger, Richard Neher
date:     05/02/2013
content:  Fits a single patient and produces likelihood plots based on one set of parameters.
          from fit_patients import fit_patient_LH is needed.
'''

#!/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/bin/python
import numpy as np
import matplotlib
matplotlib.use('PDF')
import pylab as plt
import sys
sys.path.append('src/')
sys.path.append('model_fit/')
import cPickle as pickle
#import ctl_fit
import os
import argparse
import ctlutils

params = {'backend': 'ps',  
          'axes.labelsize': 36, 
          'axes.titlesize': 36,
          'text.fontsize': 36,
'font.sans-serif': 'Helvetica',
'legend.fontsize': 25,
'xtick.labelsize': 30,
'ytick.labelsize': 30,
'lines.linewidth': 3,
'text.usetex': True}
plt.rcParams.update(params)

col=['k', 'r', 'b', 'g', 'c', 'm','y']
mar=['o','x','v','s','d','^']

patient = sys.argv[1]

sample_sizes = np.loadtxt('gt_data/'+patient+'_sample_sizes_alleles.txt')[:,1:]
allele_file = open('gt_data/'+patient+'_alleles.txt', 'r')
alleles = allele_file.readline().strip().split()[1:]
A= np.loadtxt('gt_data/'+patient+'_alleles.txt')
afreqs=A[:,1:]
tp = A[:,0]

fig=plt.figure(figsize=(8,8))
margins = [0.18, 0.14, 0.77, 0.8]
axes = plt.Axes(fig, margins) 
fig.add_axes(axes)
for ai, allele in enumerate(alleles):
    plt.plot(tp, afreqs[:, ai]/sample_sizes[:,ai], marker=mar[ai%len(mar)],markeredgecolor=col[ai%len(col)+1],markeredgewidth=3,ms=12,linestyle='--', lw=3,label=alleles[ai],c=col[ai%len(col)+1])

ax=plt.gca()
ax.set_ylim([-.1,1.1])
plt.axhline(y=0,c='k', lw=2)
plt.axhline(y=1, c='k', lw=2)
plt.legend(loc=4)
#plt.title(patient + ' allele frequencies')
if patient == 'CH40':
    plt.xticks(np.arange(0,200,40))
elif patient == 'CH58':
    plt.xticks(np.arange(0,90,20))
elif patient == 'CH77':
    plt.xticks(np.arange(0,40,10))
plt.xlabel(r'time ($\textrm{days}$)')
plt.ylabel('allele frequencies')
plt.subplots_adjust(bottom=0.15)
fig.text(0.02,0.88,'B', fontsize=66)
plt.show()

plt.savefig('figures_manuscript/' + patient + '_alleles.pdf')
#plt.close()