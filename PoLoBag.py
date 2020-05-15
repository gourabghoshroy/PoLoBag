#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# PoLoBag algorithm for d=2 polynomial features

import numpy as np
from scipy import stats
from sklearn import linear_model
import warnings
from sklearn.exceptions import ConvergenceWarning


#Input expression and output sorted edgelist file names
infilename = 'input/A_matrix.txt'
outfilename = 'A_result_network.txt'


#Algorithm parameters
n12 = 0.5    # controls number of linear features in each bootstrap sample
n22 = 3.5    # controls number of nonlinear features in each bootstrap sample
nM = 0.5     # controls bootstrap sample size
nB = 500     # total number of bootstrap samples in the ensemble
alpha = 0.1  # the Lasso regularization parameter 


#For repeatability
#np.random.seed(30)

#Ignore Lasso convergence warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)

#Algorithm

#Read input expression file
D = {}
genes = []
with open(infilename, 'r') as f:
    geneCount = -1
    for line in f:
        geneCount += 1
        if geneCount == 0:   #header
            continue
        lineData = line.rstrip().split("\t")   #tab separated input file
        genes.append(lineData[0])
        D[lineData[0]] = stats.zscore(np.array(lineData[1:]).astype(np.float))
genes = np.unique(genes)


regs = genes                                   #Potential regulators

edges = []
w = []
for t in genes:  #Regression problem for each target gene
    yt = np.transpose(D[t][np.newaxis]) 
    Xt = yt[:,[]]
    for reg in regs:
        if not reg == t:  # no auto-regulation
            Xt = np.hstack((Xt,np.transpose(D[reg][np.newaxis])))
    ntR = Xt.shape[1]
    if ntR == 0:
        continue
    valm = np.median(np.hstack((yt,Xt)))    
    tau = 3
    yt = yt - tau*valm
    Xt = Xt - tau*valm
    nlin = int(n12*np.sqrt(ntR))
    nnlin = int(n22*np.sqrt(ntR))
    wt = np.zeros((ntR,1))                     # Weight
    wtM = np.zeros((ntR,1))                    # Weight magnitude
    wtS = np.zeros((ntR,1))                    # Weight sign
    stw = np.zeros((ntR,1)) + 1e-10            # How many times a feature was chosen in a sample
    m = Xt.shape[0]
    for b in range(nB):   #For each sample
        rowindex = np.random.choice(m,int(nM*m)) 
        ytb = yt[rowindex,:]
        XtbF = Xt[rowindex,:]
        # Separate idtb for linear and nonlinear features
        idtblin = []           
        idtbnlin = []
        if nlin > 0:
            idtblin = np.random.choice(ntR,nlin,replace=False)  
        Xtb = XtbF[:,idtblin]   # linear features
        if nnlin > 0:
            allindex = [x for x in range(ntR) if x not in idtblin]
            index = np.random.choice(allindex,2*nnlin,replace=False)
            Xtb21F = XtbF[:,index[:nnlin]]
            Xtb22F = XtbF[:,index[nnlin:]]
            vala = Xtb21F*Xtb22F
            vals1 = np.sign(Xtb21F)
            vals12 = np.sign(vala)
            atbk = 0.5*np.sqrt(np.abs(vala))
            ftb21 = vals1*(1+vals12)*atbk         
            ftb22 = vals1*(1-vals12)*atbk         
            Xtb = np.hstack((Xtb,ftb21[:,:int(0.5*nnlin)],ftb22[:,int(0.5*nnlin):]))    # nonlinear features 
            for l in range(int(0.5*nnlin)):
                idtbnlin.append(str(index[l]) + '#' + str(index[l+nnlin]))
            for l in range(int(0.5*nnlin),nnlin):
                idtbnlin.append('-' + str(index[l]) + '#' + str(index[l+nnlin]))        # an indicator for category
                        
        clf = linear_model.Lasso(alpha=alpha,fit_intercept=False)  # Lasso 
        clf.fit(Xtb, ytb)
        wtb = clf.coef_
        if len(idtblin) > 0:                            # for linear features
            wtM[idtblin,0] += np.abs(wtb[0:len(idtblin)])   
            wtS[idtblin,0] += wtb[0:len(idtblin)]        
            stw[idtblin,0] += 1                                   
        for l in range(len(idtbnlin)):                   # for nonlinear features
            indi = idtbnlin[l].split("#")
            valai = np.sqrt(np.abs(wtb[l+len(idtblin)]))
            valsi = np.sign(wtb[l+len(idtblin)])*valai
            if indi[0][0] == '-':                       
                wtM[int(indi[0][1:]),0] += valai
                wtS[int(indi[0][1:]),0] += valsi
                stw[int(indi[0][1:]),0] += 1
                wtM[int(indi[1]),0] += valai
                wtS[int(indi[1]),0] += -valsi
                stw[int(indi[1]),0] += 1
            else:                                      
                wtM[int(indi[0]),0] += valai
                wtS[int(indi[0]),0] += valsi
                stw[int(indi[0]),0] += 1
                wtM[int(indi[1]),0] += valai
                wtS[int(indi[1]),0] += valsi
                stw[int(indi[1]),0] += 1
    
    wt = np.sign(wtS)*wtM/stw    # Compute weights as average
    j = 0
    for reg in regs:
        if not reg == t:
            val = wt[j,0]
            if (abs(val) > 0.0):
                edges.append(reg+"\t"+t)
                w.append(val)
            j += 1

sortindex = np.flip(np.argsort(np.abs(w))[np.newaxis],1).astype(np.int)   # Sort by absolute value of edge weights
sedgevals = [w[s] for s in sortindex[0,:]]
sedgevals = sedgevals/abs(sedgevals[0])
sedges = [edges[s] for s in sortindex[0,:]]              
ofile = open(outfilename, 'w')
for s in range(len(sedges)):
    print(sedges[s] + "\t" + str(sedgevals[s]),file=ofile)                # write to output edge list file
        
ofile.close()
