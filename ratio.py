#!/bin/python
import numpy as np
import pandas as pd
import matplotlib.pyplot as pypl
from matplotlib import rcParams
from BootStrap3 import BootStrap
import scipy.optimize as syopt
# from formatting import err_brackets
import csv
import os
import time as tm
import json
from formatting import err_brackets

def ReadEvxptdump(file, par, boots, confs=0, times=64, number=0, bin=1):
    # par should be 0 or 1 depending on whether the real or imaginary part is chosen (0=R,1=I)
    # number defines which result you want, if there are multiple in the same file
    f=open(file)
    FF=[]
    for line in f:
        strpln=line.rstrip()
        if len(strpln) > 0:
            if strpln[0]=="+" and strpln[1]=="E" and strpln[2]=="N":
                tmp = f.readline()
                tmp = f.readline()
                tmp = f.readline().split()
                times = int(tmp[5])
            # if len(strpln)>3:
            #     print(strpln[4])
            if strpln[0]=="+" and strpln[1]=="R" and strpln[2]=="P" and int(strpln[4:6])==number:
                tmp = f.readline().split()
                while tmp[0] != "nmeas":
                    tmp = f.readline().split()
                confs = int(tmp[2])
                G = np.zeros(shape=(confs,times,2))
            if strpln[0]=="+" and strpln[1]=="R" and strpln[2]=="D" and int(strpln[4:6])==number:
                for iff in range(confs):
                    for nt in range(times):
                        tmp=f.readline().split()
                        G[iff,nt,0] = tmp[1]
                        G[iff,nt,1] = tmp[2]
    f.close()
    for j in range(times):
        FF.append(BootStrap(boots,68))
        FF[-1].Import(G[:,j,par], bin=bin)
        FF[-1].Stats()
    return np.array(FF)

def linear(x,a):
    return a[0]+a[1]*x

def oneexp(x,a):
    exp = (a[0]+a[1]*1j)*np.exp((a[2]+a[3]*1j)*x)
    return [np.real(exp), np.imag(real)]

def oneexpcos(x,a):
    exp = a[0]*np.exp(-a[1]*x)*np.cos(a[2]*x)
    return exp

def chisqfn(p,fnc,x,y,cminv):
    r=fnc(x,p)-y
    return np.matmul(r, np.matmul(cminv,r.T))

def minimizerBS(args):
    func,p0,fitfnc,x,data,cv,bounds = args
    reslist=[]
    for iboot in range(len(data)):
        res=syopt.minimize(func,p0,args=(fitfnc,x,data[iboot],cv),method='L-BFGS-B',bounds=bounds,options={'disp': False})
        reslist.append(res.x)
        #print(res.fun)
    return reslist


Nboot=100
lmb=1e-4
print('lambda = ', lmb)
#up-quark
#trev=0
corr0 = ReadEvxptdump('dump.res', 0, Nboot, number=0) + 1j*ReadEvxptdump('dump.res', 1, Nboot, number=0) #Unperturbed
corr1 = ReadEvxptdump('dump.res', 0, Nboot, number=1) + 1j*ReadEvxptdump('dump.res', 1, Nboot, number=1) #lambda=1e-4
#trev=1
corr2 = ReadEvxptdump('dump.res', 0, Nboot, number=2) + 1j*ReadEvxptdump('dump.res', 1, Nboot, number=2) #Unperturbed
corr3 = ReadEvxptdump('dump.res', 0, Nboot, number=3) + 1j*ReadEvxptdump('dump.res', 1, Nboot, number=3) #lambda=1e-4

oddratio=[]
ratiophase=[]
realratio=[]
for i in range(len(corr0)):
    oddratio.append((corr1[i]*corr2[i]*(corr0[i]*corr3[i])**(-1))**(1/2))
    oddratio[-1].Stats()
    ratiophase.append(BootStrap(oddratio[-1].nboot, oddratio[-1].confidence))
    ratiophase[-1].values = np.angle(oddratio[-1].values)/(-1)
    ratiophase[-1].Stats()
    realratio.append(BootStrap(oddratio[-1].nboot, oddratio[-1].confidence))
    realratio[-1].values = np.real(oddratio[-1].values)/(-1)
    realratio[-1].Stats()


#Fitting to the phase
p0 = [1.0,1.0]
x = np.arange(3,16)
data = np.array([y.values for y in ratiophase])[x].T
avgdata = np.array([y.Avg for y in ratiophase])[x].T
yerr = np.std(data,axis=0)
cv   = np.diag(yerr**(-2))
bounds = [(-1e3,1e3),(-1e3,1e3)]
args = (chisqfn,p0,linear,x,data,cv,bounds)
result = minimizerBS(args) #result[1] = Im{E_odd}
print(np.average(result,axis=0))
print(np.average(result,axis=0)/lmb)

# #Fitting to the real ratio
# p0 = [1.0,0.001,0.0001]
# x = np.arange(3,17)
# data = np.array([y.values for y in realratio])[x].T
# avgdata = np.array([y.Avg for y in realratio])[x].T
# yerr = np.std(data,axis=0)
# cv   = np.diag(yerr**(-2))
# bounds = [(-1e2,1e2),(-1e3,1e3),(-1e-1,1e-2)]
# args = (chisqfn,p0,oneexpcos,x,data,cv,bounds)
# resultreal = minimizerBS(args) #result[1] = Im{E_odd}
# print(np.average(resultreal,axis=0))


imenergy = [] #odd energy shift ( E(lmb)-E(-lmb) )
for i, time in enumerate(ratiophase[:-1]):
    imenergy.append( -1*(ratiophase[i]-ratiophase[i+1])*(lmb)**(-1) )
    imenergy[-1].Stats()

#Plots
ratiophaseavg = [i.Avg for i in ratiophase]
ratiophaseerr = [i.Std for i in ratiophase]
pypl.figure()
pypl.errorbar(np.arange(64),ratiophaseavg,ratiophaseerr,fmt=',', linewidth=1.0)
funcvals = [linear(x,pars) for pars in result]
avgvals = np.average(funcvals,axis=0)
errorvals = np.std(funcvals,axis=0)
pypl.plot(x,avgvals)
pypl.fill_between(x,avgvals-errorvals, avgvals+errorvals, alpha=0.3, linewidth=0)
#pypl.ylim([-0.0001,0.022])
pypl.ylim([-0.022,0.0001])
pypl.axhline(0,linewidth=0.2,color='k')
pypl.savefig('Phase.pdf')
# pypl.show()

# # Real ratio
# realratioavg = [i.Avg for i in realratio]
# realratioerr = [i.Std for i in realratio]
# pypl.figure()
# pypl.errorbar(np.arange(64),realratioavg,realratioerr,fmt=',', linewidth=1.0)
# funcvals = [oneexpcos(x,pars) for pars in resultreal]
# avgvals = np.average(funcvals,axis=0)
# print(avgvals)
# errorvals = np.std(funcvals,axis=0)
# pypl.plot(x,avgvals)
# pypl.fill_between(x,avgvals-errorvals, avgvals+errorvals, alpha=0.3, linewidth=0)
# pypl.ylim(-.01,0.01)
# pypl.savefig('realratio.pdf')
# # pypl.show()



imenergyavg = [i.Avg for i in imenergy]
imenergyerr = [i.Std for i in imenergy]
pypl.figure()
pypl.errorbar(np.arange(64-1),imenergyavg,imenergyerr,fmt=',', capsize=3)
averageE = np.average(result,axis=0)[1]/(lmb)
pypl.plot(x,len(x)*[averageE], linewidth=0.7)
errorrange = np.std(result,axis=0)[1]/(lmb)
pypl.fill_between(x,len(x)*[averageE-errorrange], len(x)*[averageE+errorrange], alpha=0.3, linewidth=0)
pypl.ylim(-3e0,3e0)
pypl.xlim(-1,30)
pypl.axhline(0,linewidth=0.2,color='k')
#textstr='\n'.join( [ paramnames[i]+' = '+err_brackets(float(ratioavg[i]), float(ratiostd[i]), form='Sci', texify=True) for i in range(len(ratioavg)) ] )
pypl.title('E='+err_brackets(averageE,errorrange))
pypl.savefig('Energy.pdf')
# pypl.show()




#/home/mischa/Documents/PhD/lattice_results/Feyn-Hell_kp122130kp121756/clover_nf2p1_feyn-hell/b5p65kp122130kp121756c2p4800_g8g2_smsnk30-48x96/nucleon/kp122130kp122130/p+0+0+0/q+0+0+0/rel/dump/quark1/dump.res

# ratiophase1=[]
# for i,j,k in zip(corr1,corr2,corr0):
#     num=i
#     den=j
#     rationum = np.real(num.values)*np.real(den.values)+np.imag(num.values)*np.imag(den.values) + 1j*(np.imag(num.values)*np.real(den.values)-np.real(num.values)*np.imag(den.values))
#     ratiophase1.append(BootStrap(k.nboot, k.confidence))
#     ratiophase1[-1].values = np.arctan2(np.real(rationum), np.imag(rationum))/(-2)
#     ratiophase1[-1].Stats()

# ratiophase1avg = [i.Avg for i in ratiophase1]
# ratiophase1err = [i.Std for i in ratiophase1]
# pypl.figure()
# pypl.errorbar(np.arange(64),ratiophase1avg,ratiophase1err,fmt='.')
# pypl.plot(x, linear(x,result))
# pypl.ylim(-0.001-0.5*np.pi,0.001-0.5*np.pi)
# #pypl.ylim(0.38,0.4)
# pypl.axhline(0,linewidth=0.2,color='k')
# pypl.show()

