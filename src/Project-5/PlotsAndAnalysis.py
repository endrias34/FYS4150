import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def get_data(filename):
    data = np.loadtxt(filename, dtype = np.int8)
    return data
    
#black = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/black.txt')
#black = np.mean(black,1)
#white = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/white.txt')
#white = np.mean(white,1)
#equal = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/equal.txt')
#equal = np.mean(equal,1)
#x=np.linspace(0,100000,1001)
#plt.plot(x,black,x,white,x,equal)
#plt.xlabel('Monte Carlo cycles elapsed', fontsize=16)
#plt.ylabel('Mean Magnetization', fontsize=16)
#plt.savefig('figures/mag.png',bbox_inches = 'tight')

#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/Result1000Agents100000McCycles610.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,1000])
#fig.set_size_inches(10, 10)
#plt.title('Opinion evolution simulation', fontsize=22)
#plt.xlabel('Spin number', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
#ax.set_xticklabels([0,2000,4000,6000,8000,10000], fontsize=18)
#ax.set_yticklabels([100000,80000,60000,40000,20000,0], fontsize=18)
#plt.savefig('figures/large.png', dpi=500)
#plt.show()


#def get_data(filename):
#    data = np.loadtxt(filename, dtype = np.int16)
#    return data
#
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult1000Agents100000McCycles125.txt')
#MData=MData[0:1000]/1000
#MData0=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult10Agents1000McCycles661.txt')
#MData0=MData0[0:1000]/10
#MData2=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult100Agents10000McCycles95.txt')
#MData2=MData2[0:1000]/100
#MData3=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult4Agents400McCycles152.txt')
#MData3=MData3[0:1000]/4
#
#
#
#plt.hist([MData,MData2,MData0,MData3], label=["N = 1000","N = 100","N = 10","N = 4"])
#plt.legend(loc='upper left', fontsize=16)
#plt.xlabel('Final Mean Magnetization', fontsize=16)
#plt.ylabel('Occurence per 1.000', fontsize=16)
#plt.savefig('figures/histogram.png')
#plt.show()
#
#
#
#white=np.count_nonzero(MData == 1)
#black=np.count_nonzero(MData == -1)
#equal=np.count_nonzero(MData == 0)
#unresolved=1000-(black+white+equal)
#
#whitefrac=white/(1000-unresolved)
#blackfrac=black/(1000-unresolved)
#equalfrac=equal/(1000-unresolved)

#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/correlate.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,1000])
#fig.set_size_inches(10, 10)
#plt.title('Opinion evolution simulation', fontsize=22)
#plt.xlabel('Spin number', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
#ax.set_xticklabels([0,200,400,600,800,1000], fontsize=18)
#ax.set_yticklabels([1000,800,600,400,200,0], fontsize=18)
#plt.savefig('figures/correlate.png', dpi=500, bbox_inches = 'tight')
#plt.show()
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/correlate.txt')
#
#m=np.mean(data,1)
#x=np.linspace(0,10000,10001)
#plt.plot(x,m)
#plt.xlabel('Monte Carlo cycles elapsed', fontsize=16)
#plt.ylabel('Mean Magnetization', fontsize=16)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.savefig('figures/corrmag.png',bbox_inches = 'tight')
##
#mcorr=np.correlate(m, m, mode='full')
#plt.plot(mcorr[10000:20001] / float(mcorr.max()))
##plt.plot(x,mcorr[1000:2001])
#plt.xlabel('$\Delta t$', fontsize=16)
#plt.ylabel('Correlation', fontsize=16)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.savefig('figures/corrCurv.png',bbox_inches = 'tight')
#
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/bigcorrelate.txt')
##data = data[0:-1,:]
##data[data < 0] = 0
##fig, ax = plt.subplots(1,1)
##img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,10000,0,10000])
##fig.set_size_inches(10, 10)
##plt.title('Opinion evolution simulation', fontsize=22)
##plt.xlabel('Spin number', fontsize=22)
##plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
##ax.set_xticklabels([0,2000,4000,6000,8000,10000], fontsize=18)
##ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
##plt.savefig('figures/bigcorrelate.png', dpi=500, bbox_inches = 'tight')
#
##
#tauList=[]
#for j in range(len(data[1,:])):
#    counter=0
#    tau=1
#    for i in range(len(data[:,1])-1):
#        if data[i,j]*data[i+1,j]>0:
#            tau += 1
#        else:
#            counter +=1
#            if counter>0:
#                tauList.append(tau)
#                tau = 1
##        if i == len(data[:,1])-2:
##            tauList.append(tau)
##            tau = 1
##            print(i)
##            print('\n')
#tauList=np.array(tauList)            
#tauDist=np.zeros(np.max(tauList))
#for i in range(np.max(tauList)):
#    tauDist[i]=np.count_nonzero(tauList == i+1)#/np.sum(tauList)
#
#pTau=tauDist/np.sum(tauDist)
#tau=np.linspace(1,len(tauDist),len(tauDist))  
#
#
#
#end=600
#
#m, c = np.polyfit(np.log(tau[0:end]),np.ma.log(tauDist[0:end]), 1) # fit log(y) = m*log(x) + c
#y_fit = np.exp(m*np.log(tau) + c) # calculate the fitted values of y 
#
#plt.loglog(tau,tauDist,'.r')
#plt.plot(tau[0:end], y_fit[0:end])
#
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)
#plt.xlabel('$\\tau$', fontsize=18)
#plt.ylabel('P($\\tau$)', fontsize=18)
#
#plt.savefig('figures/power.png',bbox_inches = 'tight')
#
#m, c = np.polyfit(np.log(tau[0:end]),np.ma.log(pTau[0:end]), 1) # fit log(y) = m*log(x) + c
#y_fit = np.exp(m*np.log(tau) + c) # calculate the fitted values of y 
#
#plt.loglog(tau,pTau,'.r')
#plt.plot(tau[0:end], y_fit[0:end])
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)
#plt.xlabel('$\\tau$', fontsize=18)
#plt.ylabel('P($\\tau$)', fontsize=18)
#plt.savefig('figures/power2.png',bbox_inches = 'tight')
    
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MccResult0.710000cB100Agents10000McCycles857.txt')
#MData=MData[0:1000]/100
#MData0=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MccResult0.690000cB100Agents10000McCycles865.txt')
#MData0=MData0[0:1000]/100
#MData2=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult100Agents10000McCycles95.txt')
#MData2=MData2[0:1000]/100
#MData3=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/MResult4Agents400McCycles152.txt')
#MData3=MData3[0:1000]/4

#plt.hist([MData,MData0])#,MData0,MData3], label=["N = 1000","N = 100","N = 10","N = 4"])
#plt.legend(loc='upper left', fontsize=16)
#plt.xlabel('Final Mean Magnetization', fontsize=16)
#plt.ylabel('Occurence per 1.000', fontsize=16)
#plt.savefig('figures/histogram2.png')

#
#dist=np.zeros((9,3))
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/1.txt')
#MData=MData[0:1000]/100
#dist[0,0]=np.count_nonzero(MData == -1)/1000
#dist[0,1]=np.count_nonzero(MData ==  0)/1000
#dist[0,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/2.txt')
#MData=MData[0:1000]/100
#dist[1,0]=np.count_nonzero(MData == -1)/1000
#dist[1,1]=np.count_nonzero(MData ==  0)/1000
#dist[1,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/3.txt')
#MData=MData[0:1000]/100
#dist[2,0]=np.count_nonzero(MData == -1)/1000
#dist[2,1]=np.count_nonzero(MData ==  0)/1000
#dist[2,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/4.txt')
#MData=MData[0:1000]/100
#dist[3,0]=np.count_nonzero(MData == -1)/1000
#dist[3,1]=np.count_nonzero(MData ==  0)/1000
#dist[3,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/5.txt')
#MData=MData[0:1000]/100
#dist[4,0]=np.count_nonzero(MData == -1)/1000
#dist[4,1]=np.count_nonzero(MData ==  0)/1000
#dist[4,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/6.txt')
#MData=MData[0:1000]/100
#dist[5,0]=np.count_nonzero(MData == -1)/1000
#dist[5,1]=np.count_nonzero(MData ==  0)/1000
#dist[5,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/7.txt')
#MData=MData[0:1000]/100
#dist[6,0]=np.count_nonzero(MData == -1)/1000
#dist[6,1]=np.count_nonzero(MData ==  0)/1000
#dist[6,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/8.txt')
#MData=MData[0:1000]/100
#dist[7,0]=np.count_nonzero(MData == -1)/1000
#dist[7,1]=np.count_nonzero(MData ==  0)/1000
#dist[7,2]=np.count_nonzero(MData ==  1)/1000
#MData=get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/concentration/9.txt')
#MData=MData[0:1000]/100
#dist[8,0]=np.count_nonzero(MData == -1)/1000
#dist[8,1]=np.count_nonzero(MData ==  0)/1000
#dist[8,2]=np.count_nonzero(MData ==  1)/1000
#
#
#plt.plot(np.linspace(0.1,0.9,9),dist[:,0], '.b',np.linspace(0.1,0.9,9),dist[:,1], '*r',np.linspace(0.1,0.9,9),dist[:,2], '+g')
#plt.legend(["AAAA","ABAB","BBBB"],loc='upper center', fontsize=14)#, frameon=False)
#plt.xlabel('Concentration initated as B', fontsize=16)
#plt.ylabel('Fraction of 1000 outcomes', fontsize=16)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.savefig('figures/CBdist.png',bbox_inches = 'tight')



#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/Result1000Agents10000McCycles350.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,10000])
#ax.set_aspect(aspect =0.2)
##fig.set_size_inches(10, 30)
#plt.title('$p=10^{-5}$', fontsize=22)
#plt.xlabel('Spin #', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
##ax.set_xticklabels([0, *, 400, *,800 ], fontsize=18)
#ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
#plt.savefig('figures/p10-5.png', dpi=500, bbox_inches = 'tight')
#plt.show()
#
#    
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/Result1000Agents10000McCycles75.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,10000])
#ax.set_aspect(aspect =0.2)
##fig.set_size_inches(10, 30)
#plt.title('$p=10^{-6}$', fontsize=22)
#plt.xlabel('Spin #', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
##ax.set_xticklabels([0, *, 400, *,800 ], fontsize=18)
#ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
#plt.savefig('figures/p10-6.png', dpi=500, bbox_inches = 'tight')
#plt.show()
#    
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/Result1000Agents10000McCycles143.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,10000])
#ax.set_aspect(aspect =0.2)
##fig.set_size_inches(10, 30)
#plt.title('$p=3*10^{-6}$', fontsize=22)
#plt.xlabel('Spin #', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
##ax.set_xticklabels([0, *, 400, *,800 ], fontsize=18)
#ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
#plt.savefig('figures/p310-4.png', dpi=500, bbox_inches = 'tight')
#plt.show()
#
#
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/p-4.txt')
#
#tauList=[]
#for j in range(len(data[1,:])):
#    counter=0
#    tau=1
#    for i in range(len(data[:,1])-1):
#        if data[i,j]*data[i+1,j]>0:
#            tau += 1
#        else:
#            counter +=1
#            if counter>0:
#                tauList.append(tau)
#                tau = 1
##        if i == len(data[:,1])-2:
##            tauList.append(tau)
##            tau = 1
##            print(i)
##            print('\n')
#tauList=np.array(tauList)            
#tauDist=np.zeros(np.max(tauList))
#for i in range(np.max(tauList)):
#    tauDist[i]=np.count_nonzero(tauList == i+1)#/np.sum(tauList)
#
#pTau4=tauDist/np.sum(tauDist)
#
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/p-3.txt')
#
#tauList=[]
#for j in range(len(data[1,:])):
#    counter=0
#    tau=1
#    for i in range(len(data[:,1])-1):
#        if data[i,j]*data[i+1,j]>0:
#            tau += 1
#        else:
#            counter +=1
#            if counter>0:
#                tauList.append(tau)
#                tau = 1
##        if i == len(data[:,1])-2:
##            tauList.append(tau)
##            tau = 1
##            print(i)
##            print('\n')
#tauList=np.array(tauList)            
#tauDist=np.zeros(np.max(tauList))
#for i in range(np.max(tauList)):
#    tauDist[i]=np.count_nonzero(tauList == i+1)#/np.sum(tauList)
#
#pTau3=tauDist/np.sum(tauDist)
#
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/p-2.txt')
#
#tauList=[]
#for j in range(len(data[1,:])):
#    counter=0
#    tau=1
#    for i in range(len(data[:,1])-1):
#        if data[i,j]*data[i+1,j]>0:
#            tau += 1
#        else:
#            counter +=1
#            if counter>0:
#                tauList.append(tau)
#                tau = 1
##        if i == len(data[:,1])-2:
##            tauList.append(tau)
##            tau = 1
##            print(i)
##            print('\n')
#tauList=np.array(tauList)            
#tauDist=np.zeros(np.max(tauList))
#for i in range(np.max(tauList)):
#    tauDist[i]=np.count_nonzero(tauList == i+1)#/np.sum(tauList)
#
#pTau2=tauDist/np.sum(tauDist)
#
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/p-1.txt')
#
#tauList=[]
#for j in range(len(data[1,:])):
#    counter=0
#    tau=1
#    for i in range(len(data[:,1])-1):
#        if data[i,j]*data[i+1,j]>0:
#            tau += 1
#        else:
#            counter +=1
#            if counter>0:
#                tauList.append(tau)
#                tau = 1
##        if i == len(data[:,1])-2:
##            tauList.append(tau)
##            tau = 1
##            print(i)
##            print('\n')
#tauList=np.array(tauList)            
#tauDist=np.zeros(np.max(tauList))
#for i in range(np.max(tauList)):
#    tauDist[i]=np.count_nonzero(tauList == i+1)#/np.sum(tauList)
#
#pTau1=tauDist/np.sum(tauDist)
#
#
#tau=np.linspace(1,10000,10000)  
#end=300
##
##plt.loglog(tau[0:80],pTau1[0:80],'.r')
##plt.loglog(tau[0:391],pTau2[0:391],'.b')
##plt.loglog(tau[0:1988],pTau3[0:1988],'.g')
#plt.plot(tau[0:9329],pTau4[0:9329],'.y',tau[0:1988],pTau3[0:1988],'.g',tau[0:391],pTau2[0:391],'.b',tau[0:80],pTau1[0:80],'.r')
#
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)
#plt.xlabel('$\\tau$', fontsize=18)
#plt.ylabel('P($\\tau$)', fontsize=18)
#plt.legend(["$p=10^{-4}$","$p=10^{-3}$","$p=10^{-2}$","$p=10^{-1}$"],loc='upper right', fontsize=14)#, frameon=False)
#
#plt.savefig('figures/power4.png',bbox_inches = 'tight')

#
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,1000,0,10000])
#ax.set_aspect(aspect =0.2)
##fig.set_size_inches(10, 30)
#plt.title('$p=3*10^{-6}$', fontsize=22)
#plt.xlabel('Spin #', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
##ax.set_xticklabels([0, *, 400, *,800 ], fontsize=18)
#ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
#    
#data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/pResult10000Agents10000McCycles458.txt')
#data = data[0:-1,:]
#data[data < 0] = 0
#fig, ax = plt.subplots(1,1)
#img = ax.imshow(data,cmap=plt.cm.gray, extent=[0,10000,0,10000])
#fig.set_size_inches(10, 10)
#plt.title('Opinion evolution simulation', fontsize=22)
#plt.xlabel('Spin number', fontsize=22)
#plt.ylabel('Monte Carlo cycles elapsed', fontsize=22)
#ax.set_xticklabels([0,2000,4000,6000,8000,10000], fontsize=18)
#ax.set_yticklabels([10000,8000,6000,4000,2000,0], fontsize=18)
#plt.savefig('Figures/periodic1.png', dpi=500, bbox_inches = 'tight')

#
data = get_data('C:/Users/jbrod/Dropbox/FAM/Phd/Fys4150/IsingSocial/code/Results/pResult100000Agents10000000McCycles971.txt')
lmdList=np.zeros((len(data[:,1]),len(data[1,:])))
for i in range(len(data[:,1])):

#    lmd=1
#    for j in range(len(data[1,:])-1):
#        if data[i,j-1]+data[i,j] == data[i,j]+data[i,j+1]:
#            lmd += 1
#        else:
#            if lmd > 1:
#                lmd += 1
#            lmdList[i,j]=lmd
#            lmd = 1
##        if j == len(data[:,1])-2:
##            if data[i,j-1]+data[i,j] == data[i,j]+data[i,0]:
##                if lmd > 1:
##                    lmd += 1
##                data[i,0]+=lmd
##            else:
##                if lmd > 1:
##                    lmd += 1
##                lmdList[i,j]=lmd
##                lmd = 1
    
    color = 3
    counter = 0
    lmd = 1
    index = 0
    while counter < len(data[1,:]) - 1:
        oldcolor = color
        if data[i,counter] + data[i,counter+1] == -2:
            color = -1
            if color == oldcolor:
                lmd += 1
            elif counter > 0:
                lmdList[i,counter]=lmd
                lmd = 1
        if data[i,counter] + data[i,counter+1] == 0:
            color = 0
            if color == oldcolor:
                lmd += 1
            elif counter > 0:
                lmdList[i,counter]=lmd
                lmd = 1
        if data[i,counter] + data[i,counter+1] == 2:
            color = 1
            if color == oldcolor:
                lmd += 1
            elif counter > 0:
                lmdList[i,counter]=lmd
                lmd = 1
        counter += 1
        if  counter == len(data[1,:])-1:
            while lmdList[i,index] == 0:
                index += 1
                if lmdList[i,index] > 0:              
                    lmdList[i,index] += lmd
                    lmd = 1
            






lmdDist=np.mean(np.ma.masked_equal(lmdList[0:1000], 0), axis=1)           
plt.plot(np.linspace(1000,1000000,1000),lmdDist)
plt.title('Cluster size evolution', fontsize=16)
plt.xlabel('Monte Carlo cycles elapsed', fontsize=16)
plt.ylabel('Mean cluster size', fontsize=16)
plt.xticks(length=200000, fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('Figures/meanCluster2.png', dpi=500, bbox_inches = 'tight')


#row0=lmdList[500,:]
#row500=lmdList[500,:]
#row1000=lmdList[999,:]
#plt.hist([row500[row500>0]])#,row0[row0>0],row1000[row1000>0]], bins =100)
#
#row0=np.ma.masked_equal(lmdList[0,:], 0)
#row500=np.ma.masked_equal(lmdList[500,:], 0)
#row1000=np.ma.masked_equal(lmdList[999,:], 0)




#plt.hist(np.ma.masked_equal(lmdList[500,:], 0))
#lmdDist=np.zeros(np.max(lmdList))
#for i in range(np.max(lmdList)):
#    lmdDist[i]=np.count_nonzero(lmdList == i+1)#/np.sum(lmdList)
#
#plmd1=lmdDist/np.sum(lmdDist)


#lmd=np.linspace(1,10000,10000)  
