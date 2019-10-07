import sys
import numpy as np
import copy

from Network import networkgen
from Probselection import Probtarget

N = 100
E = 200

b=10
c = int(sys.argv[1])*0.1
R = b-c
S = -c
T = b
P = 0
d = int(sys.argv[2])*0.001
tau = -2+int(sys.argv[3])*0.2
p = int(sys.argv[4])*0.25
q = int(sys.argv[5])*0.25
TP = 0
FN = 0
FP = 0
TN = 0

m = 0.0001
t0 = 10000
Tt = 100000000

state = [0]*N
payoff = [0.0]*N
prosperity = [0.0]*N

network = networkgen(N,E)
Neigh = [[] for i in range(N)]
for i in range(N):
    for j in range(N):
        if network[i][j]==1:
            Neigh[i].append(j)

Pi = [[R,S],[T,P]]

for i in range(N):
    for j in Neigh[i]:
        payoff[i] = payoff[i] + Pi[state[i]][state[j]]
    prosperity[i] = pow(1+d,payoff[i])

t = 0
mutant = False
transition = False
transitionStart = False

transitionNum = 0
avecoop = 0.0
avedegree = 0.0
aveprosp = 0.0
avepay = 0.0

numt = 0
Avecoop = [0.0]*int(Tt/10000)
Avedegree = [0.0]*int(Tt/10000)
Aveprosp = [0.0]*int(Tt/10000)
Avepay = [0.0]*int(Tt/10000)

coop = 0.0
degree = 0.0
prosp = 0.0
pay = 0.0

Ncascadeid = []
Pcascadeid = []

NodeNcascadeid = np.zeros(N)
NodePcascadeid = np.zeros(N)

Ncascadesize = []
Pcascadesize = []

while t<Tt:
    i = np.random.randint(0,N)

    rand = np.random.rand(1)

    j = Probtarget(rand, prosperity)

    tempstatei = 0
    if t<=t0:
        tempstatei = state[j]
    else:
        rand = np.random.rand(1)
        if rand < m:
            tempstatei = 1 - state[j]
        else:
            tempstatei = state[j]

    tempneigh = []
    tempneigh1 = []

    tempneigh1.append(j)
    for k in Neigh[j]:
        tempneigh1.append(k)

    NCas = False
    PCas = False
    for k in tempneigh1:
        if state[k] == 0:
            s = np.random.normal(-0.5,np.sqrt(0.5),1)
        else:
            s = np.random.normal(0.5,np.sqrt(0.5),1)

        priv_connect = False
        if s[0]<tau:
            priv_connect = True

        pub_connect = False
        if len(Neigh[k])>avedegree:
            pub_connect = True

        choice = False
        if priv_connect and pub_connect:
            choice = True
        if not priv_connect and not pub_connect:
            choice = False
        if not priv_connect and pub_connect:
            if np.random.uniform(0, 1)<p:
                choice = True
        if priv_connect and not pub_connect:
            if np.random.uniform(0, 1)<q:
                choice = True

        if choice == True:
            if state[k] == 0:
                TP += 1
            else:
                FP += 1
                if not priv_connect:
                    NCas = True
        else:
            if state[k] == 0:
                FN += 1
                if priv_connect:
                    PCas = True
            else:
                TN += 1

        if choice and i!=k:
            tempneigh.append(k)

    if t>=t0:

        if NCas == True:
            nid = int(NodeNcascadeid[j])
            if nid == 0:
                nid = int(len(Ncascadeid)+1)
                Ncascadeid.append(nid)
                Ncascadesize.append(0)
            NodeNcascadeid[i] = nid
            Ncascadesize[nid-1] += 1
        else:
            NodeNcascadeid[i] = 0

        if PCas == True:
            pid = int(NodePcascadeid[j])
            if pid == 0:
                pid = int(len(Pcascadeid)+1)
                Pcascadeid.append(pid)
                Pcascadesize.append(0)
            NodePcascadeid[i] = pid
            Pcascadesize[pid-1] += 1
        else:
            NodePcascadeid[i] = 0

    tempneighi = copy.deepcopy(Neigh[i])
    for k in tempneighi:
        tempid = Neigh[k].index(i)
        del Neigh[k][tempid]
        payoff[k] = payoff[k]-Pi[state[k]][state[i]]
        prosperity[k] = pow(1+d,payoff[k])
        tempid = Neigh[i].index(k)
        del Neigh[i][tempid]
        network[i][k] = 0
        network[k][i] = 0

    if len(Neigh[i])!=0:
        print("System error!!!")

    if t<t0:
        TP = 0
        FN = 0
        FP = 0
        TN = 0

    state[i] = tempstatei

    payoff[i] = 0.0
    prosperity[i]=1.0
    for k in tempneigh:
        Neigh[i].append(k)
        Neigh[k].append(i)
        network[i][k] = 1
        network[k][i] = 1
        payoff[k] = payoff[k]+Pi[state[k]][state[i]]
        prosperity[k] = pow(1+d, payoff[k])
        payoff[i] = payoff[i]+Pi[state[i]][state[k]]
        prosperity[i] = pow(1+d, payoff[i])

    if sum(state)!=0 and transitionStart==False and transition== False:
        transitionStart = True

    if sum(state)==N and transitionStart == True and transition == False:
        transitionStart = False
        transition = True
        transitionNum = transitionNum + 1

    if sum(state)==0 and transitionStart == False and transition == True:
        transitionStart = True
        transition = False
        transitionNum = transitionNum + 1

    for k in range(N):
        avedegree = avedegree + len(Neigh[k])
        aveprosp = aveprosp + prosperity[k]
        avepay = avepay + payoff[k]

    avecoop = N-sum(state)
    avedegree = avedegree/N
    aveprosp = aveprosp/N
    avepay = avepay/N

    if t>=t0:
        coop += avecoop
        degree += avedegree
        prosp += aveprosp
        pay += avepay

    # if t%10000 == 0:
    #     Avecoop[numt] = avecoop
    #     Avedegree[numt] = avedegree
    #     Aveprosp[numt] = aveprosp
    #     Avepay[numt] = avepay
    #     numt += 1
    # if t%10000 == 0 and t>10000:
    #     print(Ncascadesize)
    #     print(Ncascadeid)
    #     print(Pcascadesize)
    #     print(Pcascadeid)

 #   print("time="+str(t),'q='+str(q),'theta='+str(theta),avecoop)
 #   print("time="+str(t),'q='+str(q),'theta='+str(theta),avecoop, avedegree, aveprosp,transitionNum)

    t = t+1
    avecoop = 0.0
    avedegree = 0.0
    aveprosp = 0.0
    avepay = 0.0

coop /= (Tt-t0)
degree /= (Tt-t0)
prosp /= (Tt-t0)
pay /= (Tt-t0)

maxncasc = 0
maxpcasc = 0
sumncasc = 0
sumpcasc = 0
numncasc = 0
numpcasc = 0

if len(Ncascadeid)>0:
    maxncasc = max(Ncascadesize)
    sumncasc = sum(Ncascadesize)
    numncasc = len(Ncascadesize)
if len(Pcascadeid)>0:
    maxpcasc = max(Pcascadesize)
    sumpcasc = sum(Pcascadesize)
    numpcasc = len(Pcascadesize)

print(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],coop,degree,prosp,pay,transitionNum,TP,FN,FP,TN,maxncasc,maxpcasc,sumncasc,sumpcasc,numncasc,numpcasc)

filename = 'Output'+str(sys.argv[1])+'and'+str(sys.argv[2])+'and'+str(sys.argv[3])+'and'+str(sys.argv[4])+'and'+str(sys.argv[5])
f = open(filename,'a')
f.write(str(coop)+'\t'+str(degree)+'\t'+str(prosp)+'\t'+str(pay)+'\t'+str(transitionNum)+'\n')
f.write(str(TP)+'\t'+str(FN)+'\t'+str(FP)+'\t'+str(TN)+'\t'+str(maxncasc)+'\t'+str(maxpcasc)+'\t'+str(sumncasc)+'\t'+str(sumpcasc)+'\t'+str(numncasc)+'\t'+str(numpcasc)+'\n')
f.close()

#f1 = 'Traj'+str(sys.argv[1])+str(sys.argv[2])+'.txt'
#np.savetxt(f1,[Avecoop,Avedegree,Aveprosp,Avepay])
