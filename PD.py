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
tauc = -2+int(sys.argv[3])*0.2
taud = -2+int(sys.argv[3])*0.2
TPC = 0
FNC = 0
FPC = 0
TNC = 0
TPD = 0
FND = 0
FPD = 0
TND = 0
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

    for k in tempneigh1:
        if state[k] == 0:
            sc = np.random.normal(-0.5,np.sqrt(0.5),1)
            if tempstatei == 0:
                if sc[0]<tauc and k!=i:
                    tempneigh.append(k)
                    TPC += 1
                else:
                    FNC += 1
            else:
                if sc[0]<taud and k!=i:
                    tempneigh.append(k)
                    TPD += 1
                else:
                    FND += 1
        else:
            sd = np.random.normal(0.5,np.sqrt(0.5),1)
            if tempstatei == 0:
                if sd[0]<tauc and k!=i:
                    tempneigh.append(k)
                    FPC += 1
                else:
                    TNC += 1
            else:
                if sd[0]<taud and k!=i:
                    tempneigh.append(k)
                    FPD += 1
                else:
                    TND += 1

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
        TPC = 0
        FNC = 0
        FPC = 0
        TNC = 0
        TPD = 0
        FND = 0
        FPD = 0
        TND = 0

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
#    if t%10000 == 0:
#        print('====================',t,avecoop,tauc,taud)

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

print(sys.argv[1],sys.argv[2],sys.argv[3],coop,degree,prosp,pay,transitionNum,TPC+TPD,FNC+FND,FPC+FPD,TNC+TND)
filename = 'Output'+str(sys.argv[1])+'and'+str(sys.argv[2])+'and'+str(sys.argv[3])
f = open(filename,'a')
f.write(str(coop)+'\t'+str(degree)+'\t'+str(prosp)+'\t'+str(pay)+'\t'+str(transitionNum)+'\n')
f.write(str(TPC+TPD)+'\t'+str(FNC+FND)+'\t'+str(FPC+FPD)+'\t'+str(TNC+TND)+'\n')
f.close()

#f1 = 'Traj'+str(sys.argv[1])+str(sys.argv[2])+'.txt'
#np.savetxt(f1,[Avecoop,Avedegree,Aveprosp,Avepay])
