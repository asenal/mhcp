#!/usr/bin/python 
from __future__ import division
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
random.seed(153)

class Model:
    def __init__(self,s,n):
        self.n = n
        self.s = s
    def __repr__(self):
        return self.s

R = [(n+1)*0.1 for n in range(9)]
random.shuffle(R)
D = dict()
letters = list('ABCDEFGHI')
for s in letters:
    D[s] = Model(s,R.pop())
    print D[s].s , D[s].n
letters = sorted(letters,key=lambda c: D[c].n)

#-----------------------------------
L = list()
for i,c in enumerate(letters):
    L.extend(list(c*(i+1)))
     
def get_proposal(model):
    c = random.choice(L)
    return D[c]
#-----------------------------------
rL = ['A']
N = 100000
for i in range(N):
    current = D[rL[-1]]
    p = current.n
    next = get_proposal(current)
    q = next.n
    r = q/p
    bias = L.count(str(next))/L.count(str(current))
    r /= bias
    if r >= 1:
        rL.append(str(next))
    elif r > random.random():
        rL.append(str(next))
    else:
        rL.append(str(current))
#-----------------------------------
M = int(N/10)
C = Counter(rL[M:])
for s in letters:
    print s, D[s].n, round(C[s]/D[s].n)


X = np.array(range(len(letters)))
Y = [C[c] for c in letters]
plt.xticks(X + 0.4,list(letters))
plt.bar(X,Y)
plt.savefig('example.png')
