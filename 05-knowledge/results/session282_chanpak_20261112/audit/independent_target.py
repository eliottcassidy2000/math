#!/usr/bin/env python3
from itertools import permutations,product
from fractions import Fraction
import json

def check(p,msg):
    if not p:raise RuntimeError(msg)
def tape(w):
    n=len(w)+1
    return tuple(sum(1<<j for j in range(n) if i!=j and ((i==0 and w[j-1])or(j==0 and not w[i-1])or(i and j and i<j)))for i in range(n))
def hp(T):
    n=len(T)
    return sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-1))for p in permutations(range(n)))
def full_cover_brute(T):
    n=len(T);s=0
    for p in permutations(range(n)):
        for cut in range(n+1):
            if all(T[p[i]]>>p[i+1]&1 for i in range(n-1)if i!=cut-1):s+=1
    return s
def full_three_assign(w):
    n=len(w);s=0
    for tags in product(range(3),repeat=n):
        before=[i for i in range(n)if tags[i]==0]
        after=[i for i in range(n)if tags[i]==1]
        if (not before or w[before[-1]]==0)and(not after or w[after[0]]==1):s+=2
    return s
def scalar(w):
    if not w:return 1
    return int(w[-1]==0)+int(w[0]==1)+sum(w[i]==0 and w[i+1]==1 for i in range(len(w)-1))+sum(2**(i-j-1)for i in range(len(w))for j in range(i)if w[j]==1 and w[i]==0)
controls=0
for n in range(6):
    for w in product((0,1),repeat=n):
        T=tape(w)
        check(hp(T)==scalar(w),'literal scalar')
        check(full_cover_brute(T)==full_three_assign(w),'full cover not root deletion')
        controls+=1
for m in range(1,6):
    w=(0,1)*m
    expected=Fraction(8,25)*9**m+(Fraction(4,15)*m+Fraction(178,225))*4**m+Fraction(8,9)
    check(full_three_assign(w)==expected,'F1 independent assignments')
w=(0,0,0,1,1,1)
expected=Fraction(128,1225)*729+(Fraction(148,315)+Fraction(31282,99225))*64+Fraction(128,81)
check(full_three_assign(w)==expected,'F3 independent assignments')
for m in range(1,21):
    check(9*scalar((0,1)*2*m)-25*scalar((0,0,1,1)*m)==2*m,'Jordan linear residual')
    check(scalar((0,1)*3*m)==9*scalar((0,0,0,1,1,1)*m),'scalar twins')
print(json.dumps(dict(status='PASS',literal_all_tapes_through_length=5,literal_graph_controls=controls,F1_assignments_through_m=5,F3_assignment_m=1,Jordan_controls=20),indent=2,sort_keys=True))
