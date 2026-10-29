#!/usr/bin/env python3
from itertools import product,combinations
from collections import Counter
from math import gcd
import json

def check(b,msg):
    if not b:raise RuntimeError(msg)
def jacobi(a,n):
    t=1;a%=n
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5):t=-t
        a,n=n,a
        if a%4==n%4==3:t=-t
        a%=n
    return t if n==1 else 0
s=[2]
for _ in range(8):s.append(s[-1]**2-s[-1]+1)
for a,b in combinations(s[1:],2):
    check(gcd(a,b)==1 and a%4==b%4==3,'odd coprime domain')
    check(jacobi(b,a)==1 and jacobi(a,b)==-1,'Jacobi reverse order')
# sigma tests: order (4,2,5,3,7,13); enumerate every CRT class.
args=[7,3,6,4,8,14]
counts=Counter(tuple((1-jacobi(a,q))//2 for a in args) for q in range(3,168,4) if gcd(q,168)==1)
check(len(counts)==8 and set(counts.values())=={3},'sigma CRT density')
for c in counts:
    check(c[3]==0 and c[1]^c[2]^c[4]==0 and c[0]^c[4]^c[5]==0,'sigma relations')
for q in range(3,168,4):
    if gcd(q,168)==1:check(jacobi(12,q)==jacobi(3,q),'duplicate sigma11 bit')
# Column coordinates for 3,6,2,5,30,15,10 over prime basis (2,3,5).
cols=[2,3,1,4,7,6,5]
code={tuple((z&c).bit_count()%2 for c in cols) for z in range(8)}
kernel=[w for w in product(range(2),repeat=7) if not __import__('functools').reduce(int.__xor__,(c for c,b in zip(cols,w) if b),0)]
check(Counter(map(sum,code))=={0:1,4:7},'simplex weight distribution')
check(len(kernel)==16 and min(sum(w) for w in kernel if any(w))==3,'Hamming dimension and distance')
for w in kernel:
    for c in code:check(sum(a*b for a,b in zip(w,c))%2==0,'code duality')
for i,j in combinations(range(7),2):check(set(Counter((c[i],c[j]) for c in code).values())=={2},'pairwise independent coordinates')
R=194478922284018372000
def exponent(n,b):
    e=0
    while n%b==0:n//=b;e+=1
    return e
for n,arrows in [(R,{4:2,2:5,5:3,3:7,7:4,13:2,11:1}),(R//2,{4:2,2:4,5:3,3:7,7:4,13:2,11:1})]:
    check(all(exponent(n,b)==e for b,e in arrows.items()),'composite-base exponent map')
print(json.dumps({'status':'PASS','odd_Sylvester_terms':len(s)-1,'sigma_colors':len(counts),'simplex_words':len(code),'Hamming_words':len(kernel),'composite_base_maps':'both verified'},indent=2))
