#!/usr/bin/env python3
"""Exact finite checks supporting separately proved infinite statements."""
import json
from math import isqrt

def require(x, message):
    if not x: raise RuntimeError(message)

def prime(n):
    return n >= 2 and all(n % p for p in range(2,isqrt(n)+1))

def member(n):
    return n > 0 and n % 2 == 1 and n not in (7,21) and not (n >= 151 and prime(n))

def witness(N):
    if N < 30:
        return next(a for a in range(2*N+1,3*N+1,2) if member(a) and member(a+2*N))
    return next(a for a in range(2*N+1,3*N+1,2) if a%3==0 and (a+2*N)%5==0)

def atom(n):
    return (prime(n) and n!=7) or n in (49,63,343) or (n%7==0 and prime(n//7) and n//7 not in (3,7))

small=[]
for N in range(1,10001):
    a=witness(N); b=a+2*N
    require(2*N+1<=a<b<=5*N and member(a) and member(b), f'gap {N}')
    if N<30: small.append([N,a,b])
S=[n for n in range(1,1000,2) if member(n)]
for i,a in enumerate(S):
    for b in S[i:]:
        require(member(a*b),f'closure {a},{b}')
# Independent direct irreducibility check in the full conjectured monoid.
def target(n):return n>0 and n%2==1 and n not in (7,21)
for n in range(3,20001,2):
    if not target(n):continue
    irreducible=not any(n%d==0 and target(d) and target(n//d) for d in range(3,isqrt(n)+1,2))
    require(atom(n)==irreducible, f'atom {n}')
# A CRT-free explicit arbitrary long odd interval: factorial shift.
# For 0<=i<L, (2L+3)!+2i+3 is odd composite, larger than 21.
from math import factorial
for L in range(1,101):
    A=factorial(2*L+3)+3
    for i in range(L):
        x=A+2*i
        require(x%(2*i+3)==0 and x>2*i+3 and x>21, 'interval')
# Uniform exact CRT placement across every residue of both parameters.
for N in range(1,101):
    for m in range(51,1001,2):
        candidates=[a for a in range(m-28,m+1,2) if a%3==0 and (a+2*N)%5==0]
        require(len(candidates)==1 and member(candidates[0]) and member(candidates[0]+2*N),'uniform anchor')
seeds={1,13,15,41,43}
max_depth=0
for m in range(1,100000,2):
    if not target(m): continue
    n=m; depth=0
    while n not in seeds:
        n=(n+1)//2 if n%4==1 else (n-1)//2
        require(target(n), 'inverse parent')
        depth+=1
    max_depth=max(max_depth,depth)
print(json.dumps({'status':'PASS','compact_gaps_checked':10000,'multiplicative_pairs_checked':len(S)*(len(S)+1)//2,'atom_classification_through':19999,'odd_intervals_checked_lengths':100,'small_gap_witnesses':small,'uniform_anchor_pairs_checked':47500,'inverse_tree_targets_below':100000,'inverse_tree_max_depth':max_depth,'scope':'Finite checks; infinite results proved in PROOF.md.'},indent=2))
