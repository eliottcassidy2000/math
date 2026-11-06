#!/usr/bin/env python3
from functools import cache
from fractions import Fraction
import json

def need(x,why):
    if not x: raise RuntimeError(why)
@cache
def g(u):
    if u==1:return 1
    if not u%2:return g(u//2)
    h=(u+1)//2
    return g(h)+h-h%2
@cache
def q(k):
    if not k:return 0
    h=k//2
    return q(h)+(h+h%2 if k%2 else 0)
def formula(seed,u):
    X,e,V,d=seed
    return X-V+(V-d)*u+d*g(u)+(d-e)*(1-u%2)
def direct(seed,k,length):
    X,e,V,d=seed
    for i in range(length):
        if (k>>i)&1:X,e,V=X+V-e,d,2*V-d
        else:e,V=0,2*V
    return X
seeds=[(677,27,1061,49),(759,29,927,29),(1547,101,1691,109),(821,33,1079,33)]
controls=seeds+[(1,0,2,1),(1,1,2,1)]
checks=0
for seed in controls:
    X,e,V,d=seed
    for length in range(15):
        for k in range(1<<length):
            u=k+1
            out=direct(seed,k,length)
            need(out==formula(seed,u),'word/index mismatch')
            need(u<=1+(out-X)//(V-d),'inverse bound')
            checks+=1
for u in range(1,65537):
    k=u-1
    need(g(u)>=1 and g(u)%2==1,'g positivity/parity')
    need(q(k)==k+1-k%2-g(k+1),'q/g identity')
    explicit=[1012*u+49*g(u)-384+22*(1-u%2),898*u+29*g(u)-168,1582*u+109*g(u)-144+8*(1-u%2),1046*u+33*g(u)-258]
    for seed,out in zip(seeds,explicit):need(formula(seed,u)==out and out%2==1,'expanded formula')
bound=2*sum((Fraction(1,V-d) for X,e,V,d in seeds),Fraction())
need(bound==Fraction(693609731,93988567442) and bound<Fraction(3,400),'density bound')
print(json.dumps({'status':'PASS','word_checks':checks,'maximum_word_length':14,'identity_indices':65536,'odd_relative_upper_density':str(bound),'scope':'Exact four paired O/Z orbits only; all-target and prime-track coverage remain open.'},sort_keys=True,indent=2))
