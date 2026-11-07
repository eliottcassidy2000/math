#!/usr/bin/env python3
"""Exact dependency-free Session222 shell audit; run python3 verify.py."""
from itertools import permutations
import json
from pathlib import Path

def require(condition):
    if not condition:
        raise RuntimeError("Exact verification failed")

def q(m):
    return {(a,b) for a in range(m+1) for b in range(m+1-a)}
def shell_map(a,b,m):
    if a and b: return a-1,b-1
    if a==0: return b,m+2-b
    return a-1,m+2-a

def H(n,bits):
    # Independently enumerate permutations, no recurrence assumed.
    def edge(u,v):
        x,y=max(u,v),min(u,v)
        return (u>v) != bool(bits.get((x,y),0))
    return sum(all(edge(p[i],p[i+1]) for i in range(n-1)) for p in permutations(range(1,n+1)))

def tiles(n): return [(x,y) for x in range(1,n+1) for y in range(1,x-1)]
def bits_to_q(n,bits): return {(y-1,n-x):v for (x,y),v in bits.items()}
def q_to_bits(n,bits): return {(n-b,a+1):v for (a,b),v in bits.items()}

def main():
    result={'shell_orders':[],'mode_b_exhaustive':[],'single_arc_checks':0}
    for n in range(3,101):
        m=n-3; big=q(m+2); inner={(a+1,b+1) for a,b in q(m)}
        L=big-inner; diag=big-q(m)
        f={p:shell_map(*p,m) for p in big}
        require(len(set(f.values()))==len(big) and set(f.values())==big)
        require({f[p] for p in inner}==q(m))
        require({f[p] for p in L}==diag)
        require(len(L)==len(diag)==2*n-1)
        result['shell_orders'].append(n)
    for n in range(3,6):
        cells=tiles(n)
        for mask in range(1<<len(cells)):
            bits={p:(mask>>i)&1 for i,p in enumerate(cells)}
            extended={(x+1,y+1):v for (x,y),v in bits.items()}
            require(H(n,bits)==H(n+2,extended))
        result['mode_b_exhaustive'].append({'n':n,'tilings':1<<len(cells)})
    for n in range(3,8):
        for x,y in tiles(n):
            require(H(n,{(x,y):1})==1+2**(x-y-1))
            result['single_arc_checks']+=1
    old={(3,1):1}; leg={(4,2):1}; diagonal={(5,1):1}
    require([H(3,old),H(5,leg),H(5,diagonal)]==[3,3,9])
    result['hostile']={'old_H':3,'leg_H':3,'diagonal_H':9,'leg_bit':[4,2],'diagonal_bit':[5,1]}
    result['status']='PASS'
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
