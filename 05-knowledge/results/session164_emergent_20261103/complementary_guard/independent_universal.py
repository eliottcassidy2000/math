#!/usr/bin/env python3
"""Independent direct-permutation cover audit of universal guard lift."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent / 'boolean_gates' / 'audit'))
from extended import *

def ear(a,q):return left(left(a,q),len(a))
def direct_K(a,q):
    total=0;bad=0
    for p in permutations(range(len(a))):
        split=p.index(q)+1
        if all(a[u][v] for u,v in zip(p[:split],p[1:split])) and all(a[u][v] for u,v in zip(p[split:],p[split+1:])):
            total+=1
            if split<len(a) and not a[q][p[split]]:bad+=1
    return total,bad

def main():
    graphs=ports=canonical=0
    for n in range(1,6):
        for code in range(1<<comb(n,2)):
            a=tour(n,code);h=perm_h(a);J=[];p2=sum(map(sum,direct_covers(a)))//2;graphs+=1
            for q in range(n):
                K,j=direct_K(a,q);J.append(j);need(K-h==j and j>=0,'direct bad joins')
                f=ear(a,q);m=h+2*K
                need(H(f)==m,'first ear H');need(profile(f,n)==(m,2*K,0),'first ear profile')
                if sum(a[v][q] for v in range(n))==1:need(j==profile(a,q)[2],'guard specialization')
                # Exhaust the canonical second stage through old n=4, and sample at n=5.
                if n<=4 or code%73==0:
                    c=ear(f,n);need(profile(c,n+2)==(3*m,2*m,0),'canonical lift');canonical+=1
                    if len(sccs(a))==1:need(len(sccs(f))==len(sccs(c))==1,'strong lift')
                ports+=1
            need(sum(J)==2*p2-(n-1)*h,'all-port slack')
            need(H(wrap(a))==3*h+2*p2,'wrapper H')
            need(profile(wrap(a),n+1)[2]==0,'wrapper zero defect')
            if n>=2:need(sum(J)<n*p2 and min(J)<p2,'strict wrapper improvement')
        print('Direct-cover universal audit: all labeled order',n,'PASS',flush=True)
    def C(v):h,u,d=v;return(h+u+d,h+d,h-u)
    def D(v):h,u,d=v;return(h+u+d,2*u,0)
    for m in range(1,31):
        v=(3*m,2*m,0);need(D(C(v))==tuple(3*z for z in v),'CD canonical return')
        need(D(C(D(v)))==tuple(5*z for z in v),'DCD canonical return')
    print('PASS:',graphs,'tournaments,',ports,'arbitrary ports,',canonical,'four-vertex canonical checks; singleton and empty-Q included; canonical CD/DCD algebra verified.',flush=True)

if __name__=='__main__':main()
