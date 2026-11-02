#!/usr/bin/env python3
"""Finite hostile tests of an OPEN even-order endpoint obligation."""
from itertools import combinations
from math import comb
import json,subprocess
from audit import decode,counts,need

def value(m,p):
    return 463*m[1]-256*m[2]+64*m[3]-16*m[4]+4*m[5]-2*m[6]+64*p

def main():
    rows=[];total=0
    for k in range(3):
        low=None;witness=[];count=0
        for s in subprocess.check_output(['geng','-q','8'],text=True).splitlines():
            n,e=decode(s);m,w,p=counts(n,e)
            mm=[sum((m[i-j] if 0<=i-j<len(m) else 0)*comb(k,j) for j in range(k+1)) for i in range(7)]
            f=value(mm,p+k*w);count+=1
            if low is None or f<low:low=f;witness=[s]
            elif f==low:witness.append(s)
        need(count==12346,'eight-core universe')
        rows.append({'fair_pair_components':k,'isolates_added':4-2*k,'cases':count,'minimum':low,'witnesses':witness})
        total+=count
    # Every disjoint-clique support of order12, via all integer partitions.
    def partitions(n,least=1):
        if n==0:yield ()
        for a in range(least,n+1):
            for p in partitions(n-a,a):yield (a,)+p
    clique_rows=[]
    for parts in partitions(12):
        m=[1]+[0]*6;M=W=P=0
        for size in parts:
            cm=[1]+[0]*6
            # k-matchings of K_size: choose2k vertices, then pair them.
            for k in range(1,min(6,size//2)+1):
                pairings=1
                for j in range(1,k+1):pairings*=2*j-1
                cm[k]=comb(size,2*k)*pairings
            cw=size*comb(size-1,2) if size>=3 else 0
            cp=cw*comb(size-3,2) if size>=5 else 0
            P+=cp+W*cm[1]+cw*M;M+=cm[1];W+=cw
            m=[sum(m[j]*cm[k-j] for j in range(k+1)) for k in range(7)]
        clique_rows.append({'parts':parts,'J12':value(m,P)})
    need(len(clique_rows)==77,'partition universe')
    return {'status':'finite probe complete; universal inequality OPEN',
            'scope':'All eight-vertex cores plus k=0,1,2 separate edges and enough isolated vertices to make12; all77 disjoint-clique supports of order12.',
            'core_cases':total,'core_rows':rows,'clique_minimum':min(x['J12'] for x in clique_rows),
            'clique_zero_witnesses':[x['parts'] for x in clique_rows if x['J12']==0],
            'clique_negative_witnesses':[x for x in clique_rows if x['J12']<0]}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
