"""Isomorphism-level sequences for the neutral-context metagraph through order6."""
from itertools import permutations,product
from collections import Counter
from fractions import Fraction
from context_transport import decode
from repair23 import profile,hp2
from inflation_probe import relabel
from fibre_probe import save,need

def canonical(T,r=None):
    n=len(T);free=[v for v in range(n) if v!=r]
    signatures={v:(T[v].bit_count(),int(bool(T[v]>>r&1)) if r is not None else 0) for v in free}
    cells=[[v for v in free if signatures[v]==key] for key in sorted(set(signatures.values()))]
    return min(relabel(T,sum((tuple(x) for x in order),())+(() if r is None else (r,))) for order in product(*(permutations(cell) for cell in cells)))

def run():
    levels=[]
    for n in range(2,7):
        rooted={canonical(decode(n,c),n-1) for c in range(1<<(n*(n-1)//2))};ordinary={canonical(T) for T in rooted}
        need(len(ordinary)=={2:1,3:2,4:4,5:12,6:56}[n],'A000568 control')
        transitions=Counter();records=[]
        for A in sorted(rooted):
            r=n-1;h,s,q,F=profile(A,r);full=(1<<n)-1;rev=tuple(full^(1<<v)^A[v] for v in range(n));_,e,_,G=profile(rev,r)
            P=hp2(relabel(A,list(range(n-1))))[1];delta=G-h-q+e;alpha=G+e-2*q;beta=F-P-s
            if alpha*beta+delta*P:continue
            rho=str(Fraction(-beta,P));sigma=str(Fraction(s,q));transitions[rho,sigma]+=1
            records.append(dict(rows=A,root=r,profile=(h,s,q,F),response=(alpha,delta,P,beta),input_slope=rho,output_slope=sigma))
        levels.append(dict(n=n,ordinary=len(ordinary),rooted=len(rooted),neutral=len(records),slopes=len(set(r['input_slope'] for r in records)),transitions=[dict(input=a,output=b,classes=c) for (a,b),c in sorted(transitions.items())],records=records))
        print('LEVEL',n,'ORDINARY',len(ordinary),'ROOTED',len(rooted),'NEUTRAL',len(records),'SLOPES',levels[-1]['slopes'],'TRANSITIONS',dict(transitions))
    save('context_sequences',dict(status='FINITE-EXACT',scope='Every rooted tournament up to isomorphism at n2..6, exhaustive labelled orientations with fixed root plus exact root-preserving canonicalization; ordinary counts agree with known A000568 head',levels=levels))

if __name__=='__main__':run()
