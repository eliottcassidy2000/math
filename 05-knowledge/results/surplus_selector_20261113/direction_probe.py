"""Neutral directions of exact two-observation fibres under rooted contexts."""
from fractions import Fraction
from itertools import combinations
import numpy as np
from context_transport import decode
from repair23 import profile,hp2
from inflation_probe import end_table,relabel
from fibre_probe import save,need

def run():
    count=0;singular=0;ratios={};examples={};zero_matrix=[]
    for n in range(2,7):
        for code in range(1<<(n*(n-1)//2)):
            A=decode(n,code);full=(1<<n)-1;rev=tuple(full^(1<<v)^A[v] for v in range(n))
            E=end_table(np.array(A,dtype=np.int64))
            for r in (range(n) if n<=5 else [n-1]):
                h,s,q,F=profile(A,r);e=int(E[-1,r]);G=profile(rev,r)[3]
                P=hp2(relabel(A,[v for v in range(n) if v!=r]))[1]
                delta=G-h-q+e;alpha=G+e-2*q;beta=F-P-s
                determinant=alpha*beta+delta*P
                count+=1
                if determinant==0:
                    singular+=1;ratio=Fraction(-beta,P);key=str(ratio)
                    need(alpha*ratio==delta,'nonzero neutral direction')
                    ratios[key]=ratios.get(key,0)+1
                    if key not in examples:examples[key]=dict(rows=A,root=r,profile=(h,s,q,F),outdegree=A[r].bit_count(),response=(alpha,delta,P,beta),ratio=key)
    save('direction_probe',dict(status='PASS',scope='All labelled rooted contexts of orders2..5, all labelled contexts of order6 rooted at5 (exhaustive up to rooted isomorphism); singularity is exact, no all-order classification inferred',contexts=count,singular_contexts=singular,ratios=ratios,examples=examples))
    print('CONTEXTS',count,'SINGULAR',singular,'RATIOS',ratios,'EXAMPLES',examples)

if __name__=='__main__':run()
