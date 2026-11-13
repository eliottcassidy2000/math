"""Complete <=8 entry/exit search for four specified repair directions."""
from fractions import Fraction
import sympy as sp
import numpy as np
from collections import Counter
from typed_glue import seven,glue,HERE,need,save
from context_head import ordinary_bank
from feedback_probe import all_profiles
from port_compiler import compiler,mul


def run():
    # Four actual profiles forming a basis. All attachments have order2 or3.
    basis=[((2,0),0),((2,0),1),((6,4,0),2),((2,4,1),0)]
    V=sp.Matrix([[seven(B,r)[j] for B,r in basis] for j in (0,1,3,4)])
    den=abs(int(V.det()));W=np.array(V.inv()*den,dtype=np.int64)
    need(den!=0,'realizable interpolation basis')
    targets={'1/3':(0,1,3,1),'4/9':(0,4,9,4),'6/11':(0,6,11,6),'1/2':(0,1,2,1)}
    hits={};tested=Counter();eligible=Counter();roots=0
    for (n,code),A in sorted(ordinary_bank().items()):
        if n<2:continue
        for r in range(n):
            h,s,e,q,F,G,P=seven(A,r);alpha=G+e-2*q;delta=G-h-q+e
            keys=[key for key,d in targets.items() if alpha*d[1]==delta*d[2]]
            if not keys:continue
            eligible[n]+=1
            profiles=[]
            for B,t in basis:
                T=glue(A,r,B,t);Z=all_profiles(np.array(T,dtype=np.int64))
                profiles.append(Z[:n][:,[0,1,3,4]])
            for u in range(n):
                Y=np.stack([z[u] for z in profiles],axis=1);Q=Y@W
                need(np.all(Q%den==0),'integer port coefficients');M=tuple(tuple(map(int,row)) for row in Q//den)
                roots+=1
                for key in keys:
                    d=mul(M,targets[key]);tested[n]+=1
                    if d[0] or d[1]!=d[3] or not d[2]:continue
                    out=str(Fraction(d[1],d[2]));k=key,out
                    rank=(n,u==r,code,r,u)
                    if k not in hits or rank<hits[k][0]:hits[k]=(rank,dict(input=key,output=out,n=n,rows=A,entry=r,exit=u,matrix=M,image=d))
        if n<=6:continue
    results=[v[1] for k,v in sorted(hits.items())]
    for row in results:print('PORT',row['input'],'>',row['output'],'N',row['n'],'ENTRY',row['entry'],'EXIT',row['exit'],'IMAGE',row['image'],flush=True)
    # Symbolic certificates for every minimum representative, independent of interpolation.
    for row in results:
        need(compiler(row['rows'],row['entry'],row['exit'])==row['matrix'],'independent abstract-word matrix')
    save('port_search',dict(status='FINITE-EXACT',scope='Every ordinary tournament through8, every entry root satisfying current-H neutrality for one of four directions, every output root in context. Matrices from four independent actual attachments; every reported first rule certified by exact abstract-word enumeration.',targets=targets,eligible_by_order=dict(eligible),tests_by_order=dict(tested),compiled_port_pairs=roots,rules=results))


if __name__=='__main__':run()
