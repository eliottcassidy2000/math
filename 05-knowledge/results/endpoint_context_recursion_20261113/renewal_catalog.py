"""All half-slope return matrices through8 and a literal factorization of the first large return."""
from itertools import permutations
import sympy as sp
import numpy as np
from typed_glue import HERE,need,save,seven,glue
from context_head import ordinary_bank
from feedback_probe import all_profiles
from port_compiler import compiler,mul
from inflation_probe import relabel


def run():
    basis=[((2,0),0),((2,0),1),((6,4,0),2),((2,4,1),0)]
    V=sp.Matrix([[seven(B,r)[j] for B,r in basis] for j in (0,1,3,4)])
    den=abs(int(V.det()));W=np.array(V.inv()*den,dtype=np.int64)
    d=(0,1,2,1);a=(1,0,-2,-1);b=(0,0,1,1)
    matrices={};ports=0
    for (n,code),A in sorted(ordinary_bank().items()):
        if n<2:continue
        for r in range(n):
            h,s,e,q,F,G,P=seven(A,r)
            if G+e-2*q!=2*(G-h-q+e):continue
            profiles=[all_profiles(np.array(glue(A,r,B,t),dtype=np.int64))[:n][:,[0,1,3,4]] for B,t in basis]
            for u in range(n):
                Q=np.stack([z[u] for z in profiles],axis=1)@W
                need(np.all(Q%den==0),'integral matrix')
                M=tuple(tuple(map(int,row)) for row in Q//den);image=mul(M,d)
                if not image[1] or image!=(0,image[1],2*image[1],image[1]):continue
                A0,B0=mul(M,a),mul(M,b)
                if any(z[3]+z[1]!=z[0]+z[2] for z in (A0,B0)):continue
                ports+=1
                N=((A0[0],B0[0]),(A0[0]-A0[1]+A0[3],B0[0]-B0[1]+B0[3]))
                matrices.setdefault((M,N),dict(n=n,rows=A,entry=r,exit=u,matrix=M,target_matrix=N,multiplier=image[1],determinant=N[0][0]*N[1][1]-N[0][1]*N[1][0]))
    records=list(matrices.values());unimodular=[r for r in records if abs(r['determinant'])==1]
    print('RENEWALS',ports,'MATRICES',len(records),'UNIMODULAR',len(unimodular),flush=True)
    for r in unimodular[:20]:print('UNIT',r['n'],r['target_matrix'],'SCALE',r['multiplier'],flush=True)
    # Factor the discovered order8 context as a four-vertex outer return
    # and a five-vertex inner context, retaining its original entry vertex7.
    outer=(2,4,9,3);large=(4,5,136,179,167,199,159,3);factors=[]
    for v0,v1,v3 in permutations([x for x in range(8) if x!=6],3):
        p=(v0,v1,6,v3);removed={v0,v1,6}
        if 7 in removed:continue
        if relabel(large,p)!=outer:continue
        keep=[v for v in range(8) if v not in removed]
        if not all(large[x]>>y&1 for x in removed for y in keep if y!=v3):continue
        inner=relabel(large,keep);r=keep.index(7);u=keep.index(v3)
        M=compiler(inner,r,u);MO=compiler(outer,3,2);ML=compiler(large,7,6)
        need(sp.Matrix(MO)*sp.Matrix(M)==sp.Matrix(ML),'physical context factorization matrix')
        factors.append(dict(outer_embedding=p,inner_rows=inner,inner_entry=r,inner_exit=u,matrix=M,direction_image=mul(M,d)))
    need(factors,'literal smaller factor exists')
    for f in factors:print('FACTOR',f,flush=True)
    save('renewal_catalog',dict(status='FINITE-EXACT; selected symbolic factorization',scope='Every ordinary context through8 and all entry/exit vertices, restricted to nonzero1/2 renewal preserving F+s=h+q. Distinct full matrices, not isomorphism-class counts. Target matrices need realizable profiles for inverse use.',eligible_port_occurrences=ports,distinct_matrices=len(records),unimodular_count=len(unimodular),records=records,factorizations=factors))


if __name__=='__main__':run()
