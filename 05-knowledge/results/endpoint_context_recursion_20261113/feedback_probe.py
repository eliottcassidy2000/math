"""One new vertex, one exceptional cross-edge, and every possible retained root."""
import gzip,json
from fractions import Fraction
import numpy as np
from numba import njit
from typed_glue import HERE,need,save,seven,diff,endpoint_type
from repair23 import K9,K11
from inflation_probe import end_table


@njit
def all_profiles(T):
    n=len(T);full=(1<<n)-1;D=end_table(T);H=D.sum(axis=1);H[0]=1
    R=np.zeros(n,dtype=np.int64)
    for v in range(n):R[v]=full^(1<<v)^T[v]
    E=end_table(R);out=np.zeros((n,7),dtype=np.int64)
    for r in range(n):
        out[r,0]=H[full];out[r,1]=E[full,r];out[r,2]=D[full,r];out[r,3]=H[full^(1<<r)]
        for S in range(1,1<<n):
            if S>>r&1:
                out[r,4]+=D[S,r]*H[full^S];out[r,5]+=E[S,r]*H[full^S]
        A=full^(1<<r);S=A
        while True:
            out[r,6]+=H[S]*H[A^S]
            if S==0:break
            S=(S-1)&A
    return out


def extend(T,mask):
    n=len(T)
    return tuple(T[i]|((1-((mask>>i)&1))<<n) for i in range(n))+(mask,)


def run():
    raw=json.loads(gzip.decompress((HERE/'typed_glue.json.gz').read_bytes()))
    pairs=[p for p in raw['pairs'] if p['label'] in ('2/3>1/3','5/6>2/3>1/3','5/7>4/9','7/10>6/11')]
    pairs.insert(0,dict(label='K9>K11',rows=[K9,K11],root=5))
    output=[];controls=0
    for pair in pairs:
        X,Y=pair['rows'];n=len(X);r=pair['root'];full=(1<<n)-1
        masks=set()
        if n<=7:masks.update(range(1<<n))
        else:
            for base in (0,full):
                for rootbit in (0,1):
                    m=(base&~(1<<r))|(rootbit<<r)
                    masks.add(m)
                    for v in range(n):
                        if v!=r:masks.add(m^(1<<v))
        equal_h=0;hits={}
        for mask in sorted(masks):
            A=extend(X,mask);B=extend(Y,mask)
            ZA=all_profiles(np.array(A,dtype=np.int64));ZB=all_profiles(np.array(B,dtype=np.int64))
            if n<=7 and mask in (0,full,full//3):
                for t in range(n+1):need(tuple(map(int,ZA[t]))==seven(A,t),'independent seven extraction');controls+=1
            if ZA[0,0]!=ZB[0,0]:continue
            equal_h+=1
            for t in range(n+1):
                a,b=tuple(map(int,ZA[t])),tuple(map(int,ZB[t]));d=diff(a,b)
                for mode in endpoint_type(d):
                    rho=str(Fraction(d[1 if mode=='E' else 2],d[3]));key=mode,rho,'retained' if t==r else 'changed'
                    rec=dict(mode=mode,slope=rho,root=t,old_root=r,mask=mask,rows=[A,B],profiles=[a,b],difference=d)
                    if key not in hits:hits[key]=rec
        print('FEEDBACK',pair['label'],'MASKS',len(masks),'EQUAL_H',equal_h,'RULES',[(k,v['mask'],v['root'],v['profiles'][0][0]) for k,v in hits.items()],flush=True)
        output.append(dict(label=pair['label'],source_order=n,source_root=r,source_rows=pair['rows'],masks_tested=sorted(masks),equal_H_extensions=equal_h,first_rules=list(hits.values())))
    save('feedback_probe',dict(status='FINITE-EXACT',scope='All64 common extension masks for the aligned six-core pair; for four inherited stranded pairs, all common masks constant off the old root except at most one exception. Every possible new root tested. Not all alignments or arbitrary extensions.',independent_profile_controls=controls,pairs=output))


if __name__=='__main__':run()
