"""Exact one-sided rooted amalgamation law and neutral-repair context classification."""
from itertools import combinations
import numpy as np
from repair23 import K9,K11,profile,hp2,literal_profile
from inflation_probe import end_table,relabel
from fibre_probe import save,need

def decode(n,code):
    T=[0]*n
    for k,(a,b) in enumerate(combinations(range(n),2)):
        if code>>k&1:T[a]|=1<<b
        else:T[b]|=1<<a
    return tuple(T)

def glue(A,r,B,s):
    n=len(A);others=[i for i in range(len(B)) if i!=s];mapping={s:r,**{v:n+i for i,v in enumerate(others)}}
    T=list(A)+[0]*len(others)
    for a in range(len(B)):
        for b in range(len(B)):
            if B[a]>>b&1:T[mapping[a]]|=1<<mapping[b]
    for a in range(n):
        if a!=r:
            for v in others:T[a]|=1<<mapping[v]
    return tuple(T)

def run():
    checks=0;hist={};hostile=None
    for n in range(2,5):
        for code in range(1<<(n*(n-1)//2)):
            A=decode(n,code);full=(1<<n)-1;reverse=tuple(full^(1<<v)^A[v] for v in range(n))
            ends=end_table(np.array(A,dtype=np.int64))
            for r in range(n):
                h,s,q,F=profile(A,r);e=int(ends[-1,r]);G=profile(reverse,r)[3]
                X=relabel(A,[v for v in range(n) if v!=r]);P=hp2(X)[1]
                delta=G-h-q+e;d=A[r].bit_count();need(delta>=d*(d-1),'dual weighted-cut floor')
                observed=[]
                for B in (K9,K11):
                    H,S,Q,BF=profile(B,5);T=glue(A,r,B,5);z=profile(T,r)
                    expected=((G-q)*(BF-Q)+(h-e)*Q+q*(H-S)+e*S,s*Q,q*Q,P*(BF-Q)+F*Q)
                    need(z==expected,'literal amalgamation response')
                    if n==2:need(literal_profile(T,r)==z,'permutation gluing control')
                    observed.append(z)
                diff=tuple(b-a for a,b in zip(*observed));need(diff==(-2*delta,2*s,2*q,2*(F-P)),'pure deletion transport')
                stable=diff[0]==0 and diff[0]+diff[3]-diff[1]==0
                need(stable==(d==0),'only sink contexts preserve both targets')
                if d==1:
                    u=(A[r]&-A[r]).bit_length()-1;vertices=[v for v in range(n) if v!=r];ur=vertices.index(u);xz=profile(X,ur)
                    need(F==P-xz[3] and s==xz[1],'exact one-out-neighbour response')
                hist[str(d)]=hist.get(str(d),0)+1;checks+=1
                if n==2 and d==1 and hostile is None:hostile=dict(A=A,root=r,graphs=[glue(A,r,B,5) for B in (K9,K11)],profiles=observed,next_targets=[z[0]+z[3]-z[1] for z in observed])
    save('context_transport',dict(status='PASS',universe='Every labelled rooted A of orders2..4, both six-core repairs; independent permutation checks at order2 of A',contexts=checks,outdegree_histogram=hist,hostile=hostile))
    print('PASS',checks,'contexts; HOSTILE',hostile)

if __name__=='__main__':run()
