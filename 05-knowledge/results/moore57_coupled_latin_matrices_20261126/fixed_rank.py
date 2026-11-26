from pathlib import Path
import json, sys
import numpy as np

def rank_mod(M,p=1009):
    A=np.asarray(M,dtype=np.int64).copy()%p
    m,n=A.shape;r=0
    for j in range(n):
        inds=np.flatnonzero(A[r:,j])
        if not len(inds): continue
        i=r+int(inds[0])
        A[[r,i]]=A[[i,r]]
        a=pow(int(A[r,j]),-1,p)
        A[r,j:] = A[r,j:]*a%p
        A[r+1:,j:] = (A[r+1:,j:]-A[r+1:,j:j+1]*A[r,j:])%p
        r+=1
        if r==m: break
    return r

def main(filename):
    data=iter(map(int,Path(filename).read_text().split()));m,e=next(data),next(data)
    assert m==56
    n=m*m
    D=[set() for _ in range(n)]
    for _ in range(e):
        x,y=next(data),next(data)
        D[x].add(y);D[y].add(x)
    F=[x for x,N in enumerate(D) if len(N)==m-1]
    rows=np.repeat(np.eye(m,dtype=np.int64),m,axis=1)
    cols=np.tile(np.eye(m,dtype=np.int64),(1,m))
    one=np.ones(n,dtype=np.int64)
    E=[];A=[];Gp=[];Gm=[]
    for x in F:
        u=np.zeros(n,dtype=np.int64);u[x]=1;E.append(u)
        v=np.zeros(n,dtype=np.int64);v[list(D[x])]=1;A.append(v)
        h=rows[x//m]+cols[x%m]
        w=m*m*u-m*h+one
        Xw=m*m*v+m*h-(m+1)*one
        Gp.append(Xw+8*w);Gm.append(7*w-Xw)
    complete_rows=[a for a in range(m) if all(len(D[a*m+i])==m-1 for i in range(m))]
    assert complete_rows==[0,28]
    one_result={}
    for p in [1009,1013]:
        rp,rm=rank_mod(Gp,p),rank_mod(Gm,p)
        # The sum of vectors coming from a complete row vanishes in each list,
        # giving the matching rational upper bound.
        assert rp==rm==len(F)-len(complete_rows)
        full=rank_mod(np.vstack([rows,cols[:-1],E,A]),p)
        assert full==2*m-1+rp+rm
        one_result[p]=[rp,rm,full]
    print(json.dumps({"fixed_vertices":len(F),"known_ranks":one_result,
        "undetermined":{"positive_eigenspace":1617-one_result[1009][0],
                       "negative_eigenspace":1408-one_result[1009][1],
                       "total":n-one_result[1009][2]}},indent=2))
if __name__=="__main__":main(sys.argv[1])
