#!/usr/bin/env python3
"""A bounded algebraic probe of append-block affine actions; not a completeness search."""
from itertools import product,permutations,combinations
import numpy as np


def matrices():
    states=list(product(range(3),range(3),range(2)))
    ix={s:i for i,s in enumerate(states)}
    size=2*len(states)
    matrices=[np.zeros((size,size),dtype=np.int64) for _ in range(4)]
    u=np.zeros(size,dtype=np.int64);v=np.zeros(size,dtype=np.int64)
    for block,p in enumerate(permutations(range(2))):
        off=block*len(states);u[off+ix[(0,0,0)]]=1
        for s in states:
            if 2 not in s and (s[1] !=0 or p[0]==0):v[off+ix[s]]=1
            for a in range(4):
                for gap in range(3):
                    if s[gap]==0 and gap and not (a>>p[gap-1]&1):continue
                    ns=list(s);ns[gap]=1 if gap==2 or not(a>>p[gap]&1) else 2
                    matrices[a][off+ix[s],off+ix[tuple(ns)]]+=1
    return matrices,u,v


def words(alphabet,length):
    return [w for l in range(length+1) for w in product(alphabet,repeat=l)]


def main():
    ms,u,v=matrices();ws=words(range(4),3)
    rows=[];cols=[]
    for w in ws:
        r=u.copy();c=v.copy()
        for a in w:r=r@ms[a]
        for a in reversed(w):c=ms[a]@c
        rows.append(r);cols.append(c)
    R=np.array(rows);C=np.array(cols).T
    hits=[]
    for a,b in combinations(range(4),2):
        da=R@(ms[a]-2*np.eye(len(u),dtype=np.int64))@C
        db=R@(ms[b]-2*np.eye(len(u),dtype=np.int64))@C
        ok=(da!=db)
        for w in words((a,b),3):
            rw=R.copy()
            for x in w:rw=rw@ms[x]
            ok &= ((rw@(ms[a]-2*np.eye(len(u),dtype=np.int64))@C)==da)
            ok &= ((rw@(ms[b]-2*np.eye(len(u),dtype=np.int64))@C)==db)
        for i,j in zip(*np.nonzero(ok)):
            hits.append((ws[i],ws[j],a,b,int((R@C)[i,j]),int(da[i,j]),int(db[i,j])))
    print('prefix/suffix words:',len(ws),'each; alphabet pairs: 6; middle depth <=3')
    print('distinct-offset candidate count:',len(hits))
    print('first candidates:',hits[:30])

if __name__=='__main__':main()
