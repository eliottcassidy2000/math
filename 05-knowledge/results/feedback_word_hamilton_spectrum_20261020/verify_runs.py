#!/usr/bin/env python3
from itertools import product,permutations,combinations
import numpy as np
from probe_affine import matrices, words


def need(test, reason):
    if not test:
        raise RuntimeError(reason)


def build_projectors():
    states=list(product(range(3),range(3),range(2)))
    ix={s:i for i,s in enumerate(states)}
    size=2*len(states)
    projectors=[[np.zeros((size,size),dtype=np.int64) for gap in range(3)] for a in range(4)]
    for block,p in enumerate(permutations(range(2))):
        off=block*len(states)
        for s in states:
            for a in range(4):
                for gap in range(3):
                    if s[gap]==0 and gap and not(a>>p[gap-1]&1):continue
                    ns=list(s);ns[gap]=1 if gap==2 or not(a>>p[gap]&1) else 2
                    projectors[a][gap][off+ix[s],off+ix[tuple(ns)]]+=1
    return projectors


def main():
    ms,u,v=matrices();ps=build_projectors();I=np.eye(len(u),dtype=np.int64)
    stir=[[0]*4 for _ in range(10)];stir[0][0]=1
    for r in range(1,10):
        for t in range(1,4):stir[r][t]=t*(stir[r-1][t]+stir[r-1][t-1])
    E=[]
    for a in range(4):
        need(np.array_equal(sum(ps[a]),ms[a]),('sum',a))
        for i in range(3):
            need(np.array_equal(ps[a][i]@ps[a][i],ps[a][i]),('idempotent',a,i))
            for j in range(3):need(np.array_equal(ps[a][i]@ps[a][j],ps[a][j]@ps[a][i]),('commuting',a,i,j))
        z=I.copy()
        for j in range(4):z=z@(ms[a]-j*I)
        need(not np.any(z),('annihilator',a))
        e=[np.zeros_like(I) for _ in range(4)]
        for mask in range(8):
            p=I.copy()
            for j in range(3):
                if mask>>j&1:p=p@ps[a][j]
            e[mask.bit_count()]+=p
        E.append(e)
        for r in range(10):
            expected=sum(stir[r][t]*e[t] for t in range(4))
            need(np.array_equal(np.linalg.matrix_power(ms[a],r),expected),('Stirling',a,r))
    print('four letters: 12 idempotents, 36 commutators, four degree-four annihilators: PASS')
    print('positive Stirling matrix expansion: all four letters at run lengths0..9: PASS')
    contexts=0;silent=0
    for p in words(range(4),2):
        row=u.copy()
        for a in p:row=row@ms[a]
        for s in words(range(4),2):
            col=v.copy()
            for a in reversed(s):col=ms[a]@col
            for a in range(4):
                c=[int(row@e@col) for e in E[a]]
                need(all(x>=0 for x in c),('positive coefficients',p,a,s))
                for r in range(1,10):
                    h=sum(c[t]*stir[r][t] for t in range(4))
                    need(h==int(row@np.linalg.matrix_power(ms[a],r)@col),('scalar expansion',p,a,r,s))
                    if r>3+(h.bit_length()-1):
                        need(c[2]==c[3]==0 and h==c[1],('compression',p,a,r,s))
                        silent+=1
                contexts+=1
    print('nonnegative scalar coefficient controls:',contexts,'; long silent-run checks:',silent)
    e=E[2];c=[int(u@x@v) for x in e]
    need(c==[1,3,1,0],('example',c))
    print('homogeneous letter2 coefficient vector:',c,'; H(r)=2^r+1 for r>=1')
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
