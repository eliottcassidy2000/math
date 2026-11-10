"""Exact targeted check of NEW resampling-spectrum formulas, n=3,4,5.
No inherited theorem audits. Uses direct tournament relabellings and
rational quotient transition matrices, not the Fourier proof.
"""
import itertools as it
import math
from collections import Counter
from fractions import Fraction
import sympy as s


def check(n):
    edges=list(it.combinations(range(n),2)); m=len(edges)
    index={e:i for i,e in enumerate(edges)}
    maps=[]
    for p in it.permutations(range(n)):
        maps.append([(index[tuple(sorted((p[a],p[b])))],p[a]>p[b]) for a,b in edges])
    def relabel(x,mp):
        return sum(((((x>>i)&1)^flip)<<j) for i,(j,flip) in enumerate(mp))
    canon={x:min(relabel(x,mp) for mp in maps) for x in range(1<<m)}
    reps=sorted(set(canon.values())); rid={x:i for i,x in enumerate(reps)}
    volumes=Counter(canon.values())
    ts=[1,1,1,2,4,12]
    assert len(reps)==ts[n]
    z=s.Symbol('z')
    for k in range(n):
        matrix=[]
        for x in reps:
            counts=Counter()
            for S in it.combinations(range(n),k):
                S=set(S); keep=sum(1<<i for i,(a,b) in enumerate(edges) if a in S and b in S)
                free=[i for i in range(m) if not keep>>i&1]
                for bits in range(1<<len(free)):
                    y=(x&keep)|sum(((bits>>j)&1)<<i for j,i in enumerate(free))
                    counts[rid[canon[y]]]+=1
            total=sum(counts.values())
            matrix.append([s.Rational(counts[j],total) for j in range(len(reps))])
        K=s.Matrix(matrix)
        assert all(sum(K.row(i))==1 for i in range(len(reps)))
        assert all(volumes[reps[i]]*K[i,j]==volumes[reps[j]]*K[j,i] for i in range(len(reps)) for j in range(len(reps)))
        expected=z-1
        spectrum=Counter({s.Rational(1):1})
        for v in range(1,n+1):
            multiplicity=ts[v]-ts[v-1]
            lam=s.Rational(math.comb(k,v),math.comb(n,v)) if k>=v else s.Rational(0)
            spectrum[lam]+=multiplicity
            expected*= (z-lam)**multiplicity
        assert s.expand(K.charpoly(z).as_expr()-expected)==0
        print(f'n={n}, k={k}, types={len(reps)}, spectrum='+str({str(a):b for a,b in spectrum.items() if b}))

for n in [3,4,5]: check(n)
print('PASS: direct exact quotient kernels match every predicted eigenvalue and multiplicity; stochasticity and orbit-volume reversibility hold.')
