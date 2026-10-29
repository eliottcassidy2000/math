#!/usr/bin/env python3
"""Standalone exhaustive named-color tournament refresh check; exact arithmetic.
SymPy required. No production canonicalizer or Burnside code imported.
Checks deliberately survive python -O. JSON printed to stdout.
"""
import itertools as it
import json
import math
from collections import Counter
from functools import lru_cache
import sympy as s

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

@lru_cache(None)
def edges(n):
    return tuple(it.combinations(range(n), 2))

def arc(code, n, u, v):
    if u > v:
        return 1-arc(code,n,v,u)
    return (code >> edges(n).index((u,v))) & 1

def induced(state, vertices):
    colors, code = state
    n = len(colors)
    return (tuple(colors[v] for v in vertices),
            sum(arc(code,n,vertices[i],vertices[j]) << k
                for k,(i,j) in enumerate(edges(len(vertices)))))

def run(q, N=4):
    canon, banks, sizes, aut = {}, {}, {}, {}
    for n in range(N+1):
        canon[n] = {}
        perms = tuple(it.permutations(range(n)))
        for colors in it.product(range(q), repeat=n):
            for code in range(1 << math.comb(n,2)):
                state = (colors,code)
                if state in canon[n]:
                    continue
                orbit = {induced(state,p) for p in perms}
                representative = min(orbit)
                canon[n].update({member:representative for member in orbit})
        banks[n] = sorted(set(canon[n].values()))
        sizes[n] = Counter(canon[n].values())
        aut[n] = {t:math.factorial(n)//sizes[n][t] for t in banks[n]}
    counts = [len(banks[n]) for n in range(N+1)]
    require(counts == ({1:[1,1,1,2,4],2:[1,2,4,12,48]}[q]), 'class counts')
    D,U,E,delta = {},{},{},{}
    sourcefree=[1]
    for n in range(1,N+1):
        D[n] = s.zeros(counts[n-1],counts[n])
        delta[n] = s.zeros(counts[n-1],counts[n])
        U[n-1] = s.zeros(counts[n],counts[n-1])
        free=0
        no_distinguished_source=0
        for col,t in enumerate(banks[n]):
            has_distinguished_source=False
            for v in range(n):
                small=induced(t,[w for w in range(n) if w!=v])
                row=banks[n-1].index(canon[n-1][small])
                D[n][row,col]+=1
                if all(arc(t[1],n,v,w) for w in range(n) if w!=v):
                    delta[n][row,col]+=1
                    if t[0][v]==0:
                        has_distinguished_source=True
            no_distinguished_source+=int(not has_distinguished_source)
            if not any(delta[n][:,col]):
                free+=1
        sourcefree.append(free)
        require(no_distinguished_source==counts[n]-counts[n-1], 'distinguished-source-free count equals harmonic dimension')
        require(free==counts[n]-q*counts[n-1], 'colored source-free count')
        for col,(colors,code) in enumerate(banks[n-1]):
            for color in range(q):
                for mask in range(1 << (n-1)):
                    extension=sum((arc(code,n-1,i,j) if j<n-1 else (mask>>i)&1)<<k
                                  for k,(i,j) in enumerate(edges(n)))
                    row=banks[n].index(canon[n][(colors+(color,),extension)])
                    U[n-1][row,col]+=1
        E[n-1]=U[n-1]/(q*2**(n-1))
        require(D[n].rank()==counts[n-1], 'D rational surjective')
        require(U[n-1].rank()==counts[n-1], 'U rational injective')
        for i,t in enumerate(banks[n-1]):
            for j,u in enumerate(banks[n]):
                require(U[n-1][j,i]*aut[n][u]==D[n][i,j]*aut[n-1][t], 'pointed automorphism identity')
    for n in range(N):
        prior=E[n-1]*D[n] if n else s.zeros(1)
        require(D[n+1]*E[n]-prior==s.eye(counts[n]), 'normalized Weyl identity')
        prior=U[n-1]*D[n] if n else s.zeros(1)
        require(D[n+1]*U[n]-2*prior==q*2**n*s.eye(counts[n]), 'integer attachment identity')
        prior=U[n-1]*delta[n] if n else s.zeros(1)
        require(delta[n+1]*U[n]-prior==q*s.eye(counts[n]), 'source stripping distinct Weyl pair')
    h=[counts[j]-(counts[j-1] if j else 0) for j in range(N+1)]
    spectra={}
    for n in range(1,N+1):
        K=E[n-1]*D[n]/n
        # Build every refresh directly on a fixed labeled representative.
        for col,(colors,code) in enumerate(banks[n]):
            frequencies=Counter()
            for v in range(n):
                star=[k for k,pair in enumerate(edges(n)) if v in pair]
                for color in range(q):
                    newcolors=tuple(color if w==v else colors[w] for w in range(n))
                    for mask in range(1 << (n-1)):
                        changed=code
                        for bit,k in enumerate(star):
                            changed=(changed & ~(1<<k)) | (((mask>>bit)&1)<<k)
                        frequencies[canon[n][(newcolors,changed)]]+=1
            for row,t in enumerate(banks[n]):
                require(K[row,col]==s.Rational(frequencies[t],n*q*2**(n-1)), 'direct refresh equals ED/n')
        pi=s.Matrix([s.Rational(math.factorial(n),q**n*2**math.comb(n,2)*aut[n][t]) for t in banks[n]])
        require(sum(pi)==1 and K*pi==pi, 'stationary probability')
        require(K*s.diag(*pi)==s.diag(*pi)*K.T, 'detailed balance')
        expected={1-s.Rational(j,n):h[j] for j in range(n+1) if h[j]}
        require(K.eigenvals()==expected, 'exact spectrum and multiplicities')
        for c in range(q):
            f=s.Matrix([colors.count(c)-s.Rational(n,q) for colors,code in banks[n]])
            require(K.T*f==(1-s.Rational(1,n))*f, 'centered color count eigenfunction')
        require(n*s.trace(K)==sum(counts[:n]), 'cumulative trace identity')
        spectra[n]={str(k):v for k,v in expected.items()}
    return {'q':q,'max_level':N,'counts':counts,'harmonic_dimensions':h,
            'source_free_counts':sourcefree,'spectra':spectra,
            'commutator_max_input_degree':N-1,'status':'PASS'}

if __name__=='__main__':
    print(json.dumps({'checks':[run(1),run(2)], 'arithmetic':'exact SymPy rational',
                      'sympy_version':s.__version__},indent=2,sort_keys=True))
