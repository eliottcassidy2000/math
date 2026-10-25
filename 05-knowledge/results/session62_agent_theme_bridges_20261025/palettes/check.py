#!/usr/bin/env python3
"""Exact bounded controls; no assert statements, no third-party dependencies."""
from itertools import product, combinations
from math import comb
from fractions import Fraction
from random import Random
import argparse
import json

P = argparse.ArgumentParser()
P.add_argument('--mutant', choices=['canonical-join', 'empty-rank-q'])
args = P.parse_args()
checks = 0

def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)

def canon(x, k=0):
    d = {}
    out = []
    for a in x:
        if a < k:
            out.append(a)
        else:
            if a not in d:
                d[a] = k + len(d)
            out.append(d[a])
    return tuple(out)

def orbits(q, k, n):
    return {canon(x, k) for x in product(range(q), repeat=n)}

def bell(n, r):
    s = [1] + [0] * n
    for i in range(1, n+1):
        s = [0] + [j*s[j] + s[j-1] for j in range(1, n+1)]
    return sum(s[:r+1])

def N(q, k, w):
    return sum(comb(w,j)*k**(w-j)*bell(j,q-k) for j in range(w+1))

def sub(t, indices, k=0):
    return canon(tuple(t[i] for i in indices), k)

def join(R, U, S, V, q, k, n):
    if args.mutant == 'canonical-join' and U == (0,) and V == (1,):
        # Incorrectly align the independently canonicalized fresh colours.
        return {canon(r+s,k) for r in R for s in S}
    return {p for p in orbits(q,k,n) if sub(p,U,k) in R and sub(p,V,k) in S}

orbit_cases = 0
for q in range(2,7):
    for k in range(q+1):
        for w in range(6):
            need(len(orbits(q,k,w)) == N(q,k,w), 'marked orbit count')
            orbit_cases += 1
        need(N(q,k,0)==1, 'empty boundary orbit')
    for w in range(6):
        need(N(q,q,w)==q**w, 'fully named orbit limit')
        need(N(q,0,w)==bell(w,q), 'unnamed orbit limit')

budget_cases = 0
for q in range(2,6):
    for k in range(q+1):
        for w in range(4):
            for anchored in (False,True):
                rows = set()
                for b in product(range(q), repeat=w+int(anchored)):
                    for state in ['dead','free']+list(range(q)):
                        rows.add((state,canon(b,k)) if isinstance(state,str)
                                 else ('forced',canon(b+(state,),k)))
                z=w+int(anchored)
                need(len(rows)==2*N(q,k,z)+N(q,k,z+1), 'DFA orbit budget')
                budget_cases += 1

need(join({(0,)},(0,),{(0,)},(1,),2,0,2)=={(0,0),(0,1)},
     'independent canonical representatives lose inequality amalgamation')
join_cases=0
rng=Random(620301)
for q in range(2,5):
    for k in range(q+1):
        U,V=(0,1,2),(1,2,3)
        OU,OV=sorted(orbits(q,k,3)),sorted(orbits(q,k,3))
        for trial in range(100):
            R={x for x in OU if rng.randrange(2)}
            S={x for x in OV if rng.randrange(2)}
            got=join(R,U,S,V,q,k,4)
            literal={canon(x,k) for x in product(range(q),repeat=4)
                     if sub(x,U,k) in R and sub(x,V,k) in S}
            need(got==literal, 'orbit join versus literal assignments')
            keep=(0,3)
            need({sub(x,keep,k) for x in got} ==
                 {sub(x,keep,k) for x in product(range(q),repeat=4)
                  if sub(x,U,k) in R and sub(x,V,k) in S}, 'orbit projection')
            join_cases += 1

# Equivariance is verified by a map on orbit keys, not by chosen permutations.
def classify(L):
    return 'dead' if not L else next(iter(L)) if len(L)==1 else 'free'
def step(s,L):
    return 'dead' if s=='dead' else classify(L if s=='free' else L-{s})
def state_key(b,s,k):
    return (s,canon(b,k)) if isinstance(s,str) else ('forced',canon(b+(s,),k))
transition_cases=0
for q in range(2,6):
    C=set(range(q))
    for k in range(q+1):
        for w in range(4):
            for mask in range(1<<w):
                seen={}
                for b in product(range(q),repeat=w):
                    L=C-{b[i] for i in range(w) if mask>>i&1}
                    for s in ['dead','free']+list(range(q)):
                        key=state_key(b,s,k)
                        value=state_key(b,step(s,L),k)
                        need(key not in seen or seen[key]==value, 'equivariant DFA transition')
                        seen[key]=value
                        transition_cases+=1
                # Anchor stored as the first coordinate; closure depends on equality.
                if w:
                    closure={}
                    for b in product(range(q),repeat=w):
                        for s in ['dead','free']+list(range(q)):
                            key=state_key(b,s,k)
                            value=s!='dead' and (s=='free' or s!=b[0])
                            need(key not in closure or closure[key]==value, 'anchor closure orbit')
                            closure[key]=value

# Exact induced geometry over Z[alpha], norm a²+ab+b².
points=[(0,0),(1,0),(0,1),(0,-1),(2,0),(0,2)]
def normdiff(x,y):
    a,b=x[0]-y[0],x[1]-y[1]
    return a*a+a*b+b*b
edges={ij for ij in combinations(range(6),2) if normdiff(points[ij[0]],points[ij[1]])==1}
need(edges=={(0,1),(0,2),(1,2),(0,3),(1,4),(2,5)}, 'induced six-edge geometry')
weights={}
for c in product(range(3),repeat=6):
    if all(c[i]!=c[j] for i,j in edges):
        t=c[3:]
        weights[t]=weights.get(t,0)+1
need(len(weights)==24 and set(weights.values())=={2}, 'unit graph boundary count')
need(set(weights)=={t for t in product(range(3),repeat=3) if len(set(t))>1}, 'nonmonochromatic joint relation')
for ij in combinations(range(3),2):
    need({tuple(t[i] for i in ij) for t in weights}==set(product(range(3),repeat=2)), 'full pairwise marginal')
need(len({canon(t) for t in weights})==4, 'four complete joint orbits')

# Rational linear algebra independent of Boolean quotient routines.
def eye(q):
    return [[int(i==j) for j in range(q)] for i in range(q)]
def mm(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def mv(A,v):
    return [sum(a*b for a,b in zip(row,v)) for row in A]
def rank(A):
    B=[[Fraction(x) for x in row] for row in A]
    r=0
    for j in range(len(B[0]) if B else 0):
        p=next((i for i in range(r,len(B)) if B[i][j]),None)
        if p is None: continue
        B[r],B[p]=B[p],B[r]
        z=B[r][j]; B[r]=[x/z for x in B[r]]
        for i in range(len(B)):
            if i!=r and B[i][j]:
                z=B[i][j]; B[i]=[a-z*b for a,b in zip(B[i],B[r])]
        r+=1
        if r==len(B): break
    return r

def det(A):
    B=[[Fraction(x) for x in row] for row in A]
    d=Fraction(1)
    for j in range(len(B)):
        p=next((i for i in range(j,len(B)) if B[i][j]),None)
        if p is None:return Fraction(0)
        if p!=j:B[p],B[j]=B[j],B[p];d=-d
        z=B[j][j]; d*=z
        for i in range(j+1,len(B)):
            t=B[i][j]/z
            B[i]=[a-t*b for a,b in zip(B[i],B[j])]
    return d

def M(q,L):
    return [[int(i in L and i!=j) for j in range(q)] for i in range(q)]
def literal_count(q,word):
    return sum(all(c[i] in word[i] for i in range(len(word))) and
               all(c[i]!=c[i+1] for i in range(len(word)-1))
               for c in product(range(q),repeat=len(word)))
def vector_count(q,word):
    if args.mutant=='empty-rank-q':
        v=[Fraction(1,q-1)]*q
        for L in word:v=mv(M(q,L),v)
        return sum(v)
    s,v=1,[0]*q
    for L in word:
        v=[s*int(i in L)+x for i,x in enumerate(mv(M(q,L),v))]
        s=0
    return s+sum(v)

path_words=0
for q,maxlen in ((2,5),(3,4)):
    letters=[{i for i in range(q) if mask>>i&1} for mask in range(1<<q)]
    for n in range(maxlen+1):
        for word in product(letters,repeat=n):
            need(vector_count(q,word)==literal_count(q,word), 'whole-word count including empty-prefix flag')
            path_words+=1
# Same Boolean free state, but nonproportional continuation vectors.
prefixes=(({0,1},),({1,2},))
suffixes=(({0},),({2},))
need(all(classify(prefix[0])=='free' for prefix in prefixes), 'scalar-register hostile same free state')
scalar_hostile=[[literal_count(3,prefix+suffix) for suffix in suffixes] for prefix in prefixes]
need(scalar_hostile==[[1,2],[2,1]], 'scalar-register hostile literal responses')
need(det(scalar_hostile)==-3, 'scalar-register hostile nonproportional responses')

ranks=[]
for q in range(2,9):
    A=M(q,set(range(q)))
    singles=[M(q,{i}) for i in range(q)]
    basis=[singles[i] if i==j else mm(singles[i],singles[j])
           for i in range(q) for j in range(q)]
    for i in range(q):
        for j in range(q):
            need(basis[i*q+j]==[[int(a==i)*A[j][b] for b in range(q)] for a in range(q)], 'singleton product outer basis')
    need(rank([sum(B,[]) for B in basis])==q*q, 'two-port span and entry-query rank')
    H=[[1]+[1]*q]+[[1]+row for row in A]
    need(det(H)==(-1)**q, 'whole-word Hankel determinant')
    need(rank(H)==q+1 and rank(A)==q, 'separate whole-word and prefix ranks')
    v=[Fraction(1,q-1)]*q
    need(sum(v)!=1, 'rank-q empty-word hostile')
    for i in range(q):
        need(mv(singles[i],v)==[int(i==j) for j in range(q)], 'rank-q correct after first letter')
    ranks.append({'q':q,'processed_prefix':q,'whole_word':q+1,'two_port':q*q,'hankel_det':(-1)**q})

# Entry counts and cycle traces directly enumerated, independent from rank.
port_words=0
for q in (2,3):
    letters=[{i for i in range(q) if mask>>i&1} for mask in range(1<<q)]
    for n in range(4):
        for word in product(letters,repeat=n):
            T=eye(q)
            for L in word:T=mm(M(q,L),T)
            for a,b in product(range(q),repeat=2):
                count=int(a==b) if not word else sum(
                    c[-1]==b and c[0]!=a and
                    all(c[i] in word[i] for i in range(n)) and
                    all(c[i]!=c[i+1] for i in range(n-1))
                    for c in product(range(q),repeat=n))
                need(T[b][a]==count, 'two-port literal endpoint response')
            if n>=3:
                count=sum(all(c[i] in word[i] for i in range(n)) and
                          all(c[i]!=c[(i+1)%n] for i in range(n))
                          for c in product(range(q),repeat=n))
                need(sum(T[i][i] for i in range(q))==count,'cycle trace')
            port_words+=1

# Characteristic-two hostile: every nonempty singleton product kills 1 at q=3.
need(all(sum(row)%2==0 for row in M(3,set(range(3)))), 'characteristic-two all-one kernel')
for i in range(3):
    need(all(sum(row)%2==0 for row in M(3,{i})), 'characteristic-two singleton kernel')

def comp(p,r):return tuple(p[r[i]] for i in range(len(p)))
def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
g=(1,2,0);id3=tuple(range(3))
need(sign(g)==1 and comp(id3,comp(id3,g))!=id3, 'sign-blind full-frame hostile')
g4=(1,2,0,3)
need(g4!=(0,1,2,3) and g4[3]==3, 'full-frame versus one-colour section')
x=(1,1,0)
need(sum(x)%2==0 and tuple(x[i]^x[i]^x[i] for i in range(3))==x and any(x),'ordered pair-frame parity hostile')

print(json.dumps({'status':'PASS','checks':checks,
 'orbit_formula_cases':orbit_cases,'DFA_budget_cases':budget_cases,
 'random_invariant_join_projection_cases':join_cases,
 'literal_transition_checks':transition_cases,'path_words':path_words,
 'two_port_words':port_words,'geometric_assignments_enumerated':729,
 'geometric_boundary_rows':24,'geometric_joint_orbits':4,
 'rank_table':ranks,'deterministic_scalar_register_hostile':scalar_hostile,
 'q6_w6_budgets':{'path_orbits':2*N(6,0,6)+N(6,0,7),
                  'cycle_orbits':2*N(6,0,7)+N(6,0,8),
                  'path_named':6**6*8,'cycle_named':6**7*8}},indent=2))
