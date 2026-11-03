#!/usr/bin/env python3
"""Exact bounded audit; Python standard library only, no optimizable asserts."""
import argparse
from collections import Counter
from itertools import product, permutations
from math import isqrt
import json
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


class Field:
    """F_5[z]/(monic modulus), coefficient tuples in increasing degree."""
    def __init__(self, modulus):
        self.modulus = tuple(modulus)
        self.n = len(modulus)-1
        self.elements = tuple(product(range(5), repeat=self.n))
        self.zero = (0,)*self.n
        self.one = (1,)+(0,)*(self.n-1)
    def const(self, a):
        return (a % 5,)+(0,)*(self.n-1)
    def add(self, a, b):
        return tuple((x+y)%5 for x,y in zip(a,b))
    def mul(self, a, b):
        v = [0]*(2*self.n-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                v[i+j] += x*y
        for i in range(len(v)-1,self.n-1,-1):
            for j in range(self.n):
                v[i-self.n+j] -= v[i]*self.modulus[j]
        return tuple(x%5 for x in v[:self.n])
    def power(self, a, n):
        r = self.one
        while n:
            if n & 1:
                r = self.mul(r,a)
            a = self.mul(a,a)
            n //= 2
        return r
    def validate(self):
        q = len(self.elements)
        roots = [a for a in range(5) if sum(c*a**i for i,c in enumerate(self.modulus))%5 == 0]
        if self.n in (2,3):
            need(not roots, "reducible extension modulus")
        need(all(self.power(x,q)==x for x in self.elements), "Frobenius field identity")
        need(all(self.mul(x,self.power(x,q-2))==self.one for x in self.elements if x!=self.zero), "field inverses")
        sq = Counter(self.mul(y,y) for y in self.elements)
        need(sq[self.zero]==1 and Counter(sq.values())==Counter({2:(q-1)//2,1:1}), "quadratic residue root histogram")
        return {"modulus_low_to_high":list(self.modulus), "q":q, "prime_field_roots":roots,
                "square_root_multiplicities":{"zero":1,"nonzero_square":2,"nonzero_nonsquare":0}}
    def count(self,a,b):
        sq = Counter(self.mul(y,y) for y in self.elements)
        rhs = [self.add(self.add(self.power(x,3),self.mul(self.const(a),x)),self.const(b)) for x in self.elements]
        # One point at infinity plus exact root multiplicities, never a trace recurrence.
        hist = Counter(sq[r] for r in rhs)
        n = 1+sum(sq[r] for r in rhs)
        need(n==1+hist[1]+2*hist[2], "point count histogram")
        return n, {str(k):hist[k] for k in (0,1,2)}


def mm(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def symcube(m):
    # Substitute X -> m00 X+m10 Y, Y -> m01 X+m11 Y into degree-3 monomials.
    def pm(a,b):
        out=Counter()
        for (x,y),v in a.items():
            for (u,w),z in b.items():
                out[(x+u,y+w)]+=v*z
        return out
    out=[[0]*4 for _ in range(4)]
    for j in range(4):
        poly={(0,0):1}
        for _ in range(3-j):
            poly=pm(poly,{(1,0):m[0][0],(0,1):m[1][0]})
        for _ in range(j):
            poly=pm(poly,{(1,0):m[0][1],(0,1):m[1][1]})
        for i in range(4):
            out[i][j]=poly.get((3-i,i),0)
    return out


def line_sums(g):
    return [sum(r) for r in g]+[sum(g[i][j] for i in range(3)) for j in range(3)]+[sum(g[i][i] for i in range(3)),sum(g[i][2-i] for i in range(3))]


def exhaustive_magic_count(values):
    # 8! permutations after the independently proved forced-center constraint.
    need(len(values)==9 and len(set(values))==9,"distinct arrangement control")
    if sum(values)%9:
        return 0
    center=sum(values)//9
    if center not in values:
        return 0
    count=0
    target=3*center
    for other in permutations([v for v in values if v!=center]):
        cells=other[:4]+(center,)+other[4:]
        if cells[0]+cells[1]+cells[2]!=target or cells[3]+cells[4]+cells[5]!=target:
            continue
        g=[cells[i:i+3] for i in (0,3,6)]
        count+=line_sums(g)==[target]*8
    return count


def canonical_tournaments():
    edges=list((i,j) for i in range(5) for j in range(i+1,5))
    perms=list(permutations(range(5)))
    def matrix(mask):
        a=[[0]*5 for _ in range(5)]
        for k,(i,j) in enumerate(edges):
            a[i][j]=(mask>>k)&1
            a[j][i]=1-a[i][j]
        return a
    def canon(a):
        return min(sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(edges)) for p in perms)
    reps=sorted({canon(matrix(m)) for m in range(1024)})
    spectra=[]
    dual={}
    for m in reps:
        a=matrix(m)
        dual[m]=canon([list(x) for x in zip(*a)])
        # Faddeev-LeVerrier using matrix powers and exact Newton identities.
        powers=[]; cur=[[int(i==j) for j in range(5)] for i in range(5)]
        coeff=[1]
        for k in range(1,6):
            cur=mm(cur,a); powers.append(tr(cur))
            num=-sum(coeff[k-j]*powers[j-1] for j in range(1,k+1))
            need(num%k==0,"characteristic polynomial integrality")
            coeff.append(num//k)
        spectra.append(tuple(coeff))
    fibers=sorted(Counter(spectra).values(),reverse=True)
    fixed=sum(m==dual[m] for m in reps)
    quotient=len({min(m,dual[m]) for m in reps})
    need(len(reps)==12 and len(set(spectra))==9,"tournament census")
    need(fibers==[3,2]+[1]*7 and fixed==8 and quotient==10,"tournament fiber and involution audit")
    return {"classes":len(reps),"spectrum_values":len(set(spectra)),"spectrum_fiber_partition":fibers,
            "self_converse":fixed,"converse_orbits":quotient}


def build():
    fields=[Field([0,1]),Field([2,0,1]),Field([1,2,0,1])]
    field_results=[f.validate() for f in fields]
    models=[(a,b) for a,b in product(range(5),repeat=2) if (4*a**3+27*b*b)%5]
    need(len(models)==20,"nonsingular model count")
    counts={}; rows=[]
    for a,b in models:
        data=[f.count(a,b) for f in fields]
        ns=[x[0] for x in data]; t=6-ns[0]
        need(ns[0]==1+sum((y*y-x*x*x-a*x-b)%5==0 for x,y in product(range(5),repeat=2)),"independent F5 point pairs")
        need(ns[1]==36-t*t and ns[2]==126-t**3+15*t,"extension recurrence comparison")
        counts[(a,b)]=(t,*ns)
        rows.append({"a":a,"b":b,"t":t,"N1":ns[0],"N2":ns[1],"N3":ns[2],
                     "rhs_root_histograms":{"F5":data[0][1],"F25":data[1][1],"F125":data[2][1]}})
    need(sorted({r[1] for r in counts.values()})==list(range(2,11)),"base point count range")
    need(sorted({r[0] for r in counts.values()})==list(range(-4,5)),"trace range")
    def orbit(m):
        a,b=m
        return frozenset((a*pow(u,4,5)%5,b*pow(u,6,5)%5) for u in range(1,5))
    classes=sorted({tuple(sorted(orbit(m))) for m in models})
    reps=[cl[0] for cl in classes]
    repof={m:cl[0] for cl in classes for m in cl}
    need(len(classes)==12 and set(repof)==set(models),"complete scaling quotient")
    tracefibers=Counter(counts[r][0] for r in reps)
    need(sorted(tracefibers.values(),reverse=True)==[2,2,2]+[1]*6,"elliptic trace fibers")
    twist={r:repof[(4*r[0]%5,3*r[1]%5)] for r in reps}
    need(all(twist[twist[r]]==r and twist[r]!=r for r in reps),"free quadratic twist involution")
    need(all(counts[r][0]==-counts[twist[r]][0] for r in reps),"twist trace sign")
    pairs=sorted({tuple(sorted((r,twist[r]))) for r in reps})
    need(len(pairs)==6,"twist orbit count")
    grid=[[108,148,122],[140,126,112],[130,104,144]]
    need(line_sums(grid)==[378]*8,"all eight magic lines")
    need(sorted(sum(grid,[]))==sorted({r[3] for r in counts.values()}),"N3 grid multiset")
    need(any(isqrt(x)**2!=x for x in sum(grid,[])),"not a square-of-squares grid")
    arrangements={"N3":exhaustive_magic_count(sum(grid,[])),
                  "symmetric_cube":exhaustive_magic_count([t**3-10*t for t in range(-4,5)]),
                  "literal_cube":exhaustive_magic_count([t**3 for t in range(-4,5)])}
    need(arrangements=={"N3":8,"symmetric_cube":0,"literal_cube":0},"8!-permutation arrangement counts")
    byN3={r[3]:m for m,r in sorted(counts.items(),reverse=True)}
    tracegrid=[[counts[byN3[x]][0] for x in row] for row in grid]
    need(sorted(sum(tracegrid,[]))==list(range(-4,5)),"grid trace occurrences")
    # These N1 labels are not themselves a magic arrangement: relabeling changes line predicate.
    basegrid=[[6-t for t in row] for row in tracegrid]
    need(len(set(line_sums(basegrid)))>1,"base labeling hostile")
    operators=[]
    for t in range(-4,5):
        m=[[0,-5],[1,t]]
        cube=tr(mm(mm(m,m),m)); sym=tr(symcube(m))
        need(cube==t**3-15*t and sym==t**3-10*t,"operator traces")
        need(sym-cube==5*t,"operator boundary")
        operators.append({"t":t,"power_cube_trace":cube,"symmetric_cube_trace":sym})
    n2=[36-t*t for t in range(-4,5)]
    # Integer entries require integer center equal to arithmetic mean; sum=264 yields 88/3.
    need(sum(n2)==264 and sum(n2)%9!=0,"F25 center obstruction")
    tourn=canonical_tournaments()
    inherited=json.loads((Path(__file__).parent/'sources/session143_tournaments.json').read_text())
    need(len(inherited['classes'])==tourn['classes'],"pinned tournament class census")
    need(sorted(Counter(tuple(c['characteristic_desc']) for c in inherited['classes']).values(),reverse=True)==tourn['spectrum_fiber_partition'],"pinned tournament fibers")
    need(sum(c['self_dual'] for c in inherited['classes'])==tourn['self_converse'],"pinned tournament involution")
    return {"status":"PASS: finite-exact independent audit", "field_checks":field_results,"models":rows,
            "model_count":len(models),"isomorphism_classes":[{"representative":r,"models":cl,"trace":counts[r][0]} for r,cl in zip(reps,classes)],
            "trace_fibers":dict(sorted(tracefibers.items())),"trace_fiber_partition":sorted(tracefibers.values(),reverse=True),
            "twist_pairs":pairs,"twist_fixed_points":0,"twist_orbits":6,"tournament_control":tourn,
            "magic_arrangement_counts":arrangements,"magic_grid":grid,"magic_line_sums":line_sums(grid),"grid_models":[[byN3[x] for x in row] for row in grid],
            "trace_grid":tracegrid,"base_count_grid":basegrid,"base_count_line_sums":line_sums(basegrid),
            "operator_traces":operators,"F25_values":n2,"F25_sum":sum(n2),"F25_nonmagic_reason":"center would be 88/3, not an integer entry",
            "scope":"F5 models and F5,F25,F125 point counts only; lambda-cube referent unresolved; no magic square of squares or octonion bridge"}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--certificate',type=Path); args=parser.parse_args()
    result=build()
    # Normalize tuples and integer dictionary keys before certificate comparison.
    result=json.loads(json.dumps(result,sort_keys=True))
    if args.certificate:
        supplied=json.loads(args.certificate.read_text())
        need(supplied==result,"certificate mismatch")
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
