#!/usr/bin/env python3
"""Affine-flat set compiler. Explicit tournament enumeration is audit-only."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, combinations
from math import factorial
from pathlib import Path
import argparse
import json

gates = 0


def need(ok, msg):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(msg)


class Space:
    def __init__(self, n):
        self.n = n
        self.edges = list(combinations(range(n),2))
        self.index = {e:i for i,e in enumerate(self.edges)}
        self.E = len(self.edges)
        self.rhs = 1 << self.E
        self.allbits = self.rhs-1
        self.perms = list(permutations(range(n)))
        self.ix_cache = {}

    def canonical(self, equations):
        basis = {}
        for row in equations:
            for p, b in sorted(basis.items()):
                if row >> p & 1:
                    row ^= b
            mask = row & self.allbits
            if not mask:
                if row:
                    return None
                continue
            p = (mask & -mask).bit_length()-1
            for q,b in list(basis.items()):
                if b >> p & 1:
                    basis[q] ^= row
            basis[p] = row
        return tuple(b for _,b in sorted(basis.items()))

    def meet(self, F, G):
        if F is None or G is None:
            return None
        key = (F,G) if F <= G else (G,F)
        if key not in self.ix_cache:
            self.ix_cache[key] = self.canonical(F+G)
        return self.ix_cache[key]

    def size(self, F):
        return 0 if F is None else 1 << (self.E-len(F))

    def satisfies(self, x, F):
        return F is not None and all(((x & row & self.allbits).bit_count() & 1) == (row >> self.E) for row in F)

    def mapping(self, p):
        return [(self.index[tuple(sorted((p[i],p[j])))],int(p[i]>p[j])) for i,j in self.edges]

    def act(self, F, p):
        if F is None:
            return None
        mp = self.mapping(p)
        out = []
        for row in F:
            mask, rhs = 0, row >> self.E
            for e,(dest,flip) in enumerate(mp):
                if row >> e & 1:
                    mask ^= 1 << dest
                    rhs ^= flip
            out.append(mask | (rhs << self.E))
        return self.canonical(out)

    def act_mask(self, x, p):
        return sum((((x >> e & 1)^flip) << dest) for e,(dest,flip) in enumerate(self.mapping(p)))

    def fixed(self, p):
        return self.canonical([(1 << e) ^ (1 << dest) ^ (flip << self.E)
                               for e,(dest,flip) in enumerate(self.mapping(p))])

    def path(self, p):
        eq, pathmask = [], 0
        for i,j in zip(p,p[1:]):
            e = self.index[tuple(sorted((i,j)))]
            eq.append((1 << e) | (int(i<j) << self.E))
            pathmask |= 1 << e
        return self.canonical(eq), self.allbits ^ pathmask

    def translate(self, F, mask):
        if F is None:
            return None
        # Variable matrix is unchanged, so reduced form stays reduced.
        return tuple(row ^ ((((row & mask).bit_count()) & 1) << self.E) for row in F)

    def maximal(self, flats):
        out = []
        for F in sorted(set(F for F in flats if F is not None),key=lambda F:(len(F),F)):
            if not any(self.meet(F,G) == F for G in out):
                out.append(F)
        return out

    def saturate(self, flats, group=None):
        return self.maximal(self.act(F,p) for F in flats for p in (self.perms if group is None else group))

    def step(self, flats):
        return self.maximal(self.translate(self.meet(F,L),mask)
                            for F in flats for L,mask in map(self.path,self.perms))

    def compact(self, flats):
        """Exact parallel-coset merges; no claim of a minimum flat cover."""
        flats=self.maximal(flats)
        if self.measure(self.indicator(flats))==1<<self.E:
            return [()]
        while True:
            by_direction={}
            for F in flats:
                by_direction.setdefault(tuple(r & self.allbits for r in F),[]).append(F)
            out=[];changed=False
            for group in by_direction.values():
                group.sort()
                for j in range(0,len(group)-1,2):
                    F,G=group[j:j+2]
                    different=[(f^g)>>self.E for f,g in zip(F,G)]
                    pivot=different.index(1)
                    H=self.canonical(f ^ (F[pivot] if different[i] else 0)
                                     for i,f in enumerate(F) if i!=pivot)
                    need(self.size(H)==self.size(F)+self.size(G),'parallel coset hull dimension')
                    need(self.meet(F,H)==F and self.meet(G,H)==G,'parallel coset hull inclusion')
                    out.append(H);changed=True
                if len(group)%2:out.append(group[-1])
            if not changed:return flats
            flats=self.maximal(out)

    def indicator(self, flats):
        coeff = {}
        for F in self.maximal(flats):
            delta = Counter({F:1})
            for G,c in coeff.items():
                K = self.meet(F,G)
                if K is not None:
                    delta[K] -= c
            for K,c in delta.items():
                coeff[K] = coeff.get(K,0)+c
                if not coeff[K]:
                    del coeff[K]
        return coeff

    def multiply(self, a, b):
        out = Counter()
        for F,c in a.items():
            for G,d in b.items():
                H = self.meet(F,G)
                if H is not None:
                    out[H] += c*d
        return {F:c for F,c in out.items() if c}

    def measure(self, coeff, fixed=()):
        return sum(c*self.size(self.meet(F,fixed)) for F,c in coeff.items())

    def cycle_types(self):
        groups = {}
        for p in self.perms:
            todo, cycles = set(range(self.n)), []
            while todo:
                start = min(todo); x=start; length=0
                while x in todo:
                    todo.remove(x);length+=1;x=p[x]
                cycles.append(length)
            key=tuple(sorted(cycles))
            if key not in groups:groups[key]=[p,0]
            groups[key][1]+=1
        return groups

    def orbit_count(self, coeff):
        # Input must be invariant under the full symmetric group.
        total = sum(mult*self.measure(coeff,self.fixed(p))
                    for p,mult in self.cycle_types().values())
        need(total % factorial(self.n) == 0,'Burnside integrality')
        return total//factorial(self.n)

    def module_flat(self, S):
        anchor = min(S); equations = []
        for v in range(self.n):
            if v in S:continue
            a = self.index[tuple(sorted((v,anchor)))]
            for u in S:
                if u == anchor:continue
                b = self.index[tuple(sorted((v,u)))]
                equations.append((1<<a)^(1<<b)^((int(v>anchor)^int(v>u))<<self.E))
        return self.canonical(equations)

    # The following enumerative methods are separate audit oracles.
    def members(self, flats):
        return {x for x in range(1<<self.E) if any(self.satisfies(x,F) for F in flats)}

    def direct_prime(self,x):
        def arrow(i,j):
            bit=x >> self.index[tuple(sorted((i,j)))] & 1
            return bool(bit) if i<j else not bool(bit)
        for s in range(2,self.n):
            for S in combinations(range(self.n),s):
                if all(len({arrow(v,u) for u in S})==1 for v in range(self.n) if v not in S):
                    return False
        return True

    def direct_classes(self, members):
        todo=set(members); count=0
        while todo:
            x=min(todo)
            todo.difference_update(self.act_mask(x,p) for p in self.perms)
            count+=1
        return count


def run():
    rows=[]
    for n in (3,4,5):
        S=Space(n)
        composite=S.indicator(S.module_flat(T) for k in range(2,n) for T in combinations(range(n),k))
        need(S.orbit_count({():1}) == {3:2,4:4,5:12}[n],'full-space A000568 control')
        allprime={():1}
        for F,c in composite.items():allprime[F]=allprime.get(F,0)-c
        prime_masks={x for x in range(1<<S.E) if S.direct_prime(x)}
        need(S.measure(allprime)==len(prime_masks),'module union exact cardinality')
        need(S.orbit_count(allprime)==S.direct_classes(prime_masks),'prime rank-count control')
        if n != 5:
            rows.append(dict(N=n,all_classes=S.orbit_count({():1}),prime_classes=S.orbit_count(allprime),composite_intersections=len(composite)))
            continue
        g=tuple((i+1)%n for i in range(n))
        bank=S.fixed(g)
        family=S.saturate([bank])
        for depth in range(4):
            before_compact=len(family)
            prior_members=S.members(family)
            family=S.compact(family)
            need(S.members(family)==prior_members,'coset merging preserves the whole union')
            coeff=S.indicator(family)
            members=S.members(family)
            need(S.measure(coeff)==len(members),'family union cardinality')
            for x in range(1<<S.E):
                need(sum(c for F,c in coeff.items() if S.satisfies(x,F))==int(x in members),'indicator pointwise 0/1')
            count=S.orbit_count(coeff)
            need(count==S.direct_classes(members),'ordinary unmarked quotient')
            pc=S.multiply(coeff,allprime)
            pset={x for x in members if S.direct_prime(x)}
            need(S.measure(pc)==len(pset),'restricted prime cardinality')
            need(S.orbit_count(pc)==S.direct_classes(pset),'restricted prime quotient')
            rows.append(dict(N=n,depth=depth,flats_before_compact=before_compact,flats=len(family),intersection_terms=len(coeff),
                             labelled=len(members),ordinary=count,prime_ordinary=S.orbit_count(pc),
                             prime_labelled=len(pset)))
            if depth<3:
                nxt=S.step(family)
                literal=set()
                for x in members:
                    for p in S.perms:
                        L,mask=S.path(p)
                        if S.satisfies(x,L):literal.add(x^mask)
                need(S.members(nxt)==literal,'flat native image versus literal legal moves')
                family=nxt
        # A base normalizer misses target orbits already in the full cube.
        normalizer=[p for p in S.perms if S.act(bank,p)==bank]
        full={():1}
        subgroup=sum(S.measure(full,S.fixed(p)) for p in normalizer)//len(normalizer)
        need(subgroup != S.orbit_count(full),'normalizer quotient hostile')
        rows.append(dict(N=5,normalizer_order=len(normalizer),normalizer_orbits_full_cube=subgroup,
                         ordinary_orbits_full_cube=S.orbit_count(full)))
        # Invariant full family with a named vertex: quotient only by its stabilizer.
        port_group=[p for p in S.perms if p[0]==0]
        pointed=sum(S.measure(full,S.fixed(p)) for p in port_group)//len(port_group)
        need(pointed==48,'retained named-port control')
        rows.append(dict(N=5,one_named_port_classes=pointed))
    return dict(status='FINITE-EXACT, no timing or new scalar recurrence claim',rows=rows,gates=gates)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path(__file__).with_name('flat_receipt.json'))
    args=ap.parse_args();result=run();args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
