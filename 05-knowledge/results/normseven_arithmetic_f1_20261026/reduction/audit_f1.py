#!/usr/bin/env python3
"""Independent reduction audit: Steiner designs, surviving characters, divisibility.

Imports no producer code. All gates remain active under python -O.
"""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json


def need(ok, msg):
    if not ok:
        raise ValueError(msg)


def dot(a, b):
    return (a & b).bit_count() % 2


def steiner_designs():
    triples = [frozenset(c) for c in combinations(range(7), 3)]
    pairs = {frozenset(c) for c in combinations(range(7), 2)}
    pairsets = {t: {frozenset(p) for p in combinations(t, 2)} for t in triples}
    found = []

    def recurse(left, lines):
        if not left:
            found.append(frozenset(lines))
            return
        p = min(left, key=lambda x: tuple(sorted(x)))
        for t in triples:
            if p <= t and pairsets[t] <= left:
                recurse(left - pairsets[t], lines + [t])
    recurse(pairs, [])
    need(len(found) == len(set(found)) == 30, "Steiner enumeration")
    return found


EVEN = tuple(x for x in range(128) if x.bit_count() % 2 == 0)
MARKINGS = {
    "U2": (1, 2, 4, 5, 6, 3, 7), "D2": (1, 2, 4, 6, 3, 7, 5),
    "U4": (1, 2, 4, 5, 7, 3, 6), "D4": (1, 2, 4, 7, 3, 6, 5),
}


def code(labels):
    return tuple(sum(dot(mask, lab) << i for i, lab in enumerate(labels)) for mask in range(8))


def surviving_distribution(c, d, zeros):
    zmask = sum(1 << i for i in zeros)
    full = [v for v in EVEN if not v & zmask]
    cu = [v for v in c if not v & zmask]
    cd = [v for v in d if not v & zmask]
    out = Counter()
    for twist in EVEN:
        x = len(full) if all(not dot(v, twist) for v in full) else 0
        y1 = len(cu) if all(not dot(v, twist) for v in cu) else 0
        y2 = len(cd) if all(not dot(v, twist) for v in cd) else 0
        out[x, y1 * y2] += 1
    return out


class Field:
    def __init__(self, p, extension=False):
        self.p, self.ext = p, extension
        self.q = p * p if extension else p
        self.nu = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1) if extension else 0

    def add(self, a, b):
        p = self.p
        return (a % p + b % p) % p + p * ((a // p + b // p) % p)

    def neg(self, a):
        p = self.p
        return (-a % p) + p * ((-(a // p)) % p)

    def linear_values(self, e, x, z):
        add, neg = self.add, self.neg
        return (e, add(e,x), add(e,neg(x)), add(e,z), add(e,neg(z)), add(add(e,x),z), add(add(e,neg(x)),neg(z)))

    def strata(self):
        pts = [(1, x, z) for x in range(self.q) for z in range(self.q)]
        pts += [(0, 1, z) for z in range(self.q)] + [(0, 0, 1)]
        return Counter(frozenset(i for i, v in enumerate(self.linear_values(*pt)) if v == 0) for pt in pts)


def puncture(mask, j):
    return sum(((mask >> i) & 1) << k for k, i in enumerate(i for i in range(7) if i != j))


def audit_conductor(c, d):
    checks = 0
    for j in range(7):
        cp, dp = [puncture(v,j) for v in c], [puncture(v,j) for v in d]
        e = {x ^ y: x & y for x in cp for y in dp}
        need(len(e) == 64, "punctured codes not complementary")
        for z in range(64):
            for i in range(6):
                required = max(0, max(((e[z ^ y] >> i) & 1) - (((z & y) >> i) & 1) for y in range(64)))
                need(required == 1 - ((z >> i) & 1), "conductor coefficient valuation")
                checks += 1
    return checks


def all_subspaces(n):
    current, total = {frozenset({0})}, {frozenset({0})}
    for _ in range(n):
        nxt = set()
        for s in current:
            for x in range(1 << n):
                if x not in s:
                    nxt.add(frozenset(s | {x ^ y for y in s}))
        total |= nxt
        current = nxt
    return total


def audit_arbitrary_inertia():
    cases = 0
    for n in range(1, 6):
        split = n // 2
        h = frozenset(range(1 << split))
        k = frozenset(x << split for x in range(1 << (n - split)))
        for inertia in all_subspaces(n):
            hi = {x ^ y for x in h for y in inertia}
            ki = {x ^ y for x in k for y in inertia}
            # Character annihilators supply actual rational-point counts.
            chars = [a for a in range(1 << n) if all(not dot(a,i) for i in inertia)]
            ch = [a for a in chars if all(not dot(a,x) for x in h)]
            ck = [a for a in chars if all(not dot(a,x) for x in k)]
            dim = lambda s: len(s).bit_length() - 1
            r, a, b = dim(inertia), dim(inertia & h), dim(inertia & k)
            d = r - a - b
            out = Counter()
            for t in range(1 << n):
                x = len(chars) if all(not dot(v,t) for v in chars) else 0
                y = (len(ch) if all(not dot(v,t) for v in ch) else 0) * (len(ck) if all(not dot(v,t) for v in ck) else 0)
                out[x,y] += 1
            expected = Counter({(2**(n-r),2**(n-r-d)):2**r})
            expected[0,2**(n-r-d)] += 2**(r+d)-2**r
            expected[0,0] += 2**n-2**(r+d)
            need(out == +expected, "arbitrary inertia distribution")
            need(len(hi & ki) == 2**(r+d), "arbitrary comparison kernel")
            cases += 1
    return cases


def main():
    designs = steiner_designs()
    pairs = [(a,b) for a,b in combinations(designs,2) if a.isdisjoint(b)]
    need(len(pairs) == 120, "disjoint design pairs")
    codes = {name: code(labs) for name, labs in MARKINGS.items()}
    cases = []
    for p, ext in ((3,False),(5,False),(7,False),(3,True),(5,True)):
        f = Field(p,ext)
        strata = f.strata()
        need(max(map(len,strata)) == 3, "fourfold intersection")
        triples = {s for s in strata if len(s)==3}
        expected = {frozenset(t) for t in ((0,1,2),(0,3,4),(0,5,6))}
        if p==3:
            expected |= {frozenset((1,3,6)),frozenset((2,4,5))}
        need(triples == expected, "exceptional triples")
        free = [d for d in designs if d.isdisjoint(triples)]
        freepairs = [(a,b) for a,b in pairs if a in free and b in free]
        need((len(free),len(freepairs)) == ((12,18) if p==3 else (16,32)), "free census")
        moments = []
        for root in (2,4):
            abs_sum, sq_sum = Fraction(0), Fraction(0)
            for s, multiplicity in strata.items():
                dist = surviving_distribution(codes[f'U{root}'],codes[f'D{root}'],s)
                abs_sum += multiplicity * sum(Fraction(count,64)*abs(x-y) for (x,y),count in dist.items())
                sq_sum += multiplicity * sum(Fraction(count,64)*(x-y)**2 for (x,y),count in dist.items())
                if len(s)==3:
                    expected_dist = Counter({(8,2):8,(0,2):24,(0,0):32}) if s in {frozenset((1,3,6)),frozenset((2,4,5))} else Counter({(8,1):8,(0,1):56})
                    need(dist == expected_dist, "triple twist distribution")
            need(abs_sum == Fraction(28*f.q-11,4), "absolute total")
            need(sq_sum == 112*f.q-(215 if p==3 else 251), "square total")
            moments.append({'root':root,'absolute':str(abs_sum),'square':str(sq_sum)})
        cases.append({'q':f.q,'characteristic':p,'free_kernels':len(free),'free_pairs':len(freepairs),'stratum_counts':dict(sorted(Counter({r:sum(v for s,v in strata.items() if len(s)==r) for r in range(4)}).items())),'moments':moments})
    conductor_checks = sum(audit_conductor(codes[f'U{k}'],codes[f'D{k}']) for k in (2,4))
    print(json.dumps({'status':'FINITE-EXACT independent audit','Steiner_designs':len(designs),'complementary_pairs':len(pairs),'field_cases':cases,'conductor_component_checks':conductor_checks,'arbitrary_inertia_cases':audit_arbitrary_inertia()},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
