#!/usr/bin/env python3
"""Exact finite-field cone transport, including even extension degrees.

No producer code is imported. The six graph edge lists are pinned inputs.
The output is raw character sums, never an asserted pure-motive Euler factor.
"""
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
import hashlib
import json

GATES = 0


def check(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


class Field:
    def __init__(self, p, degree):
        self.p, self.degree, self.q = p, degree, p**degree
        if degree == 1:
            self.modulus = [0, 1]
        elif degree == 2:
            d = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
            self.modulus = [(-d) % p, 0, 1]
        elif (p, degree) == (3, 3):
            self.modulus = [2, 2, 0, 1]  # x^3-x-1 has no root in F3.
        else:
            raise ValueError('undeclared field')
        self.digits = [tuple((i//p**j) % p for j in range(degree)) for i in range(self.q)]
        self.add = [[self.encode([(x+y) % p for x,y in zip(a,b)]) for b in self.digits]
                    for a in self.digits]
        self.neg = [self.encode([(-x) % p for x in a]) for a in self.digits]
        self.sub = [[self.add[a][self.neg[b]] for b in range(self.q)] for a in range(self.q)]
        self.mul = [[self.multiply_digits(a,b) for b in self.digits] for a in self.digits]
        self.inv = [0] + [next(b for b in range(1,self.q) if self.mul[a][b] == 1)
                          for a in range(1,self.q)]
        squares = {self.mul[x][x] for x in range(1,self.q)}
        self.chi = [0] + [1 if x in squares else -1 for x in range(1,self.q)]
        check(len(squares) == (self.q-1)//2, 'quadratic-character square count')
        check(all(self.power(a,self.q) == a for a in range(self.q)), 'Frobenius fixed field')
        check(all(self.power(a,self.q-1) == 1 for a in range(1,self.q)), 'field nonzero powers')
        check(all(self.chi[self.mul[a][b]] == self.chi[a]*self.chi[b]
                  for a in range(self.q) for b in range(self.q)), 'character multiplicativity')
        self.epsilon = self.chi[self.neg[1]]
        check(self.epsilon == (-1)**((self.q-1)//2), 'minus-one character')
        self.L = [sum(self.chi[self.mul[self.mul[x][self.sub[x][1]]][self.sub[x][lam]]]
                      for x in range(self.q)) for lam in range(self.q)]
        check(self.L[0] == -self.epsilon and self.L[1] == -1, 'singular Legendre fibres')
        check(all(x*x <= 4*self.q for x in self.L[2:]), 'smooth Legendre Hasse control')

    def encode(self, a):
        return sum(x*self.p**j for j,x in enumerate(a))

    def multiply_digits(self, a, b):
        p, d = self.p, self.degree
        c = [0]*(2*d-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                c[i+j] += x*y
        for k in range(2*d-2,d-1,-1):
            for j in range(d):
                c[k-d+j] -= c[k]*self.modulus[j]
        return self.encode([x % p for x in c[:d]])

    def power(self, a, n):
        answer = 1
        while n:
            if n & 1:
                answer = self.mul[answer][a]
            a = self.mul[a][a]
            n //= 2
        return answer


def cones(edges):
    E = {tuple(e) for e in edges}
    adj = {v: {u for u in range(7) if tuple(sorted((u,v))) in E} for v in range(7)}
    eligible = [v for v in range(7) if len(adj[v]) == 3]
    choices = [C for C in combinations(eligible,3)
               if all(b not in adj[a] for a,b in combinations(C,2))]
    check(bool(choices), 'three independent cubic cones')
    C = choices[0]
    return C, adj


def cone_sum(F, edges, direct=False):
    C, adj = cones(edges)
    B = [v for v in range(7) if v not in C]
    exponent = sum(u > c for c in C for u in adj[c])
    base_edges = [e for e in edges if e[0] in B and e[1] in B]
    total = 0
    for s in range(2,F.q):
        for t in range(2,F.q):
            if s == t:
                continue
            X = dict(zip(B,(0,1,s,t)))
            value = 1
            for i,j in base_edges:
                value *= F.chi[F.sub[X[j]][X[i]]]
            for c in C:
                a,b,d = sorted(adj[c])
                if direct:
                    local = 0
                    for x in range(F.q):
                        prod_value = 1
                        for u in adj[c]:
                            diff = F.sub[X[u]][x] if c < u else F.sub[x][X[u]]
                            prod_value *= F.chi[diff]
                        local += prod_value
                else:
                    delta = F.sub[X[b]][X[a]]
                    lam = F.mul[F.sub[X[d]][X[a]]][F.inv[delta]]
                    local = F.epsilon**sum(u > c for u in adj[c]) * F.chi[delta]*F.L[lam]
                value *= local
            total += value
    return total, exponent


def closed_raw(F):
    ans = [0]*4
    for s in range(2,F.q):
        for t in range(2,F.q):
            if s == t:
                continue
            ratio = F.mul[t][F.inv[s]]
            cs,cs1,ct,ct1 = F.chi[s],F.chi[F.sub[s][1]],F.chi[t],F.chi[F.sub[t][1]]
            ls,lt,lr = F.L[s],F.L[t],F.L[ratio]
            ans[0] -= cs*ct*ct1*ls*ls*lr
            ans[1] -= cs1*ct*ct1*ls*ls*lr
            ans[2] += cs*cs1*ct*ls*lt*lr
            ans[3] += cs1*ct1*ls*lt*lr
    return ans


def transitive_census(F, k):
    if F.epsilon != -1:
        return {'defined': False, 'reason': 'quadratic character is symmetric', 'count': None}
    out = [sum(1 << j for j in range(F.q) if F.chi[F.sub[j][i]] == 1) for i in range(F.q)]
    check(all(bool(out[i] >> j & 1) != bool(out[j] >> i & 1)
              for i in range(F.q) for j in range(i+1,F.q)), 'literal tournament orientation')
    @lru_cache(None)
    def count(mask, m):
        if m == 0:
            return 1
        if mask.bit_count() < m:
            return 0
        if m == 1:
            return mask.bit_count()
        total = 0
        bits = mask
        while bits:
            bit = bits & -bits
            bits -= bit
            total += count(mask & out[bit.bit_length()-1], m-1)
        return total
    return {'defined': True, 'count': count((1 << F.q)-1,k)}


def main():
    root = Path(__file__).resolve().parents[4]
    source = root/'04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/genuine37.json'
    graph_data = json.loads(source.read_text())
    ids = [2,4,5,6,9,18]
    fields = [Field(p,1) for p in (3,5,7,11,13,17,19)]
    fields += [Field(p,2) for p in (3,5,7,11,13,17,19)]
    fields += [Field(3,3)]
    prime_fields = {F.p:F for F in fields if F.degree == 1}
    rows = []
    for F in fields:
        sums, exponents = [], []
        for i in ids:
            value, exponent = cone_sum(F,graph_data[i]['edges'])
            sums.append(value); exponents.append(exponent)
            if F.q <= 27:
                check(cone_sum(F,graph_data[i]['edges'],direct=True)[0] == value,
                      'literal edge-product sum equals corrected cone formula')
        raw = closed_raw(F)
        check([-F.epsilon*raw[0],-F.epsilon*raw[1],raw[2],raw[3]]
              == [sums[0],sums[2],sums[4],sums[5]], 'closed formulas with arithmetic orientation twist')
        check(sums[1] == F.epsilon*sums[0] and sums[3] == sums[2], 'all-odd-field paired identities')
        if F.degree == 2:
            Lp = prime_fields[F.p].L
            check(all(F.L[lam] == 2*F.p-Lp[lam]**2 for lam in range(2,F.p)),
                  'independent elliptic Frobenius-square control on base-field fibres')
        rows.append({'p':F.p,'degree':F.degree,'q':F.q,'modulus':F.modulus,
                     'chi_minus_one':F.epsilon,'cone_ids':ids,'cone_sums':sums,
                     'cone_sign_exponents':exponents,'fixed_formula_values':raw,
                     'TT6':transitive_census(F,6) if F.q <= 49 else None})
    check(next(r for r in rows if r['q']==27)['TT6']['count'] == 0, 'certified TT6(P27)=0')
    check(next(r for r in rows if r['q']==9)['TT6']['count'] is None, 'undefined is not zero')
    print(json.dumps({'status':'FINITE-EXACT raw cone transport; motives and conductor OPEN',
                      'gates':GATES,'universe':'six graphs; seven prime fields, their quadratic extensions, F27',
                      'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                      'rows':rows},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
