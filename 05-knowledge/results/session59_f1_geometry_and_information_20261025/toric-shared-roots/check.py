#!/usr/bin/env python3
"""Exact, bounded checks for RESULT.md. Python standard library + SymPy.

No downloaded data. Field elements are encoded in a polynomial basis; the
quadratic-field implementation independently exercises odd extension counts.
"""
from __future__ import annotations

from itertools import product
from math import comb
from pathlib import Path
import json
import sympy as sp


class Field:
    def __init__(self, p: int, degree: int = 1):
        assert degree in (1, 2)
        self.p, self.degree, self.q = p, degree, p**degree
        if degree == 1:
            self.A = self.B = 0
        elif p == 2:
            self.A = self.B = 1  # theta^2=theta+1
        else:
            self.A = next(a for a in range(2, p)
                          if pow(a, (p-1)//2, p) == p-1)
            self.B = 0  # theta^2=A, A nonsquare

    def constant(self, x):
        return x % self.p

    def add(self, x, y):
        if self.degree == 1:
            return (x+y) % self.p
        a, b, c, d = x % self.p, x // self.p, y % self.p, y // self.p
        return (a+c) % self.p + self.p*((b+d) % self.p)

    def neg(self, x):
        if self.degree == 1:
            return (-x) % self.p
        return (-x % self.p) + self.p*((-(x//self.p)) % self.p)

    def mul(self, x, y):
        if self.degree == 1:
            return x*y % self.p
        a, b, c, d = x % self.p, x // self.p, y % self.p, y // self.p
        real = (a*c+b*d*self.A) % self.p
        imag = (a*d+b*c+b*d*self.B) % self.p
        return real+self.p*imag

    def power(self, x, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, x)
            x = self.mul(x, x)
            n //= 2
        return out

    def inv(self, x):
        assert x
        return self.power(x, self.q-2)

    def square(self, x):
        return self.mul(x, x)

    def chi(self, x):
        assert self.p != 2
        if not x:
            return 0
        y = self.power(x, (self.q-1)//2)
        assert y in (1, self.neg(1))
        return 1 if y == 1 else -1


def projective_tuples(field, ncoords):
    # The first nonzero coordinate is one, a unique field-point representative.
    for leading in range(ncoords):
        for tail in product(range(field.q), repeat=ncoords-leading-1):
            yield (0,)*leading+(1,)+tail


def raw_projective_count(field, r):
    two = field.constant(2)
    count = 0
    boundary = 0
    for point in projective_tuples(field, 2*r+1):
        c, *ab = point
        rhs = field.mul(two, field.square(c))
        if all(field.add(field.square(ab[2*j]), field.square(ab[2*j+1])) == rhs
               for j in range(r)):
            count += 1
            boundary += (c == 0)
    return count, boundary


def pair_counts(field):
    q = field.q
    square = [field.square(x) for x in range(q)]
    hist = [0]*q
    for a in range(q):
        for b in range(q):
            hist[field.add(square[a], square[b])] += 1
    two = field.constant(2)
    return [hist[field.mul(two, square[c])] for c in range(q)]


def D(r, q):
    return sum(2**s*comb(r, s)*(q-1)**(s-1) for s in range(1, r+1))


def field_checks():
    fields = [Field(p) for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]]
    fields += [Field(p, 2) for p in [2, 3, 5, 7, 11]]
    rows = []
    for field in fields:
        q = field.q
        for x in range(1, q):
            assert field.mul(x, field.inv(x)) == 1
            assert field.power(x, q-1) == 1
        counts = pair_counts(field)
        for r in range(1, 7):
            cone = sum(n**r for n in counts)
            assert (cone-1) % (q-1) == 0
            got = (cone-1)//(q-1)
            boundary = (counts[0]**r-1)//(q-1)
            if field.p != 2:
                eps = field.chi(field.neg(1))
                predicted = (q-eps)**r + (1+eps)*D(r, q)//2
                assert counts[0] == q+eps*(q-1)
                assert counts[1:] == [q-eps]*(q-1)
                assert boundary == (1+eps)*D(r, q)//2
            else:
                eps = None
                predicted = sum(q**j for j in range(r+1))
            assert got == predicted
            row = dict(p=field.p, degree=field.degree, q=q, r=r,
                       epsilon=eps, count=got, boundary=boundary,
                       independent_pair_enumeration=True)
            if ((q <= 3 and r <= 4) or (q <= 7 and r <= 2)
                    or (q == 4 and r <= 3)):
                direct, direct_boundary = raw_projective_count(field, r)
                assert direct == got and direct_boundary == boundary
                row['direct_projective_enumeration'] = True
            rows.append(row)
    return rows


def symbolic_checks():
    c, a, b, d, e = sp.symbols('c a b d e')
    i = sp.I
    u, v = (a+i*b)/(1+i), (a-i*b)/(1-i)
    w, z = (d+i*e)/(1+i), (d-i*e)/(1-i)
    assert sp.simplify(u*v-(a*a+b*b)/2) == 0
    assert sp.simplify(w*z-(d*d+e*e)/2) == 0
    assert sp.simplify(((1+i)*u+(1-i)*v)/2-a) == 0
    assert sp.simplify(((1+i)*u-(1-i)*v)/(2*i)-b) == 0
    point = {c:1, a:0, b:sp.sqrt(2), d:1, e:1}
    equations = [a*a+b*b-2*c*c, d*d+e*e-2*c*c]
    radicands = [a*a+d*d-c*c, 3*c*c-a*a-d*d,
                 c*c+a*a-d*d, c*c-a*a+d*d]
    assert [f.subs(point) for f in equations] == [0, 0]
    assert [f.subs(point) for f in radicands] == [0, 2, 0, 2]
    # Tangent vector dc=da=db=0, dd=1, de=-1.
    deriv = lambda f: (sp.diff(f,d)-sp.diff(f,e)).subs(point)
    assert [deriv(f) for f in equations] == [0, 0]
    assert [deriv(f) for f in radicands] == [2, -2, -2, 2]
    T = sp.Symbol('T')
    polynomials=[]
    for r in range(1, 9):
        plus = sp.Poly(sp.expand((T-1)**r+D(r,T)), T)
        minus = sp.Poly(sp.expand((T+1)**r), T)
        delta = plus-minus
        assert (delta.is_zero) == (r <= 2)
        if r >= 3:
            assert delta.coeff_monomial(T**(r-1)) == 2**r-2*r
        assert plus.eval(1) == 2*r
        polynomials.append(dict(r=r, split=str(plus.as_expr()),
                                nonsplit=str(minus.as_expr()),
                                difference=str(delta.as_expr())))
    return dict(coordinate_change=True,
                interior_branch_point=['1','0','sqrt(2)','1','1'],
                radicands=[0,2,0,2], tangent_radicand_values=[2,-2,-2,2],
                polynomials=polynomials)


def semigroup_checks():
    rows=[]
    for r in range(1,5):
        zero=(0,)*r
        generators={zero}
        for j in range(r):
            for sign in [-1,1]:
                vec=list(zero); vec[j]=sign
                generators.add(tuple(vec))
        sums={zero}
        for n in range(0,6):
            expected={v for v in product(range(-n,n+1), repeat=r)
                      if sum(abs(x) for x in v)<=n}
            assert sums == expected
            rows.append(dict(r=r,degree=n,distinct_monomials=len(sums)))
            sums={tuple(a+b for a,b in zip(v,g)) for v in sums for g in generators}
    return rows


def blue_assignments():
    bases=[]; plus=[]; minus=[]; full=[]
    for c,a,b,d,e in product([0,1], repeat=5):
        if not (a+b==2*c and d+e==2*c):
            continue
        bases.append((c,a,b,d,e))
        for h,hp in product([0,1], repeat=2):
            if h+c==a+d and hp+a+d==3*c:
                plus.append((c,a,b,d,e,h,hp))
        for s,t in product([0,1], repeat=2):
            if s+d==c+a and t+a==c+d:
                minus.append((c,a,b,d,e,s,t))
        for h,hp,s,t in product([0,1], repeat=4):
            if (h+c==a+d and hp+a+d==3*c
                    and s+d==c+a and t+a==c+d):
                full.append((c,a,b,d,e,h,hp,s,t))
    data=dict(base=bases,plus=plus,minus=minus,full=full)
    for rows in data.values():
        assert len(rows)==2 and all(set(row) in [{0},{1}] for row in rows)
    # Check the actual shared-root pullback, not product of marginal cardinalities.
    gluing=[p+m[5:] for p in plus for m in minus if p[:5]==m[:5]]
    assert sorted(gluing)==sorted(full)
    return data


def main():
    data=dict(status='PASS',
              scope='Bounded finite checks, not substitutes for RESULT.md proofs',
              field_rows=field_checks(), symbolic=symbolic_checks(),
              semigroup=semigroup_checks(), blue_assignments=blue_assignments())
    Path(__file__).with_name('results.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS:',len(data['field_rows']),'field/r cases;',
          sum(x.get('direct_projective_enumeration',False) for x in data['field_rows']),
          'direct projective cases;',len(data['semigroup']),'semigroup cases;')
    print('coordinate change, branch tangent, polynomial boundary, and blue gluing passed')


if __name__=='__main__':
    main()
