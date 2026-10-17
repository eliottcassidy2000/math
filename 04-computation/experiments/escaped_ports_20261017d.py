#!/usr/bin/env python3
"""Exact induced unit graphs in sqrt3 and eleventh-root escapes.

A good-prime reduction is only a sieve for edges. Every retained candidate
is tested in the characteristic-zero number field. All eight ports survive.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import numpy as np
import sympy as s
from pysat.solvers import Solver
from golden_boundary_20261017c import PORT_EDGES, patterns


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


class Field:
    def __init__(self, polynomials):
        self.polys = polynomials
        self.degrees = [len(f)-1 for f in polynomials]
        self.basis = list(product(*(range(d) for d in self.degrees)))
        self.index = {e:i for i,e in enumerate(self.basis)}
        self.n = len(self.basis)
        self.zero = (Q(0),)*self.n
        self.one = self.scalar(1)
        powers = []
        for f,d in zip(polynomials,self.degrees):
            v = [[int(i==k) for i in range(d)] for k in range(d)]
            for k in range(d,2*d-1):
                last = v[-1]
                v.append([-last[-1]*f[0]]+[last[i-1]-last[-1]*f[i] for i in range(1,d)])
            powers.append(v)
        self.table = {}
        for i,a in enumerate(self.basis):
            for j,b in enumerate(self.basis):
                z = {}
                choices = [[(k,c) for k,c in enumerate(powers[t][a[t]+b[t]]) if c] for t in range(len(a))]
                for terms in product(*choices):
                    coeff = 1
                    for k,c in terms:
                        coeff *= c
                    z[self.index[tuple(k for k,c in terms)]] = coeff
                self.table[i,j] = tuple(z.items())

    def scalar(self, x):
        return (Q(x),)+(Q(0),)*(self.n-1)

    def gen(self, i):
        e = [0]*len(self.degrees);e[i] = 1
        z = list(self.zero);z[self.index[tuple(e)]] = Q(1)
        return tuple(z)

    def add(self, a,b):return tuple(x+y for x,y in zip(a,b))
    def sub(self, a,b):return tuple(x-y for x,y in zip(a,b))
    def scale(self, a,c):return tuple(x*c for x in a)

    def mul(self, a,b):
        out = [Q(0)]*self.n
        aa = [(i,c) for i,c in enumerate(a) if c]
        bb = [(i,c) for i,c in enumerate(b) if c]
        for i,c in aa:
            for j,d in bb:
                v = c*d
                for k,e in self.table[i,j]:out[k] += v*e
        return tuple(out)

    def power(self, a,n):
        z = self.one
        while n:
            if n&1:z = self.mul(z,a)
            a = self.mul(a,a);n //= 2
        return z

    def inv(self, a):
        cols = []
        for j in range(self.n):
            e = list(self.zero);e[j] = Q(1)
            cols.append(self.mul(a,tuple(e)))
        out = tuple(Q(x) for x in s.Matrix.hstack(*(s.Matrix(c) for c in cols)).inv()[:,0])
        need(self.mul(a,out)==self.one,'exact field inverse')
        return out

    def substitution(self, gens):
        ps = [[self.power(g,k) for k in range(d)] for g,d in zip(gens,self.degrees)]
        out = []
        for e in self.basis:
            v = self.one
            for i,k in enumerate(e):v = self.mul(v,ps[i][k])
            out.append(v)
        return out

    def apply(self, a, table):
        out = self.zero
        for c,b in zip(a,table):
            if c:out = self.add(out,self.scale(b,c))
        return out


P = 1319
def fm(a,b):return ((a[0]*b[0]-a[1]*b[1])%P,(a[0]*b[1]+a[1]*b[0])%P)
def fp(a,n):
    z = (1,0)
    while n:
        if n&1:z = fm(z,a)
        a = fm(a,a);n //= 2
    return z
def root_unity(m):
    for x in range(P):
        z = fp((x,1),(P*P-1)//m)
        if fp(z,m)==(1,0) and all(fp(z,m//int(q))!=(1,0) for q in s.factorint(m)):
            return z
    raise RuntimeError('root of unity not found')


def make(mode):
    x = s.symbols('x')
    if mode=='sqrt3':
        m = 15
        polys = [[int(s.expand(s.cyclotomic_poly(m,x)).coeff(x,k)) for k in range(9)], [11,0,1],[-4,-4,2,1]]
        f = Field(polys);z,j,w = [f.gen(i) for i in range(3)]
        conj = f.substitution([f.power(z,m-1),f.scale(j,-1),w])
        escape = f.substitution([f.power(z,11),j,w])
        eta = f.add(f.one,f.power(z,5))
        roots = [root_unity(m),(0,int(s.sqrt_mod(11,P))), (1112,0)]
    else:
        m = 55
        polys = [[int(s.expand(s.cyclotomic_poly(m,x)).coeff(x,k)) for k in range(41)],[-4,-4,2,1]]
        f = Field(polys);z,w = [f.gen(i) for i in range(2)]
        j = f.zero
        for a in range(1,11):j = f.add(j,f.scale(f.power(z,5*a),int(s.legendre_symbol(a,11))))
        conj = f.substitution([f.power(z,m-1),w])
        escape = f.substitution([f.power(z,36),w])
        eta = None
        roots = [root_unity(m),(1112,0)]
    alpha = f.power(z,m//5)
    phi = f.add(f.one,f.add(alpha,f.power(alpha,4)))
    root5 = f.sub(f.scale(phi,2),f.one)
    ports = [f.zero,f.one,f.add(f.one,alpha),f.add(f.add(f.one,alpha),f.power(alpha,2)),f.scale(f.power(alpha,4),-1),
             f.scale(f.sub(f.scalar(3),f.mul(j,root5)),Q(1,8)),f.scalar(Q(3,4)),f.scale(f.add(f.scalar(3),f.mul(j,root5)),Q(1,8))]
    units = []
    for h in (f.add(w,f.scalar(3)),f.add(f.add(w,f.scalar(3)),phi)):
        units.append(f.mul(f.sub(h,j),f.inv(f.add(h,j))))
    need(f.mul(j,j)==f.scalar(-11),'quadratic Gauss sum / imaginary quadratic generator')
    need(all(f.mul(v,f.apply(v,conj))==f.one for v in units),'cubic and mixed-golden units')
    # Evaluate basis at one fixed prime. Every denominator must remain a unit.
    basis_res = []
    for e in f.basis:
        a = (1,0)
        for r,k in zip(roots,e):a = fm(a,fp(r,k))
        basis_res.append(a)
    for polynomial, r in zip(polys,roots):
        v = (0,0)
        for c in reversed(polynomial):
            a,b = fm(v,r);v = ((a+c)%P,b)
        need(v==(0,0),'all defining relations vanish at the good prime')
    def reduce(v):
        a=b=0
        for c,(x,y) in zip(v,basis_res):
            need(c.denominator%P!=0,'no denominator vanishes at the sieve prime')
            cc = c.numerator*pow(c.denominator,-1,P)%P
            a += cc*x;b += cc*y
        return a%P,b%P
    for k in range(f.n):
        e=list(f.zero);e[k]=Q(1);e=tuple(e)
        a,b=reduce(e)
        need(reduce(f.apply(e,conj))==(a,-b%P),'residue involution is exactly complex conjugation')
    return f,z,j,alpha,ports,units,eta,conj,escape,reduce


def run(mode,assembly='copies'):
    f,z,j,alpha,ports,units,eta,conj,escape,reduce = make(mode)
    if mode=='sqrt3':
        rho = f.scale(f.add(f.scalar(5),j),Q(1,6))
        a = f.add(f.one,eta)
        seed = [f.zero,f.one,eta,a,rho,f.mul(rho,eta),f.mul(rho,a)]
        multipliers = [f.mul(f.power(alpha,b),u) for b in range(5) for u in [f.one]+units]
    else:
        u,v = units
        seed = [f.zero,f.one,u,f.add(f.one,u),v,f.add(f.one,v)]
        multipliers = [f.mul(f.power(z,5*b),a) for b in range(11) for a in (f.one,alpha,units[0])]
    if assembly=='triple':
        need(mode=='sqrt3','triple assembly currently uses the seven-vertex spindle seed')
        seed = sorted({f.add(a,f.add(f.mul(units[0],b),f.mul(units[1],c)))
                       for a in seed for b in seed for c in seed})
        multipliers = [f.one]
    points = list(ports);seen = set(points)
    for p in ports:
        for a in multipliers:
            for q in seed:
                v = f.add(p,f.mul(a,q))
                if v not in seen:seen.add(v);points.append(v)
    residues = np.array([reduce(v) for v in points],dtype=np.int64)
    edges = [];candidates=0;new_edges=0
    for a in range(len(points)):
        diff = residues[a+1:]-residues[a]
        good = np.flatnonzero((diff[:,0]**2+diff[:,1]**2)%P==1)
        for off in good:
            b = a+1+int(off);candidates += 1
            v = f.sub(points[a],points[b])
            if f.mul(v,f.apply(v,conj))==f.one:
                edges.append((a,b));new_edges += v!=f.apply(v,escape)
    need([e for e in edges if e[1]<8]==sorted(tuple(sorted(e)) for e in PORT_EDGES),'the induced graph on the eight old ports is unchanged')
    # Relative degeneracy: removing an exterior vertex with <=5 total remaining
    # neighbours is an unconditional extension step for any fixed six-colouring.
    neighbours = [set() for _ in points]
    for a,b in edges:neighbours[a].add(b);neighbours[b].add(a)
    live = set(range(len(points)));order=[]
    while True:
        batch = [a for a in sorted(live) if a>=8 and len(neighbours[a]&live)<=5]
        if not batch:break
        for a in batch:live.remove(a);order.append(a)
    pats = patterns(8,PORT_EDGES);witnesses=[]
    with Solver(name='g3') as sol:
        for a in live:
            sol.add_clause([6*a+c+1 for c in range(6)])
            for c in range(6):
                for d in range(c):sol.add_clause([-6*a-c-1,-6*a-d-1])
        for a,b in edges:
            if a in live and b in live:
                for c in range(6):sol.add_clause([-6*a-c-1,-6*b-c-1])
        for pat in pats:
            ok = sol.solve(assumptions=[6*a+c+1 for a,c in enumerate(pat)])
            need(ok,'a boundary is blocked: save separately before claiming any forcing theorem')
            model = {v for v in sol.get_model() if v>0}
            colours = [-1]*len(points)
            for a in live:colours[a] = next(c for c in range(6) if 6*a+c+1 in model)
            for a in reversed(order):
                used = {colours[b] for b in neighbours[a] if colours[b]>=0}
                colours[a] = next(c for c in range(6) if c not in used)
            need(tuple(colours[:8])==pat and all(colours[a]!=colours[b] for a,b in edges),'independent literal edge/port replay')
            witnesses.append(''.join(map(str,colours)))
    payload = dict(mode=mode,assembly=assembly,field_dimension=f.n,sieve_prime=P,vertices=len(points),edges=edges,candidates=candidates,
                   escaping_edges=new_edges,relative_six_core=len(live)-8,peeling_order=order,
                   coordinates=[[str(c) for c in p] for p in points],patterns=pats,witnesses=witnesses,
                   status='FINITE-EXACT complete induced graph; all570 boundary patterns extend')
    return payload


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['sqrt3','eleven'],required=True);p.add_argument('--assembly',choices=['copies','triple'],default='copies');p.add_argument('--output',required=True);a=p.parse_args()
    result=run(a.mode,a.assembly)
    Path(a.output).write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('edges','coordinates','patterns','witnesses','peeling_order')},sort_keys=True))
