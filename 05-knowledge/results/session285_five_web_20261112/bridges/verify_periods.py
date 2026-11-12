#!/usr/bin/env python3
"""Exact period-rank controls; no numerical integration or rank tolerances."""
import json
import sympy as s
z=s.symbols('z'); I=s.I
cases=[([1],0),([I],1),([1,I],2),([1,I,1+I],3)]
out=[]
for cs, expected in cases:
    g=1+sum(c/(z-j)**2 for j,c in enumerate(cs)); h=g*g
    rows=[]
    for j,c in enumerate(cs):
        G=1+sum(d/s.Integer(j-k)**2 for k,d in enumerate(cs) if k!=j)
        D=-2*sum(d/s.Integer(j-k)**3 for k,d in enumerate(cs) if k!=j)
        r=s.expand(2*c*D); t=s.expand(j*r+2*c*G)
        if s.simplify(s.residue(h,z,j)-r)!=0: raise ValueError('h residue')
        if s.simplify(s.residue(z*h,z,j)-t)!=0: raise ValueError('zh residue')
        rows.append([s.im(t),s.im(r),-s.re(r)])
    M=s.Matrix(rows)
    if M.rank()!=expected: raise ValueError('rank')
    # Embed in a,b_real,b_imag,c,d: latter two are exactly free.
    E=M.row_join(s.zeros(len(rows),2))
    if len(E.nullspace())!=5-expected: raise ValueError('dimension')
    # Period-normalized holonomy: common factors do not affect rank or cocycle.
    mons=[]
    for row in rows:
        N=s.zeros(5)
        for k,x in enumerate(row): N[4,k]=x
        mons.append(N)
    for N in mons:
        for K in mons:
            if N*K!=s.zeros(5): raise ValueError('nilpotent product')
            if (s.eye(5)+N)*(s.eye(5)+K)!=s.eye(5)+N+K: raise ValueError('cocycle')
    out.append({'coefficients' :list(map(str,cs)),'period_matrix':[[str(x) for x in row] for row in rows],
                'rank':M.rank(),'stress_dimension':5-M.rank(),
                'determinant':str(M.det()) if M.rows==3 else None})
print(json.dumps({'status':'PASS','cases':out},indent=2,sort_keys=True))
