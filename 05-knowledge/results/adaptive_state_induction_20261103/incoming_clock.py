#!/usr/bin/env python3
"""Small exact control for the all-depth law proved in INCOMING I1-I2."""
import json
from math import gcd


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def row(v, m):
    return tuple(sum(v[i]*m[i][j] for i in range(2)) for j in range(2))


def det(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def valuation(v):
    need(any(v), 'nonzero row')
    k = 0
    while all(x % 3 == 0 for x in v):
        v = tuple(x//3 for x in v)
        k += 1
    return k


Q = ((77,152),(28,55))
U = ((3395,6688),(1232,2427))
need(det(Q) == -21 and det(U) == 49, 'determinants')
need(tuple(row(v,Q) for v in Q) == tuple(tuple(3*x for x in v) for v in U), 'Q squared = 3U')
need(all(U[i][j] == 44*Q[i][j]+7*(i == j) for i in range(2) for j in range(2)), 'Cayley Hamilton')
need(gcd(*[x for v in U for x in v]) == 1, 'Smith first divisor')
need(all((a-21*b) % 49 == 0 for a,b in U), 'row lattice inclusion; equal index proves equality')
seeds = [(31,(25,34)),(95,(75,102)),(97,(75,102))]
controls=[]
for h,beta in seeds:
    need((beta[0]-21*beta[1]) % 49 != 0, 'two-symbol return blocked')
    v=beta; levels=[]
    for j in range(13):
        val=valuation(v)
        expected=valuation(beta)+j//2 if j % 2 == 0 else valuation(row(beta,Q))+j//2
        need(val == expected, 'power control')
        levels.append(val);v=row(v,Q)
    controls.append({'seed_H':h,'bulk':beta,'next_bulk':row(beta,Q),'lattice_residue_mod49':(beta[0]-21*beta[1]) % 49,'levels_j0_through12':levels})
need(controls[0]['levels_j0_through12'] == [(j+1)//2 for j in range(13)], 'actual H31 ceil clock')
print(json.dumps({'status':'PASS','scope':'Exact two-by-two identities and finite controls; all-depth proof in INCOMING I1-I2','Q':Q,'U':U,'det_Q':det(Q),'det_U':det(U),'Smith_U':[1,49],'controls':controls},indent=2,sort_keys=True))
