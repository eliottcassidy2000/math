#!/usr/bin/env python3
"""Audit B: sizes, L/rhoS split and union of the four crit_*.json vertex sets (exact coordinates, session basis
Q(sqrt3, sqrt5, sqrt11) for x and y).  The rho S361 vertices are recognised by a nonzero sqrt5-type component
(Parts' L727 lives in Q(sqrt-3, sqrt-11); rho adds sqrt5).  Usage: python3 hn_crit_union.py crit*.json"""
import sys, json
from fractions import Fraction as Fr
U = set()
for f in sys.argv[1:]:
    d = json.load(open(f))
    V = [tuple(tuple(Fr(s) for s in part) for part in v) for v in d['vertices']]
    # session bits: 1 = sqrt3, 2 = sqrt5, 4 = sqrt11; any component with bit 2 set involves sqrt5
    rhoS = sum(1 for v in V if any(c != 0 for part in v for b, c in enumerate(part) if b & 2))
    E = d['edge_list']
    print(f"{f.split('/')[-1]}: n = {d['n']} (len {len(V)}, distinct {len(set(V))}), edges {d['edges']} (len {len(E)}), "
          f"L727 part {len(V) - rhoS}, rho S361 part {rhoS}")
    U |= set(V)
print(f"union of the vertex sets: {len(U)}")
