#!/usr/bin/env python3
"""Deterministically regenerate cubic search inputs from exact orbit data."""
from audit_basis import generate,nullspace
from pathlib import Path
root=Path(__file__).resolve().parent
m=8;reps,tg,mons,pts,basis,can=generate(m,3)
parity=[sum(((b>>(m*j))&255).bit_count()%2<<k for k,b in enumerate(basis)) for j in range(len(reps))]
ns=nullspace(parity,len(mons))
with (root/'cubic8.txt').open('w') as f:
 print(len(ns),len(reps),file=f)
 for c in ns:
  b=0
  for i,v in enumerate(basis):
   if c>>i&1:b^=v
  print(b&((1<<64)-1),b>>64,c,file=f)
 print(*can,file=f)
with (root/'cubic8.raw.txt').open('w') as f:
 print(len(basis),len(reps),file=f)
 for z in basis:
  tr=sum(((z>>(8*i))&255).bit_count()%2<<i for i in range(len(reps)))
  print(z&((1<<64)-1),z>>64,tr,file=f)
 print(*can,file=f)
print('Regenerated cubic inputs: 64 raw monomials, 56-dimensional even-output subspace.')
