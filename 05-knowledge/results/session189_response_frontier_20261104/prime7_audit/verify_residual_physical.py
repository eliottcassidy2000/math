#!/usr/bin/env python3
"""Literal singleton/pair controls for every residual singleton pattern."""
import json,hashlib
from verify_audit import P,ROOT,read,opposite,physical_rows,packed_covers,need
records=[]
for row in read('residual_mod4.json'):
 ns=tuple(2 if row['V']>>i&1 else 1 for i in range(7))
 a,c=row['query'];Q=tuple(row['Q']);physical=[physical_rows(Q,ns,a,c),physical_rows(opposite(Q),ns,a,c)]
 counts=[packed_covers(T) for T in physical]
 contrast=[(x-y)%4 for x,y in zip(*counts)]
 need(contrast==[2,0,0,0,0],'literal residual all-cover contrast')
 records.append(dict(quotient=list(Q),nontrivial_mask=row['V'],query=[a,c],sizes=list(ns),physical_rows=physical,cover_counts=counts,contrast_mod4=contrast))
print(json.dumps(dict(status='PASS',physical_graphs=12,maximum_physical_order=max(len(T) for r in records for T in r['physical_rows']),records=records,input_sha256=hashlib.sha256((P/'residual_mod4.json').read_bytes()).hexdigest()),sort_keys=True,indent=2))
