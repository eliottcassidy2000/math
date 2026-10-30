#!/usr/bin/env python3
"""Larger DP-total controls; deliberately not a raw-ring audit."""
from prime_rank import PrimeBank
from literal_ring import grouped_sum
import json
rows=[]
for p,k in ((5,4),(5,5),(7,3),(13,2)):
 B=PrimeBank(p,k);value=B.prime_total()
 if value!=grouped_sum(p,k)[0]:raise RuntimeError('DP/convolution mismatch')
 rows.append(dict(p=p,k=k,dimension=p**k-1,prime_orbit_bits=value.bit_length(),table_entries=len(B._prime)))
print(json.dumps(dict(status='PASS',scope='DP totals versus independent convolution; no raw-ring test in these larger cases',families=rows),indent=2))
