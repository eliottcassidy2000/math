#!/usr/bin/env python3
"""Checks incoming floating Boolean output against independent exact witnesses.
The exact witness proof is check.py; a small floating residual is NOT that proof.
"""
from math import isqrt
import hashlib,json
import numpy as np
N=1000000
prime=bytearray(b'\1')*(N+1);prime[:2]=b'\0\0'
for p in range(2,isqrt(N)+1):
    if prime[p]:prime[p*p:N+1:p]=b'\0'*((N-p*p)//p+1)
plist=[p for p in range(2,N+1) if prime[p]]
exact=np.zeros(N+1,dtype=bool)
for e in range(8,N+1,2):
    for p in plist:
        if p*2>=e:raise RuntimeError(('no exact witness',e))
        if prime[e-p]:exact[e]=True;break
for e in range(5,N+1,2):exact[e]=bool(prime[e-2])
P=np.frombuffer(prime,dtype=np.uint8).astype(np.float64)
size=1
while size<2*(N+1):size*=2
F=np.fft.rfft(P,size);conv=np.fft.irfft(F*F,size)[:N+1]
r=np.rint(conv).astype(np.int64)
e=np.arange(N+1);half=np.zeros(N+1,dtype=np.int64)
half[::2]=np.frombuffer(prime,dtype=np.uint8)[e[::2]//2]
distinct=r-half
wrong=np.flatnonzero((distinct>0)!=exact)
if len(wrong):raise RuntimeError(('Boolean support mismatch',wrong[:30]))
print(json.dumps({'status':'PASS','N':N,'exact_even_witnesses':len(range(8,N+1,2)),
 'support_mismatches':len(wrong),'exact_support_sha256':hashlib.sha256(exact.tobytes()).hexdigest(),
 'rounded_residual_max_diagnostic_only':float(np.max(np.abs(conv-r))),
 'proves':'This run of incoming FFT Boolean support agrees with independently exact support.',
 'does_not_prove':'No general FFT error bound; no equality of all pair multiplicities; no replay of full incoming census.'},indent=2,sort_keys=True))
