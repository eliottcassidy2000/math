#!/usr/bin/env python3
"""mac-mini-2026-07-01-S87 -- integrating the merged-metagraph structure (S81-S86) with LRC (S60-S80).
Grounds Bridge 2: the tournament spectral-twin separator (adjacency A) and the LRC floor's danger relation
composed with itself (D D^T, HYP-3571) are the SAME move -- the Gram/2nd-moment spectrum of a 0/1 relation."""
import numpy as np
r=1/14; G=1<<20; t=np.arange(G)/G
def danger(v): return (np.abs(v*t-np.round(v*t))<r)
def gram_spec(S):
    D=np.array([danger(v) for v in S],float); M=(D@D.T)/G
    return np.round(sorted(np.linalg.eigvalsh(M),reverse=True),4), round(float(np.trace(M)),4), M
for name,S in [('construction',list(range(1,13))+[182]),('AP',list(range(1,14))),('GW',[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    ev,tr,M=gram_spec(S)
    print(f"{name:14s}: danger-Gram trace={tr} top5={ev[:5]}")
print("=> same trace (1st moment), distinct spectrum (fine structure) = LRC analog of the S86 spectral twins.")
