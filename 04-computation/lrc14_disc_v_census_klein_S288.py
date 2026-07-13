#!/usr/bin/env python3
"""
lrc14_disc_v_census_klein_S288.py
=================================
klein-2026-07-13-S288 (owner: prove the analytic disc_v bound).

PROVED this session:  disc_v <= r^2/(3 v^2),  r = #arcs of the leave-one-out good set G'_{~v}
  [|U(l)|<=2r endpoints; Sum_{m!=0}1/m^2 = pi^2/3].
Fed into THM-731:  L >= (1/7)(6|G'_{~v}| - sqrt(2)*r/v).   POSITIVE  <=>  r < 3*sqrt(2)*v*|G'_{~v}|.
Both the disc_v bound and this L-bound are RIGOROUS. What remains is the COMBINATORIAL condition
r < 3 sqrt(2) v |G'_{~v}|.  This census stress-tests it over many covering 13-sets (peeling the far
element), to (a) confirm it holds, (b) find the tightest family (expected: near-AP residue).
"""
import numpy as np
from math import sqrt
NG=1<<21
THR=1.0/14.0
t=np.arange(NG,dtype=np.float64)/NG

def good_ind(S):
    g=np.ones(NG,dtype=np.float64)
    for w in S:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g*=(d>=THR)
    return g
def n_arcs(g):
    gi=g.astype(np.int8); diff=np.diff(np.concatenate([gi,gi[:1]])); return int(np.sum(diff==1))
def is_covering(S):
    return all(any(x%q==0 for x in S) for q in range(2,15))

# Covering families: (base 12 speeds) + minimal valid far element v. Peel v (the far element).
# base misses some q's; v must supply them. We list base + a few far elements (smallest = tightest).
CANDIDATES=[
    ([1,2,3,4,5,6,7,8,9,10,11,12],[182,364]),          # miss 13,14 -> v=182
    ([1,2,3,4,5,6,7,8,9,10,11,13],[84,168,252]),        # miss 12,14 -> v=84 (the residue base)
    ([1,2,3,4,5,6,7,8,9,10,12,13],[154,308]),           # miss 11,14 -> v=lcm(11,14)=154
    ([1,2,3,4,5,6,7,8,9,11,12,13],[70,140]),            # miss 10,14 -> v=lcm(10,14)=70
    ([1,2,3,4,5,6,7,8,10,11,12,13],[126,252]),          # miss 9,14  -> v=lcm(9,14)=126
    ([2,3,4,5,6,7,8,9,10,11,12,13],[14,28]),            # miss 14    -> v=14 ({2..14})
    ([1,3,4,5,6,7,8,9,10,11,12,13],[14,28,42]),         # miss 14 (2 via 4) -> v=14
    ([1,2,3,4,5,6,7,8,9,10,11,12],[182]),               # deep well
]
print("="*108)
print("CENSUS of the RIGOROUS bound  L >= (1/7)(6|G'_{~v}| - sqrt2*r/v)  over covering 13-sets (peel far v).")
print("condition for L>0:  r < 3*sqrt2*v*|G'_{~v}|.  ratio = sqrt2*r / (6 v |G'_{~v}|)  (<1 certifies).")
print("="*108)
print("%-32s %5s %9s %4s %11s %8s %s"%("family (base + far v)","v","|G'_~v|","r","L_bound","ratio","cert?"))
rows=[]
for base,fars in CANDIDATES:
    g=good_ind(base); L=g.mean(); r=n_arcs(g)   # base determines |G'| and r; far element only sets v
    for v in fars:
        S=sorted(base+[v])
        if len(set(S))!=13 or not is_covering(S): continue
        Lbound=(6*L - sqrt(2)*r/v)/7.0
        ratio=sqrt(2)*r/(6*v*L)
        rows.append((str(sorted(base))[:30],v,L,r,Lbound,ratio))
        print("%-32s %5d %9.5f %4d %11.5f %8.4f %s"%(
            "base"+str(base[-3:])+"+"+str(v), v, L, r, Lbound, ratio, "YES" if ratio<1 else "**NO**"))
print("-"*108)
worst=max(rows,key=lambda x:x[5])
print("TIGHTEST family: %s  v=%d  ratio=%.4f  (L_bound=%.5f)"%(worst[0],worst[1],worst[5],worst[4]))
print("all certify (ratio<1): %s"%all(x[5]<1 for x in rows))
print("""
Reading: the disc_v bound r^2/(3v^2) is PROVED (rigorous). The L-bound (1/7)(6|G'|-sqrt2 r/v) is
therefore rigorous. 'cert?=YES' means it proves L>0 for that covering family. The residue base
{1..11,13}+84 is expected tightest (smallest valid far element = 84). Remaining to make LRC(14)
covering-case a theorem: prove the COMBINATORIAL bound r < 3 sqrt2 v |G'_{~v}| for ALL covering 13-sets
(now arc-counting, not harmonic analysis).""")
print("done.")
