"""kps-2026-07-04: the OWNER's hint -- Fibonacci mod 7 has Pisano period 16. Explore the LRC(14)
connection: does a Fibonacci/linear-recurrence family have a PERIODIC danger structure (mod 7,
period 16) that certifies an infinite class with a FINITE check? The '7' is 1/14=1/(2*7)."""
import numpy as np
from math import gcd
from functools import reduce
def gcd_all(v): return reduce(gcd,v)
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M_of(v, grid=500000):
    v=np.array(sorted(set(int(x) for x in v)),dtype=np.float64)
    if len(v)!=13: return 9.0
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); return np.abs(np.outer(v,t)-np.round(np.outer(v,t))).min(axis=0).max()

# 1) verify Pisano period mod 7 = 16
F=[0,1]
for _ in range(60): F.append(F[-1]+F[-2])
fib_mod7=[f%7 for f in F]
print("Fib mod 7:", fib_mod7[:20], "... period 16?", fib_mod7[:16]==fib_mod7[16:32])
print("Fib mod 7 hits 0 at indices:", [i for i in range(48) if F[i]%7==0], "(period-8 zeros!)")
print("Fib mod 14 Pisano period = lcm(3,16) =", np.lcm(3,16))
print()

# 2) M of consecutive-Fibonacci families {F_k,...,F_{k+12}}, k=2.. -- periodic in k? loose?
print(f"{'k':>3} {'family (first/last)':>22} {'M':>8} {'M*14':>6} {'cov?':>5} {'k mod16':>7}")
Ms={}
for k in range(2,34):
    fam=[F[k+j] for j in range(13)]
    if gcd_all(fam)!=1: 
        g=gcd_all(fam); fam=[x//g for x in fam]
    M=M_of(fam)
    Ms[k]=M
    print(f"{k:>3} {str((fam[0],fam[-1])):>22} {M:>8.5f} {M*14:>6.3f} {str(is_cov(fam)):>5} {k%16:>7}")
print()
# 3) is M asymptotically constant? and does k mod 16 group them?
tail=[Ms[k] for k in range(18,34)]
print(f"M for large k (18..33): min={min(tail):.5f} max={max(tail):.5f} -- asymptote? all loose (M*14>1)? {min(tail)*14>1.0}")
# 4) the RESIDUE family mod 7: {F_k mod 7 : k..k+12} -- 13 residues, period 16 in k
print()
print("Residues mod 7 of {F_k..F_{k+12}} cycle with period 16 in k (=> resonance structure periodic):")
for k in [2,18]:
    print(f"  k={k}: {[F[k+j]%7 for j in range(13)]}")
