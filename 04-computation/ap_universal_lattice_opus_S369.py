# opus-2026-07-17-S369 -- MY DIRECTION C: THE AP RESONANCE LATTICE DECOUPLES.
#
# For an AP a_i = a + i*d (i = 0..k-1) with gcd(a,d) = 1:
#     sum n_i (a + i d) = a*N0 + d*N1,   N0 = sum n_i,  N1 = sum i*n_i.
# This vanishes iff a*N0 = -d*N1; since gcd(a,d)=1 we get d | N0, so
#     N0 = d*s  and  N1 = -a*s   for some s in Z.
# THE s = 0 PART is  {n : sum n_i = 0, sum i n_i = 0} -- a rank-(k-2) lattice
# that DOES NOT DEPEND ON a OR d AT ALL.  The s != 0 part needs |N0| >= d,
# hence large coordinates, hence tiny hhat products once d is large.
# PREDICTION: as d grows, delta(S) for an AP converges to a UNIVERSAL constant
# depending only on k -- the AP has a LIMITING deviation profile.  If so, the
# extremal (additively richest) corner is decidable once and for all.
from math import sin, pi, gcd
import itertools
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)
def delta(A,N):
    t=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=len(A)):
        if sum(ni*ai for ni,ai in zip(n,A))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            t+=p
    return t
def delta_universal(k,N):
    """the s=0 part: sum n_i = 0 AND sum i*n_i = 0, independent of a,d."""
    t=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=k):
        if sum(n)==0 and sum(i*ni for i,ni in enumerate(n))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            t+=p
    return t

print("(5) DOES THE AP DEVIATION CONVERGE TO A UNIVERSAL CONSTANT?")
print("    k=3 APs (a, a+d, a+2d), varying a and d:")
U3=delta_universal(3,40)
print(f"    universal s=0 value (depends only on k=3): {U3:+.6f}")
print("    a    d     AP                 delta        vs universal")
for (a,d) in [(1,1),(2,1),(1,2),(1,3),(2,3),(1,5),(3,5),(1,8),(5,8),(1,13),(4,13),(1,21),(1,34)]:
    if gcd(a,d)!=1: continue
    A=(a,a+d,a+2*d)
    dd=delta(A,40)
    print(f"    {a:3d}  {d:3d}   {str(A):18s} {dd:+.6f}     {dd-U3:+.6f}")
print()
print("    k=4 APs:")
U4=delta_universal(4,22)
print(f"    universal s=0 value (k=4): {U4:+.6f}")
for (a,d) in [(1,1),(1,2),(1,3),(1,5),(1,8),(1,13),(2,5),(3,7)]:
    A=tuple(a+i*d for i in range(4))
    dd=delta(A,22)
    print(f"    {a:3d}  {d:3d}   {str(A):18s} {dd:+.6f}     {dd-U4:+.6f}")
print()
print("  => if the residual column tends to 0 with d, the AP corner has a")
print("     LIMITING profile computable once per k, independent of the AP.")
