import sys, itertools, cmath, math
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8")
# The proof's claim: chat_T(n_j) = delta_{n_j,0}(1-|T|/7) + [n_j!=0] shat_T(n_j).
# But shat_T(n_j) = sum_{j in T} (single-sector Fourier piece). Crucially shat_T
# depends on T only through which sectors are in T -- it's a SUM over sectors in T.
# So prod_j shat_T(n_j) is NOT of the form "T contains a fixed set U".
# The proof asserts a factorization "T-sum factors through C(U)=sum_{T superset U}(-1)^|T|".
# This is only valid if shat_T(n) = sum_{s in T} f(s,n), and the product expands to
# pick one sector PER coordinate. Let me expand explicitly and test the proof's C(U).
#
# shat_T(n) = -1/(2pi i n) sum_{s in T}(e^{-2pi i n(s+1)/7}-e^{-2pi i n s/7})
#           = sum_{s in T} g(s,n),  g(s,n)= -(e^{-2pi i n(s+1)/7}-e^{-2pi i n s/7})/(2pi i n)
def g(s,n):
    a=s/7.0; b=(s+1)/7.0
    return -(cmath.exp(-2j*math.pi*n*b)-cmath.exp(-2j*math.pi*n*a))/(2j*math.pi*n)
def shat(T,n):
    return sum(g(s,n) for s in T)
SUBS=[frozenset(c) for r in range(7) for c in itertools.combinations(range(1,7),r)]

# K for an all-active relation (no zero coords): expand
# K = sum_T (-1)^|T| prod_{j} shat_T(n_j)
#   = sum_T (-1)^|T| prod_j (sum_{s_j in T} g(s_j,n_j))
#   = sum_{s_1,...,s_d in {1..6}} prod_j g(s_j,n_j) * [sum_{T superset {s_1..s_d}} (-1)^|T|]
#   = sum_{(s_j)} prod_j g(s_j,n_j) * C(U={s_1,...,s_d})
# where C(U)=sum_{T superset U}(-1)^|T| = 0 unless |U|=6.
# So for d active coords, U = set of chosen sectors has |U|<=d. C(U)=0 unless |U|=6,
# which requires d>=6 AND the s_j to hit all 6 sectors. => floor PROVED for all-active.
# VERIFY this expansion equals direct K:
def K_direct(nv):
    tot=0j
    for T in SUBS:
        p=1+0j
        for ne in nv:
            if ne==0: p*=complex(1-len(T)/7,0)
            elif ne%7==0: p*=0j
            else: p*=shat(T,ne)
        tot+=(-1)**len(T)*p
    return tot
def K_expand(active):  # all coords active
    d=len(active)
    tot=0j
    for sv in itertools.product(range(1,7),repeat=d):
        U=frozenset(sv)
        C=sum((-1)**len(T) for T in SUBS if U<=T)
        if C==0: continue
        p=1+0j
        for sj,nj in zip(sv,active): p*=g(sj,nj)
        tot+=C*p
    return tot
print("All-active relations: K_direct vs K_expand, and C(U) floor:")
for active in [[1,-1],[1,1,-1],[2,3,-1,-1],[1,1,1,1,-1,-1]]:
    kd=K_direct(active); ke=K_expand(active)
    print(f"  active={active}: K_direct={kd.real:+.3e}  K_expand={ke.real:+.3e}  match={abs(kd-ke)<1e-9}")
print()
# NOW the issue: a real LATTICE relation has ZERO coords. The proof says drop them.
# But when n_j=0, chat_T(0)=(1-|T|/7), NOT shat. Does the floor still hold?
# Expand a relation with some zero coords. The zero-coord factor (1-|T|/7) does NOT
# vanish-sum. Let's see if K=0 for a relation with support<6 BUT with zero padding.
print("Relations WITH zero coords (real lattice), support<6:")
for nv in [(1,-1,0,0,0,0,0),(1,1,-1,0,0,0,0),(0,1,2,-1,-1,0,0)]:
    sup=sum(1 for x in nv if x!=0 and x%7!=0)
    k=K_direct(nv)
    print(f"  nv={nv} support={sup}: K={k.real:+.6e}")
