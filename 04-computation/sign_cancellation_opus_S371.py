# opus-2026-07-17-S371 -- HYP-7540: WHERE DOES THE CANCELLATION LIVE?
#
# hhat(n) = sin(pi n/7)/(pi n), and sin(pi n/7) depends ONLY on n mod 14.
# Write c(r) = sin(pi r/7).  Then
#     prod_i hhat(n_i) = [prod_i c(n_i mod 14)] / (pi^k * prod_i n_i)
# so, splitting the lattice by residue vector r in (Z/14)^k,
#     delta(S) = (1/pi^k) * sum_r [prod c(r_i)] * T(r),
#          T(r) = sum over {n in Lam : n = r mod 14} of 1/prod n_i.
# ALL SIGNS SIT IN TWO PLACES: the finite coefficient table prod c(r_i)
# (14^k entries, many zero since c(r)=0 iff 7|r), and the 1/prod n_i inside
# each coset.  Measuring the two separately says WHICH one the eventual
# bound must exploit -- that is the whole point of this run.
#   FULL-ABS   = sum |prod hhat|          (S370's provable-but-loose bound)
#   COSET-ABS  = sum_r |prod c| * |T(r)|  (keeps cancellation INSIDE cosets)
#   SIGNED     = |delta|                  (keeps everything)
from math import sin, pi
import itertools
from collections import defaultdict
LAM=1.0/14
def c(r): return sin(pi*(r%14)/7)
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)

def analyse(A,N):
    k=len(A)
    signed=0.0; fullabs=0.0
    T=defaultdict(float); coef={}
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=k):
        if any(ni%7==0 for ni in n): continue
        if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
        p=1.0; prod_n=1.0; cc=1.0
        for ni in n:
            p*=hhat(ni); prod_n*=ni; cc*=c(ni)
        signed+=p; fullabs+=abs(p)
        r=tuple(ni%14 for ni in n)
        T[r]+=1.0/prod_n
        coef[r]=cc
    cosetabs=sum(abs(coef[r])*abs(T[r]) for r in T)/pi**k
    # rebuild check
    rebuilt=sum(coef[r]*T[r] for r in T)/pi**k
    return signed, rebuilt, cosetabs, fullabs, len(T)

print("(1) THE COSET DECOMPOSITION, AND WHERE CANCELLATION LIVES  (k=3)")
print("    family          signed    rebuilt    COSET-ABS   FULL-ABS   full/signed  coset/signed")
for A in [(2,3,5),(11,13,17),(31,37,41),(6,10,15),(5,7,11),(101,103,107)]:
    sg,rb,ca,fa,nc=analyse(A,45)
    print(f"    {str(A):16s} {sg:+.6f} {rb:+.6f}   {ca:.6f}   {fa:.6f}   {fa/abs(sg):8.2f}   {ca/abs(sg):8.2f}")
print()
print("    'rebuilt' must equal 'signed' -- that validates the decomposition.")
print("    If COSET-ABS is close to FULL-ABS, the cancellation is ACROSS cosets")
print("    (in the finite sign table prod c(r_i)) -- which is the exploitable case,")
print("    because that table is FINITE and explicit.")
print("    If COSET-ABS is close to SIGNED, the cancellation is INSIDE cosets")
print("    (in 1/prod n_i), which needs analytic, not combinatorial, control.")
