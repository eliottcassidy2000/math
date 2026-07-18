# opus-2026-07-17-S370 -- HYP-7530: CAN THE MULTIPLICATIVE-MINIMUM BOUND BE PROVED?
#
# THM-1080 left 1/(pi^k m7) as "the right scale, ratios 0.45-1.93" -- NOT a
# proved upper bound.  The obvious way to prove one is absolute values:
#       |delta(S)| <= sum over full-support 7-free n in Lam of prod 1/(pi|n_i|)
# BEFORE claiming that, check it CONVERGES.  Heuristic says it may not: a
# rank-(k-1) lattice has ~H^(k-1) vectors of height <= H, each contributing
# ~H^-k, giving sum_H H^(k-1) * H^-k ~ sum_H 1/H -- LOGARITHMIC DIVERGENCE.
# If so the naive bound is unprovable and the truth needs CANCELLATION.
from math import sin, pi
import itertools
LAM=1.0/14
def hhat(n):
    if n==0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)

def sums_to_N(A, Nmax):
    """(signed partial sum, absolute partial sum) at each cutoff."""
    k=len(A); out=[]
    for N in Nmax:
        sg=0.0; ab=0.0
        for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=k):
            if any(ni%7==0 for ni in n): continue
            if sum(ni*ai for ni,ai in zip(n,A))!=0: continue
            p=1.0; q=1.0
            for ni in n:
                p*=hhat(ni); q*=1.0/(pi*abs(ni))
            sg+=p; ab+=q
        out.append((N,sg,ab))
    return out

print("(1) DOES THE ABSOLUTE-VALUE BOUND CONVERGE?  k=3")
print("    family        N      signed sum (=delta)   ABSOLUTE sum")
for A in [(2,3,5),(11,13,17),(31,37,41)]:
    for (N,sg,ab) in sums_to_N(A,[10,20,40,60,80]):
        print(f"    {str(A):14s} {N:3d}    {sg:+.6f}            {ab:.6f}")
    print()
print("    If the ABSOLUTE column keeps climbing while the signed column")
print("    settles, the naive bound DIVERGES and cannot be proved this way.")
print()
print("(2) GROWTH RATE OF THE ABSOLUTE SUM (log-divergence test)")
print("    if abs-sum ~ c*log(N), successive differences at doubling N are ~constant")
for A in [(2,3,5),(11,13,17)]:
    rows=sums_to_N(A,[10,20,40,80])
    print(f"    {A}:")
    for i in range(1,len(rows)):
        print(f"      N {rows[i-1][0]:3d}->{rows[i][0]:3d}   abs-sum {rows[i-1][2]:.6f} -> {rows[i][2]:.6f}"
              f"   increment {rows[i][2]-rows[i-1][2]:+.6f}")
