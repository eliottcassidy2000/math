#!/usr/bin/env python3
"""
lrc14_riesz_verify_macmini_0615s5.py  (mac-mini-2026-06-15-S5; OPEN-Q-104 VERIFICATION)

CRITICAL CHECK of the candidate Riesz certificate ∫M·R/∫R < 1 ⟹ loose.
R(τ) = ∏_{v∈S}(1 - cos 2πvτ) ≥ 0 pointwise. Compute ∫M·R and ∫R by DIRECT grid
integration (exact, no Fourier truncation) for:
  - the extremizer core (loose, L≈0.0056)  -> must give ratio < 1 IF the certificate is real
  - the TIGHT config 14*{1..13} (L=0, covers a.e.) -> MUST give ratio ≥ 1 (else bug)
  - the dilated core {1..12}∪{14} (loose, L≈0.024)
Scan amplitude a in (1-... and a few a<1). The certificate is real ONLY IF it cleanly
separates loose (ratio<1) from tight (ratio≥1). Also report ∫R (must be > 0).
"""
import sys, math
from math import cos, pi

sys.stdout.reconfigure(line_buffering=True)

def MR_direct(S, a, Q):
    """∫M·R and ∫R by grid, R=∏(1+a cos 2π v τ), M=#{v: ||vτ||≤1/14}."""
    rad = Q//14
    sMR=0.0; sR=0.0
    for k in range(Q):
        tau=k/Q
        R=1.0
        for v in S:
            R*=(1+a*cos(2*pi*v*tau))
        if R==0.0:
            M=0  # still need M for completeness but R=0 contributes 0
        M=sum(1 for v in S if ((v*k)%Q)<=rad or ((v*k)%Q)>=Q-rad)
        sMR+=M*R; sR+=R
    return sMR/Q, sR/Q

configs = {
  "extremizer {1..13}\\{6}∪{56} (LOOSE L≈.006)": sorted(set([x for x in range(1,14) if x!=6]+[56])),
  "dilated {1..12}∪{14}        (LOOSE L≈.024)": sorted(set(list(range(1,13))+[14])),
  "TIGHT 14*{1..13}            (L=0, covers) ": [14*i for i in range(1,14)],
}

print("="*78)
print("VERIFY the Riesz certificate: R=∏(1+a cos 2πvτ)≥0, ratio=∫MR/∫R; loose<1, tight≥1?")
print("="*78)
Q=30030  # = 2*3*5*7*11*13, divisible by 14, fine grid
for name,S in configs.items():
    if len(S)!=13: print(f"  {name}: skip"); continue
    print(f"\n  {name}:  S={S}")
    for a in [-1.0, -0.99, -0.95, -0.9]:
        mr, r = MR_direct(S, a, Q)
        ratio = mr/r if r>1e-14 else float('inf')
        verdict = "ratio<1 ⟹ LOOSE-certified" if ratio<1-1e-9 else ("ratio≥1 (consistent w/ tight)" if ratio>=1-1e-9 else "?")
        print(f"     a={a:+.2f}: ∫MR={mr:.6f}, ∫R={r:.6f}, ratio={ratio:.5f}   {verdict}")
print("\n  DECISIVE: the certificate is REAL iff loose cores give ratio<1 AND the tight config")
print("  gives ratio≥1 (it must, since M≥1 a.e. there ⟹ ∫MR≥∫R for R≥0). A false <1 at tight")
print("  would expose a bug. Cross-check: ∫R must be > 0 (R≥0, not ≡0).")
print("\nDONE.")
