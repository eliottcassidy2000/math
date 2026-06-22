r"""
LRC analogue of H = product H(atoms) (kps-S31x): does safe-measure FACTOR over scale-separated
clusters? meas(safe(S1 u S2)) =?= meas(safe(S1)) * meas(safe(S2)) when S2 >> S1 (independent orbits).
If yes, the LRC reduces MULTIPLICATIVELY to irreducible single-scale ATOMS (the bounded cores), exactly
as the tournament Redei count is multiplicative over strong components.
"""
from fractions import Fraction as F
def safe_measure(S, thr=F(1,14)):
    S=sorted(set(abs(s) for s in S if s!=0)); bps={F(0),F(1)}
    for s in S:
        for m in range(0, s+1):
            for sgn in (-1,1):
                x=F(m,s)+sgn*thr/s
                if 0<=x<=1: bps.add(x)
    B=sorted(bps); tot=F(0)
    for lo,hi in zip(B,B[1:]):
        mid=(lo+hi)/2
        if all(min((s*mid)%1,1-(s*mid)%1)>=thr for s in S): tot+=hi-lo
    return tot
print("MULTIPLICATIVITY TEST: meas(safe(S1 u S2)) vs meas(safe(S1))*meas(safe(S2)), S2 scale-separated:")
tests=[([1,2,3],[1009,2011,3013]),
       ([1,2,3,4,5],[997,1999]),
       ([2,3,5],[1013,3017,5021]),
       ([1,2,3],[101,103,107])]   # less separated (resonance check)
for S1,S2 in tests:
    m1=float(safe_measure(S1)); m2=float(safe_measure(S2)); m12=float(safe_measure(S1+S2))
    sep=min(S2)/max(S1)
    print(f"  S1={S1} S2={S2} (sep~{sep:.0f}x):")
    print(f"     meas(safe S1)={m1:.5f}  meas(safe S2)={m2:.5f}  product={m1*m2:.5f}")
    print(f"     meas(safe S1uS2)={m12:.5f}  ratio actual/product={m12/(m1*m2):.4f} (->1 = INDEPENDENT/multiplicative)")
print("\n=> if ratio -> 1 for scale-separated, safe-measure is MULTIPLICATIVE over independent clusters")
print("   = the LRC's H=prod(atoms). The IRREDUCIBLE ATOMS = single-scale clusters (bounded cores); the")
print("   hardest atom = the AP (apex-7), like the tournament's forbidden K_3 atom (H=7).")
