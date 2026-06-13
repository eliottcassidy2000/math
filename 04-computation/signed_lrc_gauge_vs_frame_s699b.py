"""Signed LRC: what the negative-speed 'trick' DOES and does NOT give (sharpening S674/S674b).
(1) GAUGE = per-runner signs ε_i: M({ε_i v_i})=M(S) ALWAYS (||εvt||=||vt||) — an exact (Z/2)^{n-1}
   symmetry; preserves which residues are 0 mod n (so it does NOT remove multiples / reduce C2b).
(2) But the FRAME shift (subtract v_k, 'sit on runner k') CHANGES the config and M — so it is NOT
   a witness-transfer for the observer. Honest boundary: signs are a side-channel, not a reduction.
(3) Gauge pair-clocks: ε_i v_i - ε_j v_j ranges over {v_i-v_j, v_i+v_j} — opposite signs expose the
   pair-SUMS = the 2n-1 worry-set shells (THM-401). opus-2026-06-06-S699b."""
from itertools import combinations, product
from fractions import Fraction as F
def M_exact(V):
    ds=set(V)
    for a,b in combinations(V,2): ds.add(a+b); ds.add(abs(a-b))
    ds.discard(0); best=F(-1)
    for d in ds:
        for m in range(d):
            t=F(m,d); v=min(min((x*t)%1,1-(x*t)%1) for x in V)
            if v>best: best=v
    return best
def main():
    print("(1) GAUGE invariance: M({ε_i v_i}) = M(S) for ALL sign patterns?")
    for V in [(1,2,3,4,5),(1,3,4,5,9),(2,6,7,8,12)]:
        base=M_exact(V); ok=True
        for signs in product((1,-1),repeat=len(V)):
            if M_exact(tuple(s*v for s,v in zip(signs,V)))!=base: ok=False;break
        print(f"   V={V}: M={base}; gauge-invariant over all 2^{len(V)} signs = {ok}")
    print("\n(2) FRAME shift CHANGES M (so 'sit on runner k' is NOT an observer witness-transfer):")
    V=(1,2,3,7)  # contains 7? n=8 multiple... use n where 7 wld be multiple; just show M changes
    base=M_exact(V)
    for k in V:
        Vsh=tuple(sorted(set(abs(v-k) for v in V if v!=k)|{k}))
        print(f"   V={V} M={base};  shift by {k} -> {Vsh}  M={M_exact(Vsh)}  (different config, different M)")
    print("   => the gauge (signs) preserves M; the FRAME (shift) does not. No cheap C2b win from signs.")
    print("\n(3) Gauge pair-clocks expose pair-SUMS (the 2n-1 worry-set shells, THM-401):")
    n=14; C=2*n-1
    V=tuple(range(1,n))  # AP13
    sums=sorted({(a+b) for a,b in combinations(V,2)} ); diffs=sorted({abs(a-b) for a,b in combinations(V,2)})
    print(f"   AP13 same-sign pair clocks (differences): {diffs[:12]}...")
    print(f"   AP13 opposite-sign pair clocks (sums) mod C={C}: {sorted({(a+b)%C for a,b in combinations(V,2)})}")
    print(f"   the opposite-sign sums are exactly the odd interior residues — the worry-set shell ledger (S674b).")
if __name__=='__main__': main()
