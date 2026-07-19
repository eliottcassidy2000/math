#!/usr/bin/env python3
"""
death-star-2026-07-19-S59d -- HYP-7905: gate-completeness lemmas (numeric
verification + per-member certificates) and the D=7 rung at N=2311.

THE THREE LEMMAS (proved on paper; verified here):
  L1 (rung floor, closed form).  p=2D-1, Q=(N+1)D-1, p nmid Q, N>3D-2.  At
     a == D*p^{-1} (mod Q), every element of F_D(N) = {1..N}\\{N-1} u {D(N-1)}
     has |v a|_Q >= D.  Key: (N+1)D == 1 => D^{-1} == N+1, and p(N+1) == 1-N,
     so |v a| <= D-1 => v == -r(N-1) (|r|<=D-1), whose only representative in
     [1,N] u {x} is v = N-1 (r=-1) -- the DELETED element.  So M >= D/Q.
  L2 (small-moduli seal, N odd).  For every modulus q' <= 2N and every a, some
     element of F has distance <= 1; value <= 1/q' <= 1/(N+1).
  L3 (packing seal).  If e := |(N-1)a|_{q'} >= c := min_{u in base}|ua|_{q'},
     then the N+1 points {ua : u=0..N} are pairwise >= c apart (every index
     difference d<=N is in base or = N-1), so (N+1)c <= q' and value <= 1/(N+1).
REDUCTION.  Any candidate with value > 1/(N+1) has q' > 2N (L2) and e < c
  (L3); then family value = min(c, |xa|)/q' with |xa| = D*e (for De <= q'/2),
  and N-point packing on the deleted prefix gives e <= q'/N.  By THM-1002 the
  maximizer lives at pair-sum moduli.  Pair sums > 2N are exactly
  {x+u : u in base} u {2x}; scanning those with the scaled e-congruence
  (N-1)a == +-e covers all their divisors (non-reduced t).  At S: solutions
  need g := gcd(N-1, S) | e and e <= S/N, so if g > S/N the e-channel at S is
  EMPTY (L3 seals S entirely).  Hence
     M(F_D(N)) = max( D/Q ,  max over S, e, solutions a of min(c, De)/S ).
  ~10^5 ops instead of ~10^10: N=2311 becomes computable.

GATES: the pruned evaluator must reproduce the S59c full-evaluator values on
every F_D(N) row with N odd, p nmid Q, N > 3D-2 (members AND non-members).
"""
from fractions import Fraction as F
from math import gcd
import sys, time
sys.path.insert(0, '04-computation')

log = lambda s="": print(s, flush=True)

def inv(a, m):
    return pow(a, -1, m)

def L1_witness_check(N, D):
    """verify Lemma 1's closed-form witness directly."""
    p = 2*D - 1
    Q = (N+1)*D - 1
    if gcd(p, Q) != 1: return None
    x = D*(N-1)
    a = (D * inv(p, Q)) % Q
    fam = [v for v in range(1, N+1) if v != N-1] + [x]
    dmin = min(min((v*a) % Q, Q - (v*a) % Q) for v in fam)
    return a, dmin  # want dmin == D exactly (p and x sit at +-D)

def M_pruned(N, D):
    """exact M(F_D(N)) under L1+L2+L3 (requires N odd, p nmid Q, N > 3D-2).
    Returns (M, q, a_or_None)."""
    p = 2*D - 1
    Q = (N+1)*D - 1
    assert N % 2 == 1 and gcd(p, Q) == 1 and N > 3*D - 2
    x = D*(N-1)
    base = [v for v in range(1, N+1) if v != N-1]
    # start at the L1 floor
    aw = (D * inv(p, Q)) % Q
    bn, bd, ba = D, Q, aw
    sums = [x + u for u in base] + [2*x]
    for S in sums:
        g = gcd(N-1, S)
        E = S // N                      # e <= S/N (deleted-prefix packing)
        if g > E: continue              # channel empty; L3 seals S
        for e in range(g, E+1, g):
            if D*e*2 > S: break         # |xa| formula needs De <= S/2 (fine here)
            # value cap at this (S,e): min(c,De)/S <= De/S; skip if can't beat best
            if D*e*bd <= bn*S: continue
            # solve (N-1) a == +-e mod S: base solution a0 mod S/g, g solutions each sign
            M1 = S // g
            a0 = (e // g) * inv((N-1)//g, M1) % M1
            for sign_a in (a0 % M1, (-a0) % M1):
                for k in range(g):
                    a = sign_a + k*M1
                    if a == 0: continue
                    # c = min over base, early break below current-best threshold
                    thr_num, thr_den = bn, bd    # need value > bn/bd
                    c = S
                    dead = False
                    for u in base:
                        r = (u*a) % S
                        if r > S - r: r = S - r
                        if r < c:
                            c = r
                            if c*bd <= bn*S:     # value <= best already
                                dead = True; break
                    if dead: continue
                    xa = (x*a) % S
                    if xa > S - xa: xa = S - xa
                    val = min(c, xa)
                    if val*bd > bn*S:
                        bn, bd, ba = val, S, a
    gg = gcd(bn, bd)
    return F(bn, bd), bd, ba

def main():
    log("== HYP-7905: gate-completeness lemmas + N=2311 (death-star-S59d) ==\n")
    # --- L1 closed-form witness across the table range ---
    log("L1 closed-form witness a = D*p^{-1} mod Q (want min distance exactly D):")
    bad = 0; n_l1 = 0
    for D in (3, 4, 5, 6, 7):
        for N in range(3*D-1, 101):
            r = L1_witness_check(N, D)
            if r is None: continue
            a, dmin = r; n_l1 += 1
            if dmin != D:
                bad += 1
                log(f"   !! L1 FAIL D={D} N={N}: dmin={dmin}")
    log(f"   L1 verified on {n_l1} (D,N) pairs, failures: {bad}\n")

    # --- L2 exhaustive verification on the known member sizes ---
    log("L2 exhaustive (every q'<=2N, every a: some element at distance <=1):")
    for (N, D) in ((31, 4), (61, 4), (91, 4)):
        x = D*(N-1); fam = [v for v in range(1, N+1) if v != N-1] + [x]
        viol = 0
        for q in range(2, 2*N+1):
            for a in range(1, q//2+1):
                if min(min((v*a) % q, q-(v*a) % q) for v in fam) > 1:
                    viol += 1
        log(f"   N={N} D={D}: violations = {viol} (want 0)")
    log("")

    # --- GATE: pruned evaluator vs the S59c table (odd N, p nmid Q, N>3D-2) ---
    log("GATE: M_pruned must reproduce the full evaluator on the S59c rows:")
    from lrc_singlefar_absorption_atlas_deathstar_S59 import M_exact_wit
    mism = 0; n_gate = 0
    t0 = time.time()
    for D in (3, 4, 5, 6):
        p = 2*D-1
        for N in range(3*D-1, 101):        # odd N in table range
            if N % 2 == 0: continue
            Q = (N+1)*D - 1
            if gcd(p, Q) != 1: continue
            x = D*(N-1)
            fam = [v for v in range(1, N+1) if v != N-1] + [x]
            Mfull, qf, af = M_exact_wit(fam)
            Mp, qp, ap = M_pruned(N, D)
            n_gate += 1
            if Mfull != Mp:
                mism += 1
                log(f"   !! MISMATCH D={D} N={N}: full={Mfull}@{qf} pruned={Mp}@{qp}")
    log(f"   gate: {n_gate} rows compared, mismatches: {mism} ({time.time()-t0:.0f}s)\n")
    if mism:
        log("   GATE FAILED -- do not trust the N=2311 result"); return

    # --- per-member e-channel certificates ---
    log("Per-member exact values via the pruned (proof-backed) evaluator:")
    for (N, D, want) in ((31,4,F(4,127)), (61,4,F(4,247)), (91,4,F(4,367)), (211,6,F(6,1271))):
        t0 = time.time()
        Mp, qp, ap = M_pruned(N, D)
        tag = "OK" if Mp == want else "?!"
        log(f"   {tag} F_{D}({N}): M = {Mp} (q={qp}) [{time.time()-t0:.1f}s] (want {want})")
    log("")

    # --- THE D=7 RUNG AT N=2311 ---
    N, D = 2311, 7
    p, Q, x = 13, 2312*7 - 1, 7*2310
    log(f"== D=7 rung at N=2311: F_7(2311) = {{1..2309, 2311, {x}}} ==")
    log(f"   p=13, Q={Q} (gcd(13,Q)=1: {gcd(13,Q)==1}); window (1/2312, 2/4623), "
        f"width {float(F(2,4623)-F(1,2312)):.2e}")
    a1, d1 = L1_witness_check(N, D)
    log(f"   L1 witness: a = {a1}, min distance = {d1} (want 7) -> M >= 7/{Q}")
    t0 = time.time()
    Mp, qp, ap = M_pruned(N, D)
    log(f"   EXACT (pruned, proof-backed): M(F_7(2311)) = {Mp} at q={qp}, a={ap} "
        f"[{time.time()-t0:.0f}s]")
    rung = F(7, Q)
    inwin = F(1, 2312) < Mp < F(2, 4623)
    log(f"   rung 7/{Q}: {'ATTAINED' if Mp == rung else 'NOT attained'}; "
        f"in-window: {inwin}")
    log(f"   VERDICT: the D=7 gate at N=2311 is "
        f"{'OPEN -- third out-of-sample tower confirmation' if (Mp == rung and inwin) else 'CLOSED/OTHER -- see value above'}")

if __name__ == "__main__":
    main()
