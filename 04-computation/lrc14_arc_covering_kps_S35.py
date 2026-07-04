"""kps-S35: (1) the ARC-COVERING LEMMA rate: f_q = P(13 random danger-sets cover Z/q).
mac-mini verified f_q<1; I want the RATE (toward a rigorous 2nd-moment proof f_q=O(1/q)).
(2) the RESONANCE DICHOTOMY: does `no small pairwise resonance m_i r_i + m_j r_j == 0 mod q
(|m|<=7)` IMPLY `witness exists (danger sets don't cover)`? (contrapositive of mac-mini's
`cover => resonance`). If robust, LRC(14)-compressed reduces to the resonant families."""
import numpy as np
rng = np.random.default_rng(0)

def band_k(q): return (q-1)//14

def covers(R, q):
    """R: (S,13) residues in [1,q). returns bool (S,) : do the 13 danger-sets cover Z/q?"""
    k = band_k(q)
    a = np.arange(1, q)                     # a=0 always danger; check a in [1,q)
    # ra mod q : (S,13,q-1)
    ra = (R[:,:,None] * a[None,None,:]) % q
    m = np.minimum(ra, q-ra)
    danger = m <= k                          # (S,13,q-1)
    dangered = danger.any(axis=1)            # (S,q-1) : a is danger for some runner
    return dangered.all(axis=1)              # cover iff every a dangered

def has_pair_resonance(r, q, H=7):
    """r: length-13 residues. True if some i<j and 1<=|mi|,|mj|<=H with mi ri+mj rj==0 mod q."""
    for i in range(13):
        for j in range(i+1,13):
            for mi in range(1,H+1):
                for mj in range(-H,H+1):
                    if mj==0: continue
                    if (mi*r[i]+mj*r[j])%q==0: return True
    return False

print("=== (1) arc-covering density f_q (random nonzero residues) ===")
print(f"{'q':>4} {'f_q':>10} {'q*f_q':>8} {'main (6/7)^13':>13} {'E[#safe]/q est':>14}")
for q in [15,17,19,23,29,37,43,53,67,83,101,127]:
    S = 20000
    R = rng.integers(1, q, size=(S,13))
    cov = covers(R, q)
    fq = cov.mean()
    print(f"{q:>4} {fq:>10.5f} {q*fq:>8.2f} {(6/7)**13:>13.5f}")
print()
print("=== (2) resonance dichotomy on random residues ===")
print(f"{'q':>4} {'#cover':>7} {'#no_reson':>10} {'cover&no_reson':>14} {'no_reson=>witness?':>18}")
for q in [17,19,23,29,37,43,53,67,83,101,127]:
    S = 4000
    R = rng.integers(1, q, size=(S,13))
    cov = covers(R, q)
    ncover=0; nnores=0; viol=0
    for s in range(S):
        c = cov[s]
        nr = not has_pair_resonance(R[s], q, H=7)
        if c: ncover+=1
        if nr: nnores+=1
        if nr and c:                          # no resonance BUT covers = VIOLATION of the lemma
            viol+=1
    ok = "OK (0 viol)" if viol==0 else f"{viol} VIOLATIONS"
    print(f"{q:>4} {ncover:>7} {nnores:>10} {viol:>14} {ok:>18}")
print()
print("READING (1): if q*f_q stays bounded => f_q=O(1/q) (2nd-moment-provable arc-covering lemma).")
print("READING (2): 0 violations => THEOREM `no small pairwise resonance mod q => witness exists`;")
print("  then LRC(14)-compressed reduces to families WITH a small pairwise resonance at q*.")
