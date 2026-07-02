#!/usr/bin/env python3
"""
LOCALIZING the moment LP with the Gamma_0(N) congruence constraint -- does min m_0 become >0?
mac-mini-2026-07-01-S92, HYP-3824.  (Direct descendant of HYP-3822: the GLOBAL moment LP gave min m_0=0.)

GLOBAL LP (S91): min m_0 s.t. sum m_k=1, sum k m_k=A, sum k(k-1)m_k=Phi  ==>  min m_0 = 0 (A>1). USELESS.

LOCALIZE: partition the circle [0,1) into Q arcs I_c=[c/Q,(c+1)/Q).  Impose the PER-ARC coverage mass
  A_c = int_{I_c} C(t) dt   (the LOCAL 1st moment = the arithmetic data; Gamma_0(N)-aligned when Q | N-structure).
Within arc I_c, min |{C=0} cap I_c| given int C = A_c and C>=0 integer  =  max(0, |I_c| - A_c)  (local union).
  => LOCALIZED LP:  min m_0(Q) = sum_c max(0, 1/Q - A_c)      [a VALID lower bound on |L|, ->|L| as Q->inf]
SHARPER (local 2nd moment too):  min m_0(Q) = sum_c max(0, 1/Q - A_c + Phi_c/Cmax_c).

CLAIM tested: the GLOBAL LP (Q=1) gives 0, but LOCALIZING (Q>1, esp. Gamma_0-aligned Q=183=Phi_6(14))
makes min m_0 > 0 -- the floor is a LOCAL fact, carried by UNDER-COVERED arcs (the atom's neighborhood).
FLOOR question: over MANY covering sets, is there ALWAYS an under-covered arc? (=> a uniform floor).
"""
import numpy as np

N_ATOM = 183  # Phi_6(14)=14^2-14+1 ; the binding modulus (atom t*=14/183)

def coverage(S, r, G):
    t = (np.arange(G) + 0.5) / G
    C = np.zeros(G, dtype=np.int32)
    for v in S:
        x = (v * t) % 1.0
        C += (np.minimum(x, 1.0 - x) < r).astype(np.int32)
    return C

def localized_bound(C, Q):
    """sum_c max(0, |I_c| - A_c) [union] and the 2nd-moment-sharpened version, over Q arcs."""
    G = len(C); assert G % Q == 0, (G, Q)
    cell = G // Q
    Cc = C.reshape(Q, cell)
    A_c   = Cc.mean(axis=1) / Q            # int_{I_c} C dt  (|I_c|=1/Q, mean over cell / Q)
    Phi_c = (Cc*(Cc-1)).mean(axis=1) / Q
    Cmax_c = Cc.max(axis=1); Cmax_c = np.where(Cmax_c==0, 1, Cmax_c)
    invQ = 1.0/Q
    union = np.maximum(0.0, invQ - A_c).sum()
    pot   = np.maximum(0.0, invQ - A_c + Phi_c/Cmax_c).sum()
    n_under = int((A_c < invQ).sum())
    return union, pot, n_under

if __name__ == "__main__":
    G = 183 * 12000            # divisible by 183, 14, 6, 61, 2, 3, ... (=2196000)
    print("#"*74)
    print("# LOCALIZED moment LP: does per-arc (Gamma_0-aligned) localization make min m_0 > 0?")
    print(f"#   grid G={G}; N_atom=Phi_6(14)={N_ATOM}")
    print("#"*74)

    S_c2 = list(range(1,13)) + [182]     # {1..12,182}, |L|(r=1/14)=0.0239 (S91)
    r = 1.0/14
    C = coverage(S_c2, r, G)
    true_L = (C==0).mean()
    print(f"\nS={S_c2}  r=1/14   ACTUAL |L| = {true_L:.5f}")
    print(f"\n  {'Q':>6} {'aligned?':>10} {'union bound':>12} {'potential':>11} {'#under-cov':>11}")
    Qs = [1, 2, 3, 6, 7, 14, 61, 183, 366, 549, 2*183, 6*183, 61*183, 183*100, 183*1000]
    for Q in Qs:
        if G % Q: continue
        u, p, nu = localized_bound(C, Q)
        aligned = "Gamma_0" if (Q % 183 == 0 or 183 % Q == 0 or Q in (7,14,61)) else "generic"
        star = " <-atom" if Q==183 else ""
        print(f"  {Q:>6} {aligned:>10} {u:12.5f} {p:11.5f} {nu:11d}{star}")

    print("\n  => Q=1 (GLOBAL) gives 0 (matches S91 HYP-3822); Q>=? goes POSITIVE; ->|L| as Q->inf.")

    # ---- FLOOR question: over a sample of covering sets, is there ALWAYS an under-covered arc? ----
    print("\n" + "#"*74)
    print("# FLOOR: over many LRC14 speed sets, the localized bound (Q=183, atom-aligned).")
    print("#   min over sets of the localized bound = a proxy floor (holds for all sampled sets).")
    print("#"*74)
    rng = np.random.default_rng(12345)
    sets = {
        "{1..13}":            list(range(1,14)),
        "{1..12,182}":        list(range(1,13))+[182],
        "{1..12,182} (b)":    [1,2,3,4,5,6,7,8,9,10,11,12,183],
        "AP step2 {2..26}":   list(range(2,27,2)),
        "AP step3 {3..39}":   list(range(3,40,3)),
    }
    # a few random 13-subsets of [1,200]
    for i in range(5):
        sets[f"rand{i}"] = sorted(rng.choice(range(1,200), size=13, replace=False).tolist())
    Gs = 183*4000
    r = 1.0/14
    print(f"\n  grid={Gs}, r=1/14, Q=183 (atom-aligned)")
    print(f"  {'set':>20} {'|L|':>9} {'loc.union(Q=183)':>17} {'loc.pot':>9} {'#under':>7}")
    minbound = 1e9
    for name, S in sets.items():
        C = coverage(S, r, Gs)
        L = (C==0).mean()
        u,p,nu = localized_bound(C, 183)
        minbound = min(minbound, p)
        print(f"  {name:>20} {L:9.5f} {u:17.5f} {p:9.5f} {nu:7d}")
    print(f"\n  min localized-potential bound over sampled sets (Q=183) = {minbound:.5f}")
    print(f"  (kps floor target: inf meas >= 1/36 = {1/36:.5f}; compare.)")
