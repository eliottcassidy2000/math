"""
monad-explorer-2026-06-06-S709 — SIGNED -> UNSIGNED LRC structural transfer.

Dispatched angle: how do signed-LRC tight/shell configs transfer to the UNSIGNED LRC?

Recall (T1, gauge invariance): M({eps_i v_i}) = M({v_i}). The signed content lives in
the runner-runner pair clocks eps_i v_i - eps_j v_j, which become SUMS v_i+v_j on
bichromatic (cut) edges. A signed zero-clock <=> shell-partner v_a+v_b ≡ 0 (mod C),
C = 2n-1 (T3).

This script probes a concrete UNSIGNED consequence of carrying a shell-partner:

 (A) MIRROR-REDUNDANCY LEMMA: v_a+v_b ≡ 0 (mod C) => ||v_a t|| = ||v_b t|| for all t in
     the grid G_C = {k/C}. Then on G_C runner b is redundant (min unchanged if dropped).

 (B) On WHICH grid is M(tight) actually attained? (denominator of argmax t*.)

 (C) V*-> reduced transfer: drop one shell-partner from a shell-partner-carrying tight
     config; does grid-M survive, and is the tightness "imported" from a smaller system?

Builds on signed_lrc_worryset_signature_s699d.py (M_exact, shell_partners), THM-401,
MISTAKE-056 (first split at n=8, NOT n=14).
"""
from itertools import combinations, product
from fractions import Fraction as F
from math import gcd


def nrm(x):
    """circular distance ||x|| for Fraction x."""
    r = x % 1
    return min(r, 1 - r)


def pinch_denoms(V):
    ds = set(V)
    for a, b in combinations(V, 2):
        ds.add(a + b)
        ds.add(abs(a - b))
    ds.discard(0)
    return ds


def M_exact_arg(V):
    """Return (M, set of argmax Fractions t in [0,1))."""
    best = F(-1)
    args = []
    for d in pinch_denoms(V):
        for m in range(d):
            t = F(m, d)
            v = min(nrm(x * t) for x in V)
            if v > best:
                best = v
                args = [t]
            elif v == best:
                args.append(t)
    return best, set(args)


def shell_partners(V, C):
    return [(V[i], V[j]) for i in range(len(V)) for j in range(i + 1, len(V))
            if (V[i] + V[j]) % C == 0]


def M_grid(V, C):
    """max over t = k/C of min_i ||v_i t||."""
    best = F(-1)
    args = []
    for k in range(C):
        t = F(k, C)
        v = min(nrm(x * t) for x in V)
        if v > best:
            best = v
            args = [k]
        elif v == best:
            args.append(k)
    return best, args


def gcd_list(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return g


# -------------------------------------------------------------------------
print("=" * 78)
print("(A) MIRROR-REDUNDANCY LEMMA  (v_a+v_b≡0 mod C  =>  ||v_a t||=||v_b t|| on G_C)")
print("=" * 78)
# exhaustive check across many C and all mirror pairs
viol = 0
checked = 0
for C in range(3, 60):
    for a in range(1, C):
        b = (C - a) % C
        if b == 0 or b == a:
            continue
        for k in range(C):
            t = F(k, C)
            if nrm(a * t) != nrm(b * t):
                viol += 1
            checked += 1
print(f"  pairs (a, C-a) over C=3..59, all k: {checked} checks, {viol} violations")
print(f"  => mirror-redundancy {'HOLDS' if viol == 0 else 'FAILS'} (trivial: b≡-a => bt≡-at mod 1, ||.|| even)")

# -------------------------------------------------------------------------
print()
print("=" * 78)
print("(B) ON WHICH GRID IS M(tight) ATTAINED?  denominator of argmax t* (reduced)")
print("=" * 78)
for n in range(4, 9):
    C = 2 * n - 1
    thr = F(1, n)
    B = 2 * n
    tight = []
    for V in combinations(range(1, B + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        Mv, _ = M_exact_arg(V)
        if Mv == thr:
            tight.append(V)
    print(f" n={n} (C={C}, thr=1/{n}): {len(tight)} tight configs in [1,{B}]")
    for V in tight:
        Mv, args = M_exact_arg(V)
        denoms = sorted({t.denominator for t in args})
        on_gridC = any(t.denominator == C or C % t.denominator == 0 for t in args)
        sp = shell_partners(V, C)
        print(f"    {V}: argmax denoms={denoms}  in_G_C={on_gridC}  shell={sp if sp else '-'}")

# -------------------------------------------------------------------------
print()
print("=" * 78)
print("(C) V*-> REDUCED TRANSFER (drop one shell-partner; does grid-M survive?)")
print("=" * 78)
# n=8 corrected first-split configs (MISTAKE-056) + n=14 V*
cases = {
    8: [(1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)],
    14: [tuple(list(range(1, 12)) + [13, 24])],  # V* = AP with 12->24
}
for n, configs in cases.items():
    C = 2 * n - 1
    thr = F(1, n)
    print(f" n={n} (C={C}, thr=1/{n}):")
    for V in configs:
        sp = shell_partners(V, C)
        Mfull, argf = M_exact_arg(V)
        Mg, _ = M_grid(V, C)
        print(f"   V={V}")
        print(f"     M_full={Mfull} (=1/{n}? {Mfull==thr}); M_grid_C={Mg}; shell={sp}")
        for (a, b) in sp:
            R = tuple(x for x in V if x != b)
            MgR, _ = M_grid(R, C)
            MfR, _ = M_exact_arg(R)
            gR = gcd_list(R)
            print(f"     drop {b} (mirror of {a}): R={R}")
            print(f"        M_grid_C(R)={MgR} (== M_grid_C(V)? {MgR==Mg});  M_full(R)={MfR}; gcd(R)={gR}")
