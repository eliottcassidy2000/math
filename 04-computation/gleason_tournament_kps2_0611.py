#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-481 lab: both Gleason Type II generators are tournament-generated.
kind-pasteur-2026-06-11-S2.  Settles HYP-2409.4 (the Golay branch); tests
HYP-2410 (eQR unification, p ≡ 7 mod 8) and HYP-2411 (order-32 trichotomy).

Constructions (THM-480's TOURNAMENT GAUGE, no renormalization):
  * Paley tournament T_p (p ≡ 3 mod 4): i -> j iff chi(j - i) = +1.
  * Skew matrix S = A - A^T (entries ±1 off-diagonal), M = I + S would be the
    n×n core; the BORDERED skew-Hadamard H of order p+1: index set {inf} ∪ F_p,
    H[inf,inf] = 1, H[inf,j] = 1, H[i,inf] = -1, H[i,i] = 1, H[i,j] = chi(j-i).
    (Skew-type: H + H^T = 2I; Hadamard by the standard Paley argument.)
  * Binary code C(H) = row span over F_2 of (J - H)/2  (1 where H = -1).

Identifications (rigorous separators only; MISTAKE-067 discipline):
  * length 8: unique Type II code = ê₈ (also = RM(1,3)).
  * length 24: g₂₄ is the UNIQUE Type II [24,12] with minimum distance 8
    (equivalently no weight-4 words; Pless classification).  dim + Type II +
    d = 8 suffices.
  * length 32+: identification BY CONSTRUCTION — literal row-space EQUALITY with
    the extended QR code built independently from the generator polynomial
    (factor x^p - 1 over F_2, g(x) = prod over QR exponents; extend by overall
    parity, parity coordinate = inf).  No equivalence search, no enumerator
    matching needed (eQR₃₂ and RM(2,5) are isospectral — enumerators can't
    separate them).

Also: weight enumerators (exact, via Krawtchouk-free direct enumeration for
dim ≤ 16; for dim > 16 enumerate the dual side via MacWilliams if ever needed),
Gleason check at 8/24: the two enumerators ARE the Gleason Type II generators.
"""

import itertools, time, sys
from fractions import Fraction


# ----------------------------------------------------------- GF(2) basics

def rref(rows, ncols):
    """GF(2) row reduce a list of int bitmasks; returns (basis, pivots)."""
    basis = []
    pivots = []
    for r in rows:
        cur = r
        for b, p in zip(basis, pivots):
            if (cur >> p) & 1:
                cur ^= b
        if cur:
            p = cur.bit_length() - 1
            # reduce existing basis by new row
            nb = []
            for b2 in basis:
                if (b2 >> p) & 1:
                    b2 ^= cur
                nb.append(b2)
            basis = nb
            basis.append(cur)
            pivots = [b.bit_length() - 1 for b in basis]
            # keep sorted by pivot descending for determinism
            order = sorted(range(len(basis)), key=lambda i: -pivots[i])
            basis = [basis[i] for i in order]
            pivots = [pivots[i] for i in order]
    return basis, pivots


def span_equal(b1, b2):
    """Do two GF(2) bases span the same space? (reduce each against the other)"""
    if len(b1) != len(b2):
        return False
    def in_span(v, basis):
        cur = v
        for b in basis:
            p = b.bit_length() - 1
            if (cur >> p) & 1:
                cur ^= b
        return cur == 0
    return all(in_span(v, b2) for v in b1) and all(in_span(v, b1) for v in b2)


def weight_distribution(basis, n):
    """Exact weight distribution by enumerating the span (dim <= ~22)."""
    k = len(basis)
    dist = [0] * (n + 1)
    # Gray-code enumeration
    w = 0
    word = 0
    dist[0] += 1
    for g in range(1, 1 << k):
        low = (g & -g).bit_length() - 1  # Gray-code walk: flip one generator
        word ^= basis[low]
        dist[bin(word).count('1')] += 1
    return dist


def is_self_orthogonal(basis):
    return all(bin(a & b).count('1') % 2 == 0 for a in basis for b in basis)


def is_doubly_even(basis, n):
    # doubly-even iff all generators doubly-even and pairwise intersections even
    # (then all codewords have weight ≡ 0 mod 4)
    if not all(bin(b).count('1') % 4 == 0 for b in basis):
        return False
    return is_self_orthogonal(basis)


def min_weight(basis, n):
    dist = weight_distribution(basis, n)
    for w in range(1, n + 1):
        if dist[w]:
            return w
    return 0


# ----------------------------------------------------------- constructions

def quadratic_residues(p):
    return set((x * x) % p for x in range(1, p))


def paley_bordered_rows(p):
    """Rows of (J - H)/2 over F_2 as bitmasks; coordinates: bit p = inf,
    bits 0..p-1 = F_p. H in the THM-480 tournament gauge."""
    QR = quadratic_residues(p)
    rows = []
    # row inf: H[inf,inf]=1, H[inf,j]=1  -> (J-H)/2 = all zeros
    rows.append(0)
    for i in range(p):
        r = 0
        # H[i,inf] = -1 -> bit inf set
        r |= 1 << p
        for j in range(p):
            if j == i:
                continue  # H[i,i] = 1 -> 0
            chi = 1 if (j - i) % p in QR else -1
            if chi == -1:
                r |= 1 << j
        rows.append(r)
    return rows, p + 1


def eqr_code_rows(p):
    """Extended QR code of length p+1 from the generator polynomial over F_2.
    QR code: cyclic, generator g(x) = prod_{r in QR}(x - alpha^r) computed via
    factoring x^p - 1 over F_2: g(x) = the degree-(p-1)/2 divisor whose root
    exponent set is QR. Practical construction without finite-field roots:
    the QR code = cyclic code generated by the GF(2) idempotent
    E(x) = sum_{r in QR} x^r  (for p ≡ -1 mod 8 the idempotent of the QR code
    is sum over QR or QNR depending on convention) — to dodge convention traps
    we build BOTH cyclic codes generated by E_QR and E_QNR, take the one of
    dimension (p+1)/2, extend by parity (parity bit = inf), and ALSO return the
    other for control."""
    QR = quadratic_residues(p)
    QNR = set(range(1, p)) - QR
    out = []
    for S in (QR, QNR):
        E = 0
        for r in S:
            E |= 1 << r
        # cyclic shifts span the cyclic code generated by the idempotent
        rows = []
        for s in range(p):
            # shift E by s mod x^p - 1
            r = ((E << s) | (E >> (p - s))) & ((1 << p) - 1)
            rows.append(r)
        basis, _ = rref(rows, p)
        # extend by parity: parity bit (position p) = weight mod 2
        ext = []
        for b in basis:
            w = bin(b).count('1') % 2
            ext.append(b | (w << p))
        out.append((len(basis), ext))
    return out  # list of (dim, extended basis) for QR- and QNR-idempotents


# ----------------------------------------------------------- the lab parts

def identify(name, rows, n, expect_dim=None, full=True):
    basis, _ = rref(rows, n)
    k = len(basis)
    so = is_self_orthogonal(basis)
    de = is_doubly_even(basis, n)
    line = f"   {name}: [{n},{k}] self-orth={so} doubly-even={de}"
    dist = None
    if full and k <= 22:
        dist = weight_distribution(basis, n)
        d = next(w for w in range(1, n + 1) if dist[w])
        nz = {w: c for w, c in enumerate(dist) if c and w > 0}
        line += f" d={d} dist={nz}"
    print(line, flush=True)
    return basis, dist


def main():
    t0 = time.time()
    print("=== A. The Golay rung: C(I+S(Paley_23)) — HYP-2409.4 ===", flush=True)
    rows, n = paley_bordered_rows(23)
    b23, dist23 = identify("Paley_23 bordered, tournament gauge", rows, n)
    if dist23:
        is_golay = (len(b23) == 12 and is_doubly_even(b23, 24)
                    and is_self_orthogonal(b23)
                    and next(w for w in range(1, 25) if dist23[w]) == 8)
        print(f"   => Type II [24,12,8] (unique = g24): {is_golay}", flush=True)
        W = {w: c for w, c in enumerate(dist23) if c}
        print(f"   Gleason generator check: W == {{0:1, 8:759, 12:2576, 16:759, 24:1}}:"
              f" {W == {0:1, 8:759, 12:2576, 16:759, 24:1}}", flush=True)

    print("\n=== control: p=7 reproduces THM-480's e8 ===", flush=True)
    rows7, n7 = paley_bordered_rows(7)
    b7, dist7 = identify("Paley_7 bordered", rows7, n7)
    print(f"   W == {{0:1, 4:14, 8:1}}: "
          f"{ {w: c for w, c in enumerate(dist7) if c} == {0:1, 4:14, 8:1} }", flush=True)

    print("\n=== B. eQR unification (HYP-2410): row-space EQUALITY at p=7,23,31,47,71 ===", flush=True)
    for p in (7, 23, 31, 47, 71):
        rows, n = paley_bordered_rows(p)
        bp, _ = rref(rows, n)
        cands = eqr_code_rows(p)
        verdicts = []
        for dim, ext in cands:
            if dim == (p + 1) // 2:
                eb, _ = rref(ext, n)
                verdicts.append(span_equal(bp, eb))
        print(f"   p={p}: dim C(H)={len(bp)} (target {(p+1)//2}); "
              f"eQR row-space equality: {verdicts}", flush=True)

    print("\n=== C. order-32 trichotomy (HYP-2411) ===", flush=True)
    rows31, n31 = paley_bordered_rows(31)
    b31, dist31 = identify("Paley_31 bordered (= eQR_32 by part B)", rows31, n31)
    if dist31:
        nz = {w: c for w, c in enumerate(dist31) if c and 0 < w <= 16}
        print(f"   A4={dist31[4]} A8={dist31[8]} (d32+ has A4=120; RM(2,5)/eQR32 isospectral A4=0,A8=620)", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
