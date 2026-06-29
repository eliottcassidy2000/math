#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The number 28 and the hidden families of LRC(14): perfect=triangular(Mersenne),
the octonion/Fano apex, the 7-14-21-28 dimension tower, and the Mersenne-cap-Heegner
characterization that singles out LRC(14). (mac-mini-2026-06-29-S14)

28 = C(8,2) = arc-count of the Forcade order 8 (HYP-3546) = the SECOND even perfect number.
This script chases what is STRUCTURAL (not numerology) about 28 in the LRC(14) apex.
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def is_prime(m):
    if m < 2: return False
    i = 2
    while i * i <= m:
        if m % i == 0: return False
        i += 1
    return True


def triangular(k): return k * (k + 1) // 2


# ---------- (1) every even perfect number = T(Mersenne prime) ----------
def thread_perfect_triangular():
    print("=" * 76)
    print("(1) EVEN PERFECT NUMBER = T(Mersenne prime); 28 = T(7), apex of LRC(14)")
    print("=" * 76)
    print(f"{'p':>3} {'Mp=2^p-1':>9} {'prime?':>6} {'perfect=2^(p-1)Mp':>17} "
          f"{'=T(Mp)?':>8} {'=C(Mp+1,2)?':>11} {'LRC(2*Mp)':>9}")
    for p in range(2, 8):
        Mp = 2 ** p - 1
        perfect = 2 ** (p - 1) * Mp
        isT = (perfect == triangular(Mp))
        isC = (perfect == Mp * (Mp + 1) // 2)
        tag = is_prime(Mp)
        lrc = f"LRC({2*Mp})" if tag else "-"
        print(f"{p:>3} {Mp:>9} {('PRIME' if tag else 'comp'):>6} {perfect:>17} "
              f"{('YES' if isT else 'no'):>8} {('YES' if isC else 'no'):>11} {lrc:>9}")
    print("\n  Even perfect numbers ARE the triangular numbers of Mersenne primes:")
    print("  6=T(3), 28=T(7), 496=T(31), 8128=T(127). 28=T(7): 7=LRC(14) apex prime.")
    print("  T(Mp)=C(Mp+1,2)=arcs of the order-(Mp+1) tournament (Mp+1=2^p Forcade order).\n")


# ---------- (2) the octonion / Fano / QR apex: {1,2,4} mod 7 ----------
def paley7():
    qr = set((x * x) % 7 for x in range(1, 7))  # {1,2,4}
    arc = [[1 if (i != j and ((j - i) % 7) in qr) else 0 for j in range(7)] for i in range(7)]
    return arc, qr


def thread_octonion_apex():
    print("=" * 76)
    print("(2) The apex 7: QR{1,2,4} = Fano line = octonion rule = B_0(Paley T_7)")
    print("=" * 76)
    arc, qr = paley7()
    print(f"  QR mod 7 = {sorted(qr)}  (the Paley arc set / the 'doubling' orbit 1->2->4->1).")
    # Fano plane: lines = translates of the QR difference set {1,2,4} mod 7
    lines = sorted(set(frozenset((b + t) % 7 for b in (1, 2, 4)) for t in range(7)))
    print(f"  Fano plane lines (translates of {{1,2,4}} mod 7): {len(lines)} lines, "
          f"each |line|=3; 7 points, 7 lines (PG(2,2)).")
    # octonion multiplication: e_i e_{i+1} = e_{i+3}, indices mod 7 with the {1,2,4} (QR) pattern
    # the 7 'quaternionic triples' = the Fano lines = translates of {1,2,4}
    print(f"  Octonion imaginary units e_1..e_7: the 7 associative triples = the 7 Fano lines")
    print(f"  = translates of the QR set {{1,2,4}}. So Im(O) multiplication = Paley T_7 / QR structure.")
    out0 = [v for v in range(7) if arc[0][v]]
    print(f"  B_0(T_7) = out-nbhd of vertex 0 = {out0} = QR{{1,2,4}} (Mersenne doubling -> Fano line).")
    # the 7 sectors of LRC(14) <-> 7 Fano points / Im(O) units
    print(f"  => LRC(14)'s 7 sectors <-> 7 Fano points <-> 7 imaginary octonion units.\n")


# ---------- (3) the 7-14-21-28 dimension tower = Im(O), G2, so(7), so(8) ----------
def thread_dimension_tower():
    print("=" * 76)
    print("(3) The 7-14-21-28 tower: Im(O), G2, so(7), so(8); 28 = cut (+) cycle of an 8-tournament")
    print("=" * 76)
    rows = [
        (7,  "Im(O) / Fano pts / LRC sectors / apex prime / cut-space dim (scores) of 8-tourn"),
        (14, "dim G2 = Aut(O) / LRC(14) ORDER / 2*apex"),
        (21, "dim so(7) = C(7,2) = CYCLE-space (tiles) of an 8-tournament = T(6)"),
        (28, "dim so(8) = C(8,2) = ARC count = perfect number = T(7)"),
    ]
    for d, desc in rows:
        print(f"  {d:>3} : {desc}")
    n = 8
    cut = n - 1                 # base-path arcs = score hierarchy = cut space
    cyc = (n - 1) * (n - 2) // 2  # tiles = cycle space = C(n-1,2)
    print(f"\n  8-tournament arc space d=C(8,2)={cut+cyc}: cut(scores)={cut} (+) cycle(tiles)={cyc}.")
    print(f"  so(8) = so(7) (+) R^7  <=>  28 = 21 + 7  <=>  arcs = tiles (+) scores (cut(+)cycle).")
    print(f"  so(7) = G2 (+) (7-coset)  <=>  21 = 14 + 7  <=>  tiles = G2 (+) sectors.")
    print(f"  Tower 7,14,21,28 = 7*{{1,2,3,4}}; LRC order 14 sits at G2; arcs 28 at so(8).\n")


# ---------- (4) Mersenne AND Heegner: the family that singles out LRC(6),(14) ----------
def thread_mersenne_heegner():
    print("=" * 76)
    print("(4) Apex p both MERSENNE and HEEGNER => only {3,7} => only LRC(6), LRC(14)")
    print("=" * 76)
    heegner_primes = [3, 7, 11, 19, 43, 67, 163]       # Q(sqrt(-p)) class number 1, p odd
    mersenne_primes = [m for m in (2**k - 1 for k in range(2, 14)) if is_prime(m)]
    print(f"  Heegner primes (Q(sqrt(-p)) h=1, odd): {heegner_primes}  -- all = 3 mod 4")
    print(f"  Mersenne primes 2^k-1: {mersenne_primes}  -- all = 3 mod 4 (2^k-1 = 3 mod 4, k>=2)")
    both = sorted(set(heegner_primes) & set(mersenne_primes))
    print(f"  MERSENNE and HEEGNER: {both}  => LRC orders 2p = {[2*p for p in both]}")
    print(f"\n  For apex p in {both}: (a) p=3 mod 4 => Paley T_p exists, saddle index (p-1)/2 ODD")
    print(f"  => Borsuk-Ulam regime; (b) p Heegner => Q(sqrt(-p)) h=1 = the GENTLEST cyclotomic")
    print(f"  Fejer-Bochner SOS; (c) p Mersenne => arc count C(p+1,2) = PERFECT number = T(p).")
    print(f"  LRC(6) is PROVED; LRC(14) is the FIRST OPEN and the LAST small case with all three.")
    print(f"  saddle index (p-1)/2: p=3->1, p=7->3 (both odd); next Mersenne 31->15, 127->63 (odd).\n")


if __name__ == "__main__":
    thread_perfect_triangular()
    thread_octonion_apex()
    thread_dimension_tower()
    thread_mersenne_heegner()
    print("=" * 76)
    print("SYNTHESIS: 28 = T(7) = C(8,2) = dim so(8) = 2nd perfect number is the arc-count of")
    print("LRC(14)'s Forcade order 8; its apex 7 = Im(O)=Fano=QR{1,2,4}=B_0(T_7); the tower")
    print("7-14-21-28 = Im(O),G2,so(7),so(8) with 14=LRC order=dim G2; and apex 7 is one of only")
    print("TWO primes (3,7) that are simultaneously Mersenne (perfect arcs) and Heegner (h=1),")
    print("singling out LRC(6) and LRC(14). 14 is the octonion/G2 order; 28 its triality so(8).")
    print("=" * 76)
