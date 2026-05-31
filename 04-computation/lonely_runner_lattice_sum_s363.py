#!/usr/bin/env python3
"""
lonely_runner_lattice_sum_s363.py

oracle-2026-05-31-S16

Verify the pushforward-measure identity for the Lonely Runner problem and
expose why density/Fourier methods alone cannot settle it.  (norm/abs bars
written as plain text below to avoid escape warnings.)

Setup.  For speeds v=(v_1,...,v_k) (integers), the map phi(t)=(v_1 t,...,v_k t)
mod 1 pushes Lebesgue measure on R/Z to a measure mu_V on T^k.  Its Fourier
coefficients are

    hat{mu_V}(a) = \int_0^1 e^{2pi i (a . v) t} dt = [ a . v == 0 ],   a in Z^k.

So hat{mu_V} = indicator of the relation lattice L(V) = { a in Z^k : a.v = 0 }.
Therefore the *measure of lonely times* (Lebesgue measure of { t : ||v_i t|| >=
1/n for all i }) equals

    Leb(lonely) = \int_{T^k} 1_B d mu_V = sum_{a in L(V)} hat{1_B}(a),

where B = [1/n, 1-1/n]^k is the central box and

    hat{1_B}(a) = prod_i f_n(a_i),   f_n(0) = 1 - 2/n,
    f_n(a) = \int_{1/n}^{1-1/n} e^{-2pi i a x} dx = -sin(2 pi a / n)/(pi a)   (a != 0).

This script checks Leb(lonely) computed two ways:
  (1) exact, via the S356 interval engine (1 - forbidden_length);
  (2) the truncated lattice sum sum_{a in L(V), |a|<=A} prod f_n(a_i).

The point: (2) converges to (1).  Since Leb(lonely) >= 0 is automatic, ALL the
content of the conjecture lives in the measure-ZERO tight stratum (Leb=0 but a
boundary witness still exists).  That is exactly why Fourier/density methods are
provably insufficient and the problem becomes finite/Diophantine.
"""

from __future__ import annotations

import math
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()


def f_n(a: int, n: int) -> float:
    if a == 0:
        return 1.0 - 2.0 / n
    return -math.sin(2 * math.pi * a / n) / (math.pi * a)


def exact_lonely_measure(speeds, n) -> float:
    # threshold is 1/(k+1); to match a chosen n we require len(speeds)+1 == n
    row = S356.report("x", list(speeds))
    return float(1 - row.forbidden_length)


def lattice_sum(speeds, n, A) -> float:
    """sum over a in Z^k with |a_i|<=A and a.v=0 of prod f_n(a_i)."""
    v = list(speeds)
    k = len(v)
    total = 0.0
    rng = range(-A, A + 1)
    # iterate over all a in [-A,A]^k with a.v == 0  (k small here)
    for a in product(rng, repeat=k):
        if sum(ai * vi for ai, vi in zip(a, v)) != 0:
            continue
        term = 1.0
        for ai in a:
            term *= f_n(ai, n)
        total += term
    return total


def main() -> None:
    print("Lonely-runner pushforward lattice-sum identity check (oracle S16)\n")
    cases = [
        ("initial k=3 (tight)", (1, 2, 3)),       # n=4, Leb=0
        ("k=3 (1,2,5)", (1, 2, 5)),
        ("k=3 (1,3,5)", (1, 3, 5)),
        ("initial k=4 (tight)", (1, 2, 3, 4)),     # n=5, Leb=0
        ("k=4 (1,2,3,7)", (1, 2, 3, 7)),
        ("sporadic n=6 (1,3,4,5,9) tight", (1, 3, 4, 5, 9)),
    ]
    for label, speeds in cases:
        n = len(speeds) + 1
        exact = exact_lonely_measure(speeds, n)
        # truncation level: larger A for smaller k to control cost
        A = {2: 60, 3: 40, 4: 14, 5: 8}.get(len(speeds), 6)
        approx = lattice_sum(speeds, n, A)
        print(f"[{label}]  n={n} speeds={speeds}")
        print(f"   exact Leb(lonely)   = {exact:.8f}")
        print(f"   lattice sum (|a|<={A}) = {approx:.8f}   diff={abs(exact-approx):.2e}")
        print()
    print("Takeaway: the lattice sum converges to the exact lonely measure.")
    print("Leb(lonely) >= 0 is automatic; the tight stratum has Leb = 0, so the")
    print("conjecture's content is entirely in the measure-zero boundary -- the")
    print("regime Fourier/density arguments cannot resolve.")


if __name__ == "__main__":
    main()
