#!/usr/bin/env python3
"""Exact checks for the three-weight Gaussian nullcone degree-gap theorem.

For

    P = Z^p a(S) + b(S) + W^q c(S),   S = ZW,

the charge-zero part is indexed by a single return multiplicity k.  This
script checks the resulting one-variable factorial formula against direct
Wick expansion, samples the two strict degree-gap regimes, and records the
tournament requested by the repository methodology.

The computation is only a verification instrument.  The proof is the
uniform channel estimate recorded with the theorem/reflection.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, factorial, gcd


Poly = tuple[int, ...]


def trim(f: Poly) -> Poly:
    data = list(f)
    while len(data) > 1 and data[-1] == 0:
        data.pop()
    return tuple(data)


def add(f: Poly, g: Poly) -> Poly:
    out = [0] * max(len(f), len(g))
    for i, value in enumerate(f):
        out[i] += value
    for i, value in enumerate(g):
        out[i] += value
    return trim(tuple(out))


def scale(f: Poly, scalar: int) -> Poly:
    return trim(tuple(scalar * value for value in f))


def mul(f: Poly, g: Poly) -> Poly:
    out = [0] * (len(f) + len(g) - 1)
    for i, left in enumerate(f):
        for j, right in enumerate(g):
            out[i + j] += left * right
    return trim(tuple(out))


def power(f: Poly, exponent: int) -> Poly:
    out: Poly = (1,)
    base = f
    n = exponent
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n //= 2
    return out


def shift(f: Poly, amount: int) -> Poly:
    return (0,) * amount + f


def laplace_factorial(f: Poly) -> int:
    return sum(coefficient * factorial(degree) for degree, coefficient in enumerate(f))


def return_data(p: int, q: int, a: Poly, b: Poly, c: Poly) -> tuple[int, int, int, Poly]:
    common = gcd(p, q)
    p0, q0 = p // common, q // common
    r = p0 + q0
    h = shift(mul(power(a, q0), power(c, p0)), p * q // common)
    return p0, q0, r, h


def three_weight_moment(p: int, q: int, a: Poly, b: Poly, c: Poly, m: int) -> int:
    p0, q0, r, h = return_data(p, q, a, b, c)
    total = 0
    for k in range(m // r + 1):
        zero_count = m - r * k
        multinomial = factorial(m) // (
            factorial(q0 * k) * factorial(p0 * k) * factorial(zero_count)
        )
        radial = mul(power(h, k), power(b, zero_count))
        total += multinomial * laplace_factorial(radial)
    return total


def wick_terms(p: int, q: int, a: Poly, b: Poly, c: Poly) -> dict[tuple[int, int], int]:
    terms: dict[tuple[int, int], int] = {}

    def insert(z_degree: int, w_degree: int, coefficient: int) -> None:
        key = (z_degree, w_degree)
        terms[key] = terms.get(key, 0) + coefficient

    for radial, coefficient in enumerate(a):
        insert(p + radial, radial, coefficient)
    for radial, coefficient in enumerate(b):
        insert(radial, radial, coefficient)
    for radial, coefficient in enumerate(c):
        insert(radial, q + radial, coefficient)
    return {key: value for key, value in terms.items() if value}


def direct_wick_moment(p: int, q: int, a: Poly, b: Poly, c: Poly, m: int) -> int:
    base = wick_terms(p, q, a, b, c)
    state: dict[tuple[int, int], int] = {(0, 0): 1}
    for _ in range(m):
        nxt: dict[tuple[int, int], int] = {}
        for (z0, w0), left in state.items():
            for (z1, w1), right in base.items():
                key = (z0 + z1, w0 + w1)
                nxt[key] = nxt.get(key, 0) + left * right
        state = nxt
    return sum(value * factorial(z_degree) for (z_degree, w_degree), value in state.items() if z_degree == w_degree)


def ratio_text(value: int, reference: int) -> str:
    ratio = Fraction(value, reference)
    return f"{float(ratio):.12g}  (num_bits={abs(ratio.numerator).bit_length()}, den_bits={ratio.denominator.bit_length()})"


def tournament_fingerprint(m: int, r: int, d: int, e: int) -> dict[str, object]:
    """Tournament on return channels k, oriented by radial degree.

    Observable: D(k)-D(l), where D(k)=d*m+(e-r*d)k.
    Gauge/switch: sign of the observable; ties use the canonical k order.
    """

    vertices = list(range(m // r + 1))
    radial_degree = {k: d * m + (e - r * d) * k for k in vertices}

    def before(left: int, right: int) -> bool:
        if radial_degree[left] != radial_degree[right]:
            return radial_degree[left] > radial_degree[right]
        return left < right

    scores = {k: 0 for k in vertices}
    edges: dict[tuple[int, int], int] = {}
    for i, left in enumerate(vertices):
        for right in vertices[i + 1 :]:
            winner, loser = (left, right) if before(left, right) else (right, left)
            scores[winner] += 1
            edges[(winner, loser)] = 1

    path = sorted(vertices, key=lambda k: (-radial_degree[k], k))
    cycle_count = 0
    for i in vertices:
        for j in vertices:
            for k in vertices:
                if i < j < k:
                    cycle_count += int(
                        ((i, j) in edges and (j, k) in edges and (k, i) in edges)
                        or ((j, i) in edges and (k, j) in edges and (i, k) in edges)
                    )

    return {
        "vertices": len(vertices),
        "observable": "D(k)-D(l), D(k)=d*m+(e-r*d)k",
        "score_histogram": dict(sorted(Counter(scores.values()).items())),
        "directed_3_cycles": cycle_count,
        "scc_sizes": [1] * len(vertices),
        "hamiltonian_path": path,
        "hamiltonian_path_count": 1,
        "edge_flips_vs_k_ascending": 0 if e <= r * d else comb(len(vertices), 2),
    }


def check_example(example: dict[str, object]) -> None:
    name = str(example["name"])
    p, q = int(example["p"]), int(example["q"])
    a, b, c = example["a"], example["b"], example["c"]
    assert isinstance(a, tuple) and isinstance(b, tuple) and isinstance(c, tuple)
    p0, q0, r, h = return_data(p, q, a, b, c)
    d, e = len(trim(b)) - 1, len(trim(h)) - 1
    gap = e - r * d

    print(f"\n[{name}]")
    print(f"charges=(+{p},-{q}), primitive_counts=(q/g={q0},p/g={p0}), r={r}")
    print(f"deg(b)={d}, deg(h)={e}, e-r*d={gap}, theorem_gate={abs(gap) >= r + 1}")

    for m in range(1, 9):
        via_channels = three_weight_moment(p, q, a, b, c, m)
        via_wick = direct_wick_moment(p, q, a, b, c, m)
        assert via_channels == via_wick, (name, m, via_channels, via_wick)
    print("direct_Wick_vs_channel_formula: PASS (m=1..8)")

    if gap <= -(r + 1):
        for m in example["sample_m"]:
            assert isinstance(m, int)
            value = three_weight_moment(p, q, a, b, c, m)
            reference = laplace_factorial(power(b, m))
            print(f"m={m:3d}  M_m/L(b^m)={ratio_text(value, reference)}")
    elif gap >= r + 1:
        for n in example["sample_n"]:
            assert isinstance(n, int)
            m = r * n
            value = three_weight_moment(p, q, a, b, c, m)
            reference = (
                factorial(m)
                // (factorial(q0 * n) * factorial(p0 * n))
                * laplace_factorial(power(h, n))
            )
            print(f"n={n:3d}  M_(r*n)/(C_n L(h^n))={ratio_text(value, reference)}")
    else:
        print("resonance_band_control: theorem deliberately makes no asymptotic-ratio claim")

    print("tournament:", tournament_fingerprint(24, r, d, e))


def main() -> None:
    print("GMC(2) THREE-WEIGHT DEGREE-GAP EXACT SCOUT")
    print("moment convention: E[Z^A W^B]=A!*1[A=B]")
    print("channel identity: k copies of the primitive (+p,-q) return plus m-rk zero weights")

    examples = [
        {
            "name": "boundary_b_dominant_p1q1",
            "p": 1,
            "q": 1,
            "a": (1,),
            "b": (1, 1, 1),
            "c": (1,),
            "sample_m": (10, 20, 40, 80),
        },
        {
            "name": "boundary_b_dominant_p2q3",
            "p": 2,
            "q": 3,
            "a": (1, 1),
            "b": (1, 1, 0, 1),
            "c": (1,),
            "sample_m": (10, 20, 40),
        },
        {
            "name": "boundary_h_dominant_p1q1",
            "p": 1,
            "q": 1,
            "a": (1, 1, 0, 0, 1),
            "b": (1, 1),
            "c": (1,),
            "sample_n": (5, 10, 20, 40),
        },
        {
            "name": "resonance_band_control",
            "p": 1,
            "q": 1,
            "a": (1, 1),
            "b": (1, 1),
            "c": (1,),
        },
    ]
    for example in examples:
        check_example(example)

    print("\nTOURNAMENT METHODOLOGY")
    print("vertices considered: monomials; charges; primitive atoms; radial-deficit channels k; proof obligations")
    print("chosen vertices: radial-deficit channels k (they preserve the factorial slope used by the proof)")
    print("destroyed information: coefficient phase and cancellations inside one fixed channel")
    print("challenged assumption: first-return atoms do NOT vanish separately inside a scalar moment")
    print("tie Hamiltonian path: increasing k when D(k)=D(l); every sampled tournament is transitive")
    print("status: exact identity checks + asymptotic diagnostics; proof is external to this script")


if __name__ == "__main__":
    main()
