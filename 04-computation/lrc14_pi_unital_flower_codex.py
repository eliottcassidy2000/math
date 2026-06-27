"""Pi, flower turns, and q=3 unital ledgers around LRC14.

This scout keeps three prompts in one exact-ish ledger:

1. pi approximants: 22/7 and 31^(1/3).
2. the flower turn theta = 1/pi radians per petal.
3. q=3 geometric unital parameters v=q^3+1, k=q+1, lambda=1.

The goal is not to prove LRC14 from numerology.  The goal is to separate
which coincidences are structural enough to become proof hooks.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, pi


N = 14


def farey_count_proper(level: int) -> int:
    return sum(1 for b in range(2, level + 1) for a in range(1, b) if gcd(a, b) == 1)


def farey_packets(level: int) -> list[tuple[int, int]]:
    return [(a, b) for b in range(1, level + 1) for a in range(1, b + 1) if gcd(a, b) == 1]


def divisors(n: int) -> set[int]:
    return {d for d in range(1, n + 1) if n % d == 0}


def continued_fraction(x: float, n: int = 12) -> list[int]:
    out: list[int] = []
    for _ in range(n):
        a = int(x)
        out.append(a)
        x -= a
        if abs(x) < 1e-18:
            break
        x = 1.0 / x
    return out


def convergents(cf: list[int]) -> list[Fraction]:
    out: list[Fraction] = []
    for i in range(1, len(cf) + 1):
        value = Fraction(cf[i - 1], 1)
        for a in reversed(cf[: i - 1]):
            value = Fraction(a, 1) + Fraction(1, value)
        out.append(value)
    return out


def nearest_turn_returns(theta: float, limit: int = 120) -> list[tuple[int, int, float]]:
    """Best small circular returns for total angle k*theta near m*2*pi."""
    best: list[tuple[int, int, float]] = []
    record = 10.0
    for k in range(1, limit + 1):
        turns = k * theta / (2 * pi)
        m = round(turns)
        err = abs(k * theta - m * 2 * pi)
        if err < record:
            best.append((k, m, err))
            record = err
    return best


def unital_params(q: int) -> dict[str, int]:
    v = q**3 + 1
    k = q + 1
    lam = 1
    b = v * (v - 1) // (k * (k - 1))
    r = (v - 1) // (k - 1)
    return {"q": q, "v": v, "k": k, "lambda": lam, "blocks": b, "replication": r}


def print_pi_approximants() -> None:
    cbrt31 = 31 ** (1 / 3)
    print("[pi approximants]")
    rows = [
        ("22/7", 22 / 7, 22 / 7 - pi),
        ("31^(1/3)", cbrt31, cbrt31 - pi),
    ]
    for name, value, err in rows:
        print(f"  {name:<10} value={value:.15f} err={err:+.15e} abs={abs(err):.15e}")
    print(f"  pi^3={pi**3:.15f}; pi^3-31={pi**3 - 31:+.15e}")
    print("  31 is the nearest integer to pi^3, and 31 = 2*14 + 3.")
    print(f"  error ratio |22/7-pi|/|31^(1/3)-pi|={abs(22 / 7 - pi) / abs(cbrt31 - pi):.3f}")
    print()


def print_flower_turns() -> None:
    theta = 1 / pi
    delta = theta / (2 * pi)
    print("[flower turn theta = 1/pi radians]")
    print(f"  theta={theta:.15f} radians")
    print(f"  reciprocal-pi convergent: 7/22={7/22:.15f}, theta-7/22={theta - 7/22:+.15e}")
    print(f"  after 22 petals: 22/pi={22/pi:.15f} radians; drift from 7 radians={22/pi - 7:+.15e}")
    full_drift = (22 / pi) - (2 * pi)
    print(f"  but circle closure after 22 petals misses one full turn by {full_drift:+.15e} radians")
    print(f"  normalized turn delta=theta/(2*pi)={delta:.15f}")

    for name, x in (("theta", theta), ("delta", delta)):
        cf = continued_fraction(x, 10)
        print(f"  CF({name})={cf}")
        usable = [f for f in convergents(cf) if f.numerator != 0][:6]
        print(f"  convergents({name})={[str(f) for f in usable]}")

    print("  best full-circle returns k*theta ~= m*2*pi up to k=120:")
    for k, m, err in nearest_turn_returns(theta, 120):
        print(f"    k={k:3d}, turns={m:2d}, angular miss={err:.6e}")
    print("  readout: 22 is the denominator of theta ~= 7/22 in radian units;")
    print("           full-circle parastichy/return denominators start 19,20,59,79.")
    print()


def print_unital_ledger() -> None:
    u3 = unital_params(3)
    f3_products = {a * b for a, b in farey_packets(3)}
    f4_products = {a * b for a, b in farey_packets(4)}
    print("[q=3 unital and Farey/LRC14 checks]")
    print(
        "  q=3 unital: "
        f"points={u3['v']}, block_size={u3['k']}, blocks={u3['blocks']}, "
        f"replication={u3['replication']}, lambda={u3['lambda']}"
    )
    print(f"  proper strict Farey count |{{a/b: 0<a<b<=14}}|={farey_count_proper(14)}")
    print(f"  F_3 distinct products={sorted(f3_products)}, sum={sum(f3_products)}")
    print(f"  F_4 distinct products={sorted(f4_products)}, sum={sum(f4_products)}")
    print("  alignments:")
    print("    unital points 28 = sum(distinct products in F_4) = 2*14")
    print("    unital blocks 63 = strict proper F_14 packet count")
    print("    block size 4 = first Farey level containing 3/4 -> K_{3,4}")
    print("    replication 9 = |E(K_{3,3})|, the bipartite Kuratowski edge count")
    print("    q=3 = numerator threshold for nonplanar K_{a,b} carriers")
    print()


def print_pi_farey_packets() -> None:
    print("[pi approximants as Farey/summand/multiplicand packets]")
    packets = [
        ("22/7 reciprocal", 7, 22),
        ("3/41", 3, 41),
        ("3/4", 3, 4),
    ]
    for name, a, b in packets:
        print(
            f"  {name:<15} packet a/b={a}/{b}: sum={a+b:>3}, "
            f"product={a*b:>3}, K_{{{a},{b}}} "
            f"{'nonplanar' if min(a,b) >= 3 else 'planar'}"
        )
    print("  22/7 carries the apex numerator 7 and the flower denominator 22.")
    print("  3/41 carries the first nonplanar Farey-child numerator 3.")
    print()


def main() -> None:
    print("LRC14 PI / FLOWER / UNITAL SCOUT")
    print("=" * 78)
    print_pi_approximants()
    print_flower_turns()
    print_unital_ledger()
    print_pi_farey_packets()
    print("[proof-search readout]")
    print("  Treat design-unital as pair-uniform averaging, not as an exact tiler.")
    print("  Treat algebraic-unital as mark-preserving: maps must preserve the apex unit/identity.")
    print("  Non-unital quotients are exactly the dangerous ones: they forget whether the")
    print("  marked spectrum node stayed at 1/14 or migrated to 1/12, 2/27, 3/41, etc.")


if __name__ == "__main__":
    main()
