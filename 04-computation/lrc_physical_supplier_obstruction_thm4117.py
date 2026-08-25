#!/usr/bin/env python3
"""Fraction-exact audit of the canonical THM-4049/MISTAKE-490 row.

This is deliberately standalone: it imports no repository code and makes no
claim beyond the three displayed THM-4112 supplier shapes and exact common
dilation/primitive-normalization invariance.
"""

from fractions import Fraction
from functools import reduce
from math import gcd


U_BODY = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22)
PAIR = (1, 3)
SCALE_S = 1
SCALE_T = 2**45
ROW = tuple(sorted({SCALE_S * u for u in U_BODY} |
                   {SCALE_T * p for p in PAIR}))
DELTA = Fraction(1, 14)


def gcd_many(values):
    return reduce(gcd, values)


def primitive_normalize(values):
    values = tuple(sorted(values))
    g = gcd_many(values)
    return g, tuple(v // g for v in values)


def torus_norm(x):
    r = x % 1
    return min(r, 1 - r)


def clearance(values, theta):
    return min(torus_norm(v * theta) for v in values)


def antipodal_clearance(values, theta):
    return min(clearance(values, theta),
               clearance(values, theta + Fraction(1, 2)))


def direct_supplier(values, base, count, first, ratio=2):
    values = set(values)
    base = set(base)
    if not base <= values:
        return False, tuple(sorted(base - values)), ()
    outliers = tuple(sorted(values - base))
    ok = (len(outliers) == count and outliers[0] >= first and
          all(outliers[i + 1] >= ratio * outliers[i]
              for i in range(len(outliers) - 1)))
    return ok, (), outliers


def ap7_seam_supplier(values):
    """Recognize THM-4112 section 7 rows 2q*C union {x,y}.

    Since C contains 1,...,7, gcd(C)=1.  Thus the gcd of the even block is
    exactly 2q, so both q and C are recovered rather than guessed.
    """
    values = tuple(sorted(values))
    odd = tuple(v for v in values if v % 2)
    even = tuple(v for v in values if v % 2 == 0)
    if len(odd) != 2 or len(even) != 11:
        return False, {"odd": odd, "even_count": len(even)}
    block_gcd = gcd_many(even)
    if block_gcd % 2:
        return False, {"odd": odd, "block_gcd": block_gcd}
    q = block_gcd // 2
    core = tuple(v // block_gcd for v in even)
    base = set(range(1, 8))
    missing = tuple(sorted(base - set(core)))
    outliers = tuple(sorted(set(core) - base))
    if missing or len(outliers) != 4:
        return False, {
            "odd": odd, "block_gcd": block_gcd, "q": q,
            "core": core, "missing_ap7": missing, "outliers": outliers,
        }
    a, b, c, d = outliers
    adaptive = 84 <= a < b < c < d and (a % 2 == 0 or b % 2 == 0)
    parity_free = 47 <= a and b >= 2 * a and b < c < d
    return adaptive or parity_free, {
        "odd": odd, "block_gcd": block_gcd, "q": q,
        "core": core, "missing_ap7": (), "outliers": outliers,
        "adaptive": adaptive, "parity_free": parity_free,
    }


def assert_supplier_controls():
    ap8 = tuple(range(1, 9)) + (94, 188, 376, 752, 1504)
    d0 = (3, 4, 5, 6, 8, 10, 12, 16, 32, 64, 128, 256, 512)
    ap7_core = tuple(range(1, 8)) + (85, 86, 91, 101)
    ap7 = tuple(sorted({2 * c for c in ap7_core} | {103, 105}))
    assert direct_supplier(ap8, range(1, 9), 5, 94)[0]
    assert direct_supplier(d0, (3, 4, 5, 6, 8, 10, 12), 6, 16)[0]
    assert ap7_seam_supplier(ap7)[0]
    return ap7, ap8, d0


def main():
    assert len(ROW) == 13
    assert gcd_many(ROW) == 1
    assert SCALE_T >= max(U_BODY)
    assert gcd(SCALE_S, SCALE_T) == gcd(*PAIR) == 1
    assert sum(ROW) <= 91**12

    odd = tuple(v for v in ROW if v % 2)
    even = tuple(v for v in ROW if v % 2 == 0)
    assert odd == (1, 15)
    assert gcd_many(even) == 2
    divided = tuple(v // 2 for v in even)
    assert divided == (2, 3, 4, 5, 6, 7, 8, 9, 11, 2**44, 3 * 2**44)
    assert tuple(v % 56 for v in divided) == (2, 3, 4, 5, 6, 7, 8, 9, 11, 32, 40)

    ap8_ok, ap8_missing, _ = direct_supplier(ROW, range(1, 9), 5, 94)
    d0_ok, d0_missing, _ = direct_supplier(
        ROW, (3, 4, 5, 6, 8, 10, 12), 6, 16)
    ap7_ok, ap7_data = ap7_seam_supplier(ROW)
    assert not ap8_ok and ap8_missing == (2, 3, 5, 7)
    assert not d0_ok and d0_missing == (3, 5)
    assert not ap7_ok
    assert ap7_data["q"] == 1
    assert ap7_data["missing_ap7"] == (1,)

    # Primitive normalization is unchanged after every tested common dilation.
    # The proof is symbolic: gcd(k*S)=k*gcd(S)=k; these rows are controls.
    phase = Fraction(9, 19)
    for k in (1, 2, 3, 7, 42, 101):
        g, normalized = primitive_normalize(tuple(k * v for v in ROW))
        assert g == k and normalized == ROW
        assert clearance(tuple(k * v for v in ROW), phase / k) == clearance(ROW, phase)

    row_clearance = clearance(ROW, phase)
    anti_clearance = antipodal_clearance(ROW, phase)
    assert row_clearance == Fraction(2, 19) > DELTA
    assert anti_clearance == Fraction(1, 38) < DELTA

    firewall_times = (Fraction(1, 28), Fraction(15, 28),
                      Fraction(5, 112), Fraction(61, 112))
    firewall_clearances = tuple(clearance(ROW, x) for x in firewall_times)
    assert firewall_clearances == (
        Fraction(1, 28), Fraction(1, 28),
        Fraction(1, 56), Fraction(1, 56))

    ap7_control, ap8_control, d0_control = assert_supplier_controls()

    print("status=PASS")
    print(f"row={ROW}")
    print(f"row_gcd={gcd_many(ROW)}; row_size={len(ROW)}; row_sum={sum(ROW)}; bound={91**12}")
    print(f"typed_source=s={SCALE_S};t={SCALE_T};p,q={PAIR};U={max(U_BODY)};t_ge_U={SCALE_T >= max(U_BODY)}")
    print(f"odd_forced_tails={odd}")
    print(f"even_block={even}; even_block_gcd={gcd_many(even)}")
    print(f"forced_q={ap7_data['q']}; forced_core={ap7_data['core']}")
    print(f"forced_core_mod56={tuple(v % 56 for v in divided)}")
    print(f"AP7_seam={ap7_ok}; missing_AP7={ap7_data['missing_ap7']}")
    print(f"AP8_plus_5={ap8_ok}; missing_AP8={ap8_missing}")
    print(f"D0_plus_6={d0_ok}; missing_D0={d0_missing}")
    print(f"phase={phase}; physical_clearance={row_clearance}; threshold={DELTA}")
    print(f"same_phase_antipodal_clearance={anti_clearance}")
    print(f"firewall_clearances={firewall_clearances}")
    print("dilation_controls=(1,2,3,7,42,101); primitive_normalization_and_phase_pullback=PASS")
    print(f"positive_AP7_control={ap7_control}; recognized={ap7_seam_supplier(ap7_control)[0]}")
    print(f"positive_AP8_control={ap8_control}; recognized={direct_supplier(ap8_control, range(1, 9), 5, 94)[0]}")
    print(f"positive_D0_control={d0_control}; recognized={direct_supplier(d0_control, (3,4,5,6,8,10,12), 6, 16)[0]}")


if __name__ == "__main__":
    main()
