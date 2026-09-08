"""Exact audit of the Euler-inspired, gauge-optimized LRC leak budget.

Reproduce from the repository root:
  python 04-computation/lrc14_euler_adverse_leak_20260908.py

Universe: all R,L in {-1,0,1}^p for p=2,3,4, excluding constant R.
Independent paths: direct hinge minimization, weighted median formula,
and the vertices of the scalar-constraint dual LP. All arithmetic is rational.
Physical controls reconstruct shifted-danger intersections on the circle;
they are auxiliary densities, not canonical scalar-cover owner packets.
"""

from fractions import Fraction as Q
from itertools import product


def require(condition):
    if not condition:
        raise AssertionError("Exact audit requirement failed")


def mean(xs):
    return sum(xs, Q(0)) / len(xs)


def centered(xs):
    m = mean(xs)
    return tuple(x - m for x in xs)


def energy(xs):
    return mean(tuple(x * x for x in centered(xs)))


def budget_direct(r, leak):
    """A convex piecewise-linear minimum is at a leak-coordinate knot."""
    return min(mean(tuple(max(Q(0), -x * (y - b))
                          for x, y in zip(r, leak))) for b in set(leak))


def budget_median(r, leak):
    atoms = sorted((y, abs(x)) for x, y in zip(r, leak) if x)
    half = sum(w for _, w in atoms) / 2
    acc = Q(0)
    for b, w in atoms:
        acc += w
        if acc >= half:
            break
    mad = mean(tuple(abs(x) * abs(y - b) for x, y in zip(r, leak)))
    corr = mean(tuple(x * y for x, y in zip(r, leak)))
    return (mad - corr) / 2, b


def budget_dual(r, leak):
    """Enumerate vertices 0<=z<=1, sum r*z=0, independently of a median."""
    p = len(r)
    best = Q(0)
    for free in range(p):
        if not r[free]:
            continue
        rest = [j for j in range(p) if j != free]
        for bits in product((0, 1), repeat=p - 1):
            signed = sum((r[j] * z for j, z in zip(rest, bits)), Q(0))
            zf = -signed / r[free]
            if 0 <= zf <= 1:
                value = -(sum((r[j] * leak[j] * z
                               for j, z in zip(rest, bits)), Q(0))
                          + r[free] * leak[free] * zf) / p
                best = max(best, value)
    return best


def verify(R, leak, dual=True):
    r = centered(R)
    E = energy(R)
    require(E > 0)
    B = budget_direct(r, leak)
    B2, b = budget_median(r, leak)
    require(B == B2)
    if dual:
        require(B == budget_dual(r, leak))
    C = tuple(x + y for x, y in zip(R, leak))
    EC, EL = energy(C), energy(leak)
    require(EC * E >= (E - B) ** 2)
    require(B * B <= E * EL)
    # Constants in the leak are irrelevant to the optimized gauge.
    require(budget_direct(r, tuple(y + Q(17, 5) for y in leak)) == B)
    return E, EL, EC, B, b


def danger(center):
    """Independent physical interval model, strict endpoints irrelevant to mass."""
    radius = Q(1, 14)
    ans = []
    for shift in (-1, 0, 1):
        lo, hi = max(Q(0), center + shift - radius), min(Q(1), center + shift + radius)
        if lo < hi:
            ans.append((lo, hi))
    return ans


def overlap(A, B):
    return sum((max(Q(0), min(b, d) - max(a, c))
                for a, b in A for c, d in B), Q(0))


def main():
    count = 0
    for p in (2, 3, 4):
        rows = tuple(product(map(Q, (-1, 0, 1)), repeat=p))
        local = 0
        for R in rows:
            if len(set(R)) == 1:
                continue
            for leak in rows:
                verify(R, leak)
                local += 1
        count += local
        print(f"EXHAUSTIVE p={p} real signed profiles={local}: median=hinge=dual, bounds PASS")
    print(f"EXHAUSTIVE total={count}")

    # Exact sharpness for all these theta, and analytically for every theta>=0.
    R = (Q(-1), Q(0), Q(1))
    E = energy(R)
    for theta in (Q(0), Q(1, 4), Q(1), Q(2), Q(100)):
        leak = tuple(-theta * x for x in R)
        _, _, EC, B, _ = verify(R, leak)
        require(B == theta * E and EC * E == (E - B) ** 2)
    print("SHARPNESS theta=0,1/4,1,2,100: B=theta E, exact lower bound PASS")

    # The equality invoice is necessary for cancellation, not equivalent to it.
    witness = verify(R, (Q(1), Q(1), Q(-1)))
    require(witness[3] == witness[0] and witness[2] == Q(2, 9))
    print("HOSTILE B=E with nonconstant C=(0,1,0): E(C)=2/9 PASS")

    # Reconstruct the p=13 present-Q auxiliary density control from intervals.
    D = danger(Q(0))
    A = tuple(danger(Q(s, 13)) for s in range(13))
    R = tuple(overlap(D, a) for a in A)
    expected = tuple(Q(1, 7) if s == 0 else Q(6, 91) if s in (1, 12)
                     else Q(0) for s in range(13))
    require(R == expected)
    selected = {2, 3, 4, 5, 6}
    w = tuple(a if s in selected else D for s, a in enumerate(A))
    require(all(sum((b - a for a, b in ws), Q(0)) == Q(1, 7) for ws in w))
    C = tuple(overlap(ws, a) for ws, a in zip(w, A))
    leak = tuple(c - r for c, r in zip(C, R))
    ER, EL, EC, B, b = verify(R, leak, dual=False)
    require((ER, EL, EC, B) == (Q(2508, 1399489), Q(40, 8281),
                               Q(6018, 1399489), Q(125, 107653)))
    require(EL > ER and B < ER)
    lower = (ER - B) ** 2 / ER
    print(f"PHYSICAL p=13 all densities indicator-valued and mass=1/7: PASS")
    print(f"PHYSICAL E(R)={ER} E(L)={EL} E(C)={EC}")
    print(f"PHYSICAL B*={B} median={b} B*/E(R)={B / ER} lower={lower}")
    print("PHYSICAL old norm gate FAILS; adverse budget gate PASSES")

    constant_leak = tuple(Q(1, 7) - r for r in R)
    ER, EL, EC, B, _ = verify(R, constant_leak, dual=False)
    require(EL == ER == B and EC == 0)
    print("PHYSICAL all w_s=A_s: constant C=1/7, B*=E(R), exact cancellation PASS")
    loose_upper = 2 * ER
    wrong_bound = (loose_upper - ER) ** 2 / ER
    valid_upper_bound = max(Q(0), ER - loose_upper) ** 2 / ER
    require(wrong_bound > EC == valid_upper_bound)
    print("UPPER-BOUND HOSTILE exact cancellation, B*<=2E: square of upper-E INVALID; positive-part gate=0 PASS")
    gain = verify(R, tuple(2 * r for r in R), dual=False)
    require(gain[3] == 0 and gain[1] == 4 * gain[0] and gain[2] == 9 * gain[0])
    print("FAVORABLE L=2R: B*=0, E(L)=4E(R), E(C)=9E(R) PASS")
    print("SCOPE auxiliary real profiles only; no owner/word intertwiner or LRC(14) closure")
    print("PASS")


if __name__ == "__main__":
    main()
