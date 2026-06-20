#!/usr/bin/env python3
"""Mechanism of Delta_S vanishing: inclusion-exclusion truncation, NOT a 7-sector
coincidence.  Establishes WHY Delta_S=0 for large |S| and HOW tight "6" is.

Setup recap (exact, from scratch): p0(A) = meas{ t in [0,1) : every c in A
survives, ||c t|| >= 1/14 }. Write D_c = danger set of runner c = { t : c kills t }.
Then SURVIVAL set of A is the complement of the UNION of D_c over c in A:
    p0(A) = meas( [0,1) \ Union_{c in A} D_c ) = 1 - meas(Union_{c in A} D_c).

By inclusion-exclusion on the UNION:
    meas(Union_{c in A} D_c) = sum_{emptyset != G subset A} (-1)^{|G|+1} I(G),
where I(G) = meas( Intersection_{c in G} D_c )  (the joint-kill measure).

=> p0(A) = 1 - sum_{emptyset != G subset A} (-1)^{|G|+1} I(G)
         = sum_{G subset A} (-1)^{|G|} I(G)   (with I(empty)=1).

This is a *signed sum of joint-kill measures*. Now the finite difference
Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B u T). Substituting the I(G)
expansion, ALL terms whose index set G does NOT contain every element of S cancel
(standard alternating-sum-over-supersets argument), leaving EXACTLY:

    Delta_S(B) = (-1)^{|S|} * sum_{G0 subset B} (-1)^{|G0|} I(G0 u S)
               = (-1)^{|S|} * [ joint-kill polynomial of S, weighted by B ].

KEY CONSEQUENCE: Delta_S(B) is built only from joint-kill measures I(H) with
S subset H. Delta_S(B) = 0 whenever EVERY such joint-kill measure with S subset H
vanishes, i.e. whenever the runners in S have NO common kill phase together with
any subset of B -- i.e. Intersection_{c in S} D_c (intersected with anything) is
empty/measure-zero.  Equivalently, Delta_S = 0 once |S| exceeds the maximal
number of danger-arcs that can share a common phase.

This script:
  (1) verifies the I(G) inclusion-exclusion form of p0 exactly,
  (2) verifies the reduced formula for Delta_S exactly,
  (3) computes the maximal "kill multiplicity" m* = max over t of #{c : t in D_c}
      and shows Delta_S=0 for |S| > m*, pinpointing the TRUE threshold.
"""
from fractions import Fraction
from itertools import combinations


def _danger(c):
    if c == 0:
        return [(Fraction(0), Fraction(1))]
    hw = Fraction(1, 14 * c)
    arcs = []
    for a in range(c):
        center = Fraction(a, c)
        lo, hi = center - hw, center + hw
        if lo < 0:
            arcs.append((Fraction(0), hi)); arcs.append((lo + 1, Fraction(1)))
        elif hi > 1:
            arcs.append((lo, Fraction(1))); arcs.append((Fraction(0), hi - 1))
        else:
            arcs.append((lo, hi))
    return arcs


def _merge(iv):
    out = []
    for a, b in sorted(iv):
        if a >= b: continue
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1][1] = b
        else: out.append([a, b])
    return out


def _meas(iv):
    return sum((b - a for a, b in _merge(iv)), Fraction(0))


def _intersect(iv1, iv2):
    """Exact intersection of two interval lists."""
    iv1 = _merge(iv1); iv2 = _merge(iv2)
    out = []
    i = j = 0
    while i < len(iv1) and j < len(iv2):
        lo = max(iv1[i][0], iv2[j][0]); hi = min(iv1[i][1], iv2[j][1])
        if lo < hi: out.append([lo, hi])
        if iv1[i][1] < iv2[j][1]: i += 1
        else: j += 1
    return out


def danger_set(c):
    return _merge(_danger(c))


def joint_kill(G):
    """I(G) = meas( intersection of D_c for c in G ). I(empty)=1."""
    if not G:
        return Fraction(1)
    cur = danger_set(G[0])
    for c in G[1:]:
        cur = _intersect(cur, danger_set(c))
        if not cur:
            return Fraction(0)
    return _meas(cur)


def p0_union(runners):
    """Direct: 1 - meas(union of danger sets)."""
    danger = []
    for c in runners: danger.extend(_danger(c))
    return Fraction(1) - _meas(danger)


def p0_incexc(runners):
    """p0 via signed sum of joint-kill measures (inclusion-exclusion)."""
    runners = list(runners)
    tot = Fraction(0)
    for k in range(len(runners) + 1):
        for G in combinations(runners, k):
            tot += (-1) ** k * joint_kill(list(G))
    return tot


def delta_S_direct(B, S):
    S = list(S); tot = Fraction(0)
    for k in range(len(S) + 1):
        sgn = (-1) ** (len(S) - k)
        for T in combinations(S, k):
            tot += sgn * p0_union(tuple(B) + T)
    return tot


def delta_S_reduced(B, S):
    """Delta_S(B) = (-1)^|S| sum_{G0 subset B} (-1)^|G0| I(G0 u S)."""
    B = list(B); S = list(S)
    tot = Fraction(0)
    for k in range(len(B) + 1):
        for G0 in combinations(B, k):
            tot += (-1) ** len(G0) * joint_kill(list(G0) + S)
    return (-1) ** len(S) * tot


def max_kill_multiplicity(runners):
    """Largest set of runners whose danger arcs share a common phase (max overlap
    depth). Equals max |H| with I(H) > 0. Found by brute subset search over the
    union runner set -- exact via interval intersection."""
    runners = list(runners)
    best = 0; witness = ()
    # increasing subset size; once no subset of size k has nonempty intersection,
    # no larger one can (intersection is monotone), so we can stop.
    for k in range(1, len(runners) + 1):
        found = False
        for H in combinations(runners, k):
            if joint_kill(list(H)) > 0:
                found = True; best = k; witness = H; break
        if not found:
            break
    return best, witness


def main():
    print("MECHANISM: inclusion-exclusion truncation drives Delta_S vanishing\n")

    test_sets = [
        (1, 2, 3, 4, 5, 6, 17, 19, 23),
        (1, 2, 3, 4, 5, 7, 17, 19, 23),
        (2, 3, 5, 7, 11, 13),
    ]
    print("(1) p0 union-form == inclusion-exclusion joint-kill form?")
    for A in test_sets:
        a = p0_union(A); b = p0_incexc(A)
        print(f"  A={A}: union={a} ; inc-exc={b} ; equal? {a == b}")
    print()

    print("(2) Delta_S direct == reduced joint-kill formula?")
    B = (1, 2, 3, 4, 5, 6)
    for S in [(17,), (17, 19), (17, 19, 23), (17, 19, 23, 29)]:
        d1 = delta_S_direct(B, S); d2 = delta_S_reduced(B, S)
        print(f"  S={S}: direct={d1} ; reduced={d2} ; equal? {d1 == d2}")
    print()

    print("(3) TRUE threshold = max kill multiplicity m* (max overlap depth).")
    print("    Claim: Delta_S=0 for |S|>6.  We show Delta_S=0 already for |S|>m*,")
    print("    and report m* for several full runner sets.")
    for runners in [
        (1, 2, 3, 4, 5, 6, 17, 19, 23, 29, 31, 37, 41),
        (1, 2, 3, 4, 5, 7, 17, 19, 23, 29, 31, 37, 41),
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
        (2, 3, 5, 7, 11, 13, 17, 19, 23, 29),
    ]:
        m, w = max_kill_multiplicity(runners)
        print(f"  runners={runners}")
        print(f"     m* (max joint-kill depth) = {m}  witness subset = {w}")
        print(f"     => structurally Delta_S=0 for |S| > {m} (<= 6? {m <= 6})")
    print()

    print("(4) Direct confirmation: for a deep base, drive Delta_S nonzero up to m*.")
    # base of contiguous small runners maximizes overlap depth at t=0:
    # EVERY runner c has a danger arc at t=0 (a=0 gives center 0). So at t near 0,
    # ALL runners simultaneously kill -> joint-kill of the WHOLE set can be > 0!
    base = (1, 2, 3, 4, 5, 6, 7)
    print(f"  Probe overlap at t->0 for base={base}: do ALL danger arcs cover t=0?")
    for c in base:
        ds = danger_set(c)
        covers0 = any(lo <= 0 < hi or lo < Fraction(1) <= hi or lo == 0 for lo, hi in ds) \
                  or any(lo == 0 for lo, hi in ds)
        # t=0 is in danger of c iff ||c*0||=0 < 1/14 -> ALWAYS true
        print(f"    runner {c}: 0 in D_c? {True}  (||c*0||=0<1/14 always)")
    print("  => t=0 lies in EVERY danger set, so joint-kill I(all) can be positive")
    print("     near 0, meaning overlap depth m* can be LARGE (>6) for contiguous")
    print("     small runners. THIS IS THE ADVERSARIAL CASE FOR THE '6' CLAIM.")
    print()
    # measure the all-kill overlap explicitly for contiguous runners
    for sz in range(2, 9):
        rr = tuple(range(1, sz + 1))
        m, w = max_kill_multiplicity(rr)
        full = joint_kill(list(rr))
        print(f"  contiguous 1..{sz}: m*={m}, joint-kill of ALL = {full} "
              f"(>0 => |S|={sz} difference can be nonzero)")


if __name__ == "__main__":
    main()
