#!/usr/bin/env python3
"""Independent exact audit of the scale-two auxiliary and height-20 scouts.

This does not import either producing script.  For the finite near-AP cells it
computes the positive rational measure of G_E minus t^{-1}C, then checks an
explicit midpoint and one of its two physical lifts.  Thus its finite path is
stronger than endpoint-only noncontainment and independent of the packet
energy path used in THM-4000.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
import json
import sys


sys.stdout.reconfigure(newline="\n")
I = (F(2, 21), F(8, 63))
C = (I, (F(55, 63), F(19, 21)))
ALPHA = F(4, 63)
CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def merge(pieces):
    ans = []
    for left, right in sorted(pieces):
        if left >= right:
            continue
        if not ans or left > ans[-1][1]:
            ans.append([left, right])
        elif right > ans[-1][1]:
            ans[-1][1] = right
    return tuple((left, right) for left, right in ans)


def comb(v):
    radius = F(1, 14 * v)
    pieces = []
    for k in range(v):
        center = F(k, v)
        left, right = center - radius, center + radius
        if left < 0:
            pieces.extend(((F(0), right), (left + 1, F(1))))
        elif right > 1:
            pieces.extend(((left, F(1)), (F(0), right - 1)))
        else:
            pieces.append((left, right))
    return merge(pieces)


def covers_open(target, pieces):
    pieces = merge(pieces)
    for left, right in target:
        pos = left
        for a, b in pieces:
            if b <= pos:
                continue
            if a > left and pos == left:
                return False
            if a >= pos and pos > left and pos < right:
                return False
            pos = max(pos, b)
            if pos >= right:
                break
        if pos < right:
            return False
    return True


def clipped_gaps(target, pieces):
    left, right = target
    pieces = merge((max(left, a), min(right, b)) for a, b in pieces
                   if a < right and left < b)
    ans = []
    pos = left
    for a, b in pieces:
        if pos < a:
            ans.append((pos, a))
        pos = max(pos, b)
    if pos < right:
        ans.append((pos, right))
    return tuple(ans)


def safe_set(body, combs):
    bad = merge(piece for v in body for piece in combs[v])
    return clipped_gaps((F(0), F(1)), bad)


def strict_tail(mass, arcs):
    wall = F(arcs * arcs) * ALPHA / (3 * mass * mass * (1 - ALPHA))
    t = isqrt(wall.numerator // wall.denominator)
    while F(t * t) <= wall:
        t += 1
    return t, wall


def pullback(t):
    return tuple(((left + k) / t, (right + k) / t)
                 for k in range(t) for left, right in C)


def subtract_open(closed_pieces, open_pieces):
    """Positive-length pieces of a closed union outside an open union."""
    open_pieces = tuple(sorted(open_pieces))
    ans = []
    j = 0
    for left, right in closed_pieces:
        pos = left
        while j < len(open_pieces) and open_pieces[j][1] <= pos:
            j += 1
        k = j
        while k < len(open_pieces) and open_pieces[k][0] < right:
            a, b = open_pieces[k]
            if pos < a:
                ans.append((pos, min(a, right)))
            pos = max(pos, b)
            if pos >= right:
                break
            k += 1
        if pos < right:
            ans.append((pos, right))
    return tuple((a, b) for a, b in ans if a < b)


def clearance(v, x):
    residue = (v * x) % 1
    return min(residue, 1 - residue)


def audit_auxiliary():
    require(C[1] == (1 - I[1], 1 - I[0]), "C reflection")
    L = I[1] - I[0]
    require(L == F(2, 63), "I length")
    combs = {v: tuple(piece for piece in comb(v)
                      if any(piece[0] < right and left < piece[1]
                             for left, right in C))
             for v in range(1, 301)}

    # The analytic tooth-count bound leaves a<=25 when a<=b.  The exact
    # largest-gap reduction is deliberately recomputed here.
    require(F(4, 5 * L) == F(126, 5), "pair cutoff")
    candidates = []
    gap_table = []
    for a in range(1, 26):
        gaps = clipped_gaps(I, combs[a])
        require(gaps, f"one-comb noncover {a}")
        largest = max(right - left for left, right in gaps)
        bmax = int(F(1, 7) / largest)
        gap_table.append((a, str(largest), bmax))
        for b in range(a, bmax + 1):
            candidates.append((a, b))
            require(not covers_open(C, combs[a] + combs[b]),
                    f"pair cover {(a, b)}")
    require(len(candidates) == 20, "20 pair candidates")

    # A large finite hostile control is not part of the proof, but catches
    # sign, circle-wrap and reflection mistakes in the analytic reduction.
    require(not any(covers_open(C, combs[a] + combs[b])
                    for a in range(1, 301) for b in range(a, 301)),
            "pair hostile through 300")

    covers = []
    for c in range(1, 31):
        for a in range(1, c + 1):
            for b in range(a, c + 1):
                if covers_open(C, combs[a] + combs[b] + combs[c]):
                    covers.append((a, b, c))
    require(covers == [(8, 9, 10)], "unique triple through max 30")
    return gap_table, candidates, covers


def audit_height_twenty():
    combs = {v: comb(v) for v in range(1, 21)}
    pullbacks = {t: pullback(t) for t in range(1, 160)}
    strata = {u: [0, 0, 0, None] for u in range(11, 21)}
    total_bodies = total_cells = 0
    min_escape = None

    for body in combinations(range(1, 21), 11):
        U = body[-1]
        stat = strata[U]
        stat[0] += 1
        total_bodies += 1
        good = safe_set(body, combs)
        mass = sum((b - a for a, b in good), F(0))
        require(mass > 0, f"positive mass {body}")
        tail, wall = strict_tail(mass, len(good))
        require(F((tail - 1) ** 2) <= wall < F(tail ** 2),
                f"sharp tail {body}")
        if (tail, body) > (stat[2], stat[3] or ()):
            stat[2], stat[3] = tail, body

        first = U if U % 2 else U + 1
        for t in range(first, tail, 2):
            escaped = subtract_open(good, pullbacks[t])
            escape_mass = sum((b - a for a, b in escaped), F(0))
            require(escape_mass > 0, f"positive escape {(body, t)}")
            y = sum(escaped[0], F(0)) / 2
            require(all(clearance(v, y) >= F(1, 14) for v in body),
                    f"body midpoint {(body, t)}")
            require(not any(a < y < b for a, b in pullbacks[t]),
                    f"quotient midpoint {(body, t)}")
            lifts = []
            for k in (0, 1):
                x = (y + k) / 2
                if (clearance(t, x) >= F(1, 14)
                        and clearance(9 * t, x) >= F(1, 14)):
                    lifts.append(k)
            require(lifts, f"physical lift {(body, t)}")
            record = (escape_mass, body, t, y, tuple(lifts))
            if min_escape is None or record < min_escape:
                min_escape = record
            stat[1] += 1
            total_cells += 1

    require(total_bodies == 167960, "height-20 body count")
    require(total_cells == 574963, "height-20 finite cell count")
    require(sum(strata[u][0] for u in range(11, 20)) == 75582,
            "height-19 body count")
    require(sum(strata[u][1] for u in range(11, 20)) == 287293,
            "height-19 finite cell count")
    require(strata[19][0:2] == [43758, 172314], "max-19 stratum")
    require(strata[20][0:2] == [92378, 287670], "max-20 stratum")
    return strata, total_bodies, total_cells, min_escape


def fmt(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    gaps, candidates, triples = audit_auxiliary()
    strata, bodies, cells, minimum = audit_height_twenty()
    escape, body, t, y, lifts = minimum
    semantic = {
        "pair_gap_table": gaps,
        "pair_candidates": candidates,
        "triples_max_30": triples,
        "strata": {str(u): [s[0], s[1], s[2], list(s[3])] for u, s in strata.items()},
        "minimum_escape": [fmt(escape), list(body), t, fmt(y), list(lifts)],
    }
    digest = sha256(json.dumps(semantic, sort_keys=True,
                                separators=(",", ":")).encode("ascii")).hexdigest()
    print("LRC14_SCALE2_EXTENSION_INDEPENDENT_EXACT_AUDIT_20260824")
    print("auxiliary=no_one_or_two_comb_cover;finite_candidates=20;triple=(8,9,10)")
    print("triple_scope=unique_nondecreasing_cover_with_max_at_most_30")
    print("height19=bodies_75582;finite_cells_287293;failures_0")
    print(f"height20=bodies_{bodies};finite_cells_{cells};failures_0")
    print("finite_method=positive_exact_measure_of_G_minus_open_pullback_C_plus_explicit_physical_lift")
    print("minimum_positive_escape=" + repr((fmt(escape), body, t, fmt(y), lifts)))
    print("strata=" + repr(tuple((u, *strata[u][0:3], strata[u][3]) for u in strata)))
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
