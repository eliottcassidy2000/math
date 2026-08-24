#!/usr/bin/env python3
"""Finite-exact extension of THM-4000's fixed scale-two family to [1,20]."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
import sys


sys.stdout.reconfigure(newline="\n")
ALPHA = F(4, 63)
C = ((F(2, 21), F(8, 63)), (F(55, 63), F(19, 21)))
EXPECTED_SEMANTIC_SHA256 = "f4d01f21eaafa28929136b506df5f1ac4b94d2c93f6beb7b9e4466a080a0a4da"
CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def merge(pieces):
    out = []
    for a, b in sorted(pieces):
        if not out or a > out[-1][1]:
            out.append([a, b])
        elif b > out[-1][1]:
            out[-1][1] = b
    return tuple((a, b) for a, b in out)


def danger(v):
    r = F(1, 14*v)
    pieces = []
    for k in range(v):
        c = F(k, v)
        a, b = c-r, c+r
        if a < 0:
            pieces.extend(((F(0), b), (a+1, F(1))))
        elif b > 1:
            pieces.extend(((a, F(1)), (F(0), b-1)))
        else:
            pieces.append((a, b))
    return merge(pieces)


def safe_set(body, bads):
    bad = merge([piece for v in body for piece in bads[v]])
    good = []
    x = F(0)
    for a, b in bad:
        if x < a:
            good.append((x, a))
        x = max(x, b)
    if x < 1:
        good.append((x, F(1)))
    return tuple(good)


def tail_start(mass, arcs):
    # THM-4000: disc_t<=r^2/(3t^2), while containment forces
    # disc_t>=((1-alpha)/alpha)m^2.  The inequality here is strict.
    wall = F(arcs*arcs) * ALPHA / (3*mass*mass*(1-ALPHA))
    t = isqrt(wall.numerator // wall.denominator)
    while F(t*t) <= wall:
        t += 1
    return t


def preimage(t):
    return tuple(((a+k)/t, (b+k)/t) for k in range(t) for a, b in C)


def in_open_pieces(x, pieces):
    return any(a < x < b for a, b in pieces)


def find_witness(good, obstruction):
    j = 0
    for a, b in good:
        while j < len(obstruction) and obstruction[j][1] <= a:
            j += 1
        if j == len(obstruction):
            return a
        c, d = obstruction[j]
        if c < a and b < d:
            continue
        if not (c < a < d):
            return a
        # The left endpoint lies in the current open obstruction component,
        # but the whole closed arc does not; its right endpoint d lies in G.
        require(d <= b, "linear witness endpoint order")
        return d
    return None


def clearance(v, x):
    r = (v*x) % 1
    return min(r, 1-r)


def fmt(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    require(sum((b-a for a, b in C), F(0)) == ALPHA, "obstruction mass")
    require(all(b-a == F(2, 63) for a, b in C), "obstruction widths")
    bads = {v: danger(v) for v in range(1, 21)}
    pre = {t: preimage(t) for t in range(1, 150)}
    transcript = sha256()
    strata = {u: {"bodies":0, "cells":0, "max_tail":0, "max_tail_body":None}
              for u in range(11, 21)}
    total_bodies = total_cells = 0
    boundary_controls = []

    for body in combinations(range(1, 21), 11):
        U = body[-1]
        stat = strata[U]
        stat["bodies"] += 1
        total_bodies += 1
        good = safe_set(body, bads)
        mass = sum((b-a for a, b in good), F(0))
        require(mass > 0, f"positive body mass {body}")
        arcs = len(good)
        tail = tail_start(mass, arcs)
        if (tail, body) > (stat["max_tail"], stat["max_tail_body"] or ()):
            stat["max_tail"], stat["max_tail_body"] = tail, body
        transcript.update((repr(body)+"|"+fmt(mass)+"|"+str(arcs)+"|"+str(tail)+"\n").encode("ascii"))

        first_t = U if U % 2 else U+1
        for t in range(first_t, tail, 2):
            obstruction = pre[t]
            y = find_witness(good, obstruction)
            require(y is not None, f"finite containment failure {(body,t)}")
            require(all(clearance(v, y) >= F(1,14) for v in body),
                    f"body witness {(body,t)}")
            require(not in_open_pieces((t*y)%1, C), f"quotient witness {(body,t)}")
            lift = None
            for k in (0, 1):
                x = (y+k)/2
                if (clearance(t, x) >= F(1,14)
                    and clearance(9*t, x) >= F(1,14)):
                    lift = k
                    break
            require(lift is not None, f"physical lift {(body,t)}")
            if stat["cells"] == 0:
                x = (y+lift)/2
                require(all(clearance(2*v, x) >= F(1,14) for v in body),
                        f"stratum physical-body control {(body,t)}")
            transcript.update((str(t)+"|"+fmt(y)+"|"+str(lift)+"\n").encode("ascii"))
            stat["cells"] += 1
            total_cells += 1

        # One exact boundary control per body: the analytic tail inequality
        # is strict at tail and fails to be asserted at tail-1.
        wall = F(arcs*arcs) * ALPHA / (3*mass*mass*(1-ALPHA))
        require(F(tail*tail) > wall, f"strict tail {body}")
        require(F((tail-1)*(tail-1)) <= wall, f"minimal tail {body}")

    require(total_bodies == 167960, "full [1,20] body count")
    require(total_cells == 574963, "finite cell count")
    require(strata[20]["bodies"] == 92378 and strata[20]["cells"] == 287670,
            "new max-20 stratum")
    require(max(stat["max_tail"] for stat in strata.values()) < len(pre),
            "preimage range covers every finite tail")
    for U in range(11, 21):
        stat = strata[U]
        boundary_controls.append((U, stat["bodies"], stat["cells"],
                                  stat["max_tail"], stat["max_tail_body"]))
        transcript.update(("S|"+repr(boundary_controls[-1])+"\n").encode("ascii"))
    digest = transcript.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("LRC14_SCALE2_NEAR_AP_EXTENSION_TO_20_20260824")
    print("scope=E_any_11_subset_of_1..20;odd_t>=max(E);row=2E_union_{t,9t}")
    print(f"bodies={total_bodies};finite_exact_cells={total_cells};containment_failures=0")
    print("finite_engine=exact_closed_body_arcs_versus_open_t_pullback_of_C_(1,9);physical_lift_checked")
    print("infinite_tail=THM4000_disc_t_le_r_squared_over_3t_squared_strict_noncontainment")
    print("strata=(maxE,bodies,finite_cells,max_tail,max_tail_body)")
    for item in boundary_controls:
        print(repr(item))
    print("conclusion=every_displayed_13_speed_row_has_clearance_at_least_1/14")
    print("nonconsequence=arbitrary_eleven_speed_body_and_LRC14_remain_open")
    print(f"semantic_sha256={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
