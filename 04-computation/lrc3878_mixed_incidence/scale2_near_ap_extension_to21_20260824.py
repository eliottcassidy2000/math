#!/usr/bin/env python3
"""Finite-exact max(E)=21 extension of THM-4002's scale-two family."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import scale2_near_ap_extension_to20_20260824 as base


EXPECTED_BASE_SCRIPT_SHA256 = "f46f49b5ee8ddc48293f9390f15b61dfa191494e321e2413fc5410fd2a705334"
EXPECTED_SEMANTIC_SHA256 = "e7c6e8702fbedb9abf4a40029fc155ab109faf8d2cbaaf9b6ecdc4651e37ce94"
CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def fmt(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main():
    source = Path(__file__).with_name("scale2_near_ap_extension_to20_20260824.py")
    require(sha256(source.read_bytes()).hexdigest() == EXPECTED_BASE_SCRIPT_SHA256,
            "frozen exact-engine dependency")
    require(sum((b-a for a, b in base.C), F(0)) == base.ALPHA,
            "obstruction mass")

    bads = {v: base.danger(v) for v in range(1, 22)}
    pre = {t: base.preimage(t) for t in range(1, 150)}
    transcript = sha256()
    bodies = cells = 0
    max_tail = 0
    max_tail_body = None
    minimum_escape = None

    # These are exactly the eleven-subsets of [1,21] not already present in
    # the frozen [1,20] theorem universe: choose ten entries below 21 and adjoin 21.
    for lower in combinations(range(1, 21), 10):
        body = lower + (21,)
        bodies += 1
        good = base.safe_set(body, bads)
        mass = sum((b-a for a, b in good), F(0))
        require(mass > 0, f"positive body mass {body}")
        arcs = len(good)
        tail = base.tail_start(mass, arcs)
        if (tail, body) > (max_tail, max_tail_body or ()):
            max_tail, max_tail_body = tail, body
        transcript.update(
            (repr(body)+"|"+fmt(mass)+"|"+str(arcs)+"|"+str(tail)+"\n").encode("ascii")
        )

        for t in range(21, tail, 2):
            y = base.find_witness(good, pre[t])
            require(y is not None, f"finite containment failure {(body,t)}")
            require(all(base.clearance(v, y) >= F(1, 14) for v in body),
                    f"body witness {(body,t)}")
            require(not base.in_open_pieces((t*y) % 1, base.C),
                    f"quotient witness {(body,t)}")
            lift = None
            for k in (0, 1):
                x = (y+k)/2
                if (base.clearance(t, x) >= F(1, 14)
                        and base.clearance(9*t, x) >= F(1, 14)):
                    lift = k
                    break
            require(lift is not None, f"physical lift {(body,t)}")
            x = (y+lift)/2
            require(all(base.clearance(2*v, x) >= F(1, 14) for v in body),
                    f"physical body lift {(body,t)}")
            escape = min(
                *(base.clearance(2*v, x)-F(1, 14) for v in body),
                base.clearance(t, x)-F(1, 14),
                base.clearance(9*t, x)-F(1, 14),
            )
            require(escape >= 0, f"nonnegative clearance {(body,t)}")
            if minimum_escape is None or escape < minimum_escape:
                minimum_escape = escape
            transcript.update((str(t)+"|"+fmt(y)+"|"+str(lift)+"\n").encode("ascii"))
            cells += 1

        wall = F(arcs*arcs) * base.ALPHA / (3*mass*mass*(1-base.ALPHA))
        require(F(tail*tail) > wall, f"strict tail {body}")
        require(F((tail-1)*(tail-1)) <= wall, f"minimal tail {body}")

    require(bodies == 184756, "max-21 body count")
    require(cells == 781184, "max-21 finite cell count")
    require(max_tail == 80, "max-21 largest strict BV tail")
    require(max_tail_body == (1, 3, 4, 11, 13, 15, 16, 17, 18, 19, 21),
            "max-21 extremal tail body")
    require(max_tail < len(pre), "preimage range covers every finite tail")
    require(167960+bodies == 352716, "combined [1,21] body count")
    require(574963+cells == 1356147, "combined [1,21] finite cell count")

    transcript.update(
        ("S|"+repr((21, bodies, cells, max_tail, max_tail_body))+"\n").encode("ascii")
    )
    digest = transcript.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("LRC14_SCALE2_NEAR_AP_EXTENSION_TO_21_20260824")
    print("new_scope=max(E)=21;E_subset_1..21;size(E)=11;odd_t>=21")
    print(f"new_bodies={bodies};new_finite_exact_cells={cells};containment_failures=0")
    print("combined_scope=E_any_11_subset_of_1..21;odd_t>=max(E)")
    print("combined_bodies=352716;combined_finite_exact_cells=1356147")
    print(f"max_tail={max_tail};max_tail_body={max_tail_body}")
    print(f"minimum_sampled_clearance_surplus={fmt(minimum_escape)}")
    print("finite_engine=exact_closed_body_arcs_vs_open_pullback_plus_both_physical_lifts")
    print("infinite_tail=THM4002_strict_BV_noncontainment")
    print("conclusion=every_displayed_13_speed_row_has_clearance_at_least_1/14")
    print("nonconsequence=arbitrary_eleven_speed_body_and_LRC14_remain_open")
    print(f"semantic_sha256={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
