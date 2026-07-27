#!/usr/bin/env python3
"""Exact finite support referee for THM-2586.

Rebuild the fine THM-2584 tensor and check that every nonzero displacement
and every owner-clock cell has a positive theta=0 diagonal-rail edge.  The
two candidates are (v,t)=(0,0) and (6,12); in w=t/2 coordinates both have
w=v.  Their exact zero sets are disjoint.

The branch-freezing and all-sufficiently-late digit-cylinder continuation
are symbolic consequences of finite Perron expansion and THM-2583; they are
not finite-clock extrapolations in this script.
"""

from fractions import Fraction

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host
import lrc14_b_r5_theta_target_tensor_thm2584 as tensor


P = 13
QMOD = 7
RPACK = 13**2
DEPTH = 13**5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    print("== THM-2586: depth-five arrival-to-later-root diagonal ==")
    print("rebuild THM-2584 fine K tensor by both exact routes")
    E = base.build_set(base.PAT_E, base.ZELL)
    QB = base.build_set(host.PAT_QB, base.ZELL)
    _, _, K = tensor.exact_tensor(E, QB)
    denominator = RPACK * DEPTH * DEPTH * base.T_DEN

    zero_00 = []
    zero_6m1 = []
    choices = []
    for s in range(1, P):
        for ell in range(QMOD):
            mass_00 = K[s][0][0][ell]
            mass_6m1 = K[s][6][12][ell]
            if mass_00 == 0:
                zero_00.append((s, ell))
            if mass_6m1 == 0:
                zero_6m1.append((s, ell))
            require(mass_00 > 0 or mass_6m1 > 0,
                    "both theta-zero diagonal rails vanished")
            if mass_00 > 0:
                choices.append((s, ell, 0, 0, mass_00))
            else:
                choices.append((s, ell, 6, 12, mass_6m1))

    require(zero_00 == [(7, 4), (7, 5), (7, 6)],
            "(v,t)=(0,0) zero set changed")
    require(zero_6m1 == [(6, 4), (6, 5), (6, 6)],
            "(v,t)=(6,12) zero set changed")
    require(len(choices) == 84, "choice universe is not 12x7")
    for _, _, v, t, mass in choices:
        require((t - 2 * v) % P == 0, "selected edge has theta!=0")
        require((7 * t - v) % P == 0, "selected edge has w!=v")
        require(mass > 0, "selected diagonal mass is not positive")

    selected_00 = sum(v == 0 for _, _, v, _, _ in choices)
    selected_6m1 = len(choices) - selected_00
    minimum = min(mass for _, _, _, _, mass in choices)
    print("zero set K_ell(s,0,0):", zero_00)
    print("zero set K_ell(s,6,12):", zero_6m1)
    print("zero sets disjoint; all 84 (s!=0,ell) cells have a theta=0 edge")
    print("deterministic choice counts (0,0)/(6,12):",
          selected_00, "/", selected_6m1)
    print("minimum selected mass:", Fraction(minimum, denominator))
    print("for every selection: theta=t-2v=0 and w=7t=v: PASS")

    # The temporal identity used after freezing one positive Perron branch:
    # X=(x+a)/13^5 implies T^5 X=x, hence T^(N+5)X=T^N x.
    for a in (0, 1, DEPTH // 2, DEPTH - 1):
        for numerator in (1, 7, 25):
            x = Fraction(numerator, 26)
            if x >= 1:
                continue
            X = (x + a) / DEPTH
            require((DEPTH * X) % 1 == x,
                    "T^5 arrival-coordinate identity failed")
            for N in (1, 2, 5):
                require((P**(N + 5) * X) % 1 == (P**N * x) % 1,
                        "later-time identity failed")
    print("X=(x+a)/13^5 gives T^(N+5)X=T^N x: exact controls PASS")
    print("scope: support theorem plus exact time map; branch selection and")
    print("all-late cylinder continuation are symbolic, not finite extrapolation")
    print("no old-head/co-shift/relation-residue identification; all checks passed")


if __name__ == "__main__":
    main()
