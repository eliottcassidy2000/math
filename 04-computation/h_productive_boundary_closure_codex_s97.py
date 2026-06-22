#!/usr/bin/env python3
"""
h_productive_boundary_closure_codex_s97.py

Arithmetic proof ledger for the productive closure of the H-spectrum.

The previous even-graph-minor audit showed that arbitrary graph-minor moves are
not H-preserving.  The productive closure is the admissible one:

  * strong-component condensation, with H(T)=prod_i H(C_i);
  * OCF lower bounds inside a single strong component.

This script checks the exact small-H ledger for {7,21}.  It also records the
guardrail from the concurrent S33 correction: the forbidden set is not the
ideal 7Z, because later direct strong cores can realize values divisible by 7.
"""

from __future__ import annotations

from math import ceil


def odd_factorizations(n: int, min_factor: int = 3) -> list[tuple[int, ...]]:
    """Nondecreasing multiplicative factorizations into odd factors >= 3."""
    out: list[tuple[int, ...]] = []

    def rec(remaining: int, start: int, acc: list[int]) -> None:
        if remaining == 1:
            if acc:
                out.append(tuple(acc))
            return
        for f in range(start, remaining + 1, 2):
            if remaining % f == 0:
                acc.append(f)
                rec(remaining // f, f, acc)
                acc.pop()

    rec(n, min_factor, [])
    return out


def alpha1_lower_bound_strong(n: int) -> int:
    """Moon pancyclic lower bound for odd directed cycles in strong tournaments."""
    if n < 3:
        return 0
    return (n - 2) + sum(ceil(n / length) for length in range(5, n + 1, 2))


def h_lower_bound_strong(n: int) -> int:
    return 1 + 2 * alpha1_lower_bound_strong(n)


def route_status_for_7_21(route: tuple[int, ...]) -> str:
    if 7 in route:
        return "closed: contains forbidden factor 7 (THM-200)"
    if route == (21,):
        return "closed: single strong H=21, base <=8 plus Moon bound >=9 (THM-115)"
    return "not a {7,21} closure route"


def main() -> None:
    print("PRODUCTIVE H-CLOSURE LEDGER")
    print("admissible closure = strong components + OCF/Moon boundary bounds")
    print()

    print("Moon boundary lower bound for a single strong core:")
    print("  alpha1 >= (n-2) + sum_{odd L, 5<=L<=n} ceil(n/L)")
    print("  H >= 1 + 2*alpha1")
    for n in range(3, 16):
        alpha = alpha1_lower_bound_strong(n)
        h_lb = h_lower_bound_strong(n)
        marker = "  <-- exceeds 21" if n >= 9 and h_lb > 21 else ""
        print(f"  n={n:2d}: alpha1_lb={alpha:2d}, H_lb={h_lb:2d}{marker}")
    print()
    print("For n>=13 the c3 term alone gives alpha1>=n-2>=11, so H>=23.")
    print("The only large-core arithmetic checks needed for H=21 are n=9..12,")
    print("where the displayed pancyclic bound already gives alpha1>10.")
    print()

    for target in (7, 21):
        print(f"Target H={target}: factor routes through strong components")
        routes = odd_factorizations(target)
        for route in routes:
            print(f"  {route}: {route_status_for_7_21(route)}")
        print()

    print("Exact proof ledger:")
    print("  H=7:")
    print("    - prime, so any product realization needs a strong component with H=7.")
    print("    - THM-200 proves no tournament has H=7.")
    print("  H=21:")
    print("    - nontrivial product route is 3*7, closed by THM-200.")
    print("    - remaining route is one strong component with H=21.")
    print("    - base n<=8 closed by THM-079 Part G.")
    print("    - n>=9 closed by Moon pancyclic odd-cycle bound above.")
    print()

    print("Boundary-function interpretation:")
    print("  Strong-component condensation is the admissible boundary approach.")
    print("  It has a factor ledger, so limits do not depend on arbitrary graph paths.")
    print("  Even-graph minors are curvilinear/wild approaches: they preserve parity")
    print("  but not H, and can move H nonmonotonically.")
    print()

    print("Guardrail against the false 7Z ideal:")
    reported_strong_n7_div7 = [35, 49, 133, 147, 175]
    print("  Concurrent S33 reports strong n=7 values divisible by 7:")
    print(f"    {reported_strong_n7_div7}")
    print("  Therefore the permanent forbidden set is not the multiplicative ideal 7Z.")
    print("  The correct closure is finite low-boundary obstruction plus later direct")
    print("  strong-core re-entry.")


if __name__ == "__main__":
    main()
