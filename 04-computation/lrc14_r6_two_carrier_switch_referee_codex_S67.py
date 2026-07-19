#!/usr/bin/env python3
"""Exact referee for THM-1151's two-carrier switch classification.

Only integer arithmetic is used.  A carrier ``o`` is safe from source ``s``
at switch ``c`` exactly when

    dist(c*o/s, 14 Z) >= 1.

After clearing the positive denominator ``s``, this becomes a remainder
comparison modulo ``14*s``.  The exhaustive box is a replay, not the proof;
the script also checks the branchwise constructive certificate used in the
written proof on every row in that box.
"""

from __future__ import annotations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def carrier_safe(c: int, other: int, source: int) -> bool:
    """Return dist(c*other/source, 14 Z) >= 1 exactly."""
    modulus = 14 * source
    remainder = (c * other) % modulus
    return min(remainder, modulus - remainder) >= source


def switch_safe(M: int, u: int, v: int, source: int, c: int) -> bool:
    require(1 <= M <= 12, "core maximum out of range")
    require(M < u < v, "carriers must be ordered above the core")
    require(source in (u, v), "source is not a carrier")
    require(1 <= c <= source // M, "switch exceeds the core cap")
    return carrier_safe(c, u, source) and carrier_safe(c, v, source)


def is_exception(M: int, u: int, v: int) -> bool:
    return (
        8 <= M <= 12
        and u == M + 1
        and 13 * M + 14 <= v <= 15 * M - 1
    )


def constructive_switch(M: int, u: int, v: int) -> tuple[int, int, str] | None:
    """The exact certificate from the proof of THM-1151."""
    if switch_safe(M, u, v, u, 1):
        return (u, 1, "least carrier, c=1")

    # Failure of c=1 forces v/u into (14m-1,14m+1), m>=1.
    q = (v + u - 1) // u
    require(v > 13 * u, "failed least-carrier switch did not enter m>=1 band")

    if q % 14 != 0:
        require(q <= v // M, "ceil(v/u) violates the core cap")
        require(switch_safe(M, u, v, v, q), "ceil(v/u) certificate failed")
        return (v, q, "large carrier, ceil(v/u)")

    if q + 1 <= v // M:
        require(switch_safe(M, u, v, v, q + 1), "ceil(v/u)+1 certificate failed")
        return (v, q + 1, "large carrier, ceil(v/u)+1")

    require(is_exception(M, u, v), "uncertified row is not in the claimed list")
    return None


def brute_has_switch(M: int, u: int, v: int) -> bool:
    for source in (u, v):
        for c in range(1, source // M + 1):
            if switch_safe(M, u, v, source, c):
                return True
    return False


def main() -> None:
    # This box contains every exceptional pair and extends far into all proof
    # branches (including several 14m bands).  The theorem itself is analytic.
    u_max = 200
    v_max = 1000
    rows = 0
    failures: list[tuple[int, int, int]] = []
    branch_counts: dict[str, int] = {}

    for M in range(1, 13):
        for u in range(M + 1, u_max + 1):
            for v in range(u + 1, v_max + 1):
                rows += 1
                certificate = constructive_switch(M, u, v)
                brute = brute_has_switch(M, u, v)
                predicted = not is_exception(M, u, v)
                require(brute == predicted, f"classification mismatch at {(M, u, v)}")
                require((certificate is not None) == predicted, f"certificate mismatch at {(M, u, v)}")
                if certificate is None:
                    failures.append((M, u, v))
                else:
                    branch_counts[certificate[2]] = branch_counts.get(certificate[2], 0) + 1

    expected = [
        (M, M + 1, v)
        for M in range(8, 13)
        for v in range(13 * M + 14, 15 * M)
    ]
    require(failures == expected, "exception list mismatch")
    require(len(failures) == 30, "exception count mismatch")

    # Directly audit the two no-switch mechanisms in every exceptional row.
    for M, u, v in failures:
        require(u // M == 1, "small-source cap should be one")
        require(v // M == 14, "large-source cap should be fourteen")
        require(not switch_safe(M, u, v, u, 1), "small-source switch unexpectedly works")
        for c in range(1, 14):
            require(not carrier_safe(c, u, v), "c<=13 should miss the smaller carrier")
        require(not carrier_safe(14, v, v), "c=14 should fail at the source itself")

    print("THM-1151 exact two-carrier switch classification")
    print(f"exhaustive replay box: 1<=M<=12, M<u<={u_max}, u<v<={v_max}")
    print(f"rows checked: {rows}")
    print("constructive branch histogram:")
    for label in sorted(branch_counts):
        print(f"  {label}: {branch_counts[label]}")
    print(f"switch failures: {len(failures)}")
    print("failure formula: 8<=M<=12, u=M+1, 13M+14<=v<=15M-1")
    for M in range(8, 13):
        values = [v for m, _, v in failures if m == M]
        print(f"  M={M}: u={M + 1}, v={values[0]}..{values[-1]} ({len(values)} rows)")

    print("symbolic branch ledger:")
    print("  c=1 failure => v/u in (14m-1,14m+1), m>=1")
    print("  ceil(v/u) not 0 mod 14 => use source v and c=ceil(v/u)")
    print("  ceil(v/u)=0 mod 14 and cap allows => use c=ceil(v/u)+1")
    print("  cap can fail only at m=1, M>=8, u=M+1, v<15M")
    print("  on that band: source-u cap=1; source-v cap=14; every switch is obstructed")

    print("Tournament Analysis:")
    print("  order tournament on two carrier vertices: scores=(0,1)")
    print("  directed cycles=0; SCC sizes=(1,1); Hamiltonian paths=1")
    print("  challenged vertices: runners, carriers, ratios, switches, proof obligations")
    print("  faithful object: source-labelled switch-obstruction digraph with the core cap")
    print("  naked two-vertex tournament destroys all band and cap information")
    print("VERDICT: exact replay and constructive certificate agree; precisely 30 pairs remain")


if __name__ == "__main__":
    main()
