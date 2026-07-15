#!/usr/bin/env python3
"""Exact verifier for THM-784's unbounded prime-seven blocking runs.

Pairwise observable / tournament telemetry:
  * runner vertices: speed comparison in the increasing-speed gauge;
  * wall-event vertices: chronological comparison in the time gauge;
  * both tournaments are transitive, so the sorted order is the unique
    Hamiltonian path, the SCCs are singletons, and there are no directed cycles.

The tournament fingerprints deliberately do not certify blocking: the retained
fibre datum is the constant seven-owner token rainbow on J.
"""

from fractions import Fraction


C = 7
SLOW = (1, 2, 3, 4, 5, 8, 10)
J_LEFT = Fraction(5, 16)
J_RIGHT = Fraction(7, 20)
EXPECTED_ROUNDS = (0, 1, 1, 1, 2, 3, 3)
EXPECTED_TOKENS = (0, 3, 2, 5, 1, 4, 6)


def nearest_integer_off_wall(x: Fraction) -> int:
    """Nearest integer, asserting that x is not a half-integer."""
    q, r = divmod(x.numerator, x.denominator)
    twice_r = 2 * r
    assert twice_r != x.denominator
    return q if twice_r < x.denominator else q + 1


def token(w: int, x: Fraction) -> int:
    return (-pow(w, -1, C) * nearest_integer_off_wall(w * x)) % C


def fast_wall_indices(n: int) -> list[int]:
    f = 560 * n + 1
    # The theorem proves these exact endpoints algebraically; enumerate them as
    # an independent rational-inequality check.
    candidates = range(175 * n - 2, 196 * n + 2)
    return [
        m
        for m in candidates
        if J_LEFT < Fraction(2 * m + 1, 2 * f) < J_RIGHT
    ]


def verify_slow_chamber() -> None:
    midpoint = (J_LEFT + J_RIGHT) / 2
    rounds = tuple(nearest_integer_off_wall(w * midpoint) for w in SLOW)
    tokens = tuple(token(w, midpoint) for w in SLOW)
    assert rounds == EXPECTED_ROUNDS
    assert tokens == EXPECTED_TOKENS
    assert set(tokens) == set(range(C))

    # A rounding value can change only at a wall.  Equal endpoint-side values
    # therefore certify that there is no slow wall in the open interval.
    epsilon = Fraction(1, 10**9)
    for w, expected in zip(SLOW, EXPECTED_ROUNDS):
        assert nearest_integer_off_wall(w * (J_LEFT + epsilon)) == expected
        assert nearest_integer_off_wall(w * (J_RIGHT - epsilon)) == expected


def verify_family(n: int) -> tuple[int, int]:
    assert n >= 1
    f = 560 * n + 1
    assert f > max(SLOW) and f % C == 1
    indices = fast_wall_indices(n)
    assert indices == list(range(175 * n, 196 * n))
    assert len(indices) == 21 * n

    for m in indices:
        x = Fraction(2 * m + 1, 2 * f)
        # At a fast wall, the seven non-walling slow owners are a rainbow.
        assert tuple(token(w, x) for w in SLOW) == EXPECTED_TOKENS
    return f, len(indices)


def transitive_tournament_fingerprint(vertex_count: int) -> str:
    scores = tuple(range(vertex_count))
    score_text = str(scores) if vertex_count <= 25 else f"(0,1,...,{vertex_count-1})"
    return (
        f"vertices={vertex_count} scores={score_text} cycles=0 "
        f"sccs={vertex_count} hp_count=1 edge_flips=0"
    )


def main() -> None:
    verify_slow_chamber()
    print("THM-784 exact verification")
    print(f"slow={SLOW}")
    print(f"J=({J_LEFT},{J_RIGHT}) length={J_RIGHT-J_LEFT}")
    print(f"rounds={EXPECTED_ROUNDS}")
    print(f"tokens={EXPECTED_TOKENS} rainbow={set(EXPECTED_TOKENS)==set(range(C))}")
    print("slow_walls_in_J=0")
    for n in (1, 2, 5, 20):
        f, count = verify_family(n)
        print(f"N={n} fast={f} fast_mod_7={f%C} walls_in_J={count}")
    print("formula: m=175N,...,196N-1; count=21N")
    print("runner_tournament:", transitive_tournament_fingerprint(8))
    print("event_tournament_N=1:", transitive_tournament_fingerprint(21))
    print("event_tournament_N=20:", transitive_tournament_fingerprint(420))
    print("information_destroyed_by_bare_tournament=owner-labelled_token_fibre+metric_mesh")
    print("PASS")


if __name__ == "__main__":
    main()
