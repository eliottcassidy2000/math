"""Overlap families between the additive and square triangular towers.

This addendum to HYP-2453 focuses on the user's follow-up:

* The additive tower A partitions all positive integers into square shells.
* The square tower B only occupies alternating triangular shells.
* Find/predict the cases where B equations or B sides appear inside A
  equations/sides, including the examples 10+11+12=13+14 and 21+22+23+24.

The output gives closed classifiers rather than only finite lists.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


@dataclass(frozen=True)
class Interval:
    lo: int
    hi: int

    @property
    def size(self) -> int:
        return self.hi - self.lo + 1

    def contains(self, other: "Interval") -> bool:
        return self.lo <= other.lo and other.hi <= self.hi

    def pad(self, other: "Interval") -> tuple[int, int]:
        assert self.contains(other)
        return (other.lo - self.lo, self.hi - other.hi)

    def __str__(self) -> str:
        return f"[{self.lo},{self.hi}]"


def tri(k: int) -> int:
    return k * (k + 1) // 2


def A_shell(n: int) -> Interval:
    return Interval(n * n, n * n + 2 * n)


def A_L(n: int) -> Interval:
    return Interval(n * n, n * n + n)


def A_R(n: int) -> Interval:
    return Interval(n * n + n + 1, n * n + 2 * n)


def B_shell(m: int) -> Interval:
    start = tri(2 * m)
    return Interval(start, start + 2 * m)


def B_L(m: int) -> Interval:
    start = tri(2 * m)
    return Interval(start, start + m)


def B_R(m: int) -> Interval:
    start = tri(2 * m)
    return Interval(start + m + 1, start + 2 * m)


def square_shell_index(x: int) -> int:
    n = isqrt(x)
    while (n + 1) * (n + 1) <= x:
        n += 1
    while n * n > x:
        n -= 1
    return n


def classify_interval_in_A_side(interval: Interval) -> tuple[str, int, tuple[int, int] | None]:
    """Classify one interval relative to the additive tower.

    Returns (state, n, padding), where state is one of:
      in_AL, in_AR, crosses_A_mid, crosses_square_boundary.
    The index n is the A shell containing interval.lo.
    """

    n = square_shell_index(interval.lo)
    shell = A_shell(n)
    if not shell.contains(interval):
        return ("crosses_square_boundary", n, None)
    left = A_L(n)
    right = A_R(n)
    if left.contains(interval):
        return ("in_AL", n, left.pad(interval))
    if right.contains(interval):
        return ("in_AR", n, right.pad(interval))
    return ("crosses_A_mid", n, None)


def pell_side_aligned(limit: int) -> list[tuple[int, int, int, int]]:
    """B_m.L in A_n.L and B_m.R in A_n.R with same split.

    Equivalent to T_n = 2*T_m, or x^2 - 2y^2 = -1 with
    x=2n+1 and y=2m+1.
    """

    out: list[tuple[int, int, int, int]] = []
    n = m = 0
    while len(out) < limit:
        n, m = 3 * n + 4 * m + 3, 2 * n + 3 * m + 2
        pad = n - m
        out.append((m, n, 2 * m + 1, pad))
    return out


def exact_side_equalities(search_limit: int) -> list[str]:
    hits: list[str] = []
    for m in range(1, search_limit + 1):
        b_sides = {"B_L": B_L(m), "B_R": B_R(m)}
        for n in range(1, search_limit + 2):
            a_sides = {"A_L": A_L(n), "A_R": A_R(n)}
            for b_name, b_int in b_sides.items():
                for a_name, a_int in a_sides.items():
                    if b_int == a_int:
                        hits.append(f"{b_name}({m}) = {a_name}({n}) = {b_int}, size={b_int.size}")
    return hits


def whole_B_shell_containment(limit: int) -> list[str]:
    rows: list[str] = []
    for m in range(1, limit + 1):
        b = B_shell(m)
        n = square_shell_index(b.lo)
        a = A_shell(n)
        if a.contains(b):
            split_offset = (n * n + n) - (tri(2 * m) + m)
            rows.append(
                f"m={m:4d}, n={n:4d}: B{m} {b} inside A{n} {a}; "
                f"size={b.size}/{a.size}, pad={a.pad(b)}, split_offset={split_offset}"
            )
    return rows


def side_containment_rows(limit: int, side: str) -> list[str]:
    rows: list[str] = []
    interval_fn = B_L if side == "L" else B_R
    for m in range(1, limit + 1):
        interval = interval_fn(m)
        state, n, pad = classify_interval_in_A_side(interval)
        if state.startswith("in_"):
            rows.append(
                f"B{m}.{side} {interval} -> A{n}.{state[-2:]} "
                f"size={interval.size}, pad={pad}"
            )
    return rows


def side_state_word(limit: int) -> list[str]:
    rows: list[str] = []
    code = {
        "in_AL": "L",
        "in_AR": "R",
        "crosses_A_mid": "M",
        "crosses_square_boundary": "S",
    }
    for start in range(1, limit + 1, 40):
        chunk: list[str] = []
        for m in range(start, min(limit, start + 39) + 1):
            l_state, _, _ = classify_interval_in_A_side(B_L(m))
            r_state, _, _ = classify_interval_in_A_side(B_R(m))
            chunk.append(code[l_state] + code[r_state])
        rows.append(f"m={start:3d}..{start+len(chunk)-1:3d}: {' '.join(chunk)}")
    return rows


def main() -> None:
    print("TRIANGULAR TOWER OVERLAP FAMILIES")
    print("=================================")
    print("Addendum to HYP-2453/T797/OPEN-Q-075.")
    print()
    print("COVERAGE")
    print("========")
    print("A_n covers every positive integer exactly once by square shells:")
    print("  A_n = [n^2, (n+1)^2-1].")
    print("B_m covers only alternating triangular shells:")
    print("  B_m = [T_{2m}, T_{2m+1}-1].")
    print("The skipped gaps are:")
    for m in range(0, 7):
        gap = Interval(tri(2 * m + 1), tri(2 * m + 2) - 1)
        print(f"  gap after B_{m}: {gap}, size={gap.size}")
    print()

    print("FAMILY 1: WHOLE B EQUATION SIDE-ALIGNED INSIDE WHOLE A EQUATION")
    print("================================================================")
    print("Condition: B_m.L subset A_n.L and B_m.R subset A_n.R.")
    print("Equivalent: T_n = 2*T_m, i.e. x^2 - 2y^2 = -1 with x=2n+1, y=2m+1.")
    print("Recurrence from (n,m)=(0,0): n'=3n+4m+3, m'=2n+3m+2.")
    print("B_m is the middle of A_n with symmetric outside padding n-m.")
    for m, n, size, pad in pell_side_aligned(8):
        print(
            f"  m={m:5d}, n={n:5d}: B size={size:5d}, A size={2*n+1:5d}, "
            f"left/right outside padding={pad}"
        )
    print()

    print("FAMILY 2: EXACT WHOLE-SIDE EQUALITY")
    print("===================================")
    print("Length+start equations show the only positive whole-side equality is:")
    for row in exact_side_equalities(200):
        print(f"  {row}")
    print("Proof sketch: B_L(m)=A_R(n) forces n=m+1 and m^2-2m-3=0, so m=3;")
    print("the other side-pair equations have no positive solution.")
    print()

    print("FAMILY 3: WHOLE B SHELL INSIDE ONE A SQUARE SHELL")
    print("=================================================")
    print("Classifier: n=floor(sqrt(T_{2m})); containment iff T_{2m+1}-1 <= (n+1)^2-1.")
    print("split_offset = A_mid - B_mid; 0 means the Family-1 left/right split alignment.")
    for row in whole_B_shell_containment(140)[:45]:
        print(f"  {row}")
    print("  ...")
    print()

    print("FAMILY 4: SINGLE B SIDE INSIDE A SIDE")
    print("=====================================")
    print("For any B side I=[u,v], let n=floor(sqrt(u)).")
    print("  if v>(n+1)^2-1: I crosses an A square boundary;")
    print("  elif v<=n^2+n: I is inside A_n.L;")
    print("  elif u>=n^2+n+1: I is inside A_n.R;")
    print("  else: I crosses the A left/right midpoint.")
    print()
    print("First B_L side containments:")
    for row in side_containment_rows(160, "L")[:35]:
        print(f"  {row}")
    print("  ...")
    print("First B_R side containments:")
    for row in side_containment_rows(160, "R")[:35]:
        print(f"  {row}")
    print("  ...")
    print()

    print("SIDE STATE WORD")
    print("===============")
    print("Each token is BL-state then BR-state: L=in A_L, R=in A_R, M=crosses A midpoint, S=crosses square boundary.")
    for row in side_state_word(160):
        print(f"  {row}")
    print()

    print("TAKEAWAYS")
    print("=========")
    print("1. The 10..14 example is the first nontrivial Pell-side-aligned equation containment.")
    print("2. The 21..24 example is the unique exact whole-side equality, B_3.L=A_4.R.")
    print("3. All looser side matches are predicted by a one-line floor-sqrt classifier.")
    print("4. The side-containment word is Sturmian/Beatty-like because T_{2m} drifts through square shells with slope sqrt(2).")


if __name__ == "__main__":
    main()
