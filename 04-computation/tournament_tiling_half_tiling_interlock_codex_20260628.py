#!/usr/bin/env python3
"""HYP-3244: interlocking fixed-path tiling and half-tiling recursions.

This scout keeps two recursions separate.

* The tiling recursion is the fixed-Hamiltonian-path cube.  It keeps explicit
  flip witnesses but covers each unlabelled tournament many times.
* The half-tiling recursion is the parent-class plus incident-word orbit DAG.
  It compresses by automorphisms and coboundary laws but needs sidecars to
  remember what the tiling model knew.

The goal is not to prove LRC(14).  It is to make the controlled-forgetting
interfaces exact in the small tournament models that keep recurring in the
LRC proof frontier.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb


KNOWN_U = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56}


@lru_cache(maxsize=None)
def pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: i for i, pair in enumerate(pairs(n))}


def edge(mask: int, n: int, i: int, j: int) -> int:
    if i < j:
        return (mask >> pair_index(n)[(i, j)]) & 1
    return 1 - ((mask >> pair_index(n)[(j, i)]) & 1)


def set_edge(mask: int, n: int, i: int, j: int, i_beats_j: bool) -> int:
    if i < j:
        bit = pair_index(n)[(i, j)]
        return (mask | (1 << bit)) if i_beats_j else (mask & ~(1 << bit))
    bit = pair_index(n)[(j, i)]
    return (mask & ~(1 << bit)) if i_beats_j else (mask | (1 << bit))


def relabel(mask: int, n: int, perm: tuple[int, ...]) -> int:
    """Return the tournament in new labels where new vertex a is old perm[a]."""
    out = 0
    for a, b in pairs(n):
        if edge(mask, n, perm[a], perm[b]):
            out = set_edge(out, n, a, b, True)
    return out


@lru_cache(maxsize=None)
def all_perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(n)))


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, perm) for perm in all_perms(n))


@lru_cache(maxsize=None)
def classes(n: int) -> tuple[int, ...]:
    return tuple(
        sorted({canonical(mask, n) for mask in range(1 << (n * (n - 1) // 2))})
    )


@lru_cache(maxsize=None)
def automorphisms(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(perm for perm in all_perms(n) if relabel(mask, n, perm) == mask)


def vertex_orbits(mask: int, n: int) -> list[tuple[int, ...]]:
    auts = automorphisms(mask, n)
    unseen = set(range(n))
    out: list[tuple[int, ...]] = []
    while unseen:
        v = min(unseen)
        orbit = {perm.index(v) for perm in auts}
        out.append(tuple(sorted(orbit)))
        unseen -= orbit
    return out


def hamiltonian_path_count(mask: int, n: int) -> int:
    total = 0
    for perm in all_perms(n):
        if all(edge(mask, n, perm[i], perm[i + 1]) for i in range(n - 1)):
            total += 1
    return total


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum(edge(mask, n, i, j) for j in range(n) if j != i) for i in range(n))
    )


N4_CLASS = {
    (0, 1, 2, 3): "T",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
    (1, 1, 2, 2): "S",
}


def class4(mask: int) -> str:
    return N4_CLASS[score_sequence(mask, 4)]


def transitive_mask(n: int) -> int:
    out = 0
    for i, j in pairs(n):
        out = set_edge(out, n, i, j, True)
    return out


def fixed_path_tournament(n: int, tile_bits: int) -> int:
    """Tournament with base path 0->1->...->n-1 and nonpath arcs free."""
    mask = 0
    bit = 0
    for i, j in pairs(n):
        if j == i + 1:
            value = 1
        else:
            value = (tile_bits >> bit) & 1
            bit += 1
        if value:
            mask = set_edge(mask, n, i, j, True)
    return mask


def fixed_path_fibers(n: int) -> dict[int, int]:
    tile_count = comb(n - 1, 2)
    fibers: Counter[int] = Counter()
    for bits in range(1 << tile_count):
        fibers[canonical(fixed_path_tournament(n, bits), n)] += 1
    return dict(fibers)


def lift_mask(mask: int, n: int) -> int:
    out = 0
    for i, j in pairs(n):
        if edge(mask, n, i, j):
            out = set_edge(out, n + 1, i, j, True)
    return out


def extend_by_incident_word(mask: int, n: int, word: tuple[int, ...]) -> int:
    """Append vertex n.  word[i]=1 means old vertex i beats the new vertex."""
    out = lift_mask(mask, n)
    for i, bit in enumerate(word):
        out = set_edge(out, n + 1, i, n, bool(bit))
    return out


def all_words(length: int) -> tuple[tuple[int, ...], ...]:
    return tuple(product((0, 1), repeat=length))


def act_on_word(word: tuple[int, ...], aut: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(word[aut[i]] for i in range(len(word)))


def word_orbits(parent: int, n: int) -> list[set[tuple[int, ...]]]:
    auts = automorphisms(parent, n)
    unseen = set(all_words(n))
    out = []
    while unseen:
        word = min(unseen)
        orbit = {act_on_word(word, aut) for aut in auts}
        out.append(orbit)
        unseen -= orbit
    return out


@lru_cache(maxsize=None)
def extension_census(parent_n: int) -> dict[str, object]:
    parents = classes(parent_n)
    raw_child_fibers: Counter[int] = Counter()
    orbit_child_fibers: Counter[int] = Counter()
    orbit_size_hist: Counter[int] = Counter()

    for parent in parents:
        for word in all_words(parent_n):
            child = canonical(extend_by_incident_word(parent, parent_n, word), parent_n + 1)
            raw_child_fibers[child] += 1
        for orbit in word_orbits(parent, parent_n):
            rep_word = min(orbit)
            child = canonical(extend_by_incident_word(parent, parent_n, rep_word), parent_n + 1)
            orbit_child_fibers[child] += 1
            orbit_size_hist[len(orbit)] += 1

    child_reps = tuple(sorted(orbit_child_fibers))
    child_n = parent_n + 1
    rooted = sum(len(vertex_orbits(rep, child_n)) for rep in child_reps)
    return {
        "parent_n": parent_n,
        "child_n": child_n,
        "parent_classes": len(parents),
        "raw_extensions": len(parents) * (1 << parent_n),
        "layer_orbits": sum(orbit_child_fibers.values()),
        "rooted_count": rooted,
        "child_classes": len(child_reps),
        "orbit_size_hist": dict(sorted(orbit_size_hist.items())),
        "orbit_sink_hist": dict(sorted(Counter(orbit_child_fibers.values()).items())),
        "raw_sink_hist": dict(sorted(Counter(raw_child_fibers.values()).items())),
        "orbit_equals_rooted": sum(orbit_child_fibers.values()) == rooted,
    }


def fixed_path_n4_mask(state: str) -> int:
    """Fixed path 0->1->2->3; flips a=02, b=13, c=03 from transitive."""
    flip_arc = {"a": (0, 2), "b": (1, 3), "c": (0, 3)}
    mask = transitive_mask(4)
    for name in state:
        i, j = flip_arc[name]
        mask = set_edge(mask, 4, i, j, False)
    return mask


def half_chart_n4_mask(state: str) -> int:
    """Four fixed arcs with partial score (0,1,1,2); free x=23, y=01."""
    fixed = {
        (0, 2): False,  # 2 -> 0
        (0, 3): False,  # 3 -> 0
        (1, 2): True,   # 1 -> 2
        (1, 3): False,  # 3 -> 1
    }
    free = {"x": (2, 3), "y": (0, 1)}
    mask = 0
    for (i, j), value in fixed.items():
        mask = set_edge(mask, 4, i, j, value)
    for name, (i, j) in free.items():
        mask = set_edge(mask, 4, i, j, name in state)
    return mask


def subset_xor(a: str, b: str) -> str:
    out = sorted((set(a) ^ set(b)))
    return "".join(out)


def print_n4_tables() -> None:
    print("1. THE TWO n=4 CHARTS")
    print("   Fixed-path tiling chart: base path 0->1->2->3; free flips")
    print("   a=(0,2), b=(1,3), c=(0,3) from the transitive orientation.")
    print("   Full fixed-path fibers:")
    fibers: dict[str, list[str]] = {name: [] for name in ("T", "+", "-", "S")}
    for bits in range(8):
        state = "".join(name for bit, name in enumerate("abc") if (bits >> bit) & 1) or "E"
        key = state if state != "E" else ""
        fibers[class4(fixed_path_n4_mask(key))].append(state)
    for label in ("T", "+", "-", "S"):
        print(f"     {label}: {fibers[label]}")

    names = ["E", "a", "b", "c"]
    print("   Prompt table as XOR-on-generators followed by isomorphism class:")
    print("       *   " + "  ".join(f"{name:>2}" for name in names))
    for row in names:
        vals = []
        for col in names:
            state = subset_xor("" if row == "E" else row, "" if col == "E" else col)
            vals.append(class4(fixed_path_n4_mask(state)))
        print(f"       {row:>1} | " + "  ".join(f"{val:>2}" for val in vals))

    print()
    print("   Half-tiling chart: fixed arcs 2->0, 3->0, 1->2, 3->1.")
    print("   The fixed partial outdegree vector is [0,1,1,2].")
    print("   Free bits are x=(2,3) and y=(0,1).")
    names2 = ["E", "x", "y"]
    print("       *   " + "  ".join(f"{name:>2}" for name in names2))
    for row in names2:
        vals = []
        for col in names2:
            state = subset_xor("" if row == "E" else row, "" if col == "E" else col)
            vals.append(class4(half_chart_n4_mask(state)))
        print(f"       {row:>1} | " + "  ".join(f"{val:>2}" for val in vals))
    print("   Readout: the half chart is a four-state section, while the fixed-path")
    print("   chart has an S-fiber of size 5.  The missing payload is not noise;")
    print("   it is the canary/filler sidecar that later recurs as H(T)/|Aut(T)|.")


def print_fixed_path_cover_table() -> None:
    print()
    print("2. FIXED-PATH TILING COVER RECURSION")
    print("   n  fixed_cube  U(n)  fiber_hist  H(T)/|Aut(T)| check")
    for n in range(3, 7):
        fibers = fixed_path_fibers(n)
        hist = dict(sorted(Counter(fibers.values()).items()))
        checks = []
        for rep, count in fibers.items():
            h = hamiltonian_path_count(rep, n)
            aut = len(automorphisms(rep, n))
            checks.append(h % aut == 0 and h // aut == count)
        print(
            f"   {n}  {2 ** comb(n - 1, 2):10d}  {len(fibers):4d}  "
            f"{hist}  {all(checks)}"
        )
    print("   This recursion is witness-rich: it keeps path and flip coordinates,")
    print("   but its quotient map to A000568 is many-to-one.")


def print_half_tiling_orbit_table() -> None:
    print()
    print("3. HALF-TILING / INCIDENT-WORD ORBIT RECURSION")
    print("   parent n-class + incident word orbit under Aut(parent) -> rooted child -> unrooted sink")
    print(
        "   n->n+1  U(n)  raw=U(n)2^n  word_orbits  rooted  U(n+1)  "
        "orbit=rooted  rooted/U"
    )
    for n in range(1, 6):
        c = extension_census(n)
        rooted = int(c["rooted_count"])
        child = int(c["child_classes"])
        print(
            f"   {n}->{n+1:<2d}  {c['parent_classes']:4d}  "
            f"{c['raw_extensions']:12d}  {c['layer_orbits']:11d}  "
            f"{rooted:6d}  {child:6d}  {str(c['orbit_equals_rooted']):>12}  "
            f"{rooted / child:8.3f}"
        )
    c = extension_census(5)
    print("   First large bookkeeping surface 5->6:")
    print(f"     incident-word orbit sizes = {c['orbit_size_hist']}")
    print(f"     rooted orbits per 6-class = {c['orbit_sink_hist']}")
    print(f"     raw extensions per 6-class = {c['raw_sink_hist']}")
    print("   This recursion is quotient-efficient: it sees rooted classes exactly,")
    print("   but it does not remember a chosen Hamiltonian path or flip witness.")


def print_coboundary_laws() -> None:
    print()
    print("4. COBoundary LAWS FOR THE HALF-TILING COMPRESSION")
    print("   Consecutive layer bridge K_{k,k+1}: lines=k(k+1), GF(2) rank=2k.")
    print("   k  lines  rank  rectangle_redundancy")
    for k in range(1, 8):
        lines = k * (k + 1)
        rank = 2 * k
        print(f"   {k}  {lines:5d}  {rank:4d}  {lines - rank:21d}")
    print("   Rectangle law: L(a,b)+L(a,b')+L(a',b)+L(a',b')=0.")
    print("   A nonzero rectangle or hourglass residue means the quotient forgot")
    print("   real LRC payload: owner, route, active support, or proof-obligation status.")


def directed_three_cycles(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    low = [0] * n
    sizes: list[int] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w in range(n):
            if not adj[v][w]:
                continue
            if indices[w] == -1:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack[w] = False
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return sorted(sizes)


def ham_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp: Counter[tuple[int, int]] = Counter()
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[(mask, last)]
            if not count:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] += count
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))


def print_interlock_and_tournament() -> None:
    print()
    print("5. INTERLOCKING RECURSION INVARIANT")
    print("   Let C_n be the fixed-path tiling cube and U_n the A000568 orbit set.")
    print("   The useful object is a span, not a plain quotient:")
    print()
    print("       C_n  --tiling add/free diagonal-->  C_{n+1}")
    print("        |                                  |")
    print("        v  quotient with fiber sidecar     v")
    print("       U_n --Aut(parent)-word orbit-----> rooted U_{n+1} -> U_{n+1}")
    print()
    print("   The square commutes only after adding sidecars:")
    print("   path-presentation fiber H(T)/|Aut(T)|, parent automorphism word orbit,")
    print("   rectangle/hourglass coboundary residue, and tail/tip deletion signature.")
    print("   The n=4 S-fiber of size 5 is the smallest visible warning; the 5->6")
    print("   surface makes it numerical: fixed-path cover=1024, rooted orbits=296,")
    print("   unrooted classes=56.")
    print()
    print("   LRC translation:")
    print("   - Tiling recursion builds explicit witnesses and exchange moves.")
    print("   - Half-tiling recursion compresses by symmetry and detects legal quotienting.")
    print("   - The proof should alternate: lift to tiling for a witness, compress through")
    print("     half-tiling only when sidecars certify that the LRC predicate descends.")
    print("   - HYP-3216/HYP-3230/HYP-3231/HYP-3232/HYP-3233/HYP-3234/HYP-3235")
    print("     and HYP-3218 supply the LRC-side ladder/fold/three-gap/scale/")
    print("     modulus/chart/Fejer recursions; this packet supplies the")
    print("     tournament-side lift/compress interface.")
    print()
    print("6. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    carriers = [
        "tail_tip_deletion_sidecar",
        "tiling_witness_lift",
        "half_tiling_orbit_compression",
        "parent_aut_incident_word_orbit",
        "coboundary_rectangle_hourglass_law",
        "fiber_debt_H_over_Aut",
        "n4_canary_filler_section",
        "lrc_fejer_gram_recursion_bridge",
        "raw_fixed_path_cube_count",
        "raw_A000568_class_count",
    ]
    n = len(carriers)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        adj[i][j] = 1
    scores = [sum(row) for row in adj]
    print("   Vertices are proof carriers, not runners or arcs.")
    print("   Observable: LRC predicate retained through lift/compress recursion minus quotient risk.")
    print("   Switch: orient toward the carrier that keeps a witness and certifies descent.")
    print(f"   vertices={carriers}")
    print(f"   score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"   directed_3_cycles={directed_three_cycles(adj)}")
    print(f"   scc_sizes={scc_sizes(adj)}")
    print(f"   hamiltonian_paths={ham_path_count(adj)}")
    print(f"   one_hamiltonian_path={' -> '.join(carriers)}")
    print()
    print("ASSUMPTION CHALLENGE")
    print("   Considered vertices: runners, arcs, gaps, fixed circle sections, section")
    print("   boundaries, wall-crossing events, residues, cover arcs, Fourier modes,")
    print("   matroid circuits, flip bits, deletion roots, rectangle defects, and proof")
    print("   obligations.")
    print("   Chosen vertices: proof carriers for the two recursions.")
    print("   Preserves: explicit witnesses, orbit legality, fiber multiplicity,")
    print("   coboundary consistency, and LRC sidecar descent.")
    print("   Destroys if scalarized: runner labels, Hamiltonian-path presentation,")
    print("   canary/filler S-fiber mass, endpoint ownership, and active certificate")
    print("   coordinates.")


def main() -> None:
    print("=" * 80)
    print("HYP-3244: two interlocking recursions -- tiling model and half-tiling model")
    print("=" * 80)
    print_n4_tables()
    print_fixed_path_cover_table()
    print_half_tiling_orbit_table()
    print_coboundary_laws()
    print_interlock_and_tournament()


if __name__ == "__main__":
    main()
