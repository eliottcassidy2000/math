#!/usr/bin/env python3
"""
chessboard_weave_20261006_hp_blocking_structured.py

Structured families for HYP-9168 (beta = hall), C engine mode bnb with root symmetry breaking
(full automorphism group from networkx).  Part 3c of chessboard_weave_20261006_hp_blocking_run.py
(that run was stopped after part 3b for time):
  * all 32 circulant tournaments on Z_11 (regular: hall = 10 is forced),
  * circulant tournaments on Z_13: one representative of each of the 6 orbits of the 64
    connection sets under multiplication by Z_13^* (S -> aS is an isomorphism, S -> -S gives
    the converse, which has the same beta and hall), so all 64 are covered,
  * Paley QR_11 and QR_11 minus a vertex,
  * transitive TT_N and TT_N with one reversed arc (all choices), N = 10, 11, 12,
  * compositions h -> A => B -> h (hubs h transitive; A, B cyclic triangles / regular 5-, 7-
    tournaments / transitive), N <= 14, node limit 3e6 per instance (ABORTED is reported).
"""
import importlib.util
import itertools
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("run", os.path.join(HERE, "chessboard_weave_20261006_hp_blocking_run.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def batch(name, lines, limit=None):
    t = time.time()
    res = R.engine("bnb", lines, [str(limit)] if limit else [])
    eq = sum(1 for ln in res if ln.endswith(" EQ"))
    ab = [ln for ln in res if "ABORTED" in ln]
    ce = [ln for ln in res if "COUNTEREXAMPLE" in ln or "WITNESS_FAIL" in ln]
    ob = sum(1 for ln in res if R.field(ln, "obstruction_ok") == "1")
    halls = {}
    for ln in res:
        h = int(R.field(ln, "hall"))
        halls[h] = halls.get(h, 0) + 1
    print(f"{name}: {len(lines)} tournaments; beta == hall proved: {eq}; aborted: {len(ab)}; counterexamples: {len(ce)}; "
          f"obstruction verified: {ob}/{len(lines)}; hall histogram {dict(sorted(halls.items()))}; {time.time() - t:.1f}s")
    for ln in ab + ce:
        print("   " + ln.split()[0] + " " + " ".join(x for x in ln.split()[1:] if not x.startswith("A:")))
    sys.stdout.flush()


def main():
    R.build()
    print(__doc__.strip())
    print()
    lines = [R.line_with_auts(11, R.circulant(11, S)) for S in itertools.product(*[(s, -s) for s in range(1, 6)])]
    batch("all 32 circulants on Z_11", lines)
    reps, seen = [], set()
    for S in itertools.product(*[(s % 13, -s % 13) for s in range(1, 7)]):
        fs = frozenset(S)
        if fs in seen:
            continue
        orbit = {frozenset((a * x) % 13 for x in fs) for a in range(1, 13)}
        seen |= orbit
        reps.append(sorted(fs))
    print(f"Z_13: {len(seen)} connection sets in {len(reps)} multiplier orbits; representatives {reps}")
    batch("circulants on Z_13 (orbit representatives)", [R.line_with_auts(13, R.circulant(13, S)) for S in reps])
    qr11 = R.circulant(11, [1, 3, 4, 5, 9])
    qr11mv = {(u - 1, v - 1) for (u, v) in qr11 if u != 0 and v != 0}
    batch("Paley QR_11", [R.line_with_auts(11, qr11)])
    batch("QR_11 minus a vertex", [R.line_with_auts(10, qr11mv)])
    lines = []
    for N in (10, 11, 12):
        TT = R.transitive(N)
        lines.append(R.to_string(N, TT))
        for (i, j) in sorted(TT):
            lines.append(R.to_string(N, (TT - {(i, j)}) | {(j, i)}))
    batch("TT_N and TT_N with one reversed arc (all), N = 10, 11, 12", lines)
    C3 = R.circulant(3, [1])
    blocks = {"C3": (3, C3), "R5": (5, R.circulant(5, [1, 2])), "R7": (7, R.circulant(7, [1, 2, 4])),
              "TT1": (1, set()), "TT2": (2, {(0, 1)}), "TT3": (3, R.transitive(3))}
    lines = []
    for A in ("C3", "R5", "R7", "TT2", "TT3"):
        for B in ("C3", "R5", "R7", "TT1", "TT2", "TT3"):
            for hubs in (1, 2, 3):
                nA, aA = blocks[A]
                nB, aB = blocks[B]
                N = hubs + nA + nB
                if not (7 <= N <= 14):
                    continue
                arcs = {(i, j) for i in range(hubs) for j in range(i + 1, hubs)}
                arcs |= {(u + hubs, v + hubs) for (u, v) in aA}
                arcs |= {(u + hubs + nA, v + hubs + nA) for (u, v) in aB}
                for h in range(hubs):
                    arcs |= {(h, a) for a in range(hubs, hubs + nA)}
                    arcs |= {(b, h) for b in range(hubs + nA, N)}
                arcs |= {(a, b) for a in range(hubs, hubs + nA) for b in range(hubs + nA, N)}
                lines.append(R.line_with_auts(N, arcs))
    batch("compositions h -> A => B -> h (1-3 hubs; A, B in C3, R5, R7, TT1-3; N = 7..14)", lines, limit=3000000)


if __name__ == "__main__":
    main()
