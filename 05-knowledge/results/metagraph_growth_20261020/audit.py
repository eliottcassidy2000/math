"""Replay both independent research lanes and check small hostile witnesses.

This file does not overwrite the retained lane outputs. Newline normalization
is explicit so a Windows checkout and a Git LF blob have the same receipt.
"""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise RuntimeError(message)


def normalized(data):
    return data.replace(b"\r\n", b"\n")


def replay(job):
    name, optimized = job
    command = [sys.executable, "-X", "utf8", "-B"]
    if optimized:
        command.append("-O")
    command.append(str(HERE / (name + ".py")))
    result = subprocess.run(command, capture_output=True, check=True)
    need(not result.stderr, f"{name}: unexpected stderr")
    return name, optimized, normalized(result.stdout)


def from_bits(n, bits):
    rows = [0] * n
    index = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if bits >> index & 1 else (j, i)
            rows[u] |= 1 << v
            index += 1
    return tuple(rows)


def canon(rows):
    n = len(rows)
    return min(tuple(sum(((rows[u] >> v) & 1) << j
                         for j, v in enumerate(order)) for u in order)
               for order in permutations(range(n)))


def converse(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ row ^ (1 << i) for i, row in enumerate(rows))


def flip(rows, u, v):
    rows = list(rows)
    rows[u] ^= 1 << v
    rows[v] ^= 1 << u
    return tuple(rows)


def path_tiling(n, forward_free):
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if j == i + 1 or forward_free else (j, i)
            rows[u] |= 1 << v
    return tuple(rows)


def small_audit():
    # Independent all-labelled n<=4 canon census, with full permutations.
    counts = []
    for n in range(1, 5):
        reps = {canon(from_bits(n, b)) for b in range(1 << (n * (n - 1) // 2))}
        counts.append(len(reps))
    need(counts == [1, 1, 2, 4], "small full-permutation census")

    # The old explainer forbade triangles with two SC/NS crossing edges.
    t = from_bits(4, 0)  # all larger vertices beat smaller vertices
    a = flip(t, 2, 0)
    b = flip(t, 3, 0)
    need(len({canon(t), canon(a), canon(b)}) == 3, "three distinct classes")
    need(canon(flip(a, 3, 2)) == canon(b), "third triangle edge")
    self_converse = [canon(x) == canon(converse(x)) for x in (t, a, b)]
    need(self_converse == [True, False, True], "two-rib triangle")

    # Gap graph isomorphism loses the ordered suspension, already at n=4.
    gap_forward, gap_backward = from_bits(3, 7), from_bits(3, 0)
    need(canon(gap_forward) == canon(gap_backward), "same gap class")
    lifted_forward, lifted_backward = path_tiling(4, True), path_tiling(4, False)
    need(canon(lifted_forward) != canon(lifted_backward), "different suspension classes")
    need(canon(lifted_forward) != canon(converse(lifted_backward)),
         "different merged suspension classes")
    return counts, self_converse


def main():
    jobs = [(name, optimized) for name in ("enumeration", "growth_check")
            for optimized in (False, True)]
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(replay, jobs))
    outputs = {(name, optimized): data for name, optimized, data in results}
    print("Exact packet audit; no lane output is rewritten.")
    for name in ("enumeration", "growth_check"):
        ordinary = outputs[(name, False)]
        need(ordinary == outputs[(name, True)], f"{name}: normal/-O mismatch")
        need(ordinary == normalized((HERE / (name + ".out")).read_bytes()),
             f"{name}: retained stdout mismatch")
        print(name, "normal = optimized = retained stdout;",
              "LF SHA256", sha256(ordinary).hexdigest())
    counts, sc = small_audit()
    print("Independent full-labelled/full-permutation counts n1..4:", counts)
    print("Minimal n4 two-rib triangle self-converse flags:", sc)
    print("Gap-graph hostile: two isomorphic transitive gap triples suspend to different merged n4 nodes.")
    print("PASS: scoped algebra/finite controls; no all-size compressed class recurrence asserted.")


if __name__ == "__main__":
    main()
