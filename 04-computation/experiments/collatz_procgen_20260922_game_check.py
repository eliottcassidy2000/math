#!/usr/bin/env python3
"""Independent checks for Althofer's 3n+-1 game (= Conway's Beans-Don't-Talk).

Game: position = positive integer n (odd in Althofer's version; Guy/OEIS also allow even starts).
A move replaces n by oddpart(3n+1) or oddpart(3n-1); the player who moves to 1 wins.

This script shares NO code with the C programs.  It provides
  (A) an independent capped retrograde with Steinhaus remoteness (FIFO BFS by remoteness, per-node
      child counters), run on ALL odd n <= C (multiples of 3 included as ordinary nodes);
  (B) comparison of (A) with the 2-bit table dumped by collatz_procgen_20260922_game_dense.c at the
      SAME cap (both must be the least fixpoint of the same rules on the same finite graph);
  (C) comparison of remoteness layers with OEIS A005694-A005698 (Guy 1986, Beans-Don't-Talk),
      including even starts, and a complete check of the remoteness-2 list A005698 up to its last term
      via the Jacobsthal characterization;
  (D) a checker for proof certificates written by collatz_procgen_20260922_game_cert.c.

usage:
  game_check.py bfs C [dumpfile]          # (A)+(B)+(C)
  game_check.py cert file [file ...]      # (D)
  game_check.py certcmp certfile blocks_output   # (D') compare sampled block-run values
"""
import sys
from collections import deque


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


def children(n):
    return oddpart(3 * n + 1), oddpart(3 * n - 1)


# ---------------------------------------------------------------- (A) BFS retrograde with remoteness
def bfs_remoteness(C):
    """rem[i] for odd n=2i+1 <= C: -1 unresolved; even = P (mover loses), odd = N.
    The terminal state 1 (a player has just moved to 1) has remoteness 0; the START n=1 is handled
    separately (remoteness 1)."""
    M = (C - 1) // 2 + 1
    rem = [-1] * M
    cnt = [2] * M          # children not yet known to be N (children of n>=3 are distinct)
    q = deque()
    rem[0] = 0             # terminal
    q.append(1)
    while q:
        x = q.popleft()
        r = rem[(x - 1) // 2]
        if x % 3 == 0:
            continue       # multiples of 3 are never children
        # predecessors p with a move p -> x: 3p +- 1 = x*2^k, k >= 1
        y = x
        while True:
            y *= 2
            if (y - 1) // 3 > C:
                break
            if y % 3 == 1:
                p = (y - 1) // 3      # 3p+1 = y
            else:
                p = (y + 1) // 3      # 3p-1 = y
            if p > C:
                break
            if p == 1 or p == x:
                continue
            j = (p - 1) // 2
            if rem[j] >= 0:
                continue
            if r % 2 == 0:            # x is P (or terminal): p wins by moving to x
                rem[j] = r + 1
                q.append(p)
            else:                     # x is N
                cnt[j] -= 1
                if cnt[j] == 0:       # both children N; x is the later one => max remoteness
                    rem[j] = r + 1
                    q.append(p)
    return rem


def value_of_start(n, rem, C):
    """remoteness of an arbitrary positive start n (odd or even) from the odd table; None if unknown."""
    if n == 1:
        return 1
    a, b = children(n)
    ra = 0 if a == 1 else (rem[(a - 1) // 2] if a <= C else -1)
    rb = 0 if b == 1 else (rem[(b - 1) // 2] if b <= C else -1)
    evens = [r for r in (ra, rb) if r >= 0 and r % 2 == 0]
    if evens:
        return 1 + min(evens)
    if ra >= 0 and rb >= 0:
        return 1 + max(ra, rb)
    return None


def run_bfs(C, dumpfile=None):
    if C % 2 == 0:
        C -= 1
    rem = bfs_remoteness(C)
    M = len(rem)
    res = sum(1 for i in range(1, M) if rem[i] >= 0)
    nP = sum(1 for i in range(1, M) if rem[i] >= 0 and rem[i] % 2 == 0)
    first = next((2 * i + 1 for i in range(1, M) if rem[i] < 0), None)
    print(f"(A) BFS cap C={C}: odd 3<=n<=C resolved={res} (N={res - nP}, P={nP}), "
          f"unresolved={M - 1 - res}; first unresolved odd start={first}")
    # remoteness statistics on the fully resolved prefix
    if first:
        mx = max(range(1, (first - 1) // 2), key=lambda i: rem[i])
        print(f"    max remoteness among odd n < {first}: {rem[mx]} at n={2 * mx + 1} "
              f"(capped remoteness = upper bound; exact when optimal lines stay <= C)")
    # (B) compare with the C dump
    if dumpfile:
        import struct
        with open(dumpfile, 'rb') as f:
            cap, dumpn = struct.unpack('<QQ', f.read(16))
            data = f.read()
        assert cap == C, (cap, C)
        words = len(data) // 8

        def gs(idx):
            w = int.from_bytes(data[8 * (idx >> 5): 8 * (idx >> 5) + 8], 'little')
            return (w >> ((idx & 31) * 2)) & 3
        mism = 0
        checked = 0
        for i in range(1, M):
            n = 2 * i + 1
            if n > dumpn:
                break
            if n % 3:
                s = gs(n // 3)
                cv = 'N' if s == 2 else ('P' if s == 3 else 'U')
            else:
                a, b = children(n)
                sa = 3 if a == 1 else gs(a // 3)
                sb = gs(b // 3) if b <= C else 0
                cv = 'N' if (sa == 3 or sb == 3) else ('P' if (sa == 2 and sb == 2) else 'U')
            pv = 'U' if rem[i] < 0 else ('P' if rem[i] % 2 == 0 else 'N')
            checked += 1
            if cv != pv:
                mism += 1
                if mism <= 10:
                    print("    MISMATCH", n, "C:", cv, "python:", pv)
        print(f"(B) compared with C dump ({dumpfile}, cap {cap}, n<={dumpn}): {checked} odd positions, "
              f"{mism} mismatches")
    # (C) OEIS remoteness layers (terms as listed on oeis.org, fetched 2026-09-22)
    oeis = {
        2: ("A005698", [7, 29, 57, 227, 455, 1821, 3641, 14563, 29127, 116509, 233017, 932067, 1864135,
                        7456541, 14913081, 59652323, 119304647, 477218589, 954437177, 3817748707,
                        7635497415, 30541989661, 61083979321, 244335917283, 488671834567,
                        1954687338269]),
        3: ("A005695", [2, 9, 10, 19, 37, 39, 75, 76, 77, 149, 151, 152, 155, 299, 303, 309, 597, 605,
                        607, 619, 1195, 1211, 1213, 1214, 1237, 2389, 2421, 2427, 2475, 4779, 4843,
                        4853, 4854, 4855, 4949, 9557, 9685, 9707, 9709, 9899, 19115, 19371, 19413,
                        19417, 19419, 19797]),
        4: ("A005696", [13, 25, 50, 51, 99, 101, 103, 199, 202, 403, 404, 405, 413, 797, 807, 809, 825,
                        1593, 1618, 3229, 3235, 3236, 3237, 3299, 6371, 6457, 6471, 6473, 6599, 12743,
                        12945, 25827, 25885, 25891, 26397, 50973, 51655, 51769, 51779, 51783, 52793]),
        5: ("A005697", [4, 8, 17, 33, 34, 35, 66, 67, 69, 133, 134, 135, 137, 138, 139, 265, 266, 267,
                        269, 270, 275, 277, 531, 533, 537, 539, 549, 551, 555, 1061, 1063, 1067, 1075,
                        1076, 1077, 1078, 1079, 1099, 1100, 1101, 1109, 2123, 2124, 2125, 2133, 2149,
                        2152, 2153]),
        6: ("A005694", [6, 12, 23, 45, 46, 89, 91, 92, 93, 177, 179, 183, 185, 354, 355, 359, 367, 707,
                        708, 709, 711, 717, 718, 719, 733, 739, 1415, 1417, 1433, 1435, 1437, 1438, 1465,
                        1469, 1479, 2831, 2845, 2870, 2873, 2875, 2876, 2877, 2933, 2937, 5661, 5663,
                        5667, 5689]),
    }
    for r, (anum, terms) in sorted(oeis.items()):
        top = terms[-1]
        if r == 2:
            top = min(top, C // 12)   # direct table check only where lines of length 2 stay <= C
        ours = [n for n in range(2, top + 1) if value_of_start(n, rem, C) == r]
        want = [t for t in terms if t <= top]
        print(f"(C) remoteness {r} ({anum}) all n in [2,{top}]: ours {len(ours)} terms, OEIS {len(want)} "
              f"terms, equal={ours == want}")
    # remoteness 2, complete up to the last A005698 term, via Jacobsthal characterization
    J = set()
    for k in range(1, 50):     # 3j +- 1 = 2^k  <=> j Jacobsthal (immediate win)
        for s in (1, -1):
            if (2 ** k - s) % 3 == 0:
                j = (2 ** k - s) // 3
                if j >= 1 and j % 2 == 1:
                    J.add(j)
    J1 = sorted(J)
    # rem 2 <=> neither child is 1 and both children are in J (n odd: children (3n+-1)/2^*)
    # even n: children 3n+-1 differ by 2, both Jacobsthal impossible (checked below by search)
    cand = set()
    # every n with a child j in J: 3n +- 1 = j 2^k
    for j in J1:
        for k in range(0, 45):
            for s in (1, -1):
                y = j * 2 ** k
                if (y - s) % 3 == 0:
                    n = (y - s) // 3
                    if 2 <= n <= 1954687338269:
                        cand.add(n)
    rem2 = sorted(n for n in cand if n > 1 and all(c in J and c != 1 for c in children(n)))
    print(f"(C') remoteness-2 positions n<=1954687338269 (Jacobsthal characterization): {len(rem2)} terms; "
          f"equal to A005698 data: {rem2 == oeis[2][1]}")


# ---------------------------------------------------------------- (D) certificate checker
def check_cert_against(path, blocks_output):
    """(D') validate certificate `path` and compare its target values with the 'SAMPLE n v' lines
    printed by collatz_procgen_20260922_game_blocks.c."""
    want = {}
    for line in open(blocks_output):
        if line.startswith('SAMPLE '):
            _, n, v = line.split()
            want[int(n)] = v
    val = check_cert(path, quiet=True)
    bad = [n for n in want if val.get(n) != want[n]]
    print(f"(D') {path}: {len(want)} sampled starts from {blocks_output}; certificate values agree for "
          f"{len(want) - len(bad)}, disagree for {len(bad)} {bad[:5]}")


def check_cert(path, quiet=False):
    """Certificate format (text): header lines start with '#'; then lines 'n V c' in topological order:
    V='P': both children of n must already be listed with V='N';
    V='N': c must be a child of n and either c==1 or c already listed with V='P'.
    Targets are listed in a '# targets:' header line."""
    val = {}
    targets = []
    lines = 0
    maxnode = 0
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('# targets:'):
                    targets = [int(t) for t in line.split(':', 1)[1].split()]
                continue
            parts = line.split()
            n = int(parts[0])
            v = parts[1]
            a, b = children(n)
            if n in val:
                raise SystemExit(f"{path}: duplicate node {n}")
            if v == 'N':
                c = int(parts[2])
                if c not in (a, b):
                    raise SystemExit(f"{path}: {c} is not a child of {n}")
                if not (c == 1 or val.get(c) == 'P'):
                    raise SystemExit(f"{path}: N-node {n}: child {c} not an earlier P")
            elif v == 'P':
                if a == 1 or b == 1 or val.get(a) != 'N' or val.get(b) != 'N':
                    raise SystemExit(f"{path}: P-node {n}: children {a},{b} not both earlier N")
            else:
                raise SystemExit(f"{path}: bad value {v}")
            val[n] = v
            lines += 1
            maxnode = max(maxnode, n)
    for t in targets:
        if t not in val:
            raise SystemExit(f"{path}: target {t} not certified")
    tv = ", ".join(f"{t}={val[t]}" for t in targets[:8]) + (" ..." if len(targets) > 8 else "")
    print(f"(D) {path}: VALID certificate, {lines} nodes, {len(targets)} targets, largest node {maxnode}; {tv}")
    return val


if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == 'bfs':
        run_bfs(int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
    elif len(sys.argv) == 4 and sys.argv[1] == 'certcmp':
        check_cert_against(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 3 and sys.argv[1] == 'cert':
        for p in sys.argv[2:]:
            check_cert(p)
    else:
        print(__doc__)
