#!/usr/bin/env python3
"""AUDIT: THM-5070 (6) / HYP-9640 (E), bases 4..8 to alpha-sum 30,000.

1. Chain misses come from MY coverage run (a4_cover on my own order-<= 5 matrices, out/a4_cover_bases4to8.out).
2. Attained alpha-sums: [0, 30000] minus the chain misses, plus every m <= 80 outside THM-4940 (6)'s gap table.
3. Closure under order joins, done on VALUES: (1 + b m)(1 + b m') stays attained (H_b multiplicative over joins).
4. The remaining alpha-sums in (80, 30000] are compared with the packet's list, and every packet witness
   (out/bases_4to8_witnesses.txt) is rechecked with my C brute force a1_brute (tournament test, I(Omega, b), OCF).
Usage: a9_bases.py MY_COVER_OUT PACKET_WITNESSES BRUTE_EXE
"""
import re
import subprocess
import sys

M = 30000
GAPS80 = {4: {3, 14, 27}, 5: {3, 10, 27}, 6: {3, 10, 14, 21, 46}, 7: {3, 10, 14, 29}, 8: {3, 35, 52}}  # THM-4940 (6)


def main():
    cover, wfile, brute = sys.argv[1], sys.argv[2], sys.argv[3]
    chain_miss = {}
    b = None
    for line in open(cover):
        mm = re.match(r"== base (\d+)", line)
        if mm:
            b = int(mm.group(1))
        if line.startswith("misses"):
            vals = list(map(int, line.split(":")[1].split()))
            chain_miss[b] = sorted((v - 1) // b for v in vals)
    rem_all = {}
    for b in range(4, 9):
        att = set(range(M + 1)) - set(chain_miss[b])
        att |= set(range(81)) - GAPS80[b]
        LIM = 1 + b * M
        vals = sorted(1 + b * m for m in att)
        S = set(vals)
        changed = True
        while changed:
            changed = False
            cur = sorted(S)
            for i, u in enumerate(cur):
                if u == 1:
                    continue
                if u * u > LIM:
                    break
                for v in cur[i:]:
                    p = u * v
                    if p > LIM:
                        break
                    if p not in S:
                        S.add(p)
                        changed = True
        closed = {(v - 1) // b for v in S}
        rem = [m for m in range(81, M + 1) if m not in closed]
        rem_all[b] = rem
        print(f"base {b}: my chain misses {chain_miss[b]}")
        print(f"base {b}: open after closure (80, {M}]: {rem}")
    # witnesses
    lines, meta = [], []
    for line in open(wfile):
        if "FOUND" not in line:
            continue
        b = int(line.split()[1])
        target = int(re.search(r"target=(\d+)", line).group(1))
        n = int(re.search(r"n=(\d+)", line).group(1))
        rows = line.split("rows")[1].split()
        assert len(rows) == n and all(len(r) == n for r in rows)
        A = [[int(c) for c in r] for r in rows]
        tour = all(A[i][i] == 0 for i in range(n)) and all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n)
                                                           if i != j)
        lines.append(f"{n} 0 " + "".join(rows))
        meta.append((b, target, n, tour))
    out = subprocess.run([brute], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split("\n")
    ok = 0
    covered = {b: set() for b in range(4, 9)}
    for (b, target, n, tour), res in zip(meta, out):
        v = list(map(int, res.split()))
        H, hx = v[0], {x: v[1 + 3 * (x - 2)] for x in range(2, 9)}
        good = tour and hx[b] == target and H == hx[2] and (target - 1) % b == 0
        ok += good
        if good:
            covered[b].add((target - 1) // b)
        else:
            print("WITNESS FAILS", b, target, n, tour, hx[b], H, hx[2])
    print(f"witnesses: {len(meta)} lines, {ok} verified (tournament, I(Omega,b) = target, OCF), orders "
          f"{min(t[2] for t in meta)}..{max(t[2] for t in meta)}")
    for b in range(4, 9):
        left = sorted(set(rem_all[b]) - covered[b])
        extra = sorted(covered[b] - set(rem_all[b]))
        print(f"base {b}: open values {len(rem_all[b])}, covered by verified witnesses {len(covered[b] & set(rem_all[b]))};"
              f" still open: {left}; witnesses for values not open: {extra}")


if __name__ == "__main__":
    main()
