"""aud_compare.py - audit cross-checks for opus-2026-11-08-S86 (independent code; reads the session's OUTPUT files only
as data, never its scripts).

usage: python aud_compare.py WORKDIR SESSIONDIR S224CSV > aud_compare.out
  WORKDIR    holds census_<n>.txt (n = 1..9, written by `aud_fewcycle 0 2 n cls_n.txt` from kind-pasteur's classgen
             class lists), fc80.txt (`aud_fewcycle 80 0 13`), ic1023.txt (`aud_fewcycle 1023 1 12`),
             ic8191.txt (`aud_fewcycle 8191 1 0`)
  SESSIONDIR the session folder (few_cycle_strong_44.txt, few_cycle_strong_80.txt, intercyclic_strong_le13.txt)
  S224CSV    Session224's census/birth_fibers.csv
census line: n strong a1..a7 H sc c3=.. c5=.. c7=.. c9=.. code ; inventory line: n a1..a7 H Hdp sc code
"""
import os
import sys
from collections import Counter, defaultdict


def census(work, n):
    rows = []
    for line in open(os.path.join(work, f"census_{n}.txt")):
        t = line.split()
        cl = {k: int(v) for k, v in (x.split("=") for x in t[11:15])}
        rows.append(dict(n=int(t[0]), strong=int(t[1]), a=tuple(int(x) for x in t[2:9]), H=int(t[9]), sc=int(t[10]),
                         cl=cl, code=t[15]))
    return rows


def inv(path):
    rows = []
    for line in open(path):
        t = line.split()
        rows.append(dict(n=int(t[0]), a=tuple(int(x) for x in t[1:8]), H=int(t[8]), sc=int(t[10]), code=t[11]))
    return rows


def session_fcs(path):
    out = []
    for line in open(path):
        t = line.split()
        n = int(t[0])
        a5 = [int(x) for x in t[1:6]]
        H = int(t[6])
        r = H - 1 - sum(c << (k + 1) for k, c in enumerate(a5))
        assert r % 64 == 0
        out.append((n, tuple(a5 + [r // 64]), H))
    return out


def main():
    work, sdir, s224 = sys.argv[1], sys.argv[2], sys.argv[3]
    C = {n: census(work, n) for n in range(1, 10)}
    print("1. Census (kind-pasteur classgen class lists, every quantity recomputed by aud_fewcycle census mode;")
    print("   OCF H = sum alpha_k 2^k agreed with classgen's DP path count on every class, strong and self-converse flags agreed)")
    a000568 = [1, 1, 2, 4, 12, 56, 456, 6880, 191536]
    a051337 = [1, 0, 1, 1, 6, 35, 353, 6008, 178133]
    for n in range(1, 10):
        rows = C[n]
        strong = [r for r in rows if r["strong"]]
        codes = {r["code"] for r in rows}
        ic = [r for r in strong if r["a"][1] == 0]
        print(f"   order {n}: classes {len(rows)} (A000568 {a000568[n - 1]}), distinct own canonical codes {len(codes)}, "
              f"strong {len(strong)} (A051337 {a051337[n - 1]}), strong intercyclic {len(ic)}")
    print("   Moon's sharp count: min number of l-cycles over strong classes, and min alpha_1 vs m(n):")
    for n in range(3, 10):
        strong = [r for r in C[n] if r["strong"]]
        mins = {l: min(r["cl"][f"c{l}"] for r in strong) for l in (3, 5, 7, 9) if l <= n}
        m = sum(n - l + 1 for l in range(3, n + 1, 2))
        print(f"     n = {n}: " + ", ".join(f"min c{l} = {v} (n-l+1 = {n - l + 1})" for l, v in mins.items())
              + f"; min alpha_1 = {min(r['a'][0] for r in strong)} (m(n) = {m})")

    FC = inv(os.path.join(work, "fc80.txt"))
    print()
    print("2. Own few-cycle inventory (strong, alpha_1 <= 80; level 19 EMPTY):", len(FC), "classes; per order",
          sorted(Counter(r["n"] for r in FC).items()))
    for n in range(3, 10):
        a = {r["code"] for r in C[n] if r["strong"] and r["a"][0] <= 80}
        b = {r["code"] for r in FC if r["n"] == n}
        print(f"   order {n}: census strong alpha_1 <= 80: {len(a)}; enumerator: {len(b)}; identical code sets: {a == b}")
    own = Counter((r["n"], r["a"][:6], r["H"]) for r in FC)
    for budget, fname in ((44, "few_cycle_strong_44.txt"), (80, "few_cycle_strong_80.txt")):
        ses = Counter(session_fcs(os.path.join(sdir, fname)))
        mine = Counter({k: v for k, v in own.items() if k[1][0] <= budget})
        print(f"   session {fname}: {sum(ses.values())} classes; multiset of (order, packing vector, H) equal to own "
              f"alpha_1 <= {budget} part: {ses == mine}")
    print("   classes with alpha_6 > 0:", [(r["n"], r["a"][:6], r["H"]) for r in FC if r["a"][5] > 0])

    IC = inv(os.path.join(work, "ic1023.txt"))
    ses_ic = Counter()
    for line in open(os.path.join(sdir, "intercyclic_strong_le13.txt")):
        t = line.split()
        if int(t[1]) < 1024:
            ses_ic[(int(t[0]), int(t[1]), int(t[4]))] += 1
    own_ic = Counter((r["n"], r["a"][0], r["sc"]) for r in IC)
    print()
    print("3. Own intercyclic inventory (strong, alpha_1 <= 1023; level 14 EMPTY):", len(IC), "classes; per order",
          sorted(Counter(r["n"] for r in IC).items()))
    print("   equal to the session's intercyclic_strong_le13.txt restricted to alpha_1 < 1024 as a multiset of "
          "(order, alpha_1, self-converse):", own_ic == ses_ic)

    mult = Counter(r["a"][0] for r in IC)
    singles = sorted(a for a in mult if mult[a] == 1)
    print()
    print("4. Carried twins in base 2 for the carry-free singletons a >= 8 below 1024 (strong NON-intercyclic classes")
    print("   with H = 1 + 2a, found in the census of orders 6-9):")
    allok = True
    for a in singles:
        if a < 8:
            continue
        H = 1 + 2 * a
        hits = [(n, sum(1 for r in C[n] if r["strong"] and r["a"][1] >= 1 and r["H"] == H)) for n in range(6, 10)]
        tot = sum(c for n, c in hits)
        first = next((n for n, c in hits if c), None)
        allok &= tot > 0
        ex = next((r for n in range(6, 10) for r in C[n] if r["strong"] and r["a"][1] >= 1 and r["H"] == H), None)
        print(f"   a = {a:4d}, H = {H:5d}: twins by order {hits}; least-order example alpha = {ex['a'][:4] if ex else None}")
    print("   every carry-free singleton with 8 <= a < 1024 has a strong non-intercyclic twin:", allok)

    IC2 = inv(os.path.join(work, "ic8191.txt"))
    m2 = Counter(r["a"][0] for r in IC2)
    sc2 = {}
    for r in IC2:
        sc2.setdefault(r["a"][0], []).append((r["n"], r["sc"]))
    s2 = sorted(a for a in m2 if m2[a] == 1)
    print()
    print("5. Extended intercyclic run (alpha_1 <= 8191; level 17 EMPTY):", len(IC2), "classes")
    print(f"   carry-free singletons below 8192: {len(s2)} ({len([a for a in s2 if a < 1024])} below 1024); "
          f"those in [1024, 8192): {[a for a in s2 if a >= 1024]}")
    print("   orders of all singletons:", sorted(Counter(sc2[a][0][0] for a in s2).items()),
          "; all self-converse:", all(sc2[a][0][1] == 1 for a in s2))

    print()
    print("6. Birth fibers h(s,p) for n(h) <= 7 recomputed from the census vs Session224 birth_fibers.csv:")
    first = {}
    sp = defaultdict(lambda: [0, 0])
    for n in range(1, 8):
        for r in C[n]:
            first.setdefault(r["H"], n)
            if first[r["H"]] == n:
                sp[r["H"]][0 if r["sc"] else 1] += 1
    mine = {(h, first[h], sp[h][0], sp[h][1] // 2) for h in first}
    theirs = set()
    for line in open(s224):
        if line.startswith("H,"):
            continue
        t = line.strip().split(",")
        if int(t[1]) <= 7:
            theirs.add((int(t[0]), int(t[1]), int(t[3]), int(t[4])))
    print(f"   {len(mine)} values; identical: {mine == theirs}")
    for n in range(4, 7):
        print(f"   order {n}: " + " ".join(f"{h}({s},{p})" for (h, f, s, p) in sorted(mine) if f == n))
    print("   order 7: %d values" % sum(1 for x in mine if x[1] == 7))
    h11 = [(r["strong"], r["sc"]) for r in C[6] if r["H"] == 11]
    print("   order-6 classes with H = 11 (strong, self-converse):", h11)


if __name__ == "__main__":
    main()
