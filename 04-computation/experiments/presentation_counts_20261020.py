"""presentation_counts_20261020.py -- opus-2026-10-20-S48
Presentation-count (tc) classification data from the class files written by
tiling_line_metagraph_n9_20261020.cpp (columns H, Aut, tc, bm, SC, type, minHblack, closed):
  * tc = 1 classes versus compositions of n into parts 1 and 3 (A000930) and their H = 3^j;
  * pendant nodes of the fundamental black graph (non-SC, tc = 1) versus non-palindromic {K1, C3}-stacks;
  * tc = 3 classes versus stacks over {K1, C3, T5} with exactly one T5;
  * the strong presentation-count spectrum (closed > 0 iff strong), its minimum, the minimum strong H;
  * the odd values missing from the tc spectrum and from the H spectrum (7 and 21);
  * H-descent sinks (non-pure-blue classes with no black neighbour of smaller H).
Usage: python3 presentation_counts_20261020.py [dir]   (dir holds tiling_line_metagraph_n9_20261020_classes_n{5..9}.tsv[.gz])
"""
import sys, os, gzip, collections
from math import comb

def a000930(n):
    a = [1, 1, 1]
    while len(a) <= n: a.append(a[-1] + a[-3])
    return a[n]

def palindromic_13(n):
    return a000930(n // 2) if n % 2 == 0 else a000930((n - 1) // 2) + a000930((n - 3) // 2)

def one_t5(n):
    r = n - 5
    if r < 0: return 0
    return sum(comb(r - 2 * j, j) * (r - 2 * j + 1) for j in range(0, r // 3 + 1))

def load(d, n):
    for fn in (f"tiling_line_metagraph_n9_20261020_classes_n{n}.tsv", f"tiling_line_metagraph_n9_20261020_classes_n{n}.tsv.gz"):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            op = gzip.open if p.endswith(".gz") else open
            rows = []
            with op(p, "rt") as f:
                head = f.readline().strip().split("\t")
                for line in f:
                    r = dict(zip(head, line.rstrip("\n").split("\t")))
                    for k in list(r):
                        if k != "type": r[k] = int(r[k])
                    rows.append(r)
            return rows
    return None

def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    for n in range(5, 10):
        rows = load(d, n)
        if rows is None: print(f"n={n}: no class file"); continue
        has_closed = "closed" in rows[0]; has_min = "minHblack" in rows[0]
        tc1 = [r for r in rows if r["tc"] == 1]
        pend = [r for r in tc1 if r["SC"] == 0]
        hs = collections.Counter(r["H"] for r in tc1)
        print(f"n={n}: classes {len(rows)}; tc=1 classes {len(tc1)} (A000930 = {a000930(n)}), H histogram {dict(sorted(hs.items()))}, "
              f"pendant (non-SC) {len(pend)} (= a - p = {a000930(n) - palindromic_13(n)}), palindromic {len(tc1) - len(pend)} (= p = {palindromic_13(n)})")
        tc3 = [r for r in rows if r["tc"] == 3]
        print(f"   tc=3 classes {len(tc3)} (stacks with one T5: {one_t5(n)}); (H,Aut) multiset {sorted(collections.Counter((r['H'], r['Aut']) for r in tc3).items())}")
        tcs = set(r["tc"] for r in rows); Hs = set(r["H"] for r in rows)
        print(f"   odd tc values < 60 absent: {[v for v in range(1, 60, 2) if v not in tcs]}; odd H values < 80 absent: {[v for v in range(1, 80, 2) if v not in Hs]}")
        if has_closed:
            strong = [r for r in rows if r["closed"] > 0]
            spec = collections.Counter(r["tc"] for r in strong)
            print(f"   strong classes {len(strong)}; strong tc spectrum (smallest) {sorted(spec.items())[:8]}; min strong tc {min(spec)}; min strong H {min(r['H'] for r in strong)}; strong H values <= 30 {sorted(set(r['H'] for r in strong if r['H'] <= 30))}")
            print(f"   strong classes with tc <= 9 (tc,H,Aut,SC): {sorted(collections.Counter((r['tc'], r['H'], r['Aut'], r['SC']) for r in strong if r['tc'] <= 9).items())}")
        if has_min:
            sinks = [r for r in rows if r["type"] != "PB" and (r["minHblack"] < 0 or r["minHblack"] >= r["H"])]
            print(f"   H-descent sinks {len(sinks)} (self-converse {sum(r['SC'] for r in sinks)}); sink H histogram {dict(sorted(collections.Counter(r['H'] for r in sinks).items())[:10])}")

if __name__ == "__main__":
    main()
