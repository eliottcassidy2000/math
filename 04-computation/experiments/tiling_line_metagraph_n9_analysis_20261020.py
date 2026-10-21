"""tiling_line_metagraph_n9_analysis_20261020.py -- opus-2026-10-20-S47
Class-level statistics from tiling_line_metagraph_n9_20261020_classes_n{n}.tsv (output of the C++ census):
fundamental black graph degrees (THM-4731 (3)), the handshake law, pure-blue nodes with their blue degrees,
node-type counts, and the self-line classes.  Usage: python3 ... [n]
"""
import sys, collections

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    rows = []
    with open(f"tiling_line_metagraph_n9_20261020_classes_n{n}.tsv") as f:
        head = f.readline().strip().split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append(dict(zip(head, p)))
    for r in rows:
        for k in ("class", "H", "Aut", "tc", "bm", "SC", "selfBlue", "selfBlack", "blackComponent"):
            r[k] = int(r[k])
    T = len(rows); SC = sum(r["SC"] for r in rows)
    types = collections.Counter(r["type"] for r in rows)
    print(f"n={n}: classes {T}, SC {SC}, merged {(T + SC) // 2}, types {dict(types)}")
    # fundamental black graph B: degree tc at NS merged nodes (pair counted once), (tc-bm)/2 at SC nodes
    odd_sc = [r for r in rows if r["SC"] and ((r["tc"] - r["bm"]) // 2) % 2 == 1]
    print(f"handshake: SC classes with (tc-bm)/2 odd: {len(odd_sc)}; (T-SC)/2 = {(T - SC) // 2} -> parities {len(odd_sc) % 2} vs {((T - SC) // 2) % 2}")
    degB = collections.Counter()
    for r in rows:
        if r["SC"]: degB[(r["tc"] - r["bm"]) // 2] += 1
        else: degB[r["tc"]] += 1   # each NS class; the merged node has the same degree (twin halves)
    small = sorted(degB.items())[:12]
    print(f"fundamental black degree distribution (smallest 12 values): {small}; max degree {max(degB)}")
    pb = [r for r in rows if r["type"] == "PB"]
    print("pure-blue classes (H, Aut, tc = blue degree):", [(r["H"], r["Aut"], r["tc"]) for r in pb], "total blue degree", sum(r["tc"] for r in pb))
    sl = [r for r in rows if r["selfBlack"] > 0]
    print(f"classes carrying black self-lines: {len(sl)}; black self-lines {sum(r['selfBlack'] for r in sl)}; blue self-lines {sum(r['selfBlue'] for r in rows)}")
    print("self-line class multiplicity histogram:", sorted(collections.Counter(r["selfBlack"] for r in sl).items()))
    comps = collections.Counter(r["blackComponent"] for r in rows)
    print(f"black components: {len(comps)}; sizes {sorted(comps.values(), reverse=True)[:8]}")
    # H and Aut extremes
    print("max H class:", max(rows, key=lambda r: r["H"])["H"], "; max |Aut|:", max(r["Aut"] for r in rows), "; rigid classes:", sum(r["Aut"] == 1 for r in rows))
    print("tc distribution (smallest 8):", sorted(collections.Counter(r["tc"] for r in rows).items())[:8])

if __name__ == "__main__":
    main()
