#!/usr/bin/env python3
"""
higher_cycle_copart.py -- Co-participation in k-cycles by resonance level

KEY QUESTION: The co-participation theorem (HYP-601) shows that for Paley T_p,
chord pairs with resonance level q=3 have NO 3-cycle co-participation at some
primes (e.g., p=19). But the product law sign(h_hat) = chi(ab) must still hold.

HYPOTHESIS: q=3 pairs compensate through HIGHER-order cycle co-participation
(5-cycles, 7-cycles, etc.). The OCF H = sum 2^j * alpha_j aggregates all cycle
lengths, so the sign can flip if alpha_j (j>=2) contributions overcome alpha_1.

This script computes co-participation in k-cycles (k=3,5,7) as a function of
resonance level q, to see if there's a length-level tradeoff.

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def directed_cycles(adj, n, max_length):
    """Find all directed cycles up to max_length in adjacency matrix adj.
    Returns dict: length -> list of vertex tuples (canonical form)."""
    cycles = defaultdict(set)

    for start in range(n):
        # DFS from start, tracking path
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, mask = stack.pop()
            if len(path) > max_length:
                continue
            for w in range(n):
                if not adj[v][w]:
                    continue
                if w == start and len(path) >= 3:
                    # Found a cycle
                    clen = len(path)
                    # Canonical: minimum rotation
                    min_rot = min(path[i:] + path[:i] for i in range(clen))
                    cycles[clen].add(tuple(min_rot))
                elif not (mask & (1 << w)) and w > start:
                    # Only go to vertices > start to avoid double counting
                    stack.append((w, path + [w], mask | (1 << w)))

    return cycles


def gap_chords_of_cycle(cycle, p):
    """For a directed cycle, compute the chord indices (gap values)."""
    chords = set()
    n = len(cycle)
    for i in range(n):
        gap = (cycle[(i+1) % n] - cycle[i]) % p
        chords.add(min(gap, p - gap))
    return chords


def main():
    print("=" * 70)
    print("HIGHER-CYCLE CO-PARTICIPATION BY RESONANCE LEVEL")
    print("=" * 70)

    for p in [7, 11, 19]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        QR_set = set(QR)

        print(f"\n{'='*60}")
        print(f"p = {p}, m = {m}, QR = {QR}")
        print(f"{'='*60}")

        # Build Paley adjacency matrix
        adj = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in QR:
                adj[v][(v + s) % p] = 1

        # Find directed cycles of length 3, 5 (and 7 for small p)
        max_k = 7 if p <= 11 else 5
        print(f"\nFinding directed cycles up to length {max_k}...")

        all_cycles = directed_cycles(adj, p, max_k)
        for k in sorted(all_cycles.keys()):
            print(f"  {k}-cycles: {len(all_cycles[k])}")

        # For each cycle, extract chord pairs
        # chord_copart[k][(a,b)] = # of k-cycles using both chords a and b
        for k in sorted(all_cycles.keys()):
            chord_copart = defaultdict(int)
            for cycle in all_cycles[k]:
                chords = gap_chords_of_cycle(cycle, p)
                chord_list = sorted(chords)
                for i in range(len(chord_list)):
                    for j in range(i+1, len(chord_list)):
                        chord_copart[(chord_list[i], chord_list[j])] += 1

            print(f"\n  --- {k}-cycle co-participation by chord pair ---")
            # Group by resonance level q
            q_groups = defaultdict(list)
            for a in range(1, m+1):
                for b in range(a+1, m+1):
                    res = classify_resonance(a, b, p)
                    if not res:
                        continue
                    min_q = min(qq for qq, t in res)
                    cp = chord_copart.get((a,b), 0)
                    q_groups[min_q].append((a, b, cp))

            for q in sorted(q_groups.keys()):
                items = q_groups[q]
                vals = [cp for a, b, cp in items]
                chi_q = legendre(q, p)
                # 3-cycle copart from HYP-601
                if chi_q == 1:
                    c3_copart = legendre(q+1, p) == -1
                else:
                    c3_copart = legendre(q-1, p) == 1

                avg_cp = sum(vals) / len(vals) if vals else 0
                unique_vals = sorted(set(vals))

                print(f"    q={q:>2} (chi={chi_q:+d}, 3c-copart={'Y' if c3_copart else 'N'}): "
                      f"avg={avg_cp:>6.1f}, vals={unique_vals}, "
                      f"pairs={[(a,b) for a,b,c in items]}")

        # Summary table
        print(f"\n  --- SUMMARY: copart by (q, cycle-length) ---")
        q_vals = sorted(set(min(qq for qq, t in classify_resonance(a, b, p))
                          for a in range(1, m+1) for b in range(a+1, m+1)
                          if classify_resonance(a, b, p)))

        header = f"    {'q':>3} {'chi(q)':>6}"
        for k in sorted(all_cycles.keys()):
            header += f" {'c'+str(k)+'-copart':>12}"
        print(header)

        for q in q_vals:
            chi_q = legendre(q, p)
            row = f"    {q:>3} {chi_q:>+6}"
            for k in sorted(all_cycles.keys()):
                chord_copart = defaultdict(int)
                for cycle in all_cycles[k]:
                    chords = gap_chords_of_cycle(cycle, p)
                    chord_list = sorted(chords)
                    for i in range(len(chord_list)):
                        for j in range(i+1, len(chord_list)):
                            chord_copart[(chord_list[i], chord_list[j])] += 1

                vals = []
                for a in range(1, m+1):
                    for b in range(a+1, m+1):
                        res = classify_resonance(a, b, p)
                        if not res:
                            continue
                        min_q = min(qq for qq, t in res)
                        if min_q == q:
                            vals.append(chord_copart.get((a,b), 0))

                if vals:
                    avg = sum(vals) / len(vals)
                    row += f" {avg:>12.1f}"
                else:
                    row += f" {'N/A':>12}"
            print(row)

    # PART 2: At p=19, deeper analysis of q=3 pairs
    print("\n\n" + "=" * 70)
    print("PART 2: q=3 PAIRS AT p=19 -- CYCLE LENGTH COMPENSATION?")
    print("=" * 70)

    p = 19
    m = (p - 1) // 2
    QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    QR_set = set(QR)

    adj = [[0]*p for _ in range(p)]
    for v in range(p):
        for s in QR:
            adj[v][(v + s) % p] = 1

    print(f"\np={p}, QR={QR}")
    print(f"chi(3)={legendre(3,p):+d} => q=3 pairs have D^4 with WRONG sign")

    # Find 3 and 5 cycles
    all_cycles = directed_cycles(adj, p, 5)
    for k in sorted(all_cycles.keys()):
        print(f"  {k}-cycles: {len(all_cycles[k])}")

    # Compute chord-pair copart for 3-cycles and 5-cycles
    for k in [3, 5]:
        if k not in all_cycles:
            continue
        chord_copart = defaultdict(int)
        for cycle in all_cycles[k]:
            chords = gap_chords_of_cycle(cycle, p)
            chord_list = sorted(chords)
            for i in range(len(chord_list)):
                for j in range(i+1, len(chord_list)):
                    chord_copart[(chord_list[i], chord_list[j])] += 1

        print(f"\n  {k}-cycle co-participation at p={p}:")
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(qq for qq, t in res)
                chi_q = legendre(min_q, p)
                chi_ab = legendre(a*b, p)
                cp = chord_copart.get((a,b), 0)

                # Only show first few for brevity if many
                if min_q <= 5 or k == 3:
                    print(f"    ({a:>2},{b:>2}): q={min_q:>2}, chi(q)={chi_q:+d}, "
                          f"chi(ab)={chi_ab:+d}, {k}c-copart={cp:>4}")

    # PART 3: The key question: does 5-cycle copart correlate with chi(ab)?
    print("\n\n--- PART 3: DOES 5-CYCLE COPART CORRELATE WITH chi(ab)? ---")
    if 5 in all_cycles:
        chord_copart_5 = defaultdict(int)
        for cycle in all_cycles[5]:
            chords = gap_chords_of_cycle(cycle, p)
            chord_list = sorted(chords)
            for i in range(len(chord_list)):
                for j in range(i+1, len(chord_list)):
                    chord_copart_5[(chord_list[i], chord_list[j])] += 1

        chord_copart_3 = defaultdict(int)
        for cycle in all_cycles[3]:
            chords = gap_chords_of_cycle(cycle, p)
            chord_list = sorted(chords)
            for i in range(len(chord_list)):
                for j in range(i+1, len(chord_list)):
                    chord_copart_3[(chord_list[i], chord_list[j])] += 1

        print(f"\n  q=3 pairs (chi(3)={legendre(3,p):+d}, D^4 wrong sign):")
        print(f"  {'pair':>8} {'chi(ab)':>8} {'c3-cop':>8} {'c5-cop':>8} {'sign from c5?':>14}")
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(qq for qq, t in res)
                if min_q != 3:
                    continue
                chi_ab = legendre(a*b, p)
                c3 = chord_copart_3.get((a,b), 0)
                c5 = chord_copart_5.get((a,b), 0)
                print(f"  ({a:>2},{b:>2}) {chi_ab:>+8} {c3:>8} {c5:>8}")

        print(f"\n  q=5 pairs (chi(5)={legendre(5,p):+d}):")
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(qq for qq, t in res)
                if min_q != 5:
                    continue
                chi_ab = legendre(a*b, p)
                c3 = chord_copart_3.get((a,b), 0)
                c5 = chord_copart_5.get((a,b), 0)
                print(f"  ({a:>2},{b:>2}) {chi_ab:>+8} {c3:>8} {c5:>8}")


if __name__ == '__main__':
    main()
