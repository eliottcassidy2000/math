import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
#!/usr/bin/env python3
"""
INV-011: Does the SCORE SEQUENCE determine H(T) mod 4?

Background:
  OCF gives H(T) = I(Omega(T), 2) = 1 + 2*alpha_0 + 4*alpha_1 + ...
  where alpha_k = number of independent sets of size k in the conflict graph.
  So H(T) mod 4 = (1 + 2*alpha_0) mod 4, i.e. determined by parity of alpha_0.
  alpha_0 = number of directed odd cycles (as vertex sets with a Hamiltonian cycle).

  At n=3,4: alpha_0 = c3 (only 3-cycles exist as odd cycles).
  c3 is determined by score sequence via Moon's formula: c3 = C(n,3) - sum s_i*(s_i-1)/2.
  So at n=3,4: score => c3 => alpha_0 mod 2 => H mod 4.

  At n>=5: 5-cycles contribute to alpha_0, so alpha_0 != c3 in general.
  Question: does score still determine H mod 4?

  DEEPER: does score determine alpha_1 = total odd-cycle count?
  And is c3 constant within a score class? (Yes, by Moon's formula.)

Tests:
  n=3,4,5,6: exhaustive enumeration
  n=7: sample 10000 random tournaments
  For each score class, check constancy of H mod 4, c3, and alpha_0 (total odd cycles).
"""

import sys
import os
import random
from collections import defaultdict
from itertools import combinations, permutations

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import (
    tournament_from_bits, all_tournaments, hamiltonian_path_count,
    find_odd_cycles, random_tournament
)


def score_sequence(T):
    """Return sorted score sequence (tuple of out-degrees)."""
    n = len(T)
    scores = sorted(sum(T[v]) for v in range(n))
    return tuple(scores)


def count_3cycles(T):
    """Count directed 3-cycles using Moon's formula: c3 = C(n,3) - sum s_i*(s_i-1)/2."""
    n = len(T)
    scores = [sum(T[v]) for v in range(n)]
    return (n * (n - 1) * (n - 2) // 6) - sum(s * (s - 1) // 2 for s in scores)


def count_odd_cycle_vertex_sets(T):
    """Count the number of vertex subsets that support a directed odd cycle.
    This is alpha_0 = |V(Omega(T))| = number of vertices in the conflict graph.
    Each vertex set of odd size that contains a directed Hamiltonian cycle counts once."""
    n = len(T)
    count = 0
    for size in range(3, n + 1, 2):
        for verts in combinations(range(n), size):
            # Check if there's a directed Hamiltonian cycle on this vertex set
            if _has_ham_cycle(T, verts):
                count += 1
    return count


def _has_ham_cycle(T, verts):
    """Check if the induced tournament on verts has a directed Hamiltonian cycle.
    Fix first vertex to avoid counting rotations."""
    s = list(verts)
    first = s[0]
    rest = s[1:]
    s_set = set(verts)
    # DFS: build path from first through all vertices, last must have arc to first
    stack = [(first, frozenset([first]))]
    target_len = len(s)
    while stack:
        v, visited = stack.pop()
        if len(visited) == target_len:
            if T[v][first]:
                return True
            continue
        for u in s_set:
            if u not in visited and T[v][u]:
                stack.append((u, visited | {u}))
    return False


def count_total_directed_odd_cycles(T):
    """Count the total number of directed odd cycles (not vertex sets).
    Each vertex set may support multiple directed cycles.
    This uses find_odd_cycles from tournament_lib which returns all directed cycles
    (canonical: min vertex first, eliminating rotations)."""
    cycles = find_odd_cycles(T)
    return len(cycles)


def analyze_score_class(results_by_score):
    """Analyze whether H mod 4, c3, alpha_0 are constant within each score class.
    Returns a report dict."""
    report = {}
    for score_seq, entries in sorted(results_by_score.items()):
        h_mod4_vals = set(e['H_mod4'] for e in entries)
        c3_vals = set(e['c3'] for e in entries)
        alpha0_vals = set(e['alpha0'] for e in entries)
        total_cycles_vals = set(e['total_odd_cycles'] for e in entries)

        report[score_seq] = {
            'count': len(entries),
            'H_mod4_constant': len(h_mod4_vals) == 1,
            'H_mod4_values': sorted(h_mod4_vals),
            'c3_constant': len(c3_vals) == 1,
            'c3_values': sorted(c3_vals),
            'alpha0_constant': len(alpha0_vals) == 1,
            'alpha0_values': sorted(alpha0_vals),
            'total_odd_cycles_constant': len(total_cycles_vals) == 1,
            'total_odd_cycles_values': sorted(total_cycles_vals),
        }
    return report


def main():
    random.seed(42)

    overall_h_mod4_fails = {}  # n -> list of score seqs where H mod 4 varies

    # =========================================================================
    # Exhaustive tests for n = 3, 4, 5, 6
    # =========================================================================
    for n in range(3, 7):
        print(f"\n{'='*70}")
        print(f"n = {n}: EXHAUSTIVE enumeration")
        print(f"{'='*70}")

        m = n * (n - 1) // 2
        total_tournaments = 1 << m
        print(f"Total tournaments: {total_tournaments}")

        results_by_score = defaultdict(list)

        for bits in range(total_tournaments):
            T = tournament_from_bits(n, bits)
            H = hamiltonian_path_count(T)
            c3 = count_3cycles(T)
            alpha0 = count_odd_cycle_vertex_sets(T)
            total_odd = count_total_directed_odd_cycles(T)
            ss = score_sequence(T)

            results_by_score[ss].append({
                'H': H,
                'H_mod4': H % 4,
                'c3': c3,
                'alpha0': alpha0,
                'total_odd_cycles': total_odd,
            })

        report = analyze_score_class(results_by_score)

        # Print summary
        print(f"\nScore classes: {len(report)}")
        print(f"{'Score':<30} {'Count':>6} {'c3':>10} {'alpha0':>15} {'#OddCyc':>15} {'H%4':>12} {'H%4 const?':>12}")
        print("-" * 100)

        fails = []
        for ss, info in sorted(report.items()):
            c3_str = str(info['c3_values'][0]) if info['c3_constant'] else str(info['c3_values'])
            a0_str = str(info['alpha0_values'][0]) if info['alpha0_constant'] else str(info['alpha0_values'])
            tc_str = str(info['total_odd_cycles_values'][0]) if info['total_odd_cycles_constant'] else str(info['total_odd_cycles_values'])
            hm4_str = str(info['H_mod4_values'][0]) if info['H_mod4_constant'] else str(info['H_mod4_values'])
            const_str = "YES" if info['H_mod4_constant'] else "**NO**"

            print(f"{str(ss):<30} {info['count']:>6} {c3_str:>10} {a0_str:>15} {tc_str:>15} {hm4_str:>12} {const_str:>12}")

            if not info['H_mod4_constant']:
                fails.append(ss)

        overall_h_mod4_fails[n] = fails

        if fails:
            print(f"\n*** H mod 4 VARIES within {len(fails)} score class(es)! ***")
            for ss in fails:
                entries = results_by_score[ss]
                print(f"\n  Score {ss}:")
                # Show all distinct (H, H%4, c3, alpha0, total_odd) combos
                combos = defaultdict(int)
                for e in entries:
                    key = (e['H'], e['H_mod4'], e['c3'], e['alpha0'], e['total_odd_cycles'])
                    combos[key] += 1
                print(f"    {'H':>6} {'H%4':>5} {'c3':>5} {'alpha0':>7} {'#OddCyc':>8} {'count':>6}")
                for (H, hm4, c3, a0, tc), cnt in sorted(combos.items()):
                    print(f"    {H:>6} {hm4:>5} {c3:>5} {a0:>7} {tc:>8} {cnt:>6}")
        else:
            print(f"\n  H mod 4 is CONSTANT within every score class at n={n}.")

        # Check the formula: H mod 4 = (1 + 2*c3) mod 4 ?
        formula_holds = True
        for ss, entries in results_by_score.items():
            c3_val = entries[0]['c3']
            predicted = (1 + 2 * c3_val) % 4
            for e in entries:
                if e['H_mod4'] != predicted:
                    formula_holds = False
                    break
            if not formula_holds:
                break

        if formula_holds:
            print(f"  Formula H mod 4 = (1 + 2*c3) mod 4 HOLDS at n={n}.")
        else:
            print(f"  Formula H mod 4 = (1 + 2*c3) mod 4 FAILS at n={n}.")

        # Check: alpha0 mod 2 = c3 mod 2 within each score class?
        alpha0_c3_parity = True
        for ss, entries in results_by_score.items():
            for e in entries:
                if e['alpha0'] % 2 != e['c3'] % 2:
                    alpha0_c3_parity = False
                    break
            if not alpha0_c3_parity:
                break
        print(f"  alpha0 mod 2 = c3 mod 2 for all tournaments: {alpha0_c3_parity}")

    # =========================================================================
    # Sampled test for n = 7
    # =========================================================================
    print(f"\n{'='*70}")
    print(f"n = 7: SAMPLED (10000 random tournaments)")
    print(f"{'='*70}")

    n = 7
    num_samples = 10000
    results_by_score = defaultdict(list)

    for _ in range(num_samples):
        T = random_tournament(n)
        H = hamiltonian_path_count(T)
        c3 = count_3cycles(T)
        alpha0 = count_odd_cycle_vertex_sets(T)
        total_odd = count_total_directed_odd_cycles(T)
        ss = score_sequence(T)

        results_by_score[ss].append({
            'H': H,
            'H_mod4': H % 4,
            'c3': c3,
            'alpha0': alpha0,
            'total_odd_cycles': total_odd,
        })

    report = analyze_score_class(results_by_score)

    print(f"\nScore classes seen: {len(report)}")
    print(f"{'Score':<30} {'Count':>6} {'c3':>6} {'alpha0 const?':>14} {'H%4 values':>15} {'H%4 const?':>12}")
    print("-" * 90)

    fails_n7 = []
    for ss, info in sorted(report.items()):
        c3_str = str(info['c3_values'][0]) if info['c3_constant'] else str(info['c3_values'])
        a0_const = "YES" if info['alpha0_constant'] else "NO"
        hm4_str = str(info['H_mod4_values'])
        const_str = "YES" if info['H_mod4_constant'] else "**NO**"

        print(f"{str(ss):<30} {info['count']:>6} {c3_str:>6} {a0_const:>14} {hm4_str:>15} {const_str:>12}")

        if not info['H_mod4_constant']:
            fails_n7.append(ss)

    overall_h_mod4_fails[7] = fails_n7

    if fails_n7:
        print(f"\n*** H mod 4 VARIES within {len(fails_n7)} score class(es) at n=7! ***")
        for ss in fails_n7[:5]:  # show at most 5
            entries = results_by_score[ss]
            print(f"\n  Score {ss} ({len(entries)} samples):")
            combos = defaultdict(int)
            for e in entries:
                key = (e['H'], e['H_mod4'], e['c3'], e['alpha0'])
                combos[key] += 1
            print(f"    {'H':>6} {'H%4':>5} {'c3':>5} {'alpha0':>7} {'count':>6}")
            for (H, hm4, c3, a0), cnt in sorted(combos.items())[:20]:
                print(f"    {H:>6} {hm4:>5} {c3:>5} {a0:>7} {cnt:>6}")
    else:
        print(f"\n  H mod 4 is CONSTANT within every score class seen at n=7.")

    # Also check formula
    formula_holds_n7 = True
    for ss, entries in results_by_score.items():
        c3_val = entries[0]['c3']
        predicted = (1 + 2 * c3_val) % 4
        for e in entries:
            if e['H_mod4'] != predicted:
                formula_holds_n7 = False
                break
        if not formula_holds_n7:
            break
    print(f"  Formula H mod 4 = (1 + 2*c3) mod 4: {'HOLDS' if formula_holds_n7 else 'FAILS'}")

    # alpha0 mod 2 = c3 mod 2?
    alpha0_c3_n7 = True
    for ss, entries in results_by_score.items():
        for e in entries:
            if e['alpha0'] % 2 != e['c3'] % 2:
                alpha0_c3_n7 = False
                break
        if not alpha0_c3_n7:
            break
    print(f"  alpha0 mod 2 = c3 mod 2: {'HOLDS' if alpha0_c3_n7 else 'FAILS'}")

    # =========================================================================
    # GRAND SUMMARY
    # =========================================================================
    print(f"\n{'='*70}")
    print("GRAND SUMMARY: Does score sequence determine H(T) mod 4?")
    print(f"{'='*70}")
    for n_val in sorted(overall_h_mod4_fails.keys()):
        fails = overall_h_mod4_fails[n_val]
        if fails:
            print(f"  n={n_val}: **NO** -- {len(fails)} score class(es) with varying H mod 4")
        else:
            status = "EXHAUSTIVE" if n_val <= 6 else "SAMPLED (10000)"
            print(f"  n={n_val}: YES (verified {status})")

    print("\nKey insight: H = I(Omega,2) = 1 + 2*alpha_0 + 4*alpha_1 + ...")
    print("So H mod 4 = (1 + 2*alpha_0) mod 4, determined by parity of alpha_0.")
    print("If alpha_0 mod 2 = c3 mod 2 (and c3 is score-determined), then score => H mod 4.")
    print("Done.")


if __name__ == '__main__':
    main()
