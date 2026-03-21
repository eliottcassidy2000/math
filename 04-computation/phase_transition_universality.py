#!/usr/bin/env python3
"""
phase_transition_universality.py -- kind-pasteur-2026-03-21-S12

Is the "67% depth" phase transition universal?

Napolitano: phase transition at layer 16 of 24 = 67% depth
Our work: path homology phase transition at n=8 out of tested n=3..12

Question: Is there a universal "2/3 critical fraction" in both systems?

For tournaments, the critical phenomenon is:
- n=8: beta_3 first exceeds 1, seesaw breaks, beta_4 first appears
- n=8/12 = 0.67 if the "natural range" is n=3..12

But this is just a coincidence of our testing range. Let's look deeper.

Actually, the more interesting analogy is:
- Tournament on n vertices has m = n(n-1)/2 arcs
- At what fraction of arcs does the homological structure change?

For beta_3 onset at n=6 (first non-zero): m=15, but beta_3 is tiny
For beta_3>1 onset at n=8: m=28
For beta_4 onset at n=8: m=28

The critical vertex fraction: 8/n is not well-defined without n_max.
But the critical edge fraction: C(8,2)/C(n,2) = 28/C(n,2) IS well-defined.

For n=12: 28/66 = 0.424
For n=20: 28/190 = 0.147
Not 2/3. The analogy doesn't work this way.

Let's try another angle: the Cartan decomposition ratio at critical n.

Author: kind-pasteur-2026-03-21-S12
"""

import numpy as np


def explore_critical_ratios():
    """
    Look for universal ratios at critical points in tournament theory.
    """
    print("PHASE TRANSITION UNIVERSALITY ANALYSIS")
    print("=" * 60)

    # Tournament critical thresholds
    print("\nTournament phase transitions:")
    print("-" * 40)
    transitions = [
        ("beta_1 onset", 3, "First non-trivial homology"),
        ("beta_3 onset", 6, "First 3-holes appear"),
        ("beta_3 > 1", 8, "Multiple 3-holes; seesaw breaks"),
        ("beta_4 onset", 8, "4-dimensional holes appear"),
        ("Real roots fail", 9, "I(Omega,x) can have complex roots"),
        ("Claw-free fails", 9, "Omega not claw-free"),
        ("S_{2,1,1}-free fails", 10, "Omega not S_{2,1,1}-free"),
    ]

    for name, n_crit, desc in transitions:
        m = n_crit * (n_crit - 1) // 2
        # Active/dark at critical n
        active = n_crit * (n_crit - 1) // 2
        dark = n_crit * (n_crit + 1) // 2 - 1
        ratio = dark / (active + dark + 1)  # fraction of gl(n) that is dark
        # Napolitano's 67% comparison
        print(f"  {name:25s} n={n_crit:2d}  arcs={m:3d}  "
              f"dark_frac={ratio:.3f}  desc: {desc}")

    # The Napolitano "67%" comparison
    print("\nNapolitano phase transition:")
    print(f"  Layer 16 of 24 = {16/24:.3f} = 2/3")
    print(f"  For 2-norm arch: ~65% depth")
    print(f"  For 1-norm arch: ~83% depth")

    # Is there a 2/3 ratio hiding in tournament theory?
    print("\n2/3 ratios in tournament theory:")
    print("-" * 40)

    # Transitivity of regular tournaments
    for n in [3, 5, 7, 9, 11, 13]:
        # For regular tournament: c3 = C(n,3) - n*C((n-1)/2, 2)
        k = (n-1)//2
        c3 = n*(n-1)*(n-2)//6 - n*k*(k-1)//2
        total_triples = n*(n-1)*(n-2)//6
        trans_frac = 1 - c3/total_triples
        print(f"  n={n:2d}: transitivity of regular = {trans_frac:.4f} "
              f"(c3={c3}, total={total_triples})")
        # trans_frac = 1 - c3/C(n,3)
        # For regular: c3/C(n,3) = 1 - n*C((n-1)/2, 2)/C(n,3)
        # = 1 - (n-1)(n-3)/(4(n-2))

    print("\n  Formula: transitivity = 1 - (n-1)(n-3)/(4(n-2))")
    print("  As n->inf: transitivity -> 3/4 = 0.75")
    print("  At n=7: transitivity = 6*4/(4*5) = 24/20 = ... ")
    for n in [3, 5, 7, 9, 11]:
        val = 1 - (n-1)*(n-3)/(4*(n-2))
        print(f"    n={n}: formula = {val:.4f}")

    # The 2/3 connection: score variance
    print("\n\nScore variance ratio:")
    print("-" * 40)
    for n in [3, 5, 7, 9, 11]:
        # Regular: all scores = (n-1)/2, variance = 0
        # Transitive: scores = 0,1,...,n-1, variance = (n^2-1)/12
        max_var = (n*n - 1) / 12
        # What fraction of max variance does a typical tournament have?
        # Monte Carlo estimate
        np.random.seed(42)
        vars_list = []
        for _ in range(10000):
            # Random tournament
            T = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i+1, n):
                    if np.random.random() < 0.5:
                        T[i][j] = 1
                    else:
                        T[j][i] = 1
            scores = T.sum(axis=1)
            vars_list.append(np.var(scores))
        avg_var = np.mean(vars_list)
        frac = avg_var / max_var
        print(f"  n={n:2d}: avg_var/max_var = {frac:.4f} "
              f"(avg_var={avg_var:.3f}, max_var={max_var:.3f})")

    # Key: is there a natural 2/3 in the relationship between
    # active and dark energy for attention matrices?
    print("\n\nActive energy fraction for trained vs random:")
    print("-" * 40)
    print("Random attention (our computation):")
    for n in [4, 7, 11, 16]:
        anti_frac = {4: 0.132, 7: 0.187, 11: 0.221, 16: 0.245}[n]
        dark_frac = {4: 0.683, 7: 0.723, 11: 0.727, 16: 0.723}[n]
        print(f"  n={n:2d}: active={anti_frac:.3f}, dark={dark_frac:.3f}, "
              f"dark/(active+dark)={dark_frac/(anti_frac+dark_frac):.3f}")
    print("\nFor n=4 (Napolitano): dark/(anti+dark) = 0.838")
    print("This means 84% of non-scalar energy is 'dark' for random attention.")
    print("Napolitano's 'discovery' that dark carries information is partly")
    print("explained by this baseline imbalance.")

    print("\n\nConclusion: The '2/3' is NOT a universal constant.")
    print("It appears in Napolitano because 16/24 = 2/3 for Qwen-0.5B.")
    print("It does NOT appear at a fixed fraction in tournament homology.")
    print("The universal ratio, if any, is (n+1)/(n-1) for dark/active dims.")


if __name__ == "__main__":
    explore_critical_ratios()
