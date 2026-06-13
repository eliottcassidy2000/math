#!/usr/bin/env python3
"""
cycle_space_formula_s179.py — The Cycle Space Formula: R = 2(HC - E[HC|score])
opus-2026-03-22-S179

From kind-pasteur S20af: the OCR residual R = H - E[H|score] is
EXACTLY R = 2 × (HC - E[HC|score]).

The residual IS twice the Hamiltonian cycle deviation.
This connects the cut/cycle decomposition to the H formula:
  H = E[H|score] + 2*(HC - E[HC|score])
  = (score-determined part) + 2*(cycle fluctuation)

Verify this and understand what it means.
"""

from math import comb, factorial
from collections import defaultdict
from fractions import Fraction

print("=" * 72)
print("  THE CYCLE SPACE FORMULA: R = 2(HC - E[HC|score])")
print("  opus-2026-03-22-S179")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def count_ham_cycles(adj, n):
    dp = [0]*((1<<n)*n)
    dp[(1<<0)*n+0] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(1, n) if adj[v][0])

def scores(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

# ==========================================================================
# VERIFY AT n=5
# ==========================================================================

print("VERIFICATION AT n=5")
print("-" * 60)
print()

n = 5
score_groups = defaultdict(list)

for bits in range(1 << comb(n, 2)):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    hc = count_ham_cycles(adj, n)
    ss = scores(adj, n)
    score_groups[ss].append({'h': h, 'hc': hc, 'bits': bits})

# Compute E[H|score] and E[HC|score] for each score class
print(f"  {'score':>16s}  {'E[H|s]':>7s}  {'E[HC|s]':>7s}  {'H values':>12s}  {'HC values':>12s}")
print(f"  {'─'*16}  {'─'*7}  {'─'*7}  {'─'*12}  {'─'*12}")

for ss in sorted(score_groups.keys()):
    entries = score_groups[ss]
    eh = Fraction(sum(e['h'] for e in entries), len(entries))
    ehc = Fraction(sum(e['hc'] for e in entries), len(entries))
    h_vals = sorted(set(e['h'] for e in entries))
    hc_vals = sorted(set(e['hc'] for e in entries))
    print(f"  {str(ss):>16s}  {float(eh):7.2f}  {float(ehc):7.2f}  "
          f"{str(h_vals):>12s}  {str(hc_vals):>12s}")

print()

# Now verify: R = 2*(HC - E[HC|score])
print(f"  VERIFICATION: R = H - E[H|score] vs 2*(HC - E[HC|score])")
print()

errors = 0
for ss in sorted(score_groups.keys()):
    entries = score_groups[ss]
    eh = Fraction(sum(e['h'] for e in entries), len(entries))
    ehc = Fraction(sum(e['hc'] for e in entries), len(entries))

    for e in entries:
        R_actual = Fraction(e['h']) - eh
        R_predicted = 2 * (Fraction(e['hc']) - ehc)

        if R_actual != R_predicted:
            errors += 1
            if errors <= 3:
                print(f"  MISMATCH at bits={e['bits']}: "
                      f"R_actual={float(R_actual):.4f}, R_predicted={float(R_predicted):.4f}")

total = sum(len(entries) for entries in score_groups.values())
print(f"  Verified {total} tournaments: {errors} errors")
if errors == 0:
    print(f"  *** R = 2*(HC - E[HC|score]) EXACT AT n=5 ***")
print()

# ==========================================================================
# VERIFY AT n=4 AND n=3
# ==========================================================================

for n_check in [3, 4]:
    print(f"  VERIFICATION AT n={n_check}:")
    sg = defaultdict(list)
    for bits in range(1 << comb(n_check, 2)):
        adj = tournament_adj(n_check, bits)
        h = H_dp(adj, n_check)
        hc = count_ham_cycles(adj, n_check)
        ss = scores(adj, n_check)
        sg[ss].append({'h': h, 'hc': hc})

    err = 0
    for ss, entries in sg.items():
        eh = Fraction(sum(e['h'] for e in entries), len(entries))
        ehc = Fraction(sum(e['hc'] for e in entries), len(entries))
        for e in entries:
            R_actual = Fraction(e['h']) - eh
            R_predicted = 2 * (Fraction(e['hc']) - ehc)
            if R_actual != R_predicted:
                err += 1

    total_check = sum(len(entries) for entries in sg.values())
    print(f"    {total_check} tournaments: {err} errors → "
          f"{'EXACT ✓' if err == 0 else 'FAILS ✗'}")
    print()

# ==========================================================================
# VERIFY AT n=6
# ==========================================================================

print(f"  VERIFICATION AT n=6:")
n_check = 6
sg6 = defaultdict(list)
for bits in range(1 << comb(n_check, 2)):
    adj = tournament_adj(n_check, bits)
    h = H_dp(adj, n_check)
    hc = count_ham_cycles(adj, n_check)
    ss = scores(adj, n_check)
    sg6[ss].append({'h': h, 'hc': hc})

err6 = 0
for ss, entries in sg6.items():
    eh = Fraction(sum(e['h'] for e in entries), len(entries))
    ehc = Fraction(sum(e['hc'] for e in entries), len(entries))
    for e in entries:
        R_actual = Fraction(e['h']) - eh
        R_predicted = 2 * (Fraction(e['hc']) - ehc)
        if R_actual != R_predicted:
            err6 += 1

total6 = sum(len(entries) for entries in sg6.values())
print(f"    {total6} tournaments: {err6} errors → "
      f"{'EXACT ✓' if err6 == 0 else f'FAILS ({err6} errors)'}")
print()

# ==========================================================================
# WHAT THE FORMULA MEANS
# ==========================================================================

print()
print("=" * 72)
print("  WHAT THE FORMULA MEANS")
print("=" * 72)
print()

print(f"""
  THE EXACT FORMULA: R(T) = 2 × (HC(T) - E[HC|score(T)])

  This says: the part of H that scores CAN'T predict
  is EXACTLY twice the deviation of Hamiltonian cycles
  from their score-conditioned mean.

  EQUIVALENTLY:
    H(T) = E[H|score(T)] + 2 × (HC(T) - E[HC|score(T)])

  Which simplifies to:
    H(T) = E[H|score(T)] + 2×HC(T) - 2×E[HC|score(T)]

  Since E[H|score] - 2×E[HC|score] is a CONSTANT within each score class:
    H(T) = constant(score) + 2×HC(T)

  So WITHIN each score class: H and HC differ by a constant!
  H = 2×HC + const(score)

  Since H = n×HC + L (from the H=nHC+L decomposition):
    n×HC + L = 2×HC + const
    (n-2)×HC + L = const
    L = const - (n-2)×HC

  At n=5: L = const - 3×HC

  So L and HC are LINEARLY RELATED within each score class!
  L = a - 3×HC where a depends on the score.

  THIS IS WHY:
  - Within the PoS class: HC=1→L=6, HC=2→L=3, HC=3→L=0
    L = 9 - 3×HC ✓ (a=9 for the PoS class)

  - The formula connects the three decompositions:
    H = score part + cycle part (cut/cycle)
    H = n×HC + L (cyclic/linear)
    H = 1 + 2α₁ + 4α₂ + ... (independence polynomial)

  And ALL of them are connected by:
    R = 2×(HC - E[HC|score])
    L = const - (n-2)×HC

  THE CYCLE SPACE IS THE HAMILTONIAN CYCLE COUNT.
  The even graph determines HC, which determines R, which IS the 3% residual.
""")

# Verify L = a - (n-2)×HC within each score class at n=5
print(f"  VERIFY: L = a - {n-2}×HC within each score class at n={n}:")
print()

for ss in sorted(score_groups.keys()):
    entries = score_groups[ss]
    hc_vals = sorted(set(e['hc'] for e in entries))

    if len(hc_vals) > 1:
        # Check linearity
        pairs = [(e['hc'], e['h'] - n*e['hc']) for e in entries]
        hc_l_pairs = sorted(set(pairs))

        # L = a - (n-2)*HC → a = L + (n-2)*HC
        a_vals = set(l + (n-2)*hc for hc, l in hc_l_pairs)
        print(f"  {str(ss):>16s}: HC→L = {hc_l_pairs}, a = {a_vals}")
        if len(a_vals) == 1:
            print(f"    ✓ L = {a_vals.pop()} - {n-2}×HC (linear)")
    else:
        hc = hc_vals[0]
        l = entries[0]['h'] - n*hc
        print(f"  {str(ss):>16s}: HC={hc}, L={l} (single value)")

print()

# Does R = 2*(HC - E[HC|score]) hold at n=6?
if err6 == 0:
    print(f"  *** THE FORMULA R = 2(HC - E[HC|score]) IS EXACT AT n=3,4,5,6! ***")
    print(f"  This is likely TRUE FOR ALL n.")
elif err6 > 0:
    print(f"  The formula FAILS at n=6 with {err6} errors.")
    print(f"  The residual at n=6 has contributions beyond HC.")
    # Show the error distribution
    err_dist = defaultdict(int)
    for ss, entries in sg6.items():
        eh = Fraction(sum(e['h'] for e in entries), len(entries))
        ehc = Fraction(sum(e['hc'] for e in entries), len(entries))
        for e in entries:
            R_actual = Fraction(e['h']) - eh
            R_predicted = 2 * (Fraction(e['hc']) - ehc)
            diff = R_actual - R_predicted
            err_dist[float(diff)] += 1

    print(f"  Error distribution (R_actual - R_predicted):")
    for d, count in sorted(err_dist.items()):
        if count > 0 and d != 0:
            print(f"    Δ = {d:+.4f}: {count} tournaments")

print()
print("opus-2026-03-22-S179")
