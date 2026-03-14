#!/usr/bin/env python3
"""
alpha1_parity_analysis.py — Investigate parity of α₁ (total directed odd cycle count).

α₁ = t₃ + d₅ + d₇ + ... where t₃ = # of 3-cycle vertex sets (each supporting 1 directed
3-cycle with min-vertex start), d_k = directed k-cycles similarly canonicalized.

Key question: Is there a tournament invariant determining α₁ mod 2?
Since α₁ = t₃ + d₅ (at n=5,6), we have α₁ mod 2 = (t₃ + d₅) mod 2.
"""
from itertools import permutations, combinations
from collections import defaultdict

def make_tournament(bits, n):
    """Build adjacency matrix from bit encoding."""
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def score_sequence(A, n):
    """Sorted out-degree sequence."""
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_hp(A, n):
    """Count Hamiltonian paths."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def count_directed_cycles_by_length(A, n):
    """Count directed odd cycles by length, canonicalized (start=min vertex).
    Returns dict: length -> count."""
    counts = defaultdict(int)
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            # Fix start = min vertex = verts[0]
            for p in permutations(verts[1:]):
                order = [verts[0]] + list(p)
                is_cycle = True
                for idx in range(k):
                    if A[order[idx]][order[(idx+1) % k]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    counts[k] += 1
    return counts

def count_3cycles(A, n):
    """Count 3-cycle vertex SETS (each supports exactly 1 directed 3-cycle)."""
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check if it's a 3-cycle (not transitive)
        out = [A[i][j] + A[i][k], A[j][i] + A[j][k], A[k][i] + A[k][j]]
        if sorted(out) == [0, 1, 2]:
            # transitive triple
            pass
        else:
            # 3-cycle: each vertex has out-degree 1 within the triple
            t3 += 1
    return t3

def is_self_complementary(A, n):
    """Check if T is isomorphic to its complement tournament T^op."""
    # T^op: reverse all arcs => A_op[i][j] = A[j][i]
    # Check all permutations (feasible for small n)
    for perm in permutations(range(n)):
        is_iso = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if A[i][j] != A[perm[j]][perm[i]]:
                    is_iso = False
                    break
            if not is_iso:
                break
        if is_iso:
            return True
    return False

def analyze_n(n, exhaustive_sc=True):
    """Full parity analysis at given n."""
    total_edges = n * (n - 1) // 2
    num_t = 2 ** total_edges

    print(f"\n{'='*70}")
    print(f"  n = {n}: {num_t} tournaments, {total_edges} edges")
    print(f"{'='*70}")

    # Collect data
    data = []
    for bits in range(num_t):
        A = make_tournament(bits, n)
        H = count_hp(A, n)
        ss = score_sequence(A, n)
        t3 = count_3cycles(A, n)
        cycle_counts = count_directed_cycles_by_length(A, n)
        d5 = cycle_counts.get(5, 0)
        d7 = cycle_counts.get(7, 0)
        alpha1 = sum(cycle_counts.values())  # total directed odd cycles

        data.append({
            'bits': bits, 'H': H, 'ss': ss, 't3': t3,
            'd5': d5, 'd7': d7, 'alpha1': alpha1,
            'alpha1_mod2': alpha1 % 2,
            'H_mod4': H % 4,
            'A': A,
        })

    # Basic parity stats
    even_count = sum(1 for d in data if d['alpha1_mod2'] == 0)
    odd_count = sum(1 for d in data if d['alpha1_mod2'] == 1)
    print(f"\nα₁ parity: even={even_count}, odd={odd_count} (of {num_t})")

    # Verify H mod 4 relationship: H = 1 + 2*alpha1 + 4*alpha2 => H mod 4 = 1+2*(alpha1 mod 2)
    print(f"\nVerifying H mod 4 = 1 + 2*(α₁ mod 2)...")
    h_check_ok = True
    for d in data:
        expected = 1 + 2 * d['alpha1_mod2']
        if d['H_mod4'] != expected:
            print(f"  FAIL at bits={d['bits']}: H={d['H']}, H%4={d['H_mod4']}, α₁%2={d['alpha1_mod2']}, expected {expected}")
            h_check_ok = False
    if h_check_ok:
        print(f"  VERIFIED: H mod 4 = 1 + 2*(α₁ mod 2) for all {num_t} tournaments")

    # Cross-correlate with score sequence
    print(f"\n--- α₁ parity by score sequence ---")
    ss_groups = defaultdict(lambda: {'even': 0, 'odd': 0, 'total': 0})
    for d in data:
        ss_groups[d['ss']]['total'] += 1
        if d['alpha1_mod2'] == 0:
            ss_groups[d['ss']]['even'] += 1
        else:
            ss_groups[d['ss']]['odd'] += 1

    for ss in sorted(ss_groups.keys()):
        g = ss_groups[ss]
        det = "BOTH" if g['even'] > 0 and g['odd'] > 0 else ("ALL EVEN" if g['odd'] == 0 else "ALL ODD")
        print(f"  {ss}: {g['total']} tourn, even={g['even']}, odd={g['odd']} [{det}]")

    # Cross-correlate with t₃
    print(f"\n--- α₁ parity by t₃ ---")
    t3_groups = defaultdict(lambda: {'even': 0, 'odd': 0, 'total': 0})
    for d in data:
        t3_groups[d['t3']]['total'] += 1
        if d['alpha1_mod2'] == 0:
            t3_groups[d['t3']]['even'] += 1
        else:
            t3_groups[d['t3']]['odd'] += 1

    for t3 in sorted(t3_groups.keys()):
        g = t3_groups[t3]
        det = "BOTH" if g['even'] > 0 and g['odd'] > 0 else ("ALL EVEN" if g['odd'] == 0 else "ALL ODD")
        print(f"  t₃={t3}: {g['total']} tourn, even={g['even']}, odd={g['odd']} [{det}]")

    # α₁ = t₃ + d₅ (+ d₇ for n≥7). Check d₅ mod 2
    print(f"\n--- Decomposition: α₁ = t₃ + d₅" + (" + d₇" if n >= 7 else "") + " ---")
    print(f"--- d₅ parity analysis ---")
    d5_groups = defaultdict(lambda: {'even': 0, 'odd': 0})
    for d in data:
        if d['d5'] % 2 == 0:
            d5_groups[d['t3']]['even'] += 1
        else:
            d5_groups[d['t3']]['odd'] += 1

    for t3 in sorted(d5_groups.keys()):
        g = d5_groups[t3]
        print(f"  t₃={t3}: d₅ even={g['even']}, d₅ odd={g['odd']}")

    # Check if t₃ mod 2 determines α₁ mod 2
    print(f"\n--- Does t₃ mod 2 determine α₁ mod 2? ---")
    t3_parity_groups = defaultdict(lambda: {'match': 0, 'mismatch': 0})
    for d in data:
        if d['t3'] % 2 == d['alpha1_mod2']:
            t3_parity_groups[d['t3'] % 2]['match'] += 1
        else:
            t3_parity_groups[d['t3'] % 2]['mismatch'] += 1

    for p in sorted(t3_parity_groups.keys()):
        g = t3_parity_groups[p]
        print(f"  t₃≡{p} mod 2: α₁≡t₃ (mod 2) in {g['match']}, α₁≢t₃ in {g['mismatch']}")

    total_match = sum(g['match'] for g in t3_parity_groups.values())
    total_mismatch = sum(g['mismatch'] for g in t3_parity_groups.values())
    if total_mismatch == 0:
        print(f"  => YES! α₁ ≡ t₃ (mod 2) always => d₅ is ALWAYS EVEN")
    else:
        print(f"  => NO. Match={total_match}, Mismatch={total_mismatch}")

    # Directly check: is d₅ always even?
    print(f"\n--- Is d₅ always even? ---")
    d5_odd_count = sum(1 for d in data if d['d5'] % 2 == 1)
    print(f"  d₅ odd in {d5_odd_count}/{num_t} tournaments")
    if d5_odd_count == 0:
        print(f"  => YES! d₅ is ALWAYS EVEN at n={n}")
    else:
        print(f"  => NO. d₅ can be odd.")
        # Show examples
        for d in data:
            if d['d5'] % 2 == 1:
                print(f"    bits={d['bits']}: t₃={d['t3']}, d₅={d['d5']}, α₁={d['alpha1']}")
                break

    if n >= 7:
        d7_odd_count = sum(1 for d in data if d['d7'] % 2 == 1)
        print(f"\n--- Is d₇ always even? ---")
        print(f"  d₇ odd in {d7_odd_count}/{num_t} tournaments")

    # Self-complementary analysis (only for small n due to cost)
    if exhaustive_sc and n <= 6:
        print(f"\n--- Self-complementary tournament analysis ---")
        sc_data = []
        for d in data:
            if is_self_complementary(d['A'], n):
                sc_data.append(d)
        print(f"  {len(sc_data)} self-complementary tournaments found")
        if sc_data:
            sc_even = sum(1 for d in sc_data if d['alpha1_mod2'] == 0)
            sc_odd = sum(1 for d in sc_data if d['alpha1_mod2'] == 1)
            print(f"  α₁ parity among SC: even={sc_even}, odd={sc_odd}")
            for d in sc_data[:10]:
                print(f"    bits={d['bits']}: ss={d['ss']}, t₃={d['t3']}, d₅={d['d5']}, α₁={d['alpha1']}, α₁%2={d['alpha1_mod2']}")

    # Detailed (t₃, d₅) breakdown
    print(f"\n--- (t₃, d₅) joint distribution ---")
    td_groups = defaultdict(int)
    for d in data:
        td_groups[(d['t3'], d['d5'])] += 1
    for k in sorted(td_groups.keys()):
        t3, d5 = k
        a1 = t3 + d5
        print(f"  (t₃={t3}, d₅={d5}): {td_groups[k]} tourn, α₁={a1}, α₁%2={a1%2}")

    # Check d₅ mod 2 as function of score sequence
    print(f"\n--- d₅ mod 2 by score sequence ---")
    d5_ss = defaultdict(lambda: {'even': 0, 'odd': 0})
    for d in data:
        if d['d5'] % 2 == 0:
            d5_ss[d['ss']]['even'] += 1
        else:
            d5_ss[d['ss']]['odd'] += 1
    for ss in sorted(d5_ss.keys()):
        g = d5_ss[ss]
        det = "BOTH" if g['even'] > 0 and g['odd'] > 0 else ("ALL EVEN" if g['odd'] == 0 else "ALL ODD")
        print(f"  {ss}: d₅ even={g['even']}, d₅ odd={g['odd']} [{det}]")

    # Check formula: d₅ mod 2 in terms of edges or other invariants
    # Each directed 5-cycle (v0,v1,v2,v3,v4) with start=min has edges v0->v1->v2->v3->v4->v0
    # Can we relate d₅ mod 2 to some edge count?
    print(f"\n--- Exploring d₅ mod 2 formula ---")
    # Try: d₅ mod 2 vs t₃ mod something
    # Try: d₅ mod 2 vs C(n,2) - t₃ (transitive triples)
    # At n=5: C(5,3)=10 triples. Transitive triples = 10 - t₃
    for d in data[:5]:
        trans = n*(n-1)*(n-2)//6 - d['t3']
        print(f"  bits={d['bits']}: t₃={d['t3']}, trans={trans}, d₅={d['d5']}, d₅%2={d['d5']%2}")

    # Check: is d₅ ≡ f(t₃) mod 2 for some function f?
    print(f"\n--- d₅ mod 2 as function of t₃ ---")
    t3_to_d5_parities = defaultdict(set)
    for d in data:
        t3_to_d5_parities[d['t3']].add(d['d5'] % 2)
    for t3 in sorted(t3_to_d5_parities.keys()):
        vals = t3_to_d5_parities[t3]
        if len(vals) == 1:
            print(f"  t₃={t3}: d₅ mod 2 = {vals.pop()} (DETERMINED)")
        else:
            print(f"  t₃={t3}: d₅ mod 2 ∈ {vals} (NOT determined by t₃ alone)")

    # Check: is d₅ ≡ something mod 2 involving edge sums?
    # Try quadratic invariants: sum of A[i][j]*A[j][k]*A[k][l]*A[l][m]*A[m][i]
    # That's exactly d₅ summed over all orderings...
    # Try: sum_i (deg_i choose 2) mod 2
    print(f"\n--- Trying degree-based formulas for d₅ mod 2 ---")
    # s₂ = sum_i C(deg_i, 2) = sum_i deg_i*(deg_i-1)/2
    tests_passed = {'s2': True, 'sum_deg_sq': True, 'sum_deg_cube': True}
    for d in data:
        degs = [sum(d['A'][i]) for i in range(n)]
        s2 = sum(di*(di-1)//2 for di in degs)
        sum_sq = sum(di*di for di in degs)
        sum_cube = sum(di**3 for di in degs)

        if s2 % 2 != d['d5'] % 2:
            tests_passed['s2'] = False
        if sum_sq % 2 != d['d5'] % 2:
            tests_passed['sum_deg_sq'] = False
        if sum_cube % 2 != d['d5'] % 2:
            tests_passed['sum_deg_cube'] = False

    for name, ok in tests_passed.items():
        print(f"  d₅ ≡ {name} (mod 2): {'YES' if ok else 'NO'}")

    # More targeted: try number of 4-paths, or other subgraph counts
    # The key relation: at n=5, number of directed 5-cycles = (1/5)*sum over 5-cycles
    # For a 5-vertex tournament: d₅ = (# cyclic orderings that are directed) / 1
    # Actually the number of directed 5-cycles (Hamiltonian cycles at n=5) = HC/2 when n=5
    # Wait: at n=5, a directed 5-cycle IS a Hamiltonian cycle direction
    # HC = 2*d₅ (each undirected HC gives 2 directed ones, but min-start picks 1 from each direction)
    # Actually no: each directed HC has 5 rotations. We pick start=min. That's 1 per directed HC.
    # A directed HC and its reverse are different. So d₅ at n=5 = # directed Hamiltonian cycles / 5... no.
    # With our canonicalization: d₅ = # directed HCs with start=min vertex.
    # At n=5, the 5 rotations of a directed HC all have the same vertex set, so
    # we pick exactly 1 representative per directed HC. And each undirected HC gives 2 directed HCs.

    if n == 5:
        print(f"\n--- n=5 specific: d₅ = directed Hamiltonian cycles (canonicalized) ---")
        print(f"  d₅ = # directed HCs starting at vertex 0")
        # At n=5, also check: HC (total directed Hamiltonian CYCLES) = 5 * d₅
        # And actual total directed HCs (all start vertices) = d₅ * 5... no, each appears once.
        # Our d₅ = # distinct directed HCs (up to rotation).
        # Total directed HC count (without canonicalization) = 5 * d₅ ... no,
        # actually d₅ IS the count of distinct directed cycles, each represented once.
        # Let's verify with HC counting
        for d in data[:5]:
            A = d['A']
            hc = 0
            for perm in permutations(range(n)):
                if perm[0] != 0:
                    continue  # canonicalize start
                ok = True
                for i in range(n):
                    if A[perm[i]][perm[(i+1)%n]] != 1:
                        ok = False
                        break
                if ok:
                    hc += 1
            print(f"    bits={d['bits']}: d₅={d['d5']}, HC(start=0)={hc}, match={d['d5']==hc}")

    # H mod 4 distribution
    print(f"\n--- H mod 4 distribution ---")
    hmod4 = defaultdict(int)
    for d in data:
        hmod4[d['H_mod4']] += 1
    for v in sorted(hmod4.keys()):
        print(f"  H ≡ {v} (mod 4): {hmod4[v]} tournaments")

    # Final summary
    print(f"\n{'='*70}")
    print(f"  SUMMARY for n={n}")
    print(f"{'='*70}")
    print(f"  α₁ even: {even_count}/{num_t} ({100*even_count/num_t:.1f}%)")
    print(f"  α₁ odd:  {odd_count}/{num_t} ({100*odd_count/num_t:.1f}%)")
    d5_always_even = (d5_odd_count == 0)
    print(f"  d₅ always even: {d5_always_even}")
    if d5_always_even:
        print(f"  => α₁ mod 2 = t₃ mod 2 (since d₅ ≡ 0 mod 2)")
        print(f"  => H mod 4 = 1 + 2*(t₃ mod 2)")
    print()

    return data


# ==========================================================
# MAIN
# ==========================================================
print("="*70)
print("  α₁ PARITY ANALYSIS")
print("  α₁ = total # directed odd cycles in Ω(T)")
print("="*70)

# n=5 exhaustive
data5 = analyze_n(5, exhaustive_sc=True)

# n=6 exhaustive (2^15 = 32768 tournaments — feasible but slower)
print("\n\n" + "#"*70)
print("# n=6 analysis (this may take a few minutes)")
print("#"*70)
data6 = analyze_n(6, exhaustive_sc=False)  # SC check too expensive at n=6

# n=7 sample (2^21 too large for exhaustive with cycle counting)
# Just do a small sample
print("\n\n" + "#"*70)
print("# n=7 sample (100 random tournaments)")
print("#"*70)

import random
random.seed(42)
n = 7
total_edges = n*(n-1)//2
sample_size = 100

d7_always_even = True
d7_odd_examples = []
alpha1_parities = {'even': 0, 'odd': 0}

for _ in range(sample_size):
    bits = random.randint(0, 2**total_edges - 1)
    A = make_tournament(bits, n)
    cc = count_directed_cycles_by_length(A, n)
    t3 = cc.get(3, 0)
    d5 = cc.get(5, 0)
    d7 = cc.get(7, 0)
    alpha1 = t3 + d5 + d7

    if d5 % 2 == 1:
        d7_always_even = False  # misleading var name, checking d5 at n=7
    if d7 % 2 == 1:
        d7_odd_examples.append((bits, t3, d5, d7, alpha1))

    if alpha1 % 2 == 0:
        alpha1_parities['even'] += 1
    else:
        alpha1_parities['odd'] += 1

print(f"\nn=7 sample ({sample_size} tournaments):")
print(f"  α₁ parity: even={alpha1_parities['even']}, odd={alpha1_parities['odd']}")
print(f"  d₅ always even in sample: {d7_always_even}")
print(f"  d₇ odd examples: {len(d7_odd_examples)}")
if d7_odd_examples:
    for ex in d7_odd_examples[:5]:
        bits, t3, d5, d7, a1 = ex
        print(f"    bits={bits}: t₃={t3}, d₅={d5}, d₇={d7}, α₁={a1}, α₁%2={a1%2}")

# Check if α₁ ≡ t₃ mod 2 at n=7
print(f"\n  Checking α₁ ≡ t₃ (mod 2) at n=7:")
match_count = 0
for _ in range(sample_size):
    bits = random.randint(0, 2**total_edges - 1)
    A = make_tournament(bits, n)
    cc = count_directed_cycles_by_length(A, n)
    t3 = cc.get(3, 0)
    alpha1 = sum(cc.values())
    if alpha1 % 2 == t3 % 2:
        match_count += 1

print(f"  α₁ ≡ t₃ (mod 2) in {match_count}/{sample_size} samples")
if match_count == sample_size:
    print(f"  => Consistent with α₁ ≡ t₃ (mod 2) at n=7")
    print(f"  => This would mean d₅ + d₇ ≡ 0 (mod 2) always")
else:
    print(f"  => α₁ ≢ t₃ (mod 2) in general at n=7")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
