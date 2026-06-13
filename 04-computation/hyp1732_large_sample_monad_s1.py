#!/usr/bin/env python3
"""
monad-compute-2026-06-06-S1: Large-sample verification of HYP-1732.

HYP-1732 (OPEN-Q-053, 🔴 CRITICAL):
  For any tournament T with alpha(Omega(T)) = 2, and C* chosen by the
  pair-partner construction, the count p = #{odd cycles disjoint from C*}
  satisfies  alpha2(Omega) <= p * (m - p),
  where m = total number of directed odd cycles (= alpha1 = |V(Omega)|)
  and alpha2 = number of vertex-disjoint pairs of odd cycles.

Equivalent (all proved): alpha2 <= p(m-p)  <=>  I(Omega, -1/p) <= 0
  (since for d=2, I(Omega,x) = 1 + m*x + alpha2*x^2, and
   I(-1/p) <= 0  <=>  p^2 - m*p + alpha2 <= 0  <=>  alpha2 <= p(m-p)).
We cross-check BOTH forms per test.

Prior record (opus-2026-05-22-S2): 1637 tests at n=7..11, 0 violations.
This session: faster DFS odd-cycle enumeration + targeted near-transitive
sampling (alpha=2 needs FEW cycles, so it is rare in uniform random
tournaments at large n; near-transitive perturbations harvest many more
alpha=2 cases) to push the test count far higher and extend n.

Definitions taken EXACTLY from the validated fast_hyp1732() in
hyp1732_full_investigation.py (opus-2026-05-23-S5):
  - cycles adjacent in Omega iff they share a vertex
  - max IS = vertex-disjoint pair of cycles
  - C* = a cycle in some max-IS pair but NOT in the pair S being examined
  - B = cycles disjoint from C* (non-adjacent),  p = |B|
  - bound = p*(m-p)

MISTAKE-023 discipline: count DIRECTED odd cycles, not vertex-sets.
"""
import sys, os, time, random
from itertools import combinations

# ---------------------------------------------------------------------------
# Fast directed odd-cycle enumeration (DFS along existing arcs only).
# Each directed cycle is canonicalized by rotating its minimum vertex first;
# its reverse uses opposite arcs (generally absent), so no reverse double count.
# ---------------------------------------------------------------------------
def find_odd_cycles_dfs(A, n):
    """A[i][j]=1 iff i->j. Return list of frozensets? No -- we need vertex sets
    for the conflict graph (shared vertex => adjacent), but distinct DIRECTED
    cycles are distinct Omega-vertices. Return list of (vertex_mask, vertices)."""
    cycles = []  # list of vertex bitmasks (one per DIRECTED odd cycle)
    out = [0]*n
    for i in range(n):
        bm = 0
        for j in range(n):
            if A[i][j]:
                bm |= (1 << j)
        out[i] = bm

    # For each possible minimum start vertex s, DFS over paths using only
    # vertices > s (so s stays the minimum), closing back to s with odd length.
    for s in range(n):
        higher = ~((1 << (s+1)) - 1)  # vertices > s allowed (plus we start at s)
        # path stored as list; visited mask
        def dfs(cur, visited, length, path_mask):
            # try to close
            if length >= 3 and (length % 2 == 1) and (A[cur][s]):
                cycles.append(path_mask)
            # extend to a higher-indexed neighbour not yet visited
            nxt_mask = out[cur] & higher & ~visited
            mm = nxt_mask
            while mm:
                nb = (mm & -mm).bit_length() - 1
                mm &= mm - 1
                dfs(nb, visited | (1 << nb), length + 1, path_mask | (1 << nb))
        dfs(s, (1 << s), 1, (1 << s))
    return cycles

# ---------------------------------------------------------------------------
def conflict_adj(cycle_masks):
    """adj[i] = bitmask of cycles adjacent to i (sharing >=1 vertex)."""
    m = len(cycle_masks)
    adj = [0]*m
    for i in range(m):
        ci = cycle_masks[i]
        bi = 0
        for j in range(m):
            if i != j and (ci & cycle_masks[j]):
                bi |= (1 << j)
        adj[i] = bi
    return adj

def analyze(cycle_masks):
    """Return alpha, alpha2, and (if alpha==2) the HYP-1732 test records."""
    m = len(cycle_masks)
    if m == 0:
        return {'alpha': 0, 'm': 0}
    adj = conflict_adj(cycle_masks)
    # max IS pairs = vertex-disjoint pairs
    max_is = [(i, j) for i in range(m) for j in range(i+1, m)
              if not ((adj[i] >> j) & 1)]
    alpha2 = len(max_is)
    if alpha2 == 0:
        return {'alpha': 1, 'm': m, 'alpha2': 0}
    # detect alpha >= 3 : a triple mutually non-adjacent
    alpha_ge3 = False
    for a, b in max_is:
        # common non-neighbours of a and b
        common = (~adj[a]) & (~adj[b])
        # remove a,b themselves and any index <= b handled; just check existence
        common &= ~((1 << a) | (1 << b))
        if common & ((1 << m) - 1):
            alpha_ge3 = True
            break
    if alpha_ge3:
        return {'alpha': 3, 'm': m, 'alpha2': alpha2}

    # alpha == 2: run pair-partner tests
    tests = 0
    violations = []
    minslack = None
    for s_idx, (c1, c2) in enumerate(max_is):
        S = (1 << c1) | (1 << c2)
        for t_idx, (c3, c4) in enumerate(max_is):
            if t_idx == s_idx:
                continue
            for C_star in (c3, c4):
                if (S >> C_star) & 1:
                    continue  # C* must not be in S
                # B = cycles disjoint from C*  (non-adjacent, excluding C* itself)
                nonadj = (~adj[C_star]) & ((1 << m) - 1)
                nonadj &= ~(1 << C_star)
                p = bin(nonadj).count('1')
                bound = p * (m - p)
                holds_combo = (alpha2 <= bound)
                # equivalent quadratic form: alpha2 <= p(m-p) <=> p^2 - m p + alpha2 <= 0
                quad = p*p - m*p + alpha2
                holds_quad = (quad <= 0)
                if holds_combo != holds_quad:
                    violations.append(('FORM-MISMATCH', m, alpha2, p, bound, quad))
                tests += 1
                slack = bound - alpha2
                minslack = slack if minslack is None else min(minslack, slack)
                if not holds_combo:
                    violations.append(('VIOLATION', m, alpha2, p, bound))
                break  # one C* per partner pair (matches reference logic)
    return {'alpha': 2, 'm': m, 'alpha2': alpha2,
            'tests': tests, 'violations': violations, 'minslack': minslack}

# ---------------------------------------------------------------------------
def random_uniform(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def random_near_transitive(n, p_flip, rng):
    """Transitive base i->j for i<j, each forward arc reversed w.p. p_flip."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < p_flip:
                A[j][i] = 1   # reversed
            else:
                A[i][j] = 1
    return A

# ---------------------------------------------------------------------------
def main():
    rng = random.Random(20260606)
    t0 = time.time()
    TIME_BUDGET = 1500  # leave margin under 30-min runner cap

    print("=" * 70)
    print("HYP-1732 LARGE-SAMPLE VERIFICATION (monad-compute-2026-06-06-S1)")
    print("alpha2(Omega) <= p*(m-p) for all alpha(Omega)=2 tournaments")
    print("=" * 70)

    total_alpha2_cases = 0
    total_tests = 0
    total_violations = 0
    form_mismatches = 0
    global_minslack = None
    per_n = {}

    # Phase 1: uniform random, n=7..11 (independent of prior seed=42 record)
    print("\n--- Phase 1: uniform random tournaments ---")
    # NOTE (S1 run): the original {7:200000,8:120000,...} let n=8 alone (121M
    # tests) consume the full 1500s budget, skipping all later layers. Trimmed
    # so future runs reach Phase 2. The S1 *result artifact* covered n=7,8,9.
    uniform_plan = {7: 20000, 8: 15000, 9: 20000, 10: 15000, 11: 8000}
    for n, count in uniform_plan.items():
        a2cases = ntests = nviol = 0
        minslack = None
        for _ in range(count):
            if time.time() - t0 > TIME_BUDGET:
                break
            A = random_uniform(n, rng)
            cyc = find_odd_cycles_dfs(A, n)
            r = analyze(cyc)
            if r.get('alpha') == 2:
                a2cases += 1
                ntests += r['tests']
                for v in r['violations']:
                    if v[0] == 'FORM-MISMATCH':
                        form_mismatches += 1
                    else:
                        nviol += 1
                if r['minslack'] is not None:
                    minslack = r['minslack'] if minslack is None else min(minslack, r['minslack'])
        per_n[('uniform', n)] = (a2cases, ntests, nviol, minslack)
        total_alpha2_cases += a2cases
        total_tests += ntests
        total_violations += nviol
        if minslack is not None:
            global_minslack = minslack if global_minslack is None else min(global_minslack, minslack)
        print(f"  n={n:2d}: {count:>7d} samples -> {a2cases:>6d} alpha=2 cases, "
              f"{ntests:>7d} tests, {nviol} violations, minslack={minslack}")

    # Phase 2: near-transitive sweep, n=7..13 (harvest many alpha=2 cases)
    print("\n--- Phase 2: near-transitive sampling (p_flip sweep) ---")
    nt_plan = {7: 60000, 8: 60000, 9: 50000, 10: 40000, 11: 30000, 12: 20000, 13: 12000}
    pflips = [0.06, 0.10, 0.14, 0.18, 0.24, 0.30]
    for n, count in nt_plan.items():
        a2cases = ntests = nviol = 0
        minslack = None
        per_pf = count // len(pflips)
        for pf in pflips:
            for _ in range(per_pf):
                if time.time() - t0 > TIME_BUDGET:
                    break
                A = random_near_transitive(n, pf, rng)
                cyc = find_odd_cycles_dfs(A, n)
                # near-transitive can still blow up; skip absurdly large cycle sets
                if len(cyc) > 600:
                    continue
                r = analyze(cyc)
                if r.get('alpha') == 2:
                    a2cases += 1
                    ntests += r['tests']
                    for v in r['violations']:
                        if v[0] == 'FORM-MISMATCH':
                            form_mismatches += 1
                        else:
                            nviol += 1
                    if r['minslack'] is not None:
                        minslack = r['minslack'] if minslack is None else min(minslack, r['minslack'])
        per_n[('near-trans', n)] = (a2cases, ntests, nviol, minslack)
        total_alpha2_cases += a2cases
        total_tests += ntests
        total_violations += nviol
        if minslack is not None:
            global_minslack = minslack if global_minslack is None else min(global_minslack, minslack)
        print(f"  n={n:2d}: {count:>7d} samples -> {a2cases:>6d} alpha=2 cases, "
              f"{ntests:>7d} tests, {nviol} violations, minslack={minslack}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total alpha(Omega)=2 tournaments found: {total_alpha2_cases}")
    print(f"  Total pair-partner HYP-1732 tests:      {total_tests}")
    print(f"  Total VIOLATIONS (alpha2 > p(m-p)):     {total_violations}")
    print(f"  Form-equivalence mismatches:            {form_mismatches}")
    print(f"  Global min slack (bound - alpha2):      {global_minslack}")
    print(f"  Elapsed: {time.time()-t0:.1f}s")
    # Honest coverage report: which (mode,n) layers actually produced tests.
    covered = sorted({n for (mode, n), v in per_n.items() if v[1] > 0})
    skipped = sorted({n for (mode, n), v in per_n.items() if v[1] == 0})
    print(f"  n-levels with >0 tests (actually covered): {covered}")
    print(f"  n-levels with 0 tests (no alpha=2 found OR budget-skipped): {skipped}")
    if total_violations == 0 and form_mismatches == 0:
        print(f"\n  RESULT: HYP-1732 HOLDS in all {total_tests} tests "
              f"({total_alpha2_cases} alpha=2 tournaments).")
        print(f"  Coverage by n (with >0 tests): {covered}.")
        print(f"  Both forms (combinatorial bound and quadratic I(-1/p)<=0) agree.")
        print(f"  NOTE: a layer showing 0 tests may have been BUDGET-SKIPPED, not "
              f"empty -- check elapsed vs TIME_BUDGET above.")
    else:
        print(f"\n  RESULT: PROBLEM DETECTED -- inspect violations above.")

if __name__ == '__main__':
    main()
