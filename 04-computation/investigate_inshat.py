#!/usr/bin/env python3
"""
Investigate the paper's Remark (line 1261-1264):
"Claim A is equivalent to Σ(inshat-1)/2 = Σ mu(C)"

This requires H(T)-H(T-v) = Σ(inshat-1) = 2*Σ TypeII.
But H(T)-H(T-v) = Σ TypeII + #{orphans} (my derivation).

So we need to check: does Σ(inshat-1)/2 = Σ mu(C) hold for all (T,v)?
If Claim A holds (verified) but this doesn't, the paper's remark is incorrect.

Also investigate OPEN-Q-010: per-path formula with 3-cycles + 5-cycles at n=7.

Session: opus-2026-03-05-S3
"""
import sys
sys.path.insert(0, '/home/e/Documents/claude/math/03-artifacts/code')
from tournament_lib import *
from itertools import combinations


def inshat(T, v, path):
    """Compute inshat(v, P') for a Hamiltonian path P' of T-v.
    path: list of vertices in T-v (0-indexed in T's vertex space, but v removed).
    Returns inshat value and debug info.
    """
    n_tv = len(path)  # = n-1 vertices in T-v path
    # Signature: s[j] = T(v, path[j]) = 1 if v->path[j]
    s = [T[v][path[j]] for j in range(n_tv)]

    # Boundary term: [s[0]=1] + [s[-1]=0]
    b = (1 if s[0] == 1 else 0) + (1 if s[-1] == 0 else 0)

    # Transitions: all j where s[j] != s[j+1]
    transitions = sum(1 for j in range(n_tv - 1) if s[j] != s[j + 1])
    typeI  = sum(1 for j in range(n_tv - 1) if s[j] == 0 and s[j+1] == 1)
    typeII = sum(1 for j in range(n_tv - 1) if s[j] == 1 and s[j+1] == 0)

    ins = b + transitions
    insact_val = (1 if s[0] == 1 else 0) + (1 if s[-1] == 0 else 0) + typeI

    return ins, typeI, typeII, insact_val, s


def count_orphans(T, v):
    """Count orphan Ham paths: paths P of T where P-v is not in Ham(T-v)."""
    n = len(T)
    Tv, old_labels = delete_vertex(T, v)
    tv_ham_set = set()
    for perm in __import__('itertools').permutations(range(n-1)):
        valid = True
        for i in range(n-2):
            if not Tv[perm[i]][perm[i+1]]:
                valid = False
                break
        if valid:
            tv_ham_set.add(tuple(old_labels[x] for x in perm))

    orphans = 0
    for perm in __import__('itertools').permutations(range(n)):
        # Check if it's a Ham path of T
        valid = True
        for i in range(n-1):
            if not T[perm[i]][perm[i+1]]:
                valid = False
                break
        if not valid:
            continue
        # Remove v and check if remainder is in Ham(T-v)
        remainder = tuple(x for x in perm if x != v)
        if remainder not in tv_ham_set:
            orphans += 1
    return orphans


def get_ham_paths_tv(T, v):
    """Get all Ham paths of T-v as sequences of original vertex labels."""
    n = len(T)
    Tv, old_labels = delete_vertex(T, v)
    n_tv = n - 1
    paths = []
    for perm in __import__('itertools').permutations(range(n_tv)):
        valid = all(Tv[perm[i]][perm[i+1]] for i in range(n_tv - 1))
        if valid:
            paths.append([old_labels[x] for x in perm])
    return paths


def check_inshat_equivalence(T, v):
    """
    Check whether Σ(inshat-1)/2 = Σ mu(C).
    Also check H(T)-H(T-v) = Σ TypeII + #{orphans}.
    Returns dict with all relevant quantities.
    """
    n = len(T)
    ht = hamiltonian_path_count(T)
    Tv, old_labels = delete_vertex(T, v)
    htv = hamiltonian_path_count(Tv)
    diff = ht - htv

    # Compute Σ mu(C) via Claim A
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)
    total_mu = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)
    claim_a_holds = (diff == 2 * total_mu)

    # Compute Σ(inshat-1)/2 and Σ TypeII over all P' ∈ Ham(T-v)
    ham_paths = get_ham_paths_tv(T, v)
    sum_inshat_m1_over2 = 0
    sum_typeII = 0
    sum_insact = 0
    for path in ham_paths:
        ins, ti, tii, ia, s = inshat(T, v, path)
        assert (ins - 1) % 2 == 0, f"inshat not odd: {ins}"
        sum_inshat_m1_over2 += (ins - 1) // 2
        sum_typeII += tii
        sum_insact += ia
        assert (ins - 1) // 2 == tii, f"(inshat-1)/2 != typeII: {(ins-1)//2} vs {tii}"

    # Count orphans
    num_orphans = ht - sum_insact

    # The formula: H(T)-H(T-v) = Σ TypeII + #{orphans}
    formula_holds = (diff == sum_typeII + num_orphans)

    # Paper's claimed equivalence: Σ(inshat-1)/2 = Σ mu
    paper_remark_holds = (sum_inshat_m1_over2 == total_mu)

    return {
        'ht': ht, 'htv': htv, 'diff': diff,
        'total_mu': total_mu, 'two_mu': 2*total_mu,
        'claim_a': claim_a_holds,
        'sum_inshat_m1_over2': sum_inshat_m1_over2,
        'sum_typeII': sum_typeII,
        'sum_insact': sum_insact,
        'num_orphans': num_orphans,
        'formula_holds': formula_holds,
        'paper_remark_holds': paper_remark_holds,
        'n_cycles_v': len(cycles_v),
    }


def test_specific_case():
    """Test n=4, v=4, transitive T-v, v beats only vertex 0."""
    # T on {0,1,2,3}: T-v = transitive on {0,1,2} (0->1->2, 0->2),
    # v=3 beats only vertex 0.
    T = [[0]*4 for _ in range(4)]
    # Transitive on {0,1,2}:
    T[0][1] = 1; T[1][0] = 0
    T[0][2] = 1; T[2][0] = 0
    T[1][2] = 1; T[2][1] = 0
    # v=3 beats only 0:
    T[3][0] = 1; T[0][3] = 0
    T[1][3] = 1; T[3][1] = 0  # 1 beats 3
    T[2][3] = 1; T[3][2] = 0  # 2 beats 3

    v = 3
    result = check_inshat_equivalence(T, v)
    print("=== Specific test case (n=4, transitive T-v, v beats only 0) ===")
    for k, val in result.items():
        print(f"  {k}: {val}")
    print()
    return result


def exhaustive_test(n, max_pairs=None):
    """Exhaustive test of inshat equivalence for all (T,v) at size n."""
    print(f"=== Exhaustive test n={n} ===")
    pairs = 0
    claim_a_fails = 0
    paper_remark_fails = 0
    formula_fails = 0

    for T in all_tournaments(n):
        for v in range(n):
            if max_pairs and pairs >= max_pairs:
                break
            r = check_inshat_equivalence(T, v)
            pairs += 1
            if not r['claim_a']:
                claim_a_fails += 1
            if not r['paper_remark_holds']:
                paper_remark_fails += 1
            if not r['formula_holds']:
                formula_fails += 1

    print(f"  Pairs tested: {pairs}")
    print(f"  Claim A failures: {claim_a_fails}")
    print(f"  Paper Remark [Σ(inshat-1)/2 = Σ mu] failures: {paper_remark_fails}")
    print(f"  Formula [H-Hv = ΣTypeII + #orphans] failures: {formula_fails}")
    print()
    return paper_remark_fails, formula_fails


def investigate_open_q010_n7():
    """
    OPEN-Q-010: At n=7, mu(5-cycle) = 1 always (V\{v + 4 cycle vertices} has 2 vertices).
    Test: does Σ_{P'} Σ_{TypeII at j for 3-cycles} mu(v,P'[j],P'[j+1])
          + Σ_{P'} Σ_{TypeII at j for 5-cycles} 1
          = Σ_C mu(C)?
    Specifically: does a per-path formula counting (a) TypeII positions
    corresponding to 3-cycles with their mu weights, plus (b) TypeII positions
    corresponding to 5-cycles (mu=1), sum to Σ mu?

    Note: at n=7, a 5-cycle through v uses v + 4 other vertices, leaving 2.
    Those 2 vertices form a subtournament with no odd cycles (2 vertices can't
    have a 3-cycle), so mu(5-cycle) = 1 always. ✓

    But we also need to handle orphan paths. This function investigates the
    relationship between the per-path contributions and Σ mu.
    """
    import random
    random.seed(42)
    n = 7

    print("=== OPEN-Q-010: n=7 per-path formula investigation ===")

    failures_3only = 0
    failures_3plus5 = 0
    total = 0

    for _ in range(50):
        T = random_tournament(n)
        for v in range(n):
            Tv, old_labels = delete_vertex(T, v)
            tv_cycles = find_odd_cycles(Tv)
            cache = (Tv, old_labels, tv_cycles)

            all_cyc = find_odd_cycles(T)
            cycles_v = [c for c in all_cyc if v in set(c)]
            total_mu_rhs = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)

            ht = hamiltonian_path_count(T)
            htv = hamiltonian_path_count(Tv)
            diff = ht - htv

            # LHS using only 3-cycles
            sum_lhs_3 = 0
            sum_lhs_35 = 0
            ham_paths = get_ham_paths_tv(T, v)

            for path in ham_paths:
                ins, ti, tii, ia, s = inshat(T, v, path)
                npath = len(path)
                for j in range(npath - 1):
                    if s[j] == 1 and s[j+1] == 0:
                        # TypeII position: corresponds to a 3-cycle v->path[j]->path[j+1]->v
                        a, b = path[j], path[j+1]
                        # 3-cycle exists by TypeII def (T[v][a]=1, T[a][b]=1, T[b][v]=1)
                        cyc3 = (v, a, b) if v < a and v < b else None
                        # compute mu for this 3-cycle
                        mu3 = mu(T, v, (v, a, b), _tv_cache=cache)
                        sum_lhs_3 += mu3
                        sum_lhs_35 += mu3  # same for 3-cycles in 3+5 formula

            # For 5-cycles at n=7: mu = 1 always. But how do they contribute
            # to H(T)-H(T-v)? They show up as "orphan" paths or through
            # a different mechanism. Let's just check if sum_typeII = Σmu for 3-cycles.
            # That's the per-path identity. And if we add orphan contributions:
            num_orphans = ht - sum(inshat(T, v, p)[3] for p in ham_paths)
            sum_typeII_total = sum(inshat(T, v, p)[2] for p in ham_paths)

            # Check per-path 3-cycle identity (from THM-005 applied per-path):
            # Does Σ_{P'} #{TypeII} = Σ_{3-cycles C through v} #{P': V(C)\{v} consecutive in P'}?
            # This should always hold by THM-005.

            # The key question: does diff = 2*total_mu_rhs?
            assert diff == 2 * total_mu_rhs, f"Claim A fails at n=7!"

            # Does the per-path 3-cycle sum (sum_lhs_3) = Σ mu?
            if sum_lhs_3 != total_mu_rhs:
                failures_3only += 1

            total += 1

    print(f"  Pairs: {total}")
    print(f"  Per-path 3-cycle sum ≠ Σmu (expected failures): {failures_3only}")
    print()


if __name__ == '__main__':
    # First: specific test to confirm the discrepancy
    r = test_specific_case()
    assert not r['paper_remark_holds'], "Expected paper remark to fail!"
    assert r['claim_a'], "Expected Claim A to hold!"
    assert r['formula_holds'], "Expected H-Hv = ΣTypeII + #orphans!"
    print("✓ Confirmed: paper remark fails, but Claim A holds and H-Hv = ΣTypeII + #orphans")
    print()

    # Exhaustive tests
    for n in [3, 4, 5]:
        paper_fails, formula_fails = exhaustive_test(n)
        assert formula_fails == 0, f"Formula H-Hv = ΣTypeII + #orphans fails at n={n}!"

    print("✓ Formula H(T)-H(T-v) = Σ_{P'} TypeII(P') + #{orphans} holds exhaustively n=3,4,5")
    print()

    # Count paper remark failures at n=4 (should be many)
    print("=== Summary: paper remark failures at n=4 ===")
    fails = 0
    total = 0
    for T in all_tournaments(4):
        for v in range(4):
            r = check_inshat_equivalence(T, v)
            assert r['claim_a'], "Claim A fails!"
            assert r['formula_holds'], "Formula fails!"
            if not r['paper_remark_holds']:
                fails += 1
            total += 1
    print(f"  Total pairs: {total}")
    print(f"  Paper remark [Σ(inshat-1)/2 = Σmu] failures: {fails}/{total}")
    print(f"  (These are pairs where Σ TypeII ≠ Σ mu, i.e., #{orphans} ≠ Σ mu)")
    print()

    # OPEN-Q-010
    investigate_open_q010_n7()
