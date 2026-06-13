"""Explicit silent-cut construction at C=105 (3-prime squarefree) to test the chain prediction.
monad-explorer-S708.

The A_t*B_t=0 lemma lets us VERIFY any candidate cut in O(C^2) without enumerating 2^51 cuts.
We construct, via coset-tiling backtracking, a cut eps where flipping a subgroup half H_d is silent
(full-coset condition), and then LOCAL-SEARCH for a cut where a divisor-CHAIN combined move
H_d (+) H_e (d|e) is silent -- which (if found) confirms primitive combined silent moves exist at a
3-prime squarefree modulus, refuting 'squarefree => primes independent'.
"""
import numpy as np
import math
import random

random.seed(12345)


def build_S(C):
    n1 = (C - 1) // 2
    Hf = np.arange(1, n1 + 1)
    return n1, np.sin(2 * math.pi * np.outer(Hf, Hf) / C)


def halfmask(C, d, n1):
    g = C // d
    Kd = set((g * j) % C for j in range(d))
    mask = np.zeros(n1, dtype=bool)
    for x in Kd:
        if 0 < x <= n1:
            mask[x - 1] = True
    return mask


def is_silent(eps, Dmask, S):
    """Lemma: silent iff A_t*B_t=0 for all t.  eps in {+-1}^n1, Dmask bool."""
    Phi = S @ eps
    B = S[:, Dmask] @ eps[Dmask]
    A = Phi - B
    changed = Dmask.any()
    return changed and np.all(np.abs(A * B) < 1e-6)


def construct_fullcoset_cut(C, d, n1, tries=200000):
    """Find eps where signed nonK points form full K_d-cosets (=> flipping H_d silent).
    Coset = {r, r+C/d, ..., r+(d-1)C/d}. We choose signs so chosen points tile full cosets.
    Approach: for each runner i not in K_d, point is +-i; greedily/backtracking assign so the set
    of chosen residues is a union of full cosets."""
    g = C // d
    Kd = set((g * j) % C for j in range(d))
    # cosets of K_d (excluding K_d itself) as frozensets of residues in 1..C-1
    cosets = {}
    for r in range(1, C):
        if r in Kd:
            continue
        cs = frozenset((r + g * j) % C for j in range(d))
        cosets[min(cs)] = cs
    coset_list = list(cosets.values())
    # each runner i (1..n1, i not in Kd) can hit residue i (sign +) or C-i (sign -)
    runners = [i for i in range(1, n1 + 1) if i not in Kd]
    # map residue -> runner & sign
    res_to_run = {}
    for i in runners:
        res_to_run[i] = (i, +1)
        res_to_run[C - i] = (i, -1)
    # We must pick exactly one residue per runner s.t. chosen residue-set = union of full cosets.
    # Backtracking over cosets: a coset is either fully chosen or fully unchosen; consistency = each
    # runner used exactly once.  Build conflict structure.
    # Greedy randomized: shuffle cosets, try to select a maximal set tiling all runners.
    for _ in range(2000):
        random.shuffle(coset_list)
        used_runner = {}
        chosen = []
        ok = True
        # try to cover ALL runners: pick cosets greedily; each runner must end covered exactly once
        # Represent: choose subset of cosets that partitions all runner-residues.
        # A runner i is covered if exactly one of residues {i, C-i} is in a chosen coset.
        run_cover = {i: 0 for i in runners}
        for cs in coset_list:
            # would adding cs keep run_cover <=1 for each runner?
            inc = {}
            good = True
            for res in cs:
                if res in res_to_run:
                    run, sg = res_to_run[res]
                    inc[run] = inc.get(run, 0) + 1
                else:
                    good = False; break  # residue not a runner residue (shouldn't happen)
            if not good:
                continue
            if all(run_cover[r] + inc[r] <= 1 for r in inc):
                for r in inc:
                    run_cover[r] += inc[r]
                chosen.append(cs)
        if all(run_cover[r] == 1 for r in runners):
            # build eps
            eps = np.ones(n1, dtype=np.int64)
            for cs in chosen:
                for res in cs:
                    run, sg = res_to_run[res]
                    eps[run - 1] = sg
            return eps
    return None


def main():
    C = 105
    n1, S = build_S(C)
    print(f"C={C}=3x5x7, n1={n1}, 2^(n-2)=2^{n1-1} (enumeration infeasible; verify via A_t*B_t lemma)")
    # 1) single-subgroup silent cuts (sanity: machinery works at 105)
    for d in [3, 5, 7, 15, 21, 35]:
        Dmask = halfmask(C, d, n1)
        eps = construct_fullcoset_cut(C, d, n1)
        if eps is None:
            print(f"   H_{d}: construction failed (try larger search)")
            continue
        sil = is_silent(eps, Dmask, S)
        print(f"   H_{d}: constructed full-coset cut; flip-H_{d} silent (A_t*B_t=0)? {sil}")

    # 2) chain combined move: find cut where H_d (+) H_e (d|e) is silent, via local search
    print("\n   Chain combined moves (local search for a silent cut):")
    chains = [(3, 15), (5, 15), (7, 35), (3, 21), (5, 35), (7, 21)]
    for (d, e) in chains:
        Dmask = halfmask(C, d, n1) ^ halfmask(C, e, n1)  # combined move H_d (+) H_e
        # local search: minimize sum (A_t B_t)^2
        best = None
        for restart in range(40):
            eps = np.array([random.choice([-1, 1]) for _ in range(n1)], dtype=np.int64)
            eps[0] = 1
            def viol(e_):
                Phi = S @ e_; B = S[:, Dmask] @ e_[Dmask]; A = Phi - B
                return float(np.sum((A * B) ** 2))
            cur = viol(eps)
            improved = True
            steps = 0
            while improved and steps < 4000:
                improved = False; steps += 1
                order = list(range(1, n1))
                random.shuffle(order)
                for k in order:
                    eps[k] = -eps[k]
                    v = viol(eps)
                    if v < cur - 1e-9:
                        cur = v; improved = True
                    else:
                        eps[k] = -eps[k]
                if cur < 1e-6:
                    break
            if best is None or cur < best[0]:
                best = (cur, eps.copy())
            if cur < 1e-6:
                break
        sil = is_silent(best[1], Dmask, S)
        print(f"      H_({d})(+)H_({e}) [chain {d}|{e}]: min violation={best[0]:.3e}  silent-cut found? {sil}")


if __name__ == "__main__":
    main()
