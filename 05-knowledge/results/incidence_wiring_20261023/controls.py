"""Hostile/positive controls for the all-size mechanisms; explicit checks."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import argparse, json, math
import check as core

HERE = Path(__file__).resolve().parent
need = core.need

def matrix(n, word):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        A[i][j] = (word >> k)&1; A[j][i] = 1-A[i][j]
    return A

def partitions(n):
    if not n: yield []; return
    for p in partitions(n-1):
        for i in range(len(p)): yield [B+[n-1] if j == i else B[:] for j, B in enumerate(p)]
        yield p+[[n-1]]

def hp_start(A):
    n = len(A); dp = [[0]*n for _ in range(1 << n)]; dp[1][0] = 1
    out = [sum(a << j for j, a in enumerate(row)) for row in A]
    for m in range(1, 1 << n, 2):
        for i in range(n):
            if not dp[m][i]: continue
            z = out[i] & ~m
            for j in core.bits(z): dp[m | (1 << j)][j] += dp[m][i]
    return n*sum(dp[-1])

def hp_all(A):
    # Separate recurrence: paths spanning a subset, indexed by first vertex.
    n = len(A); dp = [[0]*n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for m in range(1, 1 << n):
        for i in core.bits(m):
            rest = m ^ (1 << i)
            if rest: dp[m][i] = sum(A[i][j]*dp[rest][j] for j in core.bits(rest))
    return sum(dp[-1])

def refinement(A, marked):
    n = len(A); sig, rho, *_ = core.profiles(A)
    tags = [(r, int(i == 0 and marked)) for i, r in enumerate(rho)]
    history = []
    for _ in range(n+1):
        keys = sorted(set(tags)); colors = [keys.index(t) for t in tags]
        history.append(sorted(Counter(colors).values()))
        following = [(colors[u], tuple(sorted(colors[v] for v in range(n) if A[u][v]))) for u in range(n)]
        if len(set(following)) == len(set(tags)): break
        tags = following
    else: raise ValueError('refinement did not stabilize')
    return history

def main():
    p = argparse.ArgumentParser(); p.add_argument('--write', action='store_true'); args = p.parse_args()
    count = partcount = decoders = 0
    for n in range(1, 6):
        pp = list(partitions(n))
        for word in range(1 << (n*(n-1)//2)):
            A = matrix(n, word); count += 1
            for cells in pp:
                for C in cells:
                    if len(C)%2 == 0:
                        degrees = [sum(A[u][v] for v in C) for u in C]
                        need(len(set(degrees)) > 1, 'even cell survives its internal-degree refinement')
                partcount += 1
            s = list(map(sum, A))
            for u, v in combinations(range(n), 2):
                if s[u] != s[v]: continue
                cells = [[u, v]]+[[w] for w in range(n) if w not in (u, v)]
                E = tuple(tuple(sum(A[a][b] for a in C for b in D) for D in cells) for C in cells)
                rho = [((s[w],), ()) for w in range(n)]
                result = core.decode_doubleton(A, rho, cells, E)
                t = len(result['mixed_observer_types'])
                need(t%2 == 1, 'mixed observer number must be odd')
                decoders += 1
    parity = []
    for t in range(1, 10, 2):
        rows = []
        for x in product((0, 1), repeat=t):
            if abs(2*sum(x)-t) != 1: continue
            # Complement every bit is a swap of the two hidden vertices.
            rows.append(min(x, tuple(1-a for a in x)))
        need(len(set(rows)) == math.comb(t, (t-1)//2), 'balanced orbit count')
        parity.append([t, len(set(rows)), t-1])
    # An inconsistent observer triangle has odd XOR around its only cycle.
    hostile = {(0, 1): 1, (1, 2): 1, (0, 2): 0}
    need(not any(all(int(x[a] == x[b]) == w for (a, b), w in hostile.items()) for x in product((0, 1), repeat=3)), 'triangle parity hostile accepted')

    # All 64 simple blue graphs on four vertices, and each of three switch cuts.
    switch_cases = Counter()
    pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    for mask in range(64):
        L = [[0]*4 for _ in range(4)]
        for k, (i, j) in enumerate(combinations(range(4), 2)):
            if mask >> k & 1:
                L[i][i] += 1; L[j][j] += 1; L[i][j] -= 1; L[j][i] -= 1
        rows = [sum((L[i][j]&1) << (j-1) for j in range(1, 4)) for i in range(1, 4)]
        V, _ = core.binary_kernel(rows, 3); V = [x << 1 for x in V]
        for (a, b), (c, d) in pairings:
            alpha = sum((((x>>a)^(x>>b))&1) << i for i, x in enumerate(V))
            beta = sum((((x>>c)^(x>>d))&1) << i for i, x in enumerate(V))
            independent = bool(alpha and beta and alpha != beta)
            delta = []
            for x in V:
                row = 0
                for j, y in enumerate(V):
                    value = 0
                    for u, v, w in ((a, c, 1), (b, d, 1), (a, d, -1), (b, c, -1)):
                        value += w*(((x>>u)&1)-((x>>v)&1))*(((y>>u)&1)-((y>>v)&1))
                    row |= (value&1) << j
                delta.append(row)
            rank = core.binary_rank(delta)
            need(rank == (2 if independent else 0), 'electrical switch rank criterion')
            need(all(((row >> i)&1) == 0 for i, row in enumerate(delta)), 'degree preserving switch should have zero diagonal')
            switch_cases[rank] += 1

    # The historical p=13 additive-energy collision is a stable-refinement hostile.
    fixtures = []
    for S in ({2, 7, 8, 9, 10, 12}, {1, 2, 4, 7, 8, 10}):
        n = 13; A = [[int((j-i)%n in S) for j in range(n)] for i in range(n)]
        sig, rho, *_ = core.profiles(A)
        H = hp_start(A); need(H == hp_all(A), 'independent Hamiltonian recurrences disagree')
        need(len(set(rho)) == 1 and rho[0][0] == (6, 21, 112), 'uniform p13 profile')
        h0, h1 = refinement(A, False), refinement(A, True)
        need(h0 == [[13]] and h1[-1] == [1]*13, 'marked refinement boundary')
        diff = [sum((a-b)%13 == k for a in S for b in S) for k in range(13)]
        sums = [sum((a+b)%13 == k for a in S for b in S) for k in range(13)]
        need(diff[0] == 6 and sums[0] == 0, 'sum and difference zero fibres confused')
        need(sum(x*x for x in diff) == sum(x*x for x in sums) == 118, 'additive energy equality')
        fixtures.append({'steps': sorted(S), 'tag': rho[0][0], 'H': H, 'unrooted_cells': h0, 'source_marked_cells': h1,
                         'difference_energy_2': sum(x*x for x in diff), 'difference_energy_3': sum(x*x*x for x in diff),
                         'sum_energy_2': sum(x*x for x in sums), 'sum_energy_3': sum(x*x*x for x in sums)})
    need([x['H'] for x in fixtures] == [3683797, 3704857], 'historical H witnesses changed')
    energy_cases = 0
    for n in (3, 5, 7, 9, 11, 13):
        for choices in product((0, 1), repeat=(n-1)//2):
            S = {k if b else n-k for k, b in enumerate(choices, 1)}
            sr = [sum((a+b)%n == k for a in S for b in S) for k in range(n)]
            dr = [sum((a-b)%n == k for a in S for b in S) for k in range(n)]
            direct = sum((a+b-c)%n in S for a in S for b in S for c in S)
            need(sr[0] == 0 and dr[0] == len(S), 'skew set zero fibres')
            need(sum(x*x for x in sr[1:]) == sum(x*x for x in dr) == direct, 'historical energy repair')
            energy_cases += 1
    out = {'status': 'FINITE-EXACT controls for elementary all-size proofs', 'labelled_tournaments': count,
           'tournament_partition_pairs': partcount, 'equal_degree_doubleton_decodings': decoders,
           'balanced_assignment_counts_t_orbits_treebits': parity, 'electrical_switch_cases_by_rank': dict(switch_cases),
           'stable_refinement_hostile_and_source_repair': fixtures, 'historical_energy_repair_skew_sets': energy_cases, 'checks': core.CHECKS}
    # JSON normalization makes integer keys/tuples compare identically under -O.
    out = json.loads(json.dumps(out))
    dest = HERE/'controls.json'
    if args.write: dest.write_text(json.dumps(out, indent=2)+'\n')
    else: need(json.loads(dest.read_text()) == out, 'control receipt differs')
    print(json.dumps(out, indent=2))

if __name__ == '__main__': main()
