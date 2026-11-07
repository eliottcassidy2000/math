"""S85: finite checks of the tournament sites and the plus construction.

Sites (all on the poset of subsets of the vertex set V of a fixed tournament T, morphisms = inclusions):
  J2  (pair topology, on finite sets): a sieve on S covers iff it contains every subset of S of size <= 2.
  Jss (strong-component topology):     a sieve on S covers iff it contains every strong component of T[S].
Plus construction on a site whose covering sieves have a smallest member R_min(S):
  F^+(S) = Match(R_min(S), F) = compatible families on R_min(S).
Checks:
  (1) J2: L (linear orders) is separated, L^+ = tournaments (2^C(n,2)), image of L = transitive ones.
  (2) Jss axioms: strong components of T[S'] lie inside strong components of T[S] for S' in S (pullback stability).
  (3) Pk2 (2-coloured packings of vertex-disjoint odd cycles): |Pk2(V)| = H(T) (Grinberg-Stanley) and the sheaf
      condition Pk2(S) = prod over strong components.
  (4) Constant presheaf A (also at the empty set): A^+(S) = A (S nonempty), A^{++}(S) = A^{c(T[S])}; A^+ is not a sheaf
      when c >= 2 and |A| >= 2 - exactly two plus steps.
  (5) Q = Pk2 / (global colour swap): not separated; one plus step gives prod_i (H(C_i)+1)/2, and the fibres of
      Q -> Q^+ are the relative-colour torsors of size 2^(m-1), m = number of components carrying a nonempty packing.
"""
import itertools, random, sys
from math import comb
from collections import Counter


def random_tournament(n, rng):
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for m in range(1 << len(pairs)):
        adj = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (m >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj


def ham_paths(adj, S):
    S = list(S)
    if not S:
        return 1
    idx = {v: k for k, v in enumerate(S)}
    m = len(S)
    dp = [[0] * m for _ in range(1 << m)]
    for k in range(m):
        dp[1 << k][k] = 1
    for M in range(1, 1 << m):
        for a in range(m):
            x = dp[M][a]
            if not x:
                continue
            for b in range(m):
                if not (M >> b) & 1 and adj[S[a]][S[b]]:
                    dp[M | (1 << b)][b] += x
    return sum(dp[(1 << m) - 1])


def strong_components(adj, S):
    S = list(S)
    reach = {u: {u} for u in S}
    for u in S:
        stack = [u]
        while stack:
            x = stack.pop()
            for y in S:
                if adj[x][y] and y not in reach[u]:
                    reach[u].add(y)
                    stack.append(y)
    comps = []
    seen = set()
    for u in S:
        if u in seen:
            continue
        c = frozenset(v for v in S if v in reach[u] and u in reach[v])
        comps.append(c)
        seen |= c
    return comps


def odd_cycles(adj, n):
    """directed odd cycles as (frozenset of vertices, canonical tuple)"""
    cyc = set()
    for k in range(3, n + 1, 2):
        for verts in itertools.combinations(range(n), k):
            first = verts[0]
            for perm in itertools.permutations(verts[1:]):
                seq = (first,) + perm
                if all(adj[seq[i]][seq[(i + 1) % k]] for i in range(k)):
                    cyc.add(seq)
    return sorted(cyc)


def packings(cycles, S):
    """all sets of pairwise vertex-disjoint cycles inside S (as tuples of cycle indices)"""
    inside = [i for i, c in enumerate(cycles) if set(c) <= set(S)]
    res = []

    def rec(start, used, chosen):
        res.append(tuple(chosen))
        for t in range(start, len(inside)):
            i = inside[t]
            vs = set(cycles[i])
            if vs & used:
                continue
            rec(t + 1, used | vs, chosen + [i])
    rec(0, set(), [])
    return res


def check_J2(nmax=5):
    out = []
    for n in range(1, nmax + 1):
        V = list(range(n))
        small = [frozenset(U) for k in range(0, 3) for U in itertools.combinations(V, k)]
        # matching families of linear orders on <=2-subsets: a linear order on a pair is an ordered pair
        pairs = [U for U in small if len(U) == 2]
        lplus = 2 ** len(pairs)          # each pair independently ordered, singletons/empty have one order
        # image of L: linear orders of V give transitive tournaments; count distinct images
        images = set()
        for perm in itertools.permutations(V):
            pos = {v: i for i, v in enumerate(perm)}
            images.add(tuple(pos[a] < pos[b] for a, b in itertools.combinations(V, 2)))
        out.append((n, lplus, 2 ** comb(n, 2), len(images), len(list(itertools.permutations(V)))))
    return out


def check_strong_site(trials=300, seed=1):
    rng = random.Random(seed)
    viol = 0
    for _ in range(trials):
        n = rng.randint(2, 7)
        adj = random_tournament(n, rng)
        V = list(range(n))
        S = [v for v in V if rng.random() < 0.8]
        comps = strong_components(adj, S)
        for k in range(len(S) + 1):
            for Sp in itertools.combinations(S, k):
                for c in strong_components(adj, Sp):
                    if not any(c <= C for C in comps):
                        viol += 1
    return viol


def plus(F_sets, F_restrict, gens_of, objects):
    """generic plus construction on a site whose minimal covering sieve of S is generated by gens_of(S).
    F_sets[S] = list of elements, F_restrict(S, x, U) = restriction of x in F(S) to U (U subset S).
    Returns Fp_sets[S] = list of matching families (tuple over gens, in order), and Fp_restrict."""
    Fp = {}
    for S in objects:
        gens = gens_of(S)
        fams = []
        for choice in itertools.product(*[F_sets[g] for g in gens]):
            ok = True
            for a in range(len(gens)):
                for b in range(a + 1, len(gens)):
                    I = gens[a] & gens[b]
                    if F_restrict(gens[a], choice[a], I) != F_restrict(gens[b], choice[b], I):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                fams.append(choice)
        Fp[S] = fams

    def Fp_restrict(S, fam, U):
        gens = gens_of(S)
        gU = gens_of(U)
        # value on each generator of U: restrict the family member of the generator of S containing it
        res = []
        for g in gU:
            host = next(k for k, G in enumerate(gens) if g <= G)
            res.append(F_restrict(gens[host], fam[host], g))
        return tuple(res)
    return Fp, Fp_restrict


def check_two_step(trials=40, seed=2, A=('x', 'y')):
    rng = random.Random(seed)
    report = Counter()
    for _ in range(trials):
        n = rng.randint(3, 6)
        adj = random_tournament(n, rng)
        V = list(range(n))
        objects = [frozenset(U) for k in range(n + 1) for U in itertools.combinations(V, k)]
        gens_of = lambda S: sorted(strong_components(adj, S), key=lambda c: sorted(c)) if S else []
        F_sets = {S: list(A) for S in objects}
        F_restrict = lambda S, x, U: x
        Fp, Fp_r = plus(F_sets, F_restrict, gens_of, objects)
        Fpp, _ = plus(Fp, Fp_r, gens_of, objects)
        for S in objects:
            c = len(gens_of(S))
            p1, p2 = len(Fp[S]), len(Fpp[S])
            exp1 = 1 if not S else len(A)
            exp2 = len(A) ** c
            report['A+ ok' if p1 == exp1 else 'A+ BAD'] += 1
            report['A++ ok' if p2 == exp2 else 'A++ BAD'] += 1
            if c >= 2:
                report['two steps needed (|A+| != |A++|)' if p1 != p2 else 'NOT needed?'] += 1
    return dict(report)


def check_packings_and_Q(nmax=6, seed=3, sample7=60):
    rng = random.Random(seed)
    res = Counter()
    fib = Counter()
    tours = []
    for n in range(1, nmax + 1):
        if n <= 5:
            tours += [(n, adj) for adj in all_tournaments(n)]
        else:
            tours += [(n, random_tournament(n, rng)) for _ in range(400)]
    tours += [(7, random_tournament(7, rng)) for _ in range(sample7)]
    for n, adj in tours:
        V = list(range(n))
        cyc = odd_cycles(adj, n)
        pk = packings(cyc, V)
        h = ham_paths(adj, V)
        if sum(2 ** len(p) for p in pk) == h:
            res['GS: sum 2^k alpha_k = H'] += 1
        else:
            res['GS FAIL'] += 1
        comps = strong_components(adj, V)
        prodH = 1
        for c in comps:
            prodH *= ham_paths(adj, c)
        res['H multiplicative' if prodH == h else 'mult FAIL'] += 1
        # Q = Pk2 / global swap: elements = (packing, colouring) mod swap; Q^+ = prod_i Q(C_i)
        qsize = (h + 1) // 2
        qplus = 1
        for c in comps:
            qplus *= (ham_paths(adj, c) + 1) // 2
        # fibre structure: count colourings per packing by components carrying cycles
        fibre_ok = True
        total = 0
        for p in pk:
            if not p:
                total += 1
                continue
            m = len({next(k for k, C in enumerate(comps) if set(cyc[i]) <= C) for i in p})
            # global classes for this packing: 2^(k-1); local image classes: prod over components 2^(k_i - 1)
            k = len(p)
            ks = Counter(next(kk for kk, C in enumerate(comps) if set(cyc[i]) <= C) for i in p)
            glob = 2 ** (k - 1)
            loc = 1
            for kk in ks.values():
                loc *= 2 ** (kk - 1)
            if glob != loc * 2 ** (m - 1):
                fibre_ok = False
            fib[m] += 1
            total += glob
        res['|Q| = (H+1)/2' if total == qsize else 'Q size FAIL'] += 1
        res['fibres are 2^(m-1) torsors' if fibre_ok else 'fibre FAIL'] += 1
        if len([c for c in comps if ham_paths(adj, c) > 1]) >= 2:
            res['Q not separated (|Q| > |Q^+|)' if qsize > qplus else 'Q separated?'] += 1
    return dict(res), dict(fib)


if __name__ == '__main__':
    print('(1) J2: (n, |L^+(V)|, 2^C(n,2), |image of L|, n!) =', check_J2())
    print('(2) Jss pullback-stability violations (300 random cases):', check_strong_site())
    print('(4) constant presheaf, two plus steps:', check_two_step())
    r, f = check_packings_and_Q()
    print('(3)+(5) packings / carrier:', r)
    print('    packings by number m of components carrying cycles:', dict(sorted(f.items())))
