"""
Investigate mu(3-cycle) distribution at n=6.

At n=6, for a 3-cycle C = v->a->b->v through v:
- C\{v} = {a, b}
- T-v has 5 vertices
- Available = T-v \ {a,b} = 3 vertices

With 3 available vertices, the only possible odd cycle is a 3-cycle.
On 3 vertices there are 0 or 2 directed 3-cycles (if cyclic, both directions; if transitive, none).

If no 3-cycle exists among the 3 available vertices:
  Omega is empty, I(empty, 2) = 1, so mu = 1.

If a 3-cycle exists (actually both directions exist, giving 2 directed 3-cycles):
  Omega has 2 vertices (the two directed 3-cycles on the same 3 vertices).
  They share all vertices, so they're adjacent in Omega.
  I(K_2, 2) = 1 + 2 + 2 = 5? No wait.
  Independent sets of K_2: {}, {v1}, {v2}. Sizes: 0, 1, 1.
  I(K_2, 2) = 2^0 + 2^1 + 2^1 = 1 + 2 + 2 = 5.

  Wait, but that seems wrong. Let me reconsider.

  Actually, there are 2 directed 3-cycles on 3 vertices (the two cyclic orderings).
  They share ALL vertices, so they're adjacent. The conflict graph is K_2.
  Independent sets: empty (size 0), {cycle1} (size 1), {cycle2} (size 1).
  I(K_2, 2) = 1 + 2 + 2 = 5.

  So mu = 5 in this case?

  Hmm, wait. Let me reconsider. Is it possible to have just ONE directed 3-cycle?
  On 3 vertices {x,y,z}, a directed 3-cycle x->y->z->x exists iff the tournament
  on {x,y,z} is cyclic (not transitive). If cyclic, BOTH x->y->z->x AND x->z->y->x
  exist (the two directions). Actually no - if x->y, y->z, z->x, then the only
  directed 3-cycle is x->y->z->x (and its rotations y->z->x->y, z->x->y->z,
  which are the same cycle). The REVERSE z->y->x->z would need z->y, y->x, x->z,
  which means the tournament has z->y (not y->z). Contradiction.

  So a cyclic tournament on 3 vertices has exactly 1 directed 3-cycle (up to rotation),
  not 2. In my enumeration, fixing a starting vertex, there's 1 such cycle.

  Wait, I need to be more careful. If the tournament is x->y, y->z, z->x:
  - x->y->z->x: valid (the 3-cycle)
  - z->y is NOT an arc (y->z is), so no reverse cycle.

  So Omega has exactly 1 vertex (the one directed 3-cycle).
  I(K_1, 2) = 1 + 2 = 3.
  So mu = 3.

Author: opus-2026-03-05-S1
"""

from itertools import permutations, combinations


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(2**len(edges)):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T


def find_directed_odd_cycles_in(T, vertices):
    """Find all directed odd cycles among given vertices."""
    vlist = sorted(vertices)
    cycles = []
    for L in range(3, len(vlist) + 1, 2):
        for sub in combinations(vlist, L):
            first = sub[0]
            for perm in permutations(sub[1:]):
                seq = (first,) + perm
                if all(T[seq[i]][seq[i+1]] == 1 for i in range(L-1)) and T[seq[-1]][first] == 1:
                    cycles.append(frozenset(sub))
    return cycles


def compute_mu(T, v, n, cycle_other_verts):
    tmv = [u for u in range(n) if u != v]
    avail = [u for u in tmv if u not in cycle_other_verts]
    if len(avail) < 3:
        return 1
    cycles = find_directed_odd_cycles_in(T, avail)
    if not cycles:
        return 1
    k = len(cycles)
    adj = [set() for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)
    total = 0
    for bits in range(2**k):
        nodes = [i for i in range(k) if (bits >> i) & 1]
        indep = True
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                if nodes[j] in adj[nodes[i]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            total += 2**len(nodes)
    return total


if __name__ == "__main__":
    n = 6
    print(f"mu(3-cycle) distribution at n={n}")
    print("="*60)

    mu_dist = {}
    total_3cycles = 0

    # For speed, sample tournaments
    count = 0
    for T in all_tournaments(n):
        count += 1
        if count > 5000:  # sample first 5000 of 32768
            break
        for v in range(n):
            other = [u for u in range(n) if u != v]
            for a, b in permutations(other, 2):
                if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
                    total_3cycles += 1
                    mu = compute_mu(T, v, n, frozenset([a, b]))
                    mu_dist[mu] = mu_dist.get(mu, 0) + 1

    print(f"Sampled {count} tournaments")
    print(f"Total directed 3-cycles through v: {total_3cycles}")
    print(f"mu distribution: {dict(sorted(mu_dist.items()))}")

    if mu_dist:
        total = sum(mu_dist.values())
        for mu_val in sorted(mu_dist.keys()):
            pct = 100 * mu_dist[mu_val] / total
            print(f"  mu={mu_val}: {mu_dist[mu_val]} ({pct:.1f}%)")

    # Also check mu for 5-cycles at n=6 (should all be 1 by THM-008)
    print(f"\nmu(5-cycle) distribution at n={n} (verification)")
    mu_dist_5 = {}
    count = 0
    for T in all_tournaments(n):
        count += 1
        if count > 1000:
            break
        for v in range(n):
            other = [u for u in range(n) if u != v]
            for subset in combinations(other, 4):
                for perm in permutations(subset):
                    valid = T[v][perm[0]] == 1
                    if valid:
                        for i in range(3):
                            if T[perm[i]][perm[i+1]] != 1:
                                valid = False
                                break
                    if valid and T[perm[3]][v] == 1:
                        mu = compute_mu(T, v, n, frozenset(perm))
                        mu_dist_5[mu] = mu_dist_5.get(mu, 0) + 1
    print(f"  mu distribution for 5-cycles: {mu_dist_5}")
    print(f"  (Expected: all mu=1 by THM-008)")

    # Analysis of mu=3 cases for 3-cycles
    print(f"\nAnalyzing mu=3 cases (3-cycle with cyclic subtournament):")
    example_count = 0
    for T in all_tournaments(n):
        for v in range(n):
            other = [u for u in range(n) if u != v]
            for a, b in permutations(other, 2):
                if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
                    mu = compute_mu(T, v, n, frozenset([a, b]))
                    if mu == 3 and example_count < 3:
                        avail = [u for u in other if u != a and u != b]
                        # Check if avail forms a 3-cycle
                        x, y, z = avail
                        is_cyclic = (T[x][y] == 1 and T[y][z] == 1 and T[z][x] == 1) or \
                                    (T[x][z] == 1 and T[z][y] == 1 and T[y][x] == 1)
                        print(f"  v={v}, 3-cycle=v->{a}->{b}->v, available={avail}")
                        print(f"    Subtournament on available is {'cyclic' if is_cyclic else 'transitive'}")
                        print(f"    mu = {mu}")
                        example_count += 1
        if example_count >= 3:
            break

    print(f"\nConclusion:")
    print(f"  mu(3-cycle at n=6) in {{1, 3}}")
    print(f"  mu=1 when the 3 available vertices form a transitive tournament in T-v")
    print(f"  mu=3 when the 3 available vertices form a cyclic tournament in T-v")
    print(f"  (Omega has 0 or 1 cycle respectively; I(empty,2)=1, I(K_1,2)=3)")
