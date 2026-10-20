"""Independent finite controls for web.md; Python standard library only."""

from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, permutations


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


V = tuple(range(4))
E = tuple(combinations(V, 2))
EI = {edge: i for i, edge in enumerate(E)}
FULL = (1 << len(E)) - 1
PERMS = tuple(permutations(V))


def edge_image(g, edge):
    return tuple(sorted(g[v] for v in edge))


def act(g, mask):
    return sum(1 << EI[edge_image(g, e)] for i, e in enumerate(E) if mask >> i & 1)


def mask_of(edges):
    return sum(1 << EI[tuple(sorted(e))] for e in edges)


def opposite(mask):
    return mask_of(set(V) - set(e) for i, e in enumerate(E) if mask >> i & 1)


def degree(mask):
    return tuple(sum(bool(mask >> i & 1) for i, e in enumerate(E) if v in e) for v in V)


def face_parity(mask):
    return tuple(
        sum(bool(mask >> i & 1) for i, e in enumerate(E) if v not in e) % 2 for v in V
    )


def cut(vertices):
    vertices = set(vertices)
    return mask_of(e for e in E if (e[0] in vertices) != (e[1] in vertices))


def sign(g):
    return (-1) ** sum(g[i] > g[j] for i, j in combinations(range(len(g)), 2))


def edge_sign(g, mask):
    return (-1) ** sum(
        g[i] > g[j] for k, (i, j) in enumerate(E) if mask >> k & 1
    )


def connected(mask):
    seen = {0}
    while True:
        new = seen | {w for i, e in enumerate(E) if mask >> i & 1 and seen.intersection(e)
                      for w in e}
        if new == seen:
            return len(seen) == 4
        seen = new


def hp_count(mask):
    def arc(u, v):
        bit = mask >> EI[tuple(sorted((u, v)))] & 1
        return (u < v) != bool(bit)
    return sum(all(arc(u, v) for u, v in zip(p, p[1:])) for p in PERMS)


B = tuple(tuple((v == j) - (v == i) for v in V) for i, j in E)
DOT = tuple(tuple(sum(x * y for x, y in zip(e, f)) for f in B) for e in B)
ADJ = tuple((i, j) for i, j in combinations(range(6), 2) if DOT[i][j])
CLIQUES = tuple(c for c in combinations(range(6), 3)
                if all(DOT[i][j] for i, j in combinations(c, 2)))
STARS = {tuple(i for i, e in enumerate(E) if v in e) for v in V}
FACES = set(CLIQUES) - STARS


def star_vector(i):
    """Hodge star of oriented e_i∧e_j, with orientation (0,1,2,3)."""
    a, b = E[i]
    c, d = sorted(set(V) - {a, b})
    result = [0] * 6
    result[EI[c, d]] = sign((a, b, c, d))
    return tuple(result)


STAR = tuple(star_vector(i) for i in range(6))


def hodge(vector):
    return tuple(sum(vector[i] * STAR[i][j] for i in range(6)) for j in range(6))


def oriented_act(g, vector):
    result = [0] * 6
    for i, (a, b) in enumerate(E):
        result[EI[edge_image(g, (a, b))]] += (1 if g[a] < g[b] else -1) * vector[i]
    return tuple(result)


def twin_quotient_defect(mask):
    """Independent triple blowup, compressed by actual false-twin classes."""
    neighborhoods = [tuple(w for w in V if w != v and mask >> EI[tuple(sorted((v, w)))] & 1)
                     for v in V]
    classes = {}
    for v, neighborhood in enumerate(neighborhoods):
        classes.setdefault(neighborhood, []).append(v)
    blocks = list(classes.values())
    weights = [3 * len(block) for block in blocks]
    qedges = {(i, j) for i, j in combinations(range(len(blocks)), 2)
              if mask >> EI[tuple(sorted((blocks[i][0], blocks[j][0])))] & 1}
    autos = [g for g in permutations(range(len(blocks)))
             if all(weights[g[i]] == weights[i] for i in range(len(blocks)))
             and {tuple(sorted((g[i], g[j]))) for i, j in qedges} == qedges]
    total = 0
    for g in autos:
        value = (-1) ** sum(weights[i] * weights[j] for i, j in qedges if g[i] > g[j])
        seen = set()
        for i in range(len(blocks)):
            if i in seen:
                continue
            orbit, j = [], i
            while j not in seen:
                seen.add(j)
                orbit.append(j)
                j = g[j]
            k = weights[i] // 3
            value *= k + 1 if len(orbit) % 2 == 0 else (1 if k % 2 == 0 else 0)
        total += value
    return -Fraction(total, len(autos))


def larger_hessian_control(n):
    vertices = tuple(range(n))
    edges = tuple(combinations(vertices, 2))
    index = {edge: i for i, edge in enumerate(edges)}
    patterns = []
    for p in permutations(vertices):
        mask = required = 0
        for a, b in zip(p, p[1:]):
            bit = 1 << index[tuple(sorted((a, b)))]
            mask |= bit
            if a > b:
                required |= bit
        patterns.append((mask, required))
    universe = 1 << len(edges)
    counts = [sum(t & m == r for m, r in patterns) for t in range(universe)]
    pairs = tuple(combinations(range(len(edges)), 2))
    incidence = [tuple((v == b) - (v == a) for v in vertices) for a, b in edges]
    pair_data = []
    for i, j in pairs:
        q = -sum(a * b for a, b in zip(incidence[i], incidence[j]))
        selected = (1 << i) | (1 << j)
        block_patterns = [(m ^ selected, r & ~selected) for m, r in patterns
                          if m & selected == selected] if n == 5 and q else None
        pair_data.append((i, j, q, block_patterns))
    fibers = {}
    checked = 0
    for t in range(universe):
        signature = []
        for i, j, q, block_patterns in pair_data:
            mixed = counts[t] - counts[t ^ (1 << i)] - counts[t ^ (1 << j)] + counts[t ^ (1 << i) ^ (1 << j)]
            signature.append(mixed)
            if q:
                expected = q * (1 - 2 * (t >> i & 1)) * (1 - 2 * (t >> j & 1))
                need(mixed * expected > 0, "All-order adjacent Hessian sign")
                if block_patterns is not None:
                    block_count = sum(t & m == r for m, r in block_patterns)
                    need(mixed == expected * block_count, "Literal formal-block coefficient")
                checked += 1
        fibers.setdefault(tuple(signature), []).append(t)
    need(all(len(fiber) == 2 and fiber[0] ^ fiber[1] == universe - 1
             for fiber in fibers.values()), "Larger Hessian converse fibers")
    print(f"n={n}: all {universe} tournaments; {checked} adjacent sign checks; "
          f"{len(fibers)} Hessian fibers, each exactly converse.")
    if n == 5:
        print("n=5 formal block counts independently match every adjacent mixed difference.")


def main():
    path = mask_of(((0, 1), (1, 2), (2, 3)))
    claw = mask_of(((0, 1), (0, 2), (0, 3)))
    triangle = mask_of(((0, 1), (0, 2), (1, 2)))
    types = {tuple(sorted(degree(path))): "P4",
             tuple(sorted(degree(claw))): "claw",
             tuple(sorted(degree(triangle))): "triangle+isolate"}
    census = Counter(types[tuple(sorted(degree(m)))] for m in range(64) if m.bit_count() == 3)
    need(census == {"P4": 12, "claw": 4, "triangle+isolate": 4}, "Three-edge census")
    need(FULL ^ claw == opposite(claw), "Claw complement is its opposite face")
    need(types[tuple(sorted(degree(FULL ^ path)))] == "P4", "Path self-complementarity")
    for original in (claw, triangle):
        selected = [e for i, e in enumerate(E) if original >> i & 1]
        need(all(set(a).intersection(b) for a, b in combinations(selected, 2)),
             "Both exceptional line graphs must be K3")
    print(f"Three-edge labelled census: {dict(census)}; both claw and triangle have line graph K3.")

    cuts = {cut(vertices) for k in range(5) for vertices in combinations(V, k)}
    cycles = {m for m in range(64) if all(d % 2 == 0 for d in degree(m))}
    need(len(cuts) == len(cycles) == 8 and len(cuts & cycles) == 4, "Binary bicycle dimension")
    need({opposite(m) for m in cuts} == cycles, "Opposition must exchange cuts and cycles")
    sign_checks = 0
    for mask in range(64):
        need(face_parity(mask) == tuple((d + mask.bit_count()) % 2 for d in degree(mask)),
             "Face-degree formula")
        need(face_parity(opposite(mask)) == tuple(d % 2 for d in degree(mask)),
             "Opposition face-degree formula")
        if mask.bit_count() % 2:
            residue = FULL ^ mask ^ opposite(mask)
            need(residue in cuts & cycles, "Odd support complement/opposition switching law")
        for g in PERMS:
            need(opposite(act(g, mask)) == act(g, opposite(mask)), "Opposition centrality")
            if act(g, mask) == mask:
                need(edge_sign(g, mask) * edge_sign(g, FULL ^ mask) == sign(g),
                     "Complement character law")
                sign_checks += 1
    need(sign_checks == 264, "Unexpected graph-automorphism pair count")
    print("Binary cut/cycle spaces: sizes 8,8; intersection size 4. Odd-support law: 32/32.")
    print(f"Complement character identity: all {sign_checks} graph/automorphism pairs.")
    for name, mask in (("P4", path), ("claw", claw), ("triangle+isolate", triangle)):
        autos = [g for g in PERMS if act(g, mask) == mask]
        data = sorted(Counter((sign(g), edge_sign(g, mask)) for g in autos).items())
        print(f"{name}: (vertex sign, edge-reversal sign) multiplicities={data}.")

    trees = {m for m in range(64) if m.bit_count() == 3 and connected(m)}
    need(len(trees) == 16, "K4 tree count")
    need(sum((FULL ^ m) in trees for m in trees) == 12, "Complementary tree count")
    need(all(opposite(FULL ^ m) in trees for m in trees), "Matroid self-duality")
    print("K4 bases: 16 trees; 12 have tree complements; opposition(complement(tree)) always a tree.")

    for i in range(6):
        unit = tuple(int(j == i) for j in range(6))
        need(hodge(hodge(unit)) == unit, "Hodge square")
        for g in PERMS:
            need(hodge(oriented_act(g, unit)) ==
                 tuple(sign(g) * x for x in oriented_act(g, hodge(unit))),
                 "Hodge sign-twisted equivariance")
    for v in V:
        gradient = tuple(edge[v] for edge in B)
        cycle = hodge(gradient)
        need(all(sum(B[i][w] * cycle[i] for i in range(6)) == 0 for w in V),
             "Hodge gradient-to-cycle map")
    print("Signed Hodge star: square=1, cut-to-cycle, and sign-twisted equivariance for all S4.")

    hps = [hp_count(mask) for mask in range(64)]
    signatures = {}
    for mask in range(64):
        # Count cycles independently, without relying on a partial score count.
        cyclic = sum(all(
            sum(((u < v) != bool(mask >> EI[tuple(sorted((u, v)))] & 1))
                for v in triple if v != u) == 1 for u in triple)
            for triple in combinations(V, 3))
        need(hps[mask] == 1 + 2 * cyclic, "Four-vertex Hamiltonian-path identity")
        M = [[0] * 6 for _ in range(6)]
        s = [1 - 2 * (mask >> i & 1) for i in range(6)]
        for i, j in combinations(range(6), 2):
            mixed = hps[mask] - hps[mask ^ (1 << i)] - hps[mask ^ (1 << j)] + hps[mask ^ (1 << i) ^ (1 << j)]
            need(mixed == -2 * DOT[i][j] * s[i] * s[j], "Signed Hessian formula")
            M[i][j] = M[j][i] = mixed
        for clique in CLIQUES:
            product_value = 1
            for i, j in combinations(clique, 2):
                product_value *= M[i][j] // 2
            need(product_value == (-1 if clique in STARS else 1), "Incidence from clique product")
        recovered = [None] * 6
        recovered[0] = 1
        queue = deque([0])
        while queue:
            i = queue.popleft()
            for j in range(6):
                if i == j or DOT[i][j] == 0:
                    continue
                value = recovered[i] * M[i][j] // (-2 * DOT[i][j])
                if recovered[j] is None:
                    recovered[j] = value
                    queue.append(j)
                else:
                    need(recovered[j] == value, "Sign-product synchronization")
        need(recovered == [s[0] * value for value in s], "Recovery up to global sign")
        signature = tuple(M[i][j] for i, j in combinations(range(6), 2))
        signatures.setdefault(signature, []).append(mask)
    need(len(signatures) == 32 and all(len(fiber) == 2 and fiber[0] ^ fiber[1] == FULL
                                      for fiber in signatures.values()), "Hessian converse fibers")
    support_autos = []
    signed_incidence_autos = []
    for p in permutations(range(6)):
        if all(bool(DOT[i][j]) == bool(DOT[p[i]][p[j]]) for i, j in combinations(range(6), 2)):
            support_autos.append(p)
            if {tuple(sorted(p[i] for i in star)) for star in STARS} == STARS:
                signed_incidence_autos.append(p)
    induced = {tuple(EI[edge_image(g, e)] for e in E) for g in PERMS}
    need(len(support_autos) == 48 and set(signed_incidence_autos) == induced,
         "Octahedron incidence-duality ambiguity")
    print("Hamiltonian-path mixed Hessian: 64/64 verified; 32 fibers, each exactly a converse pair.")
    print("Eight Hessian triangles: four negative stars, four positive faces for every tournament.")
    print("Unsigned octahedron automorphisms: 48; star-preserving automorphisms: exactly induced S4 (24).")
    larger_hessian_control(3)
    larger_hessian_control(5)
    larger_hessian_control(6)

    need(twin_quotient_defect(path) == 2 and twin_quotient_defect(FULL ^ path) == 2,
         "Path defect must survive quotient complementation")
    need(twin_quotient_defect(claw) == twin_quotient_defect(triangle) == 0,
         "Other trio blowup controls")
    expanded = [set() for _ in range(12)]
    for i, j in E:
        if path >> EI[i, j] & 1:
            for a in range(3 * i, 3 * i + 3):
                for b in range(3 * j, 3 * j + 3):
                    expanded[a].add(b)
                    expanded[b].add(a)
    complementary_neighborhoods = [frozenset(set(range(12)) - {v} - expanded[v]) for v in range(12)]
    need(max(Counter(complementary_neighborhoods).values()) == 1,
         "True graph complement has no independent twin triple")
    print("Independent triple blowups: delta(P4)=2, delta(claw)=delta(triangle+isolate)=0.")
    print("Quotient-complement P4 blowup still has delta=2; true graph complement has no twin triple.")
    print("PASS")


if __name__ == "__main__":
    main()
