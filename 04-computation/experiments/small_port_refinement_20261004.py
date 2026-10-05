"""Exact triangle-chart words, port gauges, and guarded monitor transport.

The matrix and geometry controls are independent of earlier experiment code.
All checks survive -O. No computation asserts universal Collatz payment.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product


def need(ok, message):
    if not ok:
        raise ValueError(message)


def exact(x):
    need(type(x) in (int, F), "exact rational coordinate required")
    return F(x)


def matrix(rows):
    need(type(rows) in (tuple, list) and len(rows) == 3 and
         all(type(row) in (tuple, list) and len(row) == 3 for row in rows),
         "three by three matrix required")
    return tuple(tuple(exact(x) for x in row) for row in rows)


I = matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
PERMS = tuple(permutations(range(3)))


def mul(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def apply(a, x):
    return tuple(sum(a[i][k]*x[k] for k in range(3)) for i in range(3))


def determinant(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def inverse(a):
    a = matrix(a)
    d = determinant(a)
    need(d != 0, "invertible matrix required")
    # Cyclic row/column orders already supply the cofactor signs.
    return tuple(tuple((a[(j+1) % 3][(i+1) % 3]*a[(j+2) % 3][(i+2) % 3]
                        - a[(j+1) % 3][(i+2) % 3]*a[(j+2) % 3][(i+1) % 3])
                       / d
                       for j in range(3)) for i in range(3))


def permutation_matrix(p):
    need(type(p) is tuple and all(type(i) is int for i in p) and p in PERMS,
         "three-port permutation")
    return matrix(tuple(tuple(int(i == p[j]) for j in range(3)) for i in range(3)))


def chart(p):
    need(type(p) is tuple and all(type(i) is int for i in p) and p in PERMS,
         "three-port chart order")
    a, b, _ = p
    return tuple((F(int(i == a)), F(int(i in (a, b)), 2), F(1, 3))
                 for i in range(3))


CHARTS = {p: chart(p) for p in PERMS}
INVERSES = {p: inverse(CHARTS[p]) for p in PERMS}
PORTS = {permutation_matrix(p): p for p in PERMS}


def fan_chart(c):
    need(type(c) is int and c in range(3), "centroid-fan slot")
    a, b = (c+1) % 3, (c+2) % 3
    return tuple((F(int(i == a)), F(int(i == b)), F(1, 3)) for i in range(3))


FANS = tuple(fan_chart(c) for c in range(3))
FAN_INVERSES = tuple(inverse(a) for a in FANS)
BISECTIONS = (matrix(((1, F(1, 2), 0), (0, F(1, 2), 0), (0, 0, 1))),
              matrix(((0, F(1, 2), 0), (1, F(1, 2), 0), (0, 0, 1))))


def fan_word_matrix(word):
    need(type(word) is tuple, "fan word is a tuple")
    result = I
    for c in word:
        result = mul(result, fan_chart(c))
    return result


def decode_fan(rows):
    a = matrix(rows)
    need(all(sum(a[i][j] for i in range(3)) == 1 for j in range(3)),
         "fan barycentric columns")
    d = abs(determinant(a))
    need(d.numerator == 1, "fan determinant is reciprocal power of three")
    denominator, depth = d.denominator, 0
    while denominator > 1 and denominator % 3 == 0:
        denominator //= 3
        depth += 1
    need(denominator == 1, "fan determinant is reciprocal power of three")
    word = []
    for _ in range(depth):
        center = apply(a, (F(1, 3),)*3)
        minimum = min(center)
        need(center.count(minimum) == 1 and minimum > 0, "strict fan centroid minimum")
        c = center.index(minimum)
        word.append(c)
        a = mul(FAN_INVERSES[c], a)
    need(a in PORTS, "fan terminal port frame")
    return tuple(word), PORTS[a]


@lru_cache(None)
def ternary_trees(internal):
    """Leaf=None; an internal node is its ordered triple of child trees."""
    if internal == 0:
        return (None,)
    result = []
    for first in range(internal):
        for second in range(internal-first):
            third = internal-1-first-second
            for children in product(ternary_trees(first), ternary_trees(second),
                                    ternary_trees(third)):
                result.append(children)
    return tuple(result)


def tree_word(tree):
    if tree is None:
        return ""
    return "H"+tree_word(tree[0])+"G"+tree_word(tree[1])+"G"+tree_word(tree[2])


def parse_tree(word):
    need(type(word) is str and all(c in "HG" for c in word), "H/G syntax")
    balance = 0
    for c in word:
        balance += 2 if c == "H" else -1
        need(balance >= 0, "nonnegative prefix token count")
    need(balance == 0, "balanced token word")

    def read(i):
        if i == len(word) or word[i] == "G":
            return None, i
        first, i = read(i+1)
        need(i < len(word) and word[i] == "G", "first ternary separator")
        second, i = read(i+1)
        need(i < len(word) and word[i] == "G", "second ternary separator")
        third, i = read(i+1)
        return (first, second, third), i

    result, end = read(0)
    need(end == len(word), "complete ternary parse")
    return result


def leaf_addresses(tree, prefix=()):
    if tree is None:
        return (prefix,)
    return tuple(address for c in range(3)
                 for address in leaf_addresses(tree[c], prefix+(c,)))


def tree_from_leaves(addresses):
    need(type(addresses) is tuple and len(set(addresses)) == len(addresses),
         "distinct terminal addresses")
    need(all(type(w) is tuple and all(type(c) is int and c in range(3) for c in w)
             for w in addresses), "ternary addresses")

    def read(prefix):
        below = tuple(w for w in addresses if w[:len(prefix)] == prefix)
        need(below, "every internal node has all three children")
        if prefix in below:
            need(len(below) == 1, "leaf set is prefix-free")
            return None
        return tuple(read(prefix+(c,)) for c in range(3))

    return read(())


def inside_edge(point, a, b):
    i = next(i for i in range(3) if a[i] != b[i])
    t = (point[i]-a[i])/(b[i]-a[i])
    return 0 < t < 1 and all(point[j] == a[j]+t*(b[j]-a[j]) for j in range(3))


def fan_tree_controls():
    halves = BISECTIONS
    for c in range(3):
        a, b = (c+1) % 3, (c+2) % 3
        need(mul(FANS[c], halves[0]) == CHARTS[(a, b, c)] and
             mul(FANS[c], halves[1]) == CHARTS[(b, a, c)],
             "six flags refine three centroid-fan cells")
        need(determinant(FANS[c]) == F(1, 3), "oriented fan area")
    words = 0
    for depth in range(6):
        for w in product(range(3), repeat=depth):
            m = fan_word_matrix(w)
            need(decode_fan(m) == (w, (0, 1, 2)), "fan chart word decoder")
            for p in PERMS:
                need(decode_fan(mul(m, permutation_matrix(p))) == (w, p),
                     "fan terminal frame decoder")
            words += 1
    counts, leaves_checked = {}, 0
    for internal in range(6):
        trees = ternary_trees(internal)
        counts[internal] = len(trees)
        for tree in trees:
            token_word = tree_word(tree)
            need(parse_tree(token_word) == tree, "balanced token word has its full ternary tree")
            leaves = leaf_addresses(tree)
            need(tree_from_leaves(leaves) == tree, "leaf addresses recover adaptive tree")
            mats = tuple(fan_word_matrix(w) for w in leaves)
            recovered = tuple(decode_fan(m)[0] for m in mats)
            need(tree_from_leaves(recovered) == tree, "exact leaf matrices recover grammar parse")
            cells = tuple(frozenset(tuple(m[i][j] for i in range(3)) for j in range(3)) for m in mats)
            vertices = set().union(*cells)
            edges = {frozenset(e) for cell in cells for e in combinations(cell, 2)}
            need((len(vertices), len(edges), len(cells)) ==
                 (internal+3, 3*internal+3, 2*internal+1), "adaptive fan f-vector")
            need(sum(abs(determinant(m)) for m in mats) == 1, "all leaf areas sum to parent")
            for edge in edges:
                a, b = tuple(edge)
                need(not any(inside_edge(v, a, b) for v in vertices if v not in edge),
                     "no hanging boundary vertices")
            leaves_checked += len(leaves)
    need(counts == {0: 1, 1: 1, 2: 3, 3: 12, 4: 55, 5: 273}, "ternary tree counts")
    independent_words = 0
    for size in range(13):
        for letters in product("HG", repeat=size):
            balance, legal = 0, True
            for c in letters:
                balance += 2 if c == "H" else -1
                legal &= balance >= 0
            if legal and balance == 0:
                token_word = "".join(letters)
                need(tree_word(parse_tree(token_word)) == token_word, "independent balanced-word census")
                independent_words += 1
    need(independent_words == 72, "all balanced words through length twelve")
    return words, counts, leaves_checked, independent_words


@lru_cache(None)
def mixed_trees(internal):
    if internal == 0:
        return (None,)
    result = []
    for first in range(internal):
        second = internal-1-first
        result.extend(("P", children) for children in
                      product(mixed_trees(first), mixed_trees(second)))
        for second in range(internal-first):
            third = internal-1-first-second
            result.extend(("H", children) for children in
                          product(mixed_trees(first), mixed_trees(second), mixed_trees(third)))
    return tuple(result)


def mixed_leaves(tree, current=I):
    if tree is None:
        return (current,)
    kind, children = tree
    charts = FANS if kind == "H" else BISECTIONS
    return tuple(m for child, chart in zip(children, charts)
                 for m in mixed_leaves(child, mul(current, chart)))


def mixed_node_counts(tree):
    if tree is None:
        return 0, 0
    kind, children = tree
    counts = tuple(mixed_node_counts(child) for child in children)
    return int(kind == "H")+sum(x[0] for x in counts), int(kind == "P")+sum(x[1] for x in counts)


def mixed_controls():
    alphabet = dict(zip(("A0", "A1", "A2", "B0", "B1"), FANS+BISECTIONS))

    def compose_word(word):
        result = I
        for letter in word:
            result = mul(result, alphabet[letter])
        return result

    left = ("A1", "A0", "A1", "B0", "B0")
    right = ("B0", "B0", "A0", "A1", "A1")
    expected = matrix(((F(4, 9), F(7, 12), F(16, 27)),
                       (F(1, 9), F(1, 12), F(4, 27)),
                       (F(4, 9), F(1, 3), F(7, 27))))
    need(left != right and compose_word(left) == compose_word(right) == expected,
         "mixed chart histories have the same full marked matrix")
    need(determinant(expected) == F(1, 108), "collision has three fan and two binary steps")
    seen, collisions, counts = {I: ()}, [], {0: 1}
    frontier = [((), I)]
    total = 1
    for depth in range(1, 6):
        next_frontier = []
        for word, m in frontier:
            for letter, chart in alphabet.items():
                new_word, new_matrix = word+(letter,), mul(m, chart)
                if new_matrix in seen:
                    collisions.append((seen[new_matrix], new_word))
                else:
                    seen[new_matrix] = new_word
                next_frontier.append((new_word, new_matrix))
        counts[depth] = len({m for _, m in next_frontier})
        total += len(next_frontier)
        frontier = next_frontier
    need(total == 3906 and len(seen) == 3902 and len(collisions) == 4 and
         all(len(a) == len(b) == 5 for a, b in collisions),
         "complete first-collision universe")
    multiplicities = Counter(m for _, m in frontier)
    need(sorted(v for v in multiplicities.values() if v > 1) == [2]*4,
         "exactly four two-word collision fibers at depth five")
    need(counts == {0: 1, 1: 5, 2: 25, 3: 125, 4: 625, 5: 3121},
         "mixed matrix counts by depth")
    checked, tree_counts = 0, {}
    for internal in range(4):
        trees = mixed_trees(internal)
        tree_counts[internal] = len(trees)
        for tree in trees:
            h, p = mixed_node_counts(tree)
            leaves = mixed_leaves(tree)
            need(h+p == internal and len(leaves) == 2*h+p+1, "mixed tree leaf count")
            need(sum(abs(determinant(m)) for m in leaves) == 1, "mixed leaf areas cover parent")
            checked += 1
    need(tree_counts == {0: 1, 1: 2, 2: 10, 3: 66}, "mixed typed tree census")
    # Fixed-port operations already produce a hanging point after H,H,P.
    hanging_tree = ("H", (None, None, ("H", (("P", (None, None)), None, None))))
    leaves = mixed_leaves(hanging_tree)
    point = (F(1, 6), F(2, 3), F(1, 6))
    vertices = {tuple(m[i][j] for i in range(3)) for m in leaves for j in range(3)}
    untouched = frozenset(tuple(FANS[0][i][j] for i in range(3)) for j in range(3))
    need(point in vertices and point not in untouched and
         inside_edge(point, (0, 1, 0), (F(1, 3),)*3),
         "a midpoint hangs on the boundary of an untouched adjacent leaf")
    need(FANS[0] in leaves, "unrefined neighbor exists")
    return total, counts, len(collisions), checked, tree_counts, expected


def word_matrix(word):
    need(type(word) is tuple, "chart word is a tuple")
    result = I
    for p in word:
        result = mul(result, chart(p))
    return result


def decode_matrix(rows, allow_port_frame=False):
    """Return exact chart word and terminal permutation; never guess a frame."""
    a = matrix(rows)
    need(all(sum(a[i][j] for i in range(3)) == 1 for j in range(3)),
         "affine barycentric columns must sum to one")
    d = abs(determinant(a))
    need(d.numerator == 1, "word determinant must be a reciprocal power of six")
    denominator, depth = d.denominator, 0
    while denominator > 1 and denominator % 6 == 0:
        denominator //= 6
        depth += 1
    need(denominator == 1, "word determinant must be a reciprocal power of six")
    word = []
    for _ in range(depth):
        centroid = apply(a, (F(1, 3),)*3)
        need(len(set(centroid)) == 3, "child centroid has strict coordinate order")
        p = tuple(sorted(range(3), key=lambda i: -centroid[i]))
        need(centroid[p[-1]] > 0, "child centroid lies in triangle interior")
        word.append(p)
        a = mul(INVERSES[p], a)
    need(a in PORTS, "decoded remainder is not a port frame")
    need(allow_port_frame or a == I, "unadvertised terminal port permutation")
    return tuple(word), PORTS[a]


def normalize(x):
    x = tuple(exact(v) for v in x)
    need(len(x) == 3 and min(x) >= 0 and sum(x) > 0, "nonzero nonnegative triple")
    return tuple(v/sum(x) for v in x)


def point_chart(lam):
    lam = normalize(lam)
    p = tuple(sorted(range(3), key=lambda i: (-lam[i], i)))
    local = apply(INVERSES[p], lam)
    need(min(local) >= 0 and sum(local) == 1, "point is in its selected chart")
    return p, local


def stratum(x):
    a, b, c = sorted(x, reverse=True)
    need(c >= 0 and a > 0, "nonzero triangle address")
    if b == 0:
        return "vertex"
    if c == 0:
        return "equal-edge" if a == b else "unequal-edge"
    if a == c:
        return "center"
    if a == b:
        return "two-large"
    if b == c:
        return "two-small"
    return "distinct-interior"


def stats(word):
    p, q, b = 1, 1, 0
    for a in word:
        need(type(a) is int and a >= 1, "positive valuation letter")
        p, q, b = 3*p, (2**a)*q, 3*b+q
    return p, q, b


def monitor_matrix(word):
    p, q, b = stats(word)
    return matrix(((p, 0, b), (0, q, 0), (0, 0, q)))


def replay(x, word):
    need(type(x) is int and x > 0 and x % 2 == 1, "positive odd current value")
    states = [x]
    for a in word:
        raw = 3*x+1
        k = (raw & -raw).bit_length()-1
        if k != a:
            return None
        x = raw >> k
        states.append(x)
    return tuple(states)


def guard(x, word):
    p, q, b = stats(word)
    return type(x) is int and x > 0 and x % 2 == 1 and (p*x+b-q) % (2*q) == 0


def graph_controls():
    vertices = tuple(tuple(F(int(i in face), len(face)) for i in range(3))
                     for size in (1, 2, 3) for face in combinations(range(3), size))
    rank = [sum(v > 0 for v in x) for x in vertices]
    cells = {frozenset(tuple(CHARTS[p][i][j] for i in range(3)) for j in range(3))
             for p in PERMS}
    edges = {frozenset(e) for cell in cells for e in combinations(cell, 2)}
    need((len(vertices), len(edges), len(cells)) == (7, 12, 6), "barycentric mesh")
    untyped, typed = [], []
    for p in permutations(range(7)):
        image = dict(zip(vertices, (vertices[i] for i in p)))
        if {frozenset(image[v] for v in e) for e in edges} == edges:
            untyped.append(p)
            if all(rank[i] == rank[p[i]] for i in range(7)):
                typed.append(p)
    need(len(untyped) == 12 and len(typed) == 6, "bare wheel loses parent ranks")
    return len(vertices), len(edges), len(cells), len(untyped), len(typed)


def main():
    for p in PERMS:
        need(mul(CHARTS[p], INVERSES[p]) == I and mul(INVERSES[p], CHARTS[p]) == I,
             "independent exact inverse")
        need(abs(determinant(CHARTS[p])) == F(1, 6), "six equal-area chart cells")
    print("Barycentric triangle:", graph_controls(), "= vertices,edges,cells,untypedAut,typedAut")
    words, frame_checks, by_depth = [], 0, {}
    for depth in range(5):
        matrices = set()
        for w in product(PERMS, repeat=depth):
            m = word_matrix(w)
            need(decode_matrix(m) == (w, (0, 1, 2)), "lossless chart word")
            matrices.add(m)
            words.append(w)
            if depth <= 3:
                for frame in PERMS:
                    need(decode_matrix(mul(m, permutation_matrix(frame)), True) == (w, frame),
                         "explicit output port frame")
                    frame_checks += 1
        need(len(matrices) == 6**depth, "no same-depth matrix collision")
        by_depth[depth] = len(matrices)
    print("Exact chart-word counts by depth:", by_depth, "; total", len(words))
    print("Independent terminal port frames:", frame_checks)
    short = [w for w in words if len(w) <= 2]
    for u, v in product(short, repeat=2):
        need(mul(word_matrix(u), word_matrix(v)) == word_matrix(u+v),
             "concatenation composes affine charts")
    port_checks = 0
    for p, q, frame in product(PERMS, repeat=3):
        gauge = permutation_matrix(frame)
        need(mul(mul(CHARTS[p], gauge), mul(inverse(gauge), CHARTS[q]))
             == mul(CHARTS[p], CHARTS[q]), "internal frame cancels across interface")
        port_checks += 1
    flip = permutation_matrix((1, 0, 2))
    b = CHARTS[(0, 1, 2)]
    need(mul(mul(b, flip), b) != mul(b, b), "forgotten frame changes refinement")
    center = (F(1, 3),)*3
    need(all(apply(CHARTS[p], (0, 0, 1)) == center for p in PERMS),
         "boundary point alone loses six incident chart choices")
    print("Composition pairs", len(short)**2, "; compensated interface gauges", port_checks)
    print("Hostiles: uncompensated port swap changes matrix; all six chart tips give one center.")
    fan_words, tree_counts, leaf_count, balanced_words = fan_tree_controls()
    print("Centroid fan: three cells; each splits into two flags; all fan words depth0..5:", fan_words)
    print("Full adaptive ternary-tree counts:", tree_counts, "; exact leaves checked", leaf_count)
    print("Balanced H(+2)/G(-1) words through length12:", balanced_words)
    print("For m fan splits: (V,E,F)=(m+3,3m+3,2m+1); exact area and no hanging nodes.")
    mixed_total, mixed_counts, collisions, mixed_checked, tree_counts, witness = mixed_controls()
    print("Mixed alphabet A0,A1,A2,B0,B1 words length0..5:", mixed_total,
          "; distinct matrices by depth", mixed_counts)
    print("First marked-matrix collisions at depth5:", collisions,
          "; A1 A0 A1 B0 B0 = B0 B0 A0 A1 A1")
    print("Collision matrix:", witness)
    print("Mixed H/P typed trees through3 nodes:", tree_counts, "; area/leaf checks", mixed_checked)
    print("Mixed leaf count2h+p+1; midpoint(1/6,2/3,1/6) is a hanging-node hostile.")

    point_checks = 0
    strata = {}
    for degree in range(1, 13):
        census = Counter()
        for x in range(degree+1):
            for y in range(degree-x+1):
                z = degree-x-y
                lam = tuple(F(v, degree) for v in (x, y, z))
                p, local = point_chart(lam)
                need(apply(CHARTS[p], local) == lam, "all weak-order boundaries reconstruct")
                census[stratum((x, y, z))] += 1
                point_checks += 1
        strata[degree] = census
    need(set(strata[6]) == set(strata[12])-{"two-large"}, "degree six omits a stratum")
    need(all(len(strata[d]) < 7 for d in range(1, 12)) and len(strata[12]) == 7,
         "first degree displaying every continuous orbit type")
    print("Point/boundary reconstructions degrees1..12:", point_checks)
    print("Degree6 strata:", dict(sorted(strata[6].items())))
    print("Degree12 first realizes all seven strata:", dict(sorted(strata[12].items())))

    transport_checks = 0
    valuation_words = [w for length in range(1, 4) for w in product(range(1, 5), repeat=length)]
    for w in valuation_words:
        p, q, carry = stats(w)
        residue = (q-carry)*pow(p, -1, 2*q) % (2*q)
        mw = monitor_matrix(w)
        for n, lift in product((3, 7, 27, 703), (0, 1)):
            x = residue+2*q*lift
            need(guard(x, w), "exact source cylinder")
            route = replay(x, w)
            need(route is not None, "independent literal guard replay")
            y = route[-1]
            before = normalize((x, n, 1))
            after = normalize((y, n, 1))
            sigma, local = point_chart(before)
            tau, expected = point_chart(after)
            transition = mul(mul(INVERSES[tau], mw), CHARTS[sigma])
            observed = normalize(apply(transition, local))
            need(observed == expected and normalize(apply(mw, before)) == after,
                 "chart-conjugated arithmetic edge")
            decoded = apply(CHARTS[tau], observed)
            need(decoded[0]/decoded[2] == y and decoded[1]/decoded[2] == n,
                 "both current and original source survive")
            need((decoded[0] < decoded[1]) == (y < n), "source payment wall")
            transport_checks += 1
    # Compose chart interfaces at actual intermediate values, including walls.
    telescope_checks = 0
    for n, x in product((3, 7, 27), range(1, 40, 2)):
        def next_step(v):
            raw = 3*v+1
            a = (raw & -raw).bit_length()-1
            return a, raw >> a
        a, y = next_step(x)
        c, z = next_step(y)
        sigma, _ = point_chart((x, n, 1))
        tau, _ = point_chart((y, n, 1))
        rho, _ = point_chart((z, n, 1))
        first = mul(mul(INVERSES[tau], monitor_matrix((a,))), CHARTS[sigma])
        second = mul(mul(INVERSES[rho], monitor_matrix((c,))), CHARTS[tau])
        together = mul(mul(INVERSES[rho], monitor_matrix((a, c))), CHARTS[sigma])
        need(mul(second, first) == together, "transition gauges telescope exactly")
        need(guard(x, (a,)) and guard(y, (c,)) and guard(x, (a, c)),
             "sequential guard intersection")
        telescope_checks += 1
    need(not guard(27, (1, 2, 1, 2)) and guard(27, (1, 2)),
         "word repetition does not pay the stronger guard")
    need(not 5 < 3 and 5 < 7, "current value alone cannot report payment")
    need(point_chart((5, 3, 1))[0] == (0, 1, 2) and
         point_chart((5, 7, 1))[0] == (1, 0, 2), "two physically occupied open chambers")
    print("Guarded chart transport:", transport_checks, "cases from", len(valuation_words), "words")
    print("Exact telescoping/guard intersections:", telescope_checks)
    print("Payment hostile: current5 with original3 is unpaid, with original7 is paid.")
    print("Only two open monitor chambers occur when current and original are >1 and unequal.")
    print("Formal repeat12 at27 is legal once, not twice; chart transport does not repair it.")
    rejected = 0
    for thunk in (lambda: chart((0, 0, 1)), lambda: chart((False, 1, 2)),
                  lambda: matrix(((True, 0, 0), (0, 1, 0), (0, 0, 1))),
                  lambda: decode_matrix(((1, 0, 0), (0, F(1, 2), 0), (0, F(1, 2), 1))),
                  lambda: decode_matrix(mul(b, flip)), lambda: normalize((1, -1, 1)),
                  lambda: parse_tree("G"), lambda: parse_tree("H"),
                  lambda: tree_from_leaves(((0,), (1,))), lambda: decode_fan(b)):
        try:
            thunk()
        except ValueError:
            rejected += 1
        else:
            raise ValueError("malformed chart control accepted")
    print("Rejected malformed/undeclared-frame controls:", rejected)
    print("PASS: pure-chart decoders, typed mixed partitions, guarded interfaces; no universal payment claim.")


if __name__ == "__main__":
    main()
