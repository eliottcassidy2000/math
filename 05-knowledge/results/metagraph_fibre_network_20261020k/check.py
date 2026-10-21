"""Exact finite controls; no orbit/metagraph connectivity oracle is imported."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys

CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(message)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError("exact integer outside domain")
    return x


def reduced_forms(m):
    integer(m, 1)
    D = -12 * m * m
    forms = []
    for a in range(1, isqrt((-D) // 3) + 1):
        for b in range(-a, a + 1):
            numerator = b * b - D
            if numerator % (4 * a):
                continue
            c = numerator // (4 * a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    return tuple(forms)


def factor(n):
    integer(n, 1)
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return tuple(out)


def class_number(m):
    integer(m, 1)
    h = 1
    for p, e in factor(m):
        if p == 2 or p == 3:
            h *= p ** e
        else:
            chi = 1 if p % 3 == 1 else -1
            h *= p ** (e - 1) * (p - chi)
    return h


def ambiguous(form):
    a, b, c = form
    return b == 0 or abs(b) == a or a == c


def power_two_hostile(a):
    integer(a, 3)
    m = 1 << a
    b = next(b for b in (-6, -4, -2, 2, 4, 6) if (b - 3 * m) % 7 == 0)
    numerator = b * b + 12 * m * m
    if numerator % 28:
        raise ArithmeticError("integral reduced form")
    return (7, b, numerator // 28)


MARKED_CLASSES = {
    0: ((1, 0, 3),),
    1: ((1, 0, 12), (3, 0, 4)),
    2: ((1, 0, 48), (3, 0, 16), (4, 4, 13), (7, 2, 7)),
}


def marked_class(bits):
    if type(bits) is not tuple or len(bits) > 2 or any(type(b) is not int or b not in (0, 1) for b in bits):
        raise ValueError("zero, one, or two marked binary coordinates required")
    code = sum(b << i for i, b in enumerate(bits))
    return MARKED_CLASSES[len(bits)][code]


def tournament(k):
    integer(k, 2)
    n = 2 * k + 1
    return tuple(tuple((i < j) != ((i, j) in ((0, k), (k, 2 * k)) or
                                     (j, i) in ((0, k), (k, 2 * k))) if i != j else False
                       for j in range(n)) for i in range(n))


def choice(k, x):
    integer(k, 2)
    integer(x)
    top = 1 << (k - 1)
    if x > top:
        raise ValueError("choice outside source")
    return set(range(k)) if x == top else {i + 1 for i in range(k - 1) if x >> i & 1}


def word(k, x, y):
    I, K = choice(k, x), choice(k, y)
    A, B = set(range(k)), set(range(k + 1, 2 * k + 1))
    J = {2 * k - v for v in A - K}
    return tuple(sorted(I) + sorted(J) + [k] + sorted(A - I) + sorted(B - J))


def validate_word(k, w):
    rows = tournament(k)
    if type(w) is not tuple or any(type(v) is not int for v in w) or sorted(w) != list(range(2 * k + 1)):
        raise ValueError("literal source word required")
    if not all(rows[u][v] for u, v in zip(w, w[1:])):
        raise ValueError("not a Hamiltonian path of the source")
    return rows


def blue(k, w):
    rows = validate_word(k, w)
    rev = dict(zip(w, reversed(w)))
    return all(rows[u][v] != rows[rev[u]][rev[v]] for u, v in combinations(range(len(rows)), 2))


def paths(k):
    rows = tournament(k)
    answer = set()
    def visit(w, left):
        if not left:
            answer.add(w)
        for v in left:
            if not w or rows[w[-1]][v]:
                visit(w + (v,), left - {v})
    visit((), set(range(2 * k + 1)))
    return answer


def forget(k, x):
    choice(k + 1, x)
    return (1 << (k - 1)) if x == (1 << k) else x & ((1 << (k - 1)) - 1)


def delete_middle(k, w):
    integer(k, 2)
    validate_word(k + 1, w)
    survivor = tuple(v if v < k else v - 1 if v == k + 1 else v - 2
                     for v in w if v not in (k, k + 2))
    validate_word(k, survivor)
    return survivor


@dataclass(frozen=True)
class Network:
    n: int
    edges: tuple

    def __post_init__(self):
        integer(self.n, 1)
        if type(self.edges) is not tuple:
            raise ValueError("exact tuple edge list")
        last = (-1, -1)
        for edge in self.edges:
            if type(edge) is not tuple or len(edge) != 3:
                raise ValueError("edge triple")
            u, v, c = edge
            integer(u)
            integer(v)
            if not u < v < self.n or (u, v) <= last or type(c) not in (int, F) or c <= 0:
                raise ValueError("simple ordered positive conductance edges")
            last = (u, v)


def complete(weights):
    if type(weights) is not tuple or not weights:
        raise ValueError("nonempty weight tuple")
    for x in weights:
        integer(x, 1)
    return Network(len(weights), tuple((i, j, weights[i] * weights[j])
                                      for i, j in combinations(range(len(weights)), 2)))


def wire(net, projection):
    if not isinstance(net, Network) or type(projection) is not tuple or len(projection) != net.n:
        raise ValueError("source-bound projection")
    for x in projection:
        integer(x)
    q = max(projection) + 1
    if set(projection) != set(range(q)):
        raise ValueError("surjective contiguous fibre labels")
    counts, collapsed = Counter(), []
    for u, v, c in net.edges:
        a, b = sorted((projection[u], projection[v]))
        if a == b:
            collapsed.append((u, v, c))
        else:
            counts[a, b] += c
    return Network(q, tuple((a, b, c) for (a, b), c in sorted(counts.items()))), tuple(collapsed)


def energy(net, potentials):
    if type(potentials) is not tuple or len(potentials) != net.n or any(type(x) not in (int, F) for x in potentials):
        raise ValueError("exact source potential tuple")
    return sum(c * (potentials[u] - potentials[v]) ** 2 for u, v, c in net.edges)


def components(net):
    near = [set() for _ in range(net.n)]
    for u, v, _ in net.edges:
        near[u].add(v)
        near[v].add(u)
    out, unseen = [], set(range(net.n))
    while unseen:
        todo, found = [min(unseen)], set()
        while todo:
            v = todo.pop()
            if v not in found:
                found.add(v)
                todo.extend(near[v] - found)
        out.append(tuple(sorted(found)))
        unseen -= found
    return tuple(out)


def resistance(net, a, b):
    integer(a)
    integer(b)
    if max(a, b) >= net.n:
        raise ValueError("terminal outside source")
    if a == b:
        return F(0)
    comp = next(c for c in components(net) if a in c)
    if b not in comp:
        return None  # infinite resistance; never replace by a finite wired value
    labels = [v for v in comp if v != b]
    ix = {v: i for i, v in enumerate(labels)}
    n = len(labels)
    M = [[F(0) for _ in range(n + 1)] for _ in range(n)]
    M[ix[a]][n] = F(1)
    for u, v, c in net.edges:
        if u not in comp:
            continue
        for x, y in ((u, v), (v, u)):
            if x != b:
                M[ix[x]][ix[x]] += c
                if y != b:
                    M[ix[x]][ix[y]] -= c
    for j in range(n):
        p = next(i for i in range(j, n) if M[i][j])
        M[j], M[p] = M[p], M[j]
        scale = M[j][j]
        M[j] = [x / scale for x in M[j]]
        for i in range(n):
            if i != j:
                c = M[i][j]
                M[i] = [x - c * y for x, y in zip(M[i], M[j])]
    return M[ix[a]][n]


def reject(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        check(True, "typed hostile rejected")
    else:
        check(False, "typed hostile accepted")


def verify():
    global CHECKS
    CHECKS = 0
    small = {}
    for m in range(1, 65):
        forms = reduced_forms(m)
        check(len(forms) == class_number(m), "independent forms versus ring class-number formula")
        check(all(ambiguous(t) for t in forms) == (m in (1, 2, 4)), "restricted classification finite control")
        if m in (1, 2, 4, 8):
            small[str(m)] = forms
    for m in range(1, 2001):
        odd = m
        while odd % 2 == 0:
            odd //= 2
        check(odd == 1 or class_number(m) % 3 == 0, "odd-factor 3-divisibility")
    for a in range(3, 101):
        x, b, c = power_two_hostile(a)
        check(b * b - 4 * x * c == -12 * 4 ** a, "hostile exact discriminant")
        check(0 < abs(b) < x < c and gcd(gcd(x, b), c) == 1, "strict nonambiguous primitive form")
    for d in range(3):
        bits = list(product((0, 1), repeat=d))
        check(set(map(marked_class, bits)) == set(reduced_forms(1 << d)), "marked class torsor bijection")
        for b, v, w in product(bits, repeat=3):
            add = lambda x, y: tuple(a ^ c for a, c in zip(x, y))
            check(marked_class(add(add(b, v), w)) == marked_class(add(b, add(v, w))), "marked regular action")
    path_rows = []
    for k in range(2, 6):
        q, qp = (1 << (k - 1)) + 1, (1 << k) + 1
        decoded = {word(k, x, y) for x, y in product(range(q), repeat=2)}
        check(decoded == paths(k), "literal DFS equals complete pointed-product decoder")
        fibres = Counter()
        colours = Counter()
        for x, y in product(range(qp), repeat=2):
            w = word(k + 1, x, y)
            xx, yy = forget(k, x), forget(k, y)
            check(delete_middle(k, w) == word(k, xx, yy), "literal central deletion")
            check(blue(k + 1, w) == (x == y), "literal blue anti-automorphism")
            fibres[xx, yy] += 1
            colours[xx, yy, x == y] += 1
        for x, y in product(range(q), repeat=2):
            wx, wy = (1 if x == q - 1 else 2), (1 if y == q - 1 else 2)
            check(fibres[x, y] == wx * wy, "4/2/1 fibre law")
            check(colours[x, y, True] == (wx if x == y else 0), "blue lifts and collapsed black sidecar")
        path_rows.append(dict(k=k, source_order=2 * k + 3, target_order=2 * k + 1,
                              fine_paths=qp * qp, coarse_paths=q * q,
                              collapsed_unoriented_black=q - 1))
    hostile = word(3, 0, 2)
    check(not blue(3, hostile) and blue(2, delete_middle(2, hostile)), "actual black-to-blue deletion")
    network_rows = []
    for q in (3, 5, 9):
        weights = (2,) * (q - 1) + (1,)
        W = sum(weights)
        projection = tuple(x for x, w in enumerate(weights) for _ in range(w))
        fine = complete((1,) * W)
        coarse, collapsed = wire(fine, projection)
        check(coarse == complete(weights), "actual complete-edge aggregation")
        check(len(collapsed) == q - 1, "actual within-fibre black edges retained")
        # Wiring composes at the energy level, including newly collapsed edges.
        rho = tuple(x % 2 for x in range(q))
        twice, _ = wire(coarse, rho)
        once, _ = wire(fine, tuple(rho[x] for x in projection))
        check(twice == once, "source-bound wiring composition")
        for z in range(17):
            u = tuple(((x + 1) * z + x * x) % 11 - 5 for x in range(q))
            fine_u = tuple(u[x] for x in projection)
            identity = W * sum(w * a * a for w, a in zip(weights, u)) - sum(w * a for w, a in zip(weights, u)) ** 2
            check(energy(fine, fine_u) == energy(coarse, u) == identity, "exact wiring energy")
        for x, y in combinations(range(q), 2):
            check(resistance(coarse, x, y) == (F(1, weights[x]) + F(1, weights[y])) / W,
                  "independent rational nodal resistance")
        Q = tuple(tuple(F(weights[y] - (x == y), W - 1) for y in range(q)) for x in range(q))
        for a in range(W):
            for y in range(q):
                actual = F(sum(b != a and projection[b] == y for b in range(W)), W - 1)
                check(actual == Q[projection[a]][y], "literal fine walk strong lumpability with holds")
        for x, y in product(range(q), repeat=2):
            check(F(weights[x], W) * Q[x][y] == F(weights[y], W) * Q[y][x], "stationary detailed balance")
            count = sum(projection[a] == x and projection[b] == y for a, b in product(range(W), repeat=2))
            check(F(count, W * W) == F(weights[x] * weights[y], W * W), "actual uniform path pushforward")
        check(Q[0][0] == F(1, W - 1) and Q[-1][-1] == 0, "different fibre holding clocks")
        check(Q[0][1] != F(weights[1], W - weights[0]), "suppressing holds changes the transition operator")
        H = F(2 * sum(w * (w == 2) for w in weights), W)
        check(H == F(4 * (q - 1), 2 * q - 1), "conditional lift entropy in bits")
        network_rows.append(dict(q=q, weights=weights, mass=W, collapsed_edges=len(collapsed),
                                 finite_finite_R=str(resistance(coarse, 0, 1)),
                                 finite_infinity_R=str(resistance(coarse, 0, q - 1)),
                                 conditional_entropy_bits=str(H)))
    disconnected = Network(4, ((0, 1, 1), (2, 3, 1)))
    glued, _ = wire(disconnected, (0, 1, 1, 2))
    check(len(components(disconnected)) == 2 and len(components(glued)) == 1, "wiring creates false fine connectivity")
    check(resistance(disconnected, 0, 3) is None and resistance(glued, 0, 2) == 2, "wiring replaces infinite resistance by two")
    # Same coarse two-edge network; adding the internal connector is invisible to P^T L P.
    connected = Network(4, ((0, 1, 1), (1, 2, 1), (2, 3, 1)))
    check(wire(connected, (0, 1, 1, 2))[0] == glued, "lost internal connectivity coordinate")
    check(resistance(connected, 0, 3) == 3, "actual connector cost is retained only in fine network")
    for g in (F(1, 2), F(1), F(2), F(7)):
        with_connector = Network(4, ((0, 1, 1), (1, 2, g), (2, 3, 1)))
        check(resistance(with_connector, 0, 3) == 2 + 1 / g, "retained two-port response recovers connector cost")
        check(wire(with_connector, (0, 1, 1, 2))[0] == glued, "wiring erases every positive connector value")
    # Non-equitable fibre: vertices zero and one have different rates to the other fibre.
    uneven = Network(3, ((0, 2, 1),))
    check(wire(uneven, (0, 0, 1))[0] == Network(2, ((0, 1, 1),)), "aggregate network exists without lumpability")
    check(sum(c for u, v, c in uneven.edges if 0 in (u, v)) !=
          sum(c for u, v, c in uneven.edges if 1 in (u, v)), "fine states have unequal cross rates")
    reject(reduced_forms, True)
    reject(power_two_hostile, 2)
    reject(marked_class, (False, 1))
    reject(marked_class, (0.0, 1))
    reject(choice, 3, True)
    reject(delete_middle, 2, tuple(float(v) for v in hostile))
    reject(delete_middle, 3, hostile)
    reject(Network, 2, ((0, True, 1),))
    reject(Network, 2, ((0, 1, True),))
    reject(Network, 2, ((0, 1, 1.0),))
    reject(wire, connected, (0, 1, 1, True))
    reject(wire, connected, (0, 1, 1, 3))
    reject(energy, connected, (0.0, 1, 2, 3))
    reject(resistance, connected, False, 3)
    return dict(schema=1, status="PROVED elementary transfers using cited classical class-group facts; FINITE-EXACT controls",
                checks=CHECKS, arithmetic=dict(m_enumerated=[1, 64], divisibility_controls=[1, 2000],
                power_two_hostile_exponents=[3, 100], small_reduced_forms=small),
                path_fibres=path_rows, literal_black_to_blue=dict(source=hostile, target=delete_middle(2, hostile)),
                networks=network_rows, wiring_hostile=dict(fine_edges=disconnected.edges, projection=[0, 1, 1, 2],
                coarse_edges=glued.edges, fine_resistance="infinite", coarse_resistance=2, connected_fine_resistance=3),
                limits=["Choice-state complete graphs are not tournament metagraphs.",
                        "The arithmetic torsor bridge requires a marked origin and generators.",
                        "Wiring energy does not certify arbitrary fine path lifting.",
                        "No claim about completeness of the general 65-idoneal list is used."])


if __name__ == "__main__":
    payload = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    target = Path(__file__).with_name("receipt.json")
    if "--write" in sys.argv:
        temp = target.with_suffix(".tmp")
        temp.write_text(payload, encoding="utf-8", newline="\n")
        temp.replace(target)
    if "--verify" in sys.argv:
        if target.read_text(encoding="utf-8") != payload:
            raise SystemExit("saved receipt mismatch")
    if "--write" in sys.argv or "--verify" in sys.argv:
        print(f"{CHECKS} checks; LF SHA256 {hashlib.sha256(payload.encode()).hexdigest()}")
    else:
        print(payload, end="")
