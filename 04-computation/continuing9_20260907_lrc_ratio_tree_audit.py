"""Independent exact referee: ratio-tree gluing and the composite 7200 pair.

No producer import or execution. Uses integer matrix kernels, intersections
of literal periodic interval lists, and a full literal grid at EVERY wall
and chamber representative, instead of the producer's interval-floor counter.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm
from hashlib import sha256
import argparse
import json
import sys
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", newline="\n")
HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--repo", type=Path,
                    default=HERE.parent if HERE.name == "04-computation" else Path("C:/w/s0905"))
parser.add_argument("--certificate", type=Path)
args = parser.parse_args()
REPO = args.repo
CERT = args.certificate or HERE / "continuing9_20260907_lrc_ratio_tree_certificate.json"
if not CERT.exists():
    CERT = REPO / "05-knowledge/results" / CERT.name
J = json.loads(CERT.read_bytes())
gates = 0


def need(test, message):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(message)


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


producer_stem = "continuing9_20260907_lrc_ratio_tree"
producer_pins = {
    ".py": "eacacacdf4bd4e063036d82785c72bb358ebbf315894c76e0fce91646bde4187",
    ".out": "45d19166b682d4098280d7d3f0f0feec06825294d6af7d8da4e1620e9435ac9a",
    "_certificate.json": "f93d2c8b4f112027bbba23fa529742cd12ac8107e7c9bfdf57f15d6666799533",
}
for suffix, expected in producer_pins.items():
    path = HERE / (producer_stem + suffix)
    if not path.exists():
        path = REPO / ("04-computation" if suffix == ".py" else "05-knowledge/results") / path.name
    need(sha256(path.read_bytes()).hexdigest() == expected, "frozen producer pin " + suffix)


def matrix_row(n, edges):
    matrix = []
    for i, j, numerator, denominator in edges:
        row = [0]*n
        row[i], row[j] = numerator, -denominator
        matrix.append(row)
    basis = sp.Matrix(matrix).nullspace()
    need(len(basis) == 1, "connected ratio constraints have a one-dimensional kernel")
    vector = basis[0]
    scale = lcm(*(int(x.q) for x in vector))
    integers = [int(x*scale) for x in vector]
    sign = 1 if integers[0] > 0 else -1
    gg = gcd(*integers)
    answer = tuple(sign*x//gg for x in integers)
    need(all(x > 0 for x in answer) and gcd(*answer) == 1, "primitive positive kernel")
    return answer


def allpair_factors(a, b):
    gg = gcd(a, b)
    return gg, min(a, b)//gg, max(a, b)//gg


def atlas_sum(s):
    return all(int(prime) % 3 == 2 and int(exponent) <= 2
               for prime, exponent in sp.factorint(s).items())


raw = (REPO/"05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json").read_bytes()
need(sha256(raw).hexdigest() == J["inherited_sha256"] ==
     "580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d",
     "inherited clock certificate pin")
clock = next(row for row in json.loads(raw)["clocks"] if row["t"] == 7200)
weights = {tuple(k): v for k, v in clock["weights"]}
raw = (REPO/"04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json").read_bytes()
need(sha256(raw).hexdigest() == J["profile_sha256"] ==
     "935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f",
     "inherited profile pin")
inherited = json.loads(raw)
allowed = {int(k): {(int(c), tuple(v)) for c, v in rec["profiles"]}
           for k, rec in inherited["levels"].items()}


def valid(word):
    if gcd(*word) != 1:
        return False
    # Traverse masks in numerical order, independently of the producer's
    # decreasing body-cardinality combination ordering.
    partial = [0]*128
    for mask in range(1, 128):
        bit = mask & -mask
        partial[mask] = gcd(partial[mask ^ bit], word[bit.bit_length()-1])
    for mask in range(1, 127):
        c = partial[mask]
        other = sorted(gcd(c, word[i]) for i in range(7) if not mask >> i & 1)
        if (c, tuple(other)) not in allowed[len(other)]:
            return False
    return True


def excess(word):
    # Direct ceil of a rational quantity, different from the producer's
    # ceiling encoded with 7200+7d-1.
    return sum(d*(-((-7200)//(7*d))) for d in word)-7200


# General prime-potential formula on fresh finite trees and every possible
# divisor-labelled margin vector in a small complete universe.
generic = 0
for row in [(1, 4, 6), (2, 3, 5), (3, 7, 10), (4, 9, 25), (5, 6, 14), (1, 2, 1)]:
    edges = [[0, 1, row[1], row[0]], [1, 2, row[2], row[1]]]
    u = matrix_row(3, edges)
    need(u == row, "literal primitive row recovered")
    rel = [F(x, row[0]) for x in row]
    for t in range(1, 31):
        primes = sp.factorint(t)
        predicted = [1]*3
        for prime, cap in primes.items():
            potentials = []
            for xx in rel:
                num, den = xx.numerator, xx.denominator
                value = 0
                while num % prime == 0:
                    value += 1
                    num //= prime
                while den % prime == 0:
                    value -= 1
                    den //= prime
                potentials.append(value)
            base = min(potentials)
            for i, value in enumerate(potentials):
                predicted[i] *= int(prime)**min(int(cap), value-base)
        need(predicted == [gcd(t, x) for x in u], "prime potential margin formula")
        divisors = sp.divisors(t)
        for margins in product(divisors, repeat=3):
            got = len(set(u)) == 3 and tuple(predicted) == margins
            actual = len(set(u)) == 3 and margins == tuple(gcd(t, x) for x in u)
            need(got == actual, "all positive and false margin labellings")
            generic += 1
need((gcd(2, 1), gcd(2, 4), gcd(2, 6)) == (1, 2, 2), "local-margin hostile global failure")
need(all(atlas_sum(5) for _ in range(2)), "hostile uses strict atlas edges")

# Recover all named fixed-tree cases by matrix kernels. MST optimality uses
# the cycle property, not a second execution of the producer's Kruskal code.
counts = [0, 0, 0, 0]
valid_words = set()
need(len(J["fixed_tree_bank"]) == len(clock["survivors"]) == 15, "15 complete projected survivors")
for bank in J["fixed_tree_bank"]:
    word = bank["word"]
    tree = bank["tree"]
    need(valid(word) and excess(word) == bank["E"], "all inherited profiles on fixed word")
    need(len(tree) == 6 and len({tuple(sorted((i, j))) for i, j, *rest in tree}) == 6,
         "six distinct tree edges")
    adjacency = [[] for _ in range(7)]
    choices = []
    total = 0
    for i, j, e, p, q in tree:
        cost = weights[tuple(sorted((word[i], word[j])))][0]
        total += cost
        adjacency[i].append((j, cost))
        adjacency[j].append((i, cost))
        pair = []
        if (e*gcd(7200//e, p), e*gcd(7200//e, q)) == (word[i], word[j]):
            pair.append((i, j, q, p))
        if (e*gcd(7200//e, q), e*gcd(7200//e, p)) == (word[i], word[j]):
            pair.append((i, j, p, q))
        choices.append(pair)
    need(total == bank["M"], "tree attains inherited total")
    for start, finish in combinations(range(7), 2):
        todo = [(start, -1, 0)]
        found = None
        for here, prior, cost in todo:
            if here == finish:
                found = cost
                break
            todo.extend((nxt, here, max(cost, edge_cost))
                        for nxt, edge_cost in adjacency[here] if nxt != prior)
        need(found is not None and found <= weights[tuple(sorted((word[start], word[finish])))][0],
             "MST cycle criterion")
    expected = {tuple(tuple(e) for e in edges) for edges in product(*choices)}
    actual = {tuple(tuple(e) for e in row["ratios"]) for row in bank["orientations"]}
    need(expected == actual, "complete local first-attainer orientation bank")
    for record in bank["orientations"]:
        u = matrix_row(7, record["ratios"])
        need(list(u) == record["row"], "independent exact primitive realization")
        margins = tuple(gcd(7200, x) for x in u)
        good = len(set(u)) == 7 and margins == tuple(word)
        need(list(margins) == record["margins"] and good == record["valid"], "declared bank status")
        counts[0] += 1
        if good:
            counts[1] += 1
            valid_words.add(tuple(word))
        elif len(set(u)) < 7:
            counts[2] += 1
        else:
            counts[3] += 1
need(counts == [29, 20, 8, 1] and len(valid_words) == 10, "bank summary")

# The actual seven-tail control, including all non-atlas overlaps.
U = tuple(J["actual_row"])
need(gcd(*U) == 1 and len(set(U)) == 7, "actual seven-tail primitive distinct")
need(tuple(gcd(7200, x) for x in U) == tuple(J["actual_margins"]), "actual margins")
graph = []
for i, j in combinations(range(7), 2):
    D, p, q = allpair_factors(U[i], U[j])
    if p+q <= 356 and atlas_sum(p+q):
        graph.append([i, j, gcd(7200, D), p, q, D//gcd(7200, D)])
need(graph == J["actual_graph"], "entire strict graph reconstructed")
marginals = [0]*7
pairs = {(i, j): 0 for i, j in combinations(range(7), 2)}
histogram = [0]*8
for j in range(7200):
    active = []
    for i, speed in enumerate(U):
        rem = (speed*(2*j+1)) % 14400
        if 14*min(rem, 14400-rem) < 14400:
            active.append(i)
            marginals[i] += 1
    histogram[len(active)] += 1
    for pair in combinations(active, 2):
        pairs[pair] += 1
need(marginals == J["danger_marginals"] and histogram == J["danger_histogram"], "literal7200 histogram")
need([[i, j, pairs[i, j]] for i, j in pairs] == J["paircounts"], "all21 pair intersections")
need(all(pairs[i, j] == 0 for i, j, *rest in graph), "all atlas minima coexist")
need(sum(marginals)-(7200-histogram[0]) == 2036 and pairs[0, 4] == 232,
     "non-atlas overlap observed")

# Independent continuous geometry: intersect periodic danger intervals for
# each speed by a two-pointer sweep. Intervals straddling zero are retained
# as single real intervals; the translate straddling one is dropped.
n, p, q = 3600, 578, 801
lists = [[(F(14*k-1, 14*v), F(14*k+1, 14*v)) for k in range(v+1)] for v in (p, q)]
i = j = 0
intervals = []
while i < len(lists[0]) and j < len(lists[1]):
    a, b = lists[0][i]
    c, d = lists[1][j]
    lo, hi = max(a, c), min(b, d)
    if lo < hi and lo+hi < 2:
        intervals.append((lo, hi))
    if b < d:
        i += 1
    elif d < b:
        j += 1
    else:
        i += 1
        j += 1
composite = J["composite"]
L = composite["L"]
need(intervals == [(F(a, L), F(b, L)) for a, b in composite["intervals"]],
     "197 intervals independently reconstructed")
events = sorted({n*x % 1 for ab in intervals for x in ab})
need(events == [F(x, L) for x in composite["walls"]], "all344 exact phase walls")
phases = events + [(events[i]+(events[i+1] if i+1 < len(events) else events[0]+1))/2
                   for i in range(len(events))]
counts_by_phase = {}
for alpha in phases:
    den = alpha.denominator*n
    count = 0
    for j in range(n):
        location = j*alpha.denominator+alpha.numerator
        r1, r2 = p*location % den, q*location % den
        if 14*min(r1, den-r1) < den and 14*min(r2, den-r2) < den:
            count += 1
    counts_by_phase[alpha] = count
need(len(phases) == 688, "every wall and every chamber represented")
for r, count in composite["phase2_counts"]:
    need(counts_by_phase[F(r, 2*L)] == count, "literal full grid at every arrangement representative")
need(min(counts_by_phase.values()) == 57 and max(counts_by_phase.values()) == 116,
     "all-phase minimum57 and maximum116")
need(counts_by_phase[F(*composite["minimizer"])] == 57, "attained wall minimizer")
need(2*57 == 114 and 2*counts_by_phase[F(1, 2)] == 232, "uniform114 versus phase232")

# Complete inherited-profile conditional universe, retaining all126 proper
# position subsets. Duplicate margins are allowed and their positions retained.
alphabet = [d for d in inherited["levels"]["6"]["gcds"] if 7200 % d == 0]
conditional = []
attempted = 0
for rest in combinations_with_replacement(sorted(alphabet), 4):
    attempted += 1
    word = tuple(sorted((4, 18, 24)+rest))
    if valid(word):
        conditional.append([list(word), excess(word)])
conditional.sort()
need(attempted == 23751, "complete four-free-position universe")
need(len(conditional) == 1033 and conditional == J["conditional_words"], "all admissible words independently recovered")
need(sha256(canonical(conditional)).hexdigest() == J["conditional_words_sha256"], "complete word bank pin")
need(max(e for word, e in conditional) == 108 and 114-108 == 6, "strict uniform closure margin6")

# Full13 positive control; this is not an actual W=Vdec assertion.
control = J["full13_control"]
body, tail = tuple(control["body"]), tuple(control["tails"])
need(len(body) == 6 and len(tail) == 7 and gcd(*body) == 7200 and
     len(set(body+tail)) == 13 and gcd(*(body+tail)) == 1, "full13 primitive distinct control")
need(valid(tuple(sorted(gcd(7200, x) for x in tail))), "all126 projected profiles on full13 control")
edges = []
for i, j in combinations(range(7), 2):
    gg, pp, qq = allpair_factors(tail[i], tail[j])
    if pp+qq <= 356 and atlas_sum(pp+qq):
        edges.append([i, j])
need(edges == control["complement_edges"] == [[0, 1], [1, 2]], "complement graph disconnected")
need(tail[:3] == (2*578, 2*132, 2*801), "actual two-arm family")
alpha = F(*control["alpha"])
den = alpha.denominator*7200
safe = []
for j in range(7200):
    pos = j*alpha.denominator+alpha.numerator
    if all(14*min(v*pos % den, den-(v*pos % den)) >= den for v in body+tail):
        safe.append(j)
need(safe == control["safe_lifts"], "every declared full13 safe lift independently reconstructed")
need(len(safe) >= 114-excess(tuple(gcd(7200, x) for x in tail)), "literal family consequence")

print("GENERAL prime-potential normalization: complete false-margin controls", generic)
print("FIXED_TREE bank:29 orientations,20 distinct realizations on10 words,8 duplicate failures,1 valuation failure")
print("ACTUAL seven-tail: all21 intersections and7200 cells agree; six atlas credits0 at alpha1/2")
print("COMPOSITE: independent197 interval intersections;344 walls+344 chambers; every3600-cell grid enumerated")
print("COMPOSITE min57 max116; sheet credit114; declared-phase credit232; minimizer", F(*composite["minimizer"]))
print("CONDITIONAL:23751 candidates,1033 complete126-profile survivors,maximum108,uniform margin6")
print("FULL13 positive control: all safe lifts", len(safe), "; no W=Vdec inference")
print("CERTIFICATE_SHA256", sha256(CERT.read_bytes()).hexdigest())
print("PASS", gates, "always-active exact gates; raw LF")
