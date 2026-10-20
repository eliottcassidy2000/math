"""Exact run-shuffle transport checks; standard library; checks survive -O."""
from dataclasses import dataclass
from functools import reduce
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import factorial, gcd
from pathlib import Path

_prior = Path(__file__).resolve().parents[1] / "metagraph_growth_20261020" / "growth_check.py"
_spec = spec_from_file_location("growth_reference", _prior)
ref = module_from_spec(_spec)
_spec.loader.exec_module(ref)
need = ref.need


@dataclass(frozen=True)
class State:
    q: tuple
    weights: tuple
    runs: tuple
    shuffles: tuple


def multinomial(counts):
    return factorial(sum(counts)) // reduce(lambda a, b: a * factorial(b), counts, 1)


def rank_shuffle(word):
    counts = [word.count(j) for j in range(1, max(word, default=0) + 1)]
    rank = 0
    for letter in word:
        for j in range(letter - 1):
            if counts[j]:
                counts[j] -= 1
                rank += multinomial(counts)
                counts[j] += 1
        counts[letter - 1] -= 1
    return rank


def unrank_shuffle(counts, rank):
    counts = list(counts)
    need(0 <= rank < multinomial(counts), "shuffle rank range")
    word = []
    while sum(counts):
        for j in range(len(counts)):
            if not counts[j]:
                continue
            counts[j] -= 1
            size = multinomial(counts)
            if rank < size:
                word.append(j + 1)
                break
            rank -= size
            counts[j] += 1
        else:
            raise RuntimeError("unrank exhaustion")
    return tuple(word)


def owners(weights):
    return tuple((i, r) for i, size in enumerate(weights) for r in range(size))


def weighted(q, weights):
    labels = owners(weights)
    return tuple(sum((int(r < s) if i == j else ref.arc(q, i, j)) << v
                     for v, (j, s) in enumerate(labels))
                 for i, r in labels)


def encode(q, weights, path):
    labels = owners(weights)
    runs, counts = [], [0] * len(q)
    shuffles = [[0] * size for size in weights]
    last = None
    for v in path:
        i, r = labels[v]
        if i == last:
            runs[-1] = (i, runs[-1][1] + 1)
        else:
            counts[i] += 1
            runs.append((i, 1))
        shuffles[i][r] = counts[i]
        last = i
    state = State(q, weights, tuple(runs), tuple(tuple(s) for s in shuffles))
    need(decode(state) == path, "encode internal-run order")
    return state


def decode(state):
    offsets = [sum(state.weights[:i]) for i in range(len(state.q))]
    visits = [0] * len(state.q)
    path = []
    for i, size in state.runs:
        visits[i] += 1
        slots = [offsets[i] + r for r, letter in enumerate(state.shuffles[i])
                 if letter == visits[i]]
        need(len(slots) == size, "run multiplicity")
        path.extend(slots)
    need(sorted(path) == list(range(sum(state.weights))), "path bijection")
    return tuple(path)


def matrix(state):
    return ref.induced(weighted(state.q, state.weights), decode(state))


def normalize(a):
    chains = ref.clone_chains(a)
    q = ref.induced(a, tuple(c[0] for c in chains))
    weights = tuple(map(len, chains))
    concatenation = tuple(v for c in chains for v in c)
    inverse = {v: j for j, v in enumerate(concatenation)}
    path = tuple(inverse[v] for v in range(len(a)))
    state = encode(q, weights, path)
    need(matrix(state) == a, "canonical chain normalization")
    return state


def opposite(a):
    return tuple(((1 << len(a)) - 1) ^ row ^ (1 << v) for v, row in enumerate(a))


def reverse_shuffle(word):
    visits = max(word)
    return tuple(visits + 1 - letter for letter in reversed(word))


def reflect(state):
    shuffles = tuple(reverse_shuffle(s) for s in state.shuffles)
    return State(opposite(state.q), state.weights, tuple(reversed(state.runs)), shuffles)


def blue_from_carrier(state):
    mirror = {}
    for (i, size), (j, other_size) in zip(state.runs, reversed(state.runs)):
        if size != other_size or (i in mirror and mirror[i] != j):
            return False
        mirror[i] = j
    if len(set(mirror.values())) != len(state.q):
        return False
    for i, j in mirror.items():
        if state.weights[i] != state.weights[j] or mirror[j] != i:
            return False
        expected = reverse_shuffle(state.shuffles[i])
        if state.shuffles[j] != expected:
            return False
    return all(ref.arc(state.q, i, j) == ref.arc(state.q, mirror[j], mirror[i])
               for i in range(len(state.q)) for j in range(len(state.q)) if i != j)


def closed_from_carrier(state):
    first_block = state.runs[0][0]
    last_block = state.runs[-1][0]
    if first_block != last_block:
        return bool(ref.arc(state.q, last_block, first_block))
    word = state.shuffles[first_block]
    first_clone = word.index(1)
    last_visit = max(word)
    last_clone = max(r for r, letter in enumerate(word) if letter == last_visit)
    return last_clone < first_clone


def dilate(state, k):
    return State(state.q, tuple(k * w for w in state.weights),
                 tuple((i, k * size) for i, size in state.runs),
                 tuple(tuple(letter for letter in s for _ in range(k))
                       for s in state.shuffles))


def rle_lengths(word):
    lengths = []
    previous = None
    for letter in word:
        if letter == previous:
            lengths[-1] += 1
        else:
            lengths.append(1)
        previous = letter
    return tuple(lengths)


def path_scale(state):
    return reduce(gcd, (m for s in state.shuffles for m in rle_lengths(s)))


def undilate(state, k):
    need(path_scale(state) % k == 0, "nondivisible shuffle scale")
    root = State(state.q, tuple(w // k for w in state.weights),
                 tuple((i, size // k) for i, size in state.runs),
                 tuple(s[::k] for s in state.shuffles))
    need(dilate(root, k) == state, "decorated unique root")
    return root


def sparse_tau(state):
    # Converse weighted core plus a sparse Hamilton-path defect.
    path = decode(state)
    labels = owners(state.weights)
    offsets = [sum(state.weights[:i]) for i in range(len(state.q))]
    reverse_clone = {v: offsets[i] + state.weights[i] - 1 - r
                     for v, (i, r) in enumerate(labels)}
    moved = tuple(reverse_clone[v] for v in path)
    b = list(weighted(opposite(state.q), state.weights))
    for u, v in zip(moved, moved[1:]):
        b[u] ^= 1 << v
        b[v] ^= 1 << u
    return ref.induced(tuple(b), moved)


def append_clone(state, block, ordinal):
    labels = owners(state.weights)
    last_block, last_rank = labels[decode(state)[-1]]
    good = ordinal > last_rank if block == last_block else ref.arc(state.q, last_block, block)
    if not good:
        return None
    runs = list(state.runs)
    shuffles = [list(s) for s in state.shuffles]
    if block == last_block:
        letter = max(shuffles[block])
        runs[-1] = (block, runs[-1][1] + 1)
    else:
        letter = max(shuffles[block]) + 1
        runs.append((block, 1))
    shuffles[block].insert(ordinal, letter)
    weights = list(state.weights)
    weights[block] += 1
    return State(state.q, tuple(weights), tuple(runs), tuple(tuple(s) for s in shuffles))


def inversion_profile(state):
    return tuple(tuple(sum(s[u] == a and s[v] == b
                           for u in range(len(s)) for v in range(u + 1, len(s)))
                       for a in range(1, max(s) + 1) for b in range(1, a))
                 for s in state.shuffles)


def main():
    print("FINITE-EXACT: weighted-core run-shuffle transport.")
    shuffle_cases = 0
    for length in range(1, 8):
        for word in product((1, 2, 3), repeat=length):
            if set(word) != set(range(1, max(word) + 1)):
                continue
            counts = tuple(word.count(j) for j in range(1, max(word) + 1))
            need(unrank_shuffle(counts, rank_shuffle(word)) == word, "enumerative code")
            shuffle_cases += 1
    print("Multiset shuffle rank/unrank words (length <=7, alphabet <=3):", shuffle_cases)

    states = growth = dilation = root_tests = endpoint_cases = 0
    for n in range(1, 7):
        for code in range(1 << len(ref.free_pairs(n))):
            a = ref.tournament(n, code, fixed_path=True)
            state = normalize(a)
            need(blue_from_carrier(state) == (a == ref.transpose(a)), "compressed blue test")
            closure = closed_from_carrier(state)
            need(closure == bool(ref.arc(a, n - 1, 0)), "compressed endpoint closure")
            need(closed_from_carrier(reflect(state)) == closure, "reflection closure")
            need(matrix(reflect(state)) == ref.transpose(a), "reflection rule")
            need(reflect(reflect(state)) == state, "reflection involution")
            need(sparse_tau(state) == ref.tau(a), "sparse tau")
            target = normalize(sparse_tau(state))
            if n >= 3:
                need(closed_from_carrier(target) != closure, "tau flips closure")
                endpoint_cases += 1
            need(sparse_tau(target) == a, "tau involution after normalization")
            need(matrix(reflect(target)) == sparse_tau(reflect(state)), "sigma-tau commutation")
            for k in (2, 3):
                child = dilate(state, k)
                need(matrix(child) == ref.inflate(a, k), "dilation matrix")
                need(path_scale(child) == k * path_scale(state), "path scale multiplication")
                need(reflect(child) == dilate(reflect(state), k), "reflection-dilation")
                need(undilate(child, k) == state, "dilation inverse")
                if n >= 2:
                    need(closed_from_carrier(child) == closure, "dilation closure")
                dilation += 1
            scale = path_scale(state)
            for k in range(1, scale + 1):
                if scale % k == 0:
                    need(matrix(dilate(undilate(state, k), k)) == a, "scale root")
                    root_tests += 1
            for i, size in enumerate(state.weights):
                for ordinal in range(size + 1):
                    child = append_clone(state, i, ordinal)
                    if child is None:
                        continue
                    b = matrix(child)
                    need(ref.induced(b, tuple(range(n))) == a, "append preserves prefix")
                    need(ref.arc(b, n - 1, n), "append HP edge")
                    need(all(len(c) == 1 for c in ref.clone_chains(child.q)), "append core")
                    growth += 1
            states += 1
    print("Complete fixed-path cubes n=1..6: states =", states)
    print("Exact reflection/tau/renormalization and compressed blue controls =", states)
    print("Endpoint closure controls (including tau flip, n>=3) =", endpoint_cases)
    print("Dilation controls =", dilation, "; unique decorated-root controls =", root_tests)
    print("Admissible single-clone endpoint growth controls =", growth)

    c3 = ref.tournament(3, 0, fixed_path=True)
    weights = (4, 1, 2)
    runs = ((0, 2), (1, 1), (2, 2), (0, 2))
    s = State(c3, weights, runs, ((1, 2, 2, 1), (1,), (1, 1)))
    t = State(c3, weights, runs, ((2, 1, 1, 2), (1,), (1, 1)))
    need(inversion_profile(s) == inversion_profile(t), "same inversion profile")
    need(not closed_from_carrier(s) and not closed_from_carrier(t),
         "closure bit does not repair inversion profile")
    scores = [sorted(row.bit_count() for row in sparse_tau(x)) for x in (s, t)]
    need(scores == [[1, 2, 3, 3, 3, 4, 5], [1, 2, 2, 4, 4, 4, 4]],
         "inversion hostile degrees")
    merged_scores = [min(tuple(d), tuple(sorted(6 - x for x in d))) for d in scores]
    need(merged_scores[0] != merged_scores[1], "hostile survives converse")
    print("Equal inversion profiles: A-shuffles1221/2112, weights(4,1,2).")
    print("Both inversion-hostile presentations are OPEN; endpoint bit also agrees.")
    print("Tau degree multisets:", scores)
    print("Lexicographic multiset ranks:", rank_shuffle(s.shuffles[0]), rank_shuffle(t.shuffles[0]))

    for n in range(5, 101):
        transitive = ref.tournament(n, (1 << len(ref.free_pairs(n))) - 1, fixed_path=True)
        target = ref.tau(transitive)
        need(all(len(c) == 1 for c in ref.clone_chains(target)), "unbounded core growth")
    print("Unbounded-growth theorem controls: one core vertex -> N core vertices, N=5..100.")
    print("PASS: no assertions; all checks survive optimized Python.")


if __name__ == "__main__":
    main()
