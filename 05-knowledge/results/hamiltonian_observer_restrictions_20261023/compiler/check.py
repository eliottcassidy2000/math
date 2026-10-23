"""Actual hierarchical tournament paths, partial pins, and targeted arc flags."""
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
from random import Random
from fractions import Fraction
import json
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def arc(x, y):
    return next(((b-a) % 3 == 1 for a, b in zip(x, y) if a != b), False)


class Carrier:
    def __init__(self, cells, depth, levels):
        self.m, self.k = cells, depth
        self.words = tuple(product(range(3), repeat=depth))
        self.n = len(self.words)
        self.N = self.m*self.n
        self.ids = {x: i for i, x in enumerate(self.words)}
        self.levels = levels
        self.components = []
        for r in range(depth):
            parent = list(range(cells))
            def find(x):
                while parent[x] != x:
                    x = parent[x]
                return x
            for (i, j), ell in levels.items():
                if ell > r:
                    a, b = find(i), find(j)
                    parent[max(a, b)] = min(a, b)
            self.components.append(tuple(find(i) for i in range(cells)))
        self.keys = tuple((r, c, w) for r in range(depth)
                          for c in sorted(set(self.components[r]))
                          for w in product(range(3), repeat=r))
        self.key_index = {q: i for i, q in enumerate(self.keys)}
        self.owner_keys = [tuple(self.key_index[(r, self.components[r][i], x[:r])]
                                 for r in range(depth))
                           for i in range(cells) for x in self.words]
        self.P = tuple(i*self.n+j for i in reversed(range(cells)) for j in range(self.n))
        self.A = tuple(tuple(self.edge(i, j) for j in range(self.N)) for i in range(self.N))

    def edge(self, u, v):
        i, a = divmod(u, self.n)
        j, b = divmod(v, self.n)
        if i == j:
            return arc(self.words[a], self.words[b])
        if i > j:
            return not self.edge(v, u)
        ell = self.levels.get((i, j), 0)
        equal = self.words[a][:ell] == self.words[b][:ell]
        outer = 1 <= (j-i) % self.m <= (self.m-1)//2
        return equal if outer else not equal

    def frame(self, coordinates):
        return tuple(i*self.n+self.ids[tuple((x[r]+coordinates[q]) % 3
                                             for r, q in enumerate(self.owner_keys[i*self.n+j]))]
                     for i in range(self.m) for j, x in enumerate(self.words))

    def pins(self, answers, path=None):
        path = self.P if path is None else path
        forced = {}
        for position, target in answers:
            if not 0 <= position < self.N or not 0 <= target < self.N:
                return None
            owner = path[position]
            i, x = divmod(owner, self.n)
            j, y = divmod(target, self.n)
            if i != j:
                return None
            for r, q in enumerate(self.owner_keys[owner]):
                a = (self.words[y][r]-self.words[x][r]) % 3
                if q in forced and forced[q] != a:
                    return None
                forced[q] = a
        return forced

    def exposed(self, positions, path=None):
        path = self.P if path is None else path
        return set().union(*(set(self.owner_keys[path[j]]) for j in positions)) if positions else set()

    def actual_path(self, path, converse=False):
        return len(path) == self.N and len(set(path)) == self.N and all(
            self.A[v][u] if converse else self.A[u][v] for u, v in zip(path, path[1:]))

    def preserves(self, g):
        return all(self.A[i][j] == self.A[g[i]][g[j]]
                   for i in range(self.N) for j in range(self.N))


def audit_answers(C, frames, pins, query):
    forced = C.pins(pins)
    retained = [g for g in frames if all(g[C.P[pos]] == owner for pos, owner in pins)]
    expected = 0 if forced is None else 3**(len(C.keys)-len(forced))
    need(len(retained) == expected, "partial indexed answer fibre")
    if not retained:
        return False
    exposure = C.exposed(query)
    observed = Counter(tuple(g[C.P[j]] for j in query) for g in retained)
    need(len(observed) == 3**len(exposure-set(forced)), "conditional observer image count")
    need(set(observed.values()) == {3**(len(C.keys)-len(set(forced)|exposure))},
         "conditional observer uniform lift count")
    # Actual converse changes clocks and reverses endpoints, not chart coordinates.
    reverse = tuple(reversed(C.P))
    reversed_pins = [(C.N-1-pos, owner) for pos, owner in pins]
    need(C.pins(reversed_pins, reverse) == forced, "converse owner clock transport")
    return True


def small_controls():
    C = Carrier(1, 2, {})
    frames = tuple(C.frame(a) for a in product(range(3), repeat=len(C.keys)))
    need(len(set(frames)) == 81, "native nine-cell parameterization")
    full = {p for p in permutations(range(9)) if C.preserves(p)}
    need(full == set(frames), "full S9 automorphism audit")
    need(C.actual_path(C.P) and C.actual_path(tuple(reversed(C.P)), True), "literal source and converse paths")
    for g in frames:
        need(C.actual_path(tuple(g[x] for x in C.P)), "actual transported Hamiltonian path")
        need(C.actual_path(tuple(g[x] for x in reversed(C.P)), True), "actual transported converse path")
    cases = positives = 0
    for k in range(3):
        for positions in combinations(range(9), k):
            for targets in product(range(9), repeat=k):
                pins = list(zip(positions, targets))
                query = tuple(j for j in (0, 3, 6) if j not in positions)
                positives += audit_answers(C, frames, pins, query)
                cases += 1
    retained = [g for g in frames if g[0] == 0]
    need(len(retained) == 9, "one owner leaves two bottom rotations")
    query = Counter((g[3], g[6]) for g in retained)
    need(len(query) == 9 and set(query.values()) == {1}, "two future owners retain actual nine lifts")
    # Explicit non-native filter hostile: no electrical meaning is claimed.
    guarded = [g for g in retained if g[3] != 5]
    need(len(guarded) == 6, "arbitrary exact guard need not retain a power of three")
    # Nontrivial stabilizer cannot fix a spanning ordered path.
    need(all(tuple(g[x] for x in C.P) != C.P for g in retained if g != tuple(range(9))),
         "no equivariant single-path selector under residual symmetry")
    return dict(frames=81, exhaustive_pin_cases=cases, feasible_pin_cases=positives,
                one_owner_remaining_frames=9, two_future_owner_tuples=9,
                explicitly_non_native_filter_remaining=6)


def mixed_controls():
    C = Carrier(3, 2, {(0, 1): 1, (0, 2): 1, (1, 2): 2})
    need(len(C.keys) == 7, "mixed-depth coordinate dimension")
    frames = tuple(C.frame(a) for a in product(range(3), repeat=len(C.keys)))
    need(len(set(frames)) == 2187, "mixed source frames distinct")
    need(C.actual_path(C.P), "actual mixed source path")
    for g in frames:
        need(C.preserves(g), "all mixed frames preserve every actual arc")
        need(C.actual_path(tuple(g[x] for x in C.P)), "all mixed observer paths genuine")
    # Independent recovery from all 81^3 combinations of genuine internal charts.
    one = Carrier(1, 2, {})
    G = tuple(one.frame(a) for a in product(range(3), repeat=4))
    allowed = {}
    for i, j in combinations(range(3), 2):
        allowed[i, j] = {(a, b) for a, g in enumerate(G) for b, h in enumerate(G)
                          if all(C.A[i*9+u][j*9+v] == C.A[i*9+g[u]][j*9+h[v]]
                                 for u in range(9) for v in range(9))}
    recovered = {tuple(i*9+u for i, a in enumerate(choice) for u in G[a])
                 for choice in product(range(81), repeat=3)
                 if all((choice[i], choice[j]) in allowed[i, j]
                        for i, j in combinations(range(3), 2))}
    need(recovered == set(frames), "full named-cell group recovered from actual blocks")
    rng = Random(57031)
    feasible = 0
    for test in range(300):
        g = rng.choice(frames)
        positions = rng.sample(range(27), rng.randrange(7))
        pins = [(p, g[C.P[p]]) for p in positions]
        if test % 3 == 0 and pins:
            p, owner = pins[0]
            cell, leaf = divmod(owner, 9)
            pins[0] = (p, cell*9+(leaf+1) % 9)
        query = rng.sample(range(27), rng.randrange(8))
        feasible += audit_answers(C, frames, pins, query)
    # Equal prefix coordinates shared through fine edge cannot be independently pinned.
    pos1 = C.P.index(9)
    pos2 = C.P.index(18)
    need(C.pins([(pos1, 9), (pos2, 19)]) is None, "shared full-resolution owner conflict")
    need(C.pins([(pos1, 18)]) is None, "named-cell guard")
    return dict(order=27, source_frame_count=2187, coordinate_count=7,
                audited_internal_chart_triples=81**3, seeded_pin_cases=300,
                feasible_pin_cases=feasible, seed=57031)


def order63_controls():
    fine = {1, 3, 5}
    levels = {(i, j): 2 if i in fine and j in fine else 1
              for i, j in combinations(range(7), 2)}
    C = Carrier(7, 2, levels)
    need(len(C.keys) == 16 and C.N == 63, "incoming native host dimension")
    need({sum(row) for row in C.A} == {31}, "actual regular order63")
    need(C.actual_path(C.P) and C.actual_path(tuple(reversed(C.P)), True), "actual63 path and converse")
    for q in range(len(C.keys)):
        a = [0]*len(C.keys)
        a[q] = 1
        g = C.frame(a)
        need(C.preserves(g), "every independent native coordinate generator")
        need(C.actual_path(tuple(g[x] for x in C.P)), "native generator path")
    bottom = sorted(q for q, key in enumerate(C.keys) if key[0] == 1)
    ids = {q: i for i, q in enumerate(bottom)}
    slot = [ids[C.owner_keys[u][-1]] for u in C.P]
    flag_masks = tuple(sorted({(1 << a) | (1 << b) for a, b in zip(slot, slot[1:])}))
    n = len(bottom)
    adjacent = [0]*n
    for mask in flag_masks:
        if mask.bit_count() == 2:
            a, b = [i for i in range(n) if mask >> i & 1]
            adjacent[a] |= 1 << b
            adjacent[b] |= 1 << a

    @lru_cache(None)
    def matching(mask):
        if not mask:
            return 0
        first = mask & -mask
        i = first.bit_length()-1
        rest = mask ^ first
        answer = matching(rest)
        choices = adjacent[i] & rest
        while choices:
            bit = choices & -choices
            answer = max(answer, 1+matching(rest ^ bit))
            choices ^= bit
        return answer

    @lru_cache(None)
    def actual_flag_cover(mask):
        if not mask:
            return 0
        first = mask & -mask
        return 1+min(actual_flag_cover(mask & ~f) for f in flag_masks if f & first)

    histogram = Counter()
    for required in range(1 << n):
        expected = required.bit_count()-matching(required)
        need(actual_flag_cover(required) == expected, "every targeted flag cost from actual indexed arcs")
        histogram[expected] += 1
    need(n == 15 and matching((1 << n)-1) == 6 and actual_flag_cover((1 << n)-1) == 9,
         "inherited63 shared-incidence hostile")
    inherited = HERE.parents[1] / "session42_shared_lattices_square_colors_20261023/hierarchical-kernels/receipt.json"
    old = json.loads(inherited.read_text(encoding="utf-8"))["mixed_order63_hostile"]
    flags = old["arc_flag_positions"]
    answers = [(p, C.P[p]) for j in flags for p in (j, j+1)]
    forced = C.pins(answers)
    need(forced is not None and len(forced) == 16, "incoming exact flag clocks kill every native coordinate")
    need(old["chart_order_power_of_3"] == 16 and old["arc_flag_base"] == 9,
         "incoming hostile receipt retained")
    # For a targeted subset only its unresolved slots matter, even after nonzero answers.
    rng = Random(639)
    for _ in range(200):
        a = [rng.randrange(3) for _ in C.keys]
        g = C.frame(a)
        pinned = rng.sample(range(63), rng.randrange(12))
        query = rng.sample(range(63), rng.randrange(12))
        F = C.pins([(j, g[C.P[j]]) for j in pinned])
        need(F is not None, "actual nonzero partial answers")
        unobserved_slots = {slot[j] for j in query}-{slot[j] for j in pinned}
        required = sum(1 << s for s in unobserved_slots)
        known = C.exposed(query) <= set(F)
        need(known == (required == 0), "conditional full-owner determination equals slot coverage")
        need(actual_flag_cover(required) == len(unobserved_slots)-matching(required),
             "conditional targeted flag cost")
    return dict(order=63, coordinate_count=16, frame_count=3**16,
                bottom_slots=15, full_path_flags=9,
                exhaustive_target_slot_sets=1 << n, target_cost_histogram=dict(sorted(histogram.items())),
                incoming_actual_arc_clocks=flags, seeded_conditional_cases=200)


def move_controls():
    C = Carrier(3, 2, {(0, 1): 1, (0, 2): 1, (1, 2): 2})
    coordinates = tuple(product(range(3), repeat=len(C.keys)))
    frame = {a: C.frame(a) for a in coordinates}
    directed = 0
    for a, g in frame.items():
        for q, (r, component, prefix) in enumerate(C.keys):
            cell = C.components[r].index(component)
            source = cell*C.n+C.ids[prefix+(0,)*(C.k-r)]
            target_prefix = C.words[g[source] % C.n][:r]
            for delta in (1, 2):
                b = list(a)
                b[q] = (b[q]+delta) % 3
                g1 = frame[tuple(b)]
                literal = []
                for image in g:
                    i, leaf = divmod(image, C.n)
                    word = list(C.words[leaf])
                    if C.components[r][i] == component and tuple(word[:r]) == target_prefix:
                        word[r] = (word[r]+delta) % 3
                    literal.append(i*C.n+C.ids[tuple(word)])
                need(tuple(literal) == g1, "coordinate move equals actual target-prefix rotation")
                need(C.actual_path(tuple(g1[u] for u in C.P)), "every moved observer is actual Hamiltonian")
                directed += 1
    # Every subset of actual owner clocks on R9, with nonzero answers.
    C = Carrier(1, 2, {})
    all_a = tuple(product(range(3), repeat=4))
    all_g = {a: C.frame(a) for a in all_a}
    rng = Random(913)
    dimensions = Counter()
    for mask in range(1 << 9):
        source_a = rng.choice(all_a)
        source_g = all_g[source_a]
        pins = [(j, source_g[C.P[j]]) for j in range(9) if mask >> j & 1]
        F = C.pins(pins)
        free = [q for q in range(4) if q not in F]
        allowed = {a for a in all_a if all(a[q] == value for q, value in F.items())}
        start = min(allowed)
        distance = {start: 0}
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for q in free:
                for delta in (1, 2):
                    b = list(a)
                    b[q] = (b[q]+delta) % 3
                    b = tuple(b)
                    need(b in allowed and all(all_g[b][C.P[j]] == owner for j, owner in pins),
                         "actual neighbouring path preserves all pins")
                    if b not in distance:
                        distance[b] = distance[a]+1
                        queue.append(b)
        need(set(distance) == allowed, "conditional native Hamming graph connected")
        need(max(distance.values()) == len(free), "conditional native Hamming graph eccentricity")
        dimensions[len(free)] += 1
    # A literal extra path-label predicate can remove every available move.
    guarded = [a for a in all_a if all_g[a][0] == 0
               and all_g[a][3] % 3 == all_g[a][6] % 3]
    need(len(guarded) == 3 and all(sum(x != y for x, y in zip(a, b)) == 2
                                 for a, b in combinations(guarded, 2)),
         "rechecking an arbitrary guard does not preserve native-move connectivity")
    # A full rational tensor eigenbasis: each coordinate has constant and
    # two independent zero-sum columns. Check every vector at every state.
    cols = ((1, 1, 1), (1, -1, 0), (1, 1, -2))
    spectrum = {}
    for r in range(1, 5):
        states = tuple(product(range(3), repeat=r))
        eigencount = Counter()
        for types in states:
            j = sum(t != 0 for t in types)
            eigenvalue = Fraction(r-j, r)
            eigencount[str(eigenvalue)] += 1
            def value(state):
                answer = 1
                for pos in range(r):
                    answer *= cols[types[pos]][state[pos]]
                return answer
            for state in states:
                total = 0
                for pos in range(r):
                    for digit in range(3):
                        next_state = list(state)
                        next_state[pos] = digit
                        total += value(next_state)
                need(Fraction(total, 3*r) == eigenvalue*value(state),
                     "exact rational full heatbath eigenbasis")
        spectrum[str(r)] = dict(sorted(eigencount.items()))
    return dict(actual_mixed27_directed_coordinate_moves=directed,
                exhaustive_nine_owner_clock_subsets=1 << 9,
                disconnected_guard_survivors=len(guarded),
                conditional_dimension_histogram=dict(sorted(dimensions.items())),
                exact_heatbath_spectra=spectrum)


def footprint_controls():
    C = Carrier(1, 2, {})
    states = tuple(product(range(3), repeat=4))
    frames = {a: C.frame(a) for a in states}
    adjacent = {a: {b for b in states if sum(x != y for x, y in zip(a, b)) == 1}
                for a in states}
    visible_adjacent = {x: {y for y in range(9)
                           if (x//3 == y//3) != (x % 3 == y % 3)} for x in range(9)}

    def components(vertices, neighbours):
        unseen = set(vertices)
        result = set()
        while unseen:
            start = unseen.pop()
            comp = {start}
            pending = [start]
            while pending:
                v = pending.pop()
                new = neighbours[v] & unseen
                comp |= new
                unseen -= new
                pending.extend(new)
            result.add(frozenset(comp))
        return result

    def diameter(vertices, neighbours):
        answer = 0
        for start in vertices:
            distance = {start: 0}
            queue = deque([start])
            while queue:
                v = queue.popleft()
                for w in neighbours[v] & vertices:
                    if w not in distance:
                        distance[w] = distance[v]+1
                        queue.append(w)
            need(set(distance) == set(vertices), "guard component literal connectivity")
            answer = max(answer, max(distance.values()))
        return answer

    histogram = Counter()
    for mask in range(512):
        A = {x for x in range(9) if mask >> x & 1}
        retained = {a for a, g in frames.items() if g[0] in A}
        need(len(retained) == 9*len(A), "actual indexed observer footprint count")
        visible_components = components(A, visible_adjacent)
        actual_components = components(retained, adjacent)
        expected = {frozenset(a for a in retained if frames[a][0] in piece)
                    for piece in visible_components}
        need(actual_components == expected, "all observer-subset guards retain exact invisible product fibres")
        histogram[len(actual_components)] += 1
    reports = {}
    for name, A in (("connected", {0, 1}), ("disconnected", {0, 4, 8})):
        retained = {a for a, g in frames.items() if g[0] in A}
        actual_components = components(retained, adjacent)
        sizes, diameters = [], []
        for piece in actual_components:
            visible = {frames[a][0] for a in piece}
            actual_diameter = diameter(piece, adjacent)
            need(actual_diameter == diameter(visible, visible_adjacent)+2,
                 "native guarded component diameter retains invisible coordinates")
            sizes.append(len(piece))
            diameters.append(actual_diameter)
        reports[name] = dict(accepted_observer_values=sorted(A),
                             component_sizes=sorted(sizes), component_diameters=sorted(diameters))
    return dict(actual_carrier_order=9, conditional_dimension=4,
                visible_dimension=2, invisible_dimension=2,
                all_observer_subset_guards=512,
                component_histogram=dict(sorted(histogram.items())), examples=reports)


def main():
    out = dict(status="FINITE-EXACT on supplied actual native path orbits",
               regular_nine=small_controls(), mixed_twenty_seven=mixed_controls(),
               regular_sixty_three=order63_controls(), native_moves=move_controls(),
               observer_footprints=footprint_controls(), checks=CHECKS)
    text = json.dumps(out, indent=2, sort_keys=True)+"\n"
    target = HERE / "receipt.json"
    if "--write" in sys.argv:
        target.write_text(text, encoding="utf-8")
    else:
        need(json.loads(target.read_text(encoding="utf-8")) == json.loads(text), "frozen exact receipt")
    print(json.dumps(dict(status="PASS", checks=CHECKS,
                          exhaustive63targets=out["regular_sixty_three"]["exhaustive_target_slot_sets"],
                          native_orbit_sizes=[81, 2187, 3**16])))


if __name__ == "__main__":
    main()
