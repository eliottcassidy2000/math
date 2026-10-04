"""Exact triangle/star flow fibres and Heawood surface-boundary controls.

No third-party imports. All checks survive -O. The graph, rotations, linear
algebra and flow generator are reconstructed here rather than imported.
"""
from collections import Counter, deque
from itertools import combinations, permutations, product
from math import prod


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


D = (1, 2, 4)
EDGES = tuple((p, 7+t) for p in range(7) for t in range(7) if (p-t) % 7 in D)
EINDEX = {edge: i for i, edge in enumerate(EDGES)}


class Group:
    def __init__(self, q, xor=False):
        self.q, self.xor = q, xor
        self.add = (lambda a, b: a ^ b) if xor else (lambda a, b: (a+b) % q)
        self.neg = (lambda a: a) if xor else (lambda a: -a % q)

    def sub(self, a, b):
        return self.add(a, self.neg(b))


def rotations(line_sign):
    orders = {p: [7+(p-d) % 7 for d in D] for p in range(7)}
    orders.update({7+t: [(t+d) % 7 for d in D] for t in range(7)})
    darts = tuple(d for u, v in EDGES for d in ((u, v), (v, u)))
    successor = {}
    for u, v in darts:
        cycle = orders[v]
        step = 1 if v < 7 else line_sign
        successor[(u, v)] = (v, cycle[(cycle.index(u)+step) % 3])
    unseen, faces = set(darts), []
    while unseen:
        first = min(unseen)
        current, face = first, []
        while current in unseen:
            unseen.remove(current)
            face.append(current)
            current = successor[current]
        need(current == first, 'face closes at its own starting dart')
        faces.append(tuple(face))
    owner = {dart: f for f, face in enumerate(faces) for dart in face}
    dual = Counter(tuple(sorted((owner[(u, v)], owner[(v, u)]))) for u, v in EDGES)
    boundaries = []
    for face in faces:
        row = [0]*len(EDGES)
        for u, v in face:
            e = (u, v) if u < 7 else (v, u)
            row[EINDEX[e]] += 1 if u < 7 else -1
        boundaries.append(tuple(row))
    return faces, dual, boundaries


def rref(rows, modulus):
    rows = [[x % modulus for x in row] for row in rows]
    lead = 0
    pivots = []
    for col in range(len(EDGES)):
        found = next((i for i in range(lead, len(rows)) if rows[i][col]), None)
        if found is None:
            continue
        rows[lead], rows[found] = rows[found], rows[lead]
        scale = pow(rows[lead][col], -1, modulus)
        rows[lead] = [scale*x % modulus for x in rows[lead]]
        for i in range(len(rows)):
            if i != lead and rows[i][col]:
                scale = rows[i][col]
                rows[i] = [(x-scale*y) % modulus for x, y in zip(rows[i], rows[lead])]
        pivots.append(col)
        lead += 1
        if lead == len(rows):
            break
    return tuple((p, tuple(rows[i])) for i, p in enumerate(pivots))


def reduce_mod_boundaries(flow, basis, group):
    current = list(flow)
    for pivot, row in basis:
        factor = current[pivot]
        if not factor:
            continue
        if group.xor:
            current = [x ^ (factor if y else 0) for x, y in zip(current, row)]
        else:
            current = [(x-factor*y) % group.q for x, y in zip(current, row)]
    return tuple(current)


def spanning_tree():
    # Make the first two point-0 edges cotree variables, for normalized F2^3 enumeration.
    adjacency = [[] for _ in range(14)]
    for e, (u, v) in enumerate(EDGES):
        if e in (0, 1):
            continue
        adjacency[u].append((v, e))
        adjacency[v].append((u, e))
    order, parent, queue = [0], {}, deque([0])
    visited = {0}
    while queue:
        u = queue.popleft()
        for v, e in adjacency[u]:
            if v not in visited:
                visited.add(v)
                order.append(v)
                parent[v] = (u, e)
                queue.append(v)
    need(len(visited) == 14, 'connected tree after two exclusions')
    tree_edges = {e for u, e in parent.values()}
    chords = tuple(e for e in range(21) if e not in tree_edges)
    need(len(chords) == 8 and chords[:2] == (0, 1), 'eight cycle coordinates')
    return tuple(reversed(order[1:])), parent, chords


LEAF_ORDER, PARENT, CHORDS = spanning_tree()


def evaluate_flow(chord_values, group):
    values, balance = [0]*21, [0]*14
    for e, value in zip(CHORDS, chord_values):
        values[e] = value
        u, v = EDGES[e]
        balance[u] = group.add(balance[u], value)
        balance[v] = group.sub(balance[v], value)
    for v in LEAF_ORDER:
        u, edge = PARENT[v]
        value = group.neg(balance[v]) if v < 7 else balance[v]
        if value == 0:
            return None
        values[edge] = value
        balance[u] = group.add(balance[u], value if u < 7 else group.neg(value))
    need(balance[0] == 0, 'final Kirchhoff equation')
    return tuple(values)


def nowhere_zero_flows(group, normalized=False):
    exponent = 6 if normalized else 8
    for free in product(range(1, group.q), repeat=exponent):
        flow = evaluate_flow((1, 2)+free if normalized else free, group)
        if flow is not None:
            yield flow


def check_flow(values, group):
    total = [0]*14
    for (u, v), value in zip(EDGES, values):
        total[u] = group.add(total[u], value)
        total[v] = group.sub(total[v], value)
    return not any(total)


def lift_to_paley(flow, choices, group):
    # Triangle L_t is oriented t+1 -> t+2 -> t+4 -> t+1.
    values = {}
    for t in range(7):
        x, y, z = [(t+d) % 7 for d in D]
        bx, by, bz = [flow[EINDEX[(p, 7+t)]] for p in (x, y, z)]
        need(group.add(group.add(bx, by), bz) == 0, 'line-vertex conservation')
        a = choices[t]
        b, c = group.add(a, by), group.sub(a, bx)
        values[(x, y)], values[(y, z)], values[(z, x)] = a, b, c
    need(len(values) == 21, 'Fano triangles partition K7 edges')
    balances = [0]*7
    for (u, v), value in values.items():
        balances[u] = group.add(balances[u], value)
        balances[v] = group.sub(balances[v], value)
    need(not any(balances), 'lift is a K7 circulation')
    return values


def project_from_paley(values, group):
    result = [0]*21
    for t in range(7):
        x, y, z = [(t+d) % 7 for d in D]
        a, b, c = values[(x, y)], values[(y, z)], values[(z, x)]
        for point, value in ((x, group.sub(a, c)), (y, group.sub(b, a)),
                             (z, group.sub(c, b))):
            result[EINDEX[(point, 7+t)]] = value
    need(check_flow(result, group), 'projected boundary is a Heawood flow')
    return tuple(result)


def main():
    print('GEOMETRY / FLOW OBSTRUCTIONS: EXACT AUDIT')
    local_cases = 0
    for group in [Group(q) for q in range(2, 9)]+[Group(4, True), Group(8, True)]:
        for bx, by in product(range(group.q), repeat=2):
            bz = group.neg(group.add(bx, by))
            count = sum(all((a, group.add(a, by), group.sub(a, bx)))
                        for a in range(group.q))
            forbidden = {0, group.neg(by), bx}
            need(count == group.q-len(forbidden), 'all triangle boundary types')
            if bx and by and bz:
                need(count == group.q-3, 'nonzero star has q-3 lifts')
            local_cases += 1
    print('Local triangle fibres:', local_cases, 'boundary triples across Z2..Z8,F2^2,F2^3.')
    print('Full flow quotient K7 -> Heawood: dimensions15 -> 8, seven triangle circulations lost.')
    for group in (Group(2), Group(3), Group(4, True), Group(5), Group(8, True)):
        hostile = lift_to_paley((0,)*21, (1,)*7, group)
        need(all(hostile.values()) and not any(project_from_paley(hostile, group)),
             'nowhere-zero flow projects to zero')
    print('Hostile: value1 on every oriented Paley triangle is nowhere-zero, but maps to zero.')

    matching_components = Counter()
    for pi in permutations(range(7)):
        if not all((p-pi[p]) % 7 in D for p in range(7)):
            continue
        adjacency = [[] for _ in range(14)]
        for p, line in EDGES:
            if line != 7+pi[p]:
                adjacency[p].append(line)
                adjacency[line].append(p)
        seen, lengths = set(), []
        for v in range(14):
            if v in seen:
                continue
            stack, length = [v], 0
            seen.add(v)
            while stack:
                u = stack.pop()
                length += 1
                for w in adjacency[u]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            lengths.append(length)
        matching_components[tuple(sorted(lengths))] += 1
    need(matching_components == Counter({(14,): 24}), 'independent perfect-matching census')
    tait_count = sum(amount*2**len(cycles) for cycles, amount in matching_components.items())
    need(tait_count == 48, 'independent three-edge-color count')
    print('Independent q4 check:24 perfect matchings, each complementary to one14-cycle, giving48 Tait colorings.')

    embeddings = []
    for sign in (1, -1):
        faces, dual, boundaries = rotations(sign)
        need(all(a != b for a, b in dual), 'no dual loops')
        expected_faces = 7 if sign == 1 else 3
        need(len(faces) == expected_faces, 'face count')
        expected_dual = Counter({edge: 1 for edge in combinations(range(7), 2)}) if sign == 1 else Counter({edge: 7 for edge in combinations(range(3), 2)})
        need(dual == expected_dual, 'exact embedded dual')
        genus = (2-14+21-len(faces))//2
        for p in (2, 5):
            need(len(rref(boundaries, p)) == len(faces)-1, 'face-boundary rank')
        need(all(check_flow(tuple(v % 5 for v in row), Group(5)) for row in boundaries),
             'boundary of face is a flow')
        embeddings.append((genus, boundaries))
        print('Heawood embedding sign', sign, 'face lengths', [len(f) for f in faces],
              'genus', genus, 'dual', 'K7' if sign == 1 else '7*K3',
              'boundary rank', len(faces)-1, 'homology rank', 2*genus)

    need(check_flow((2,)*21, Group(6)), 'explicit nowhere-zero Z6 flow')
    integer_flow = tuple(-2 if (p-(line-7)) % 7 == 1 else 1 for p, line in EDGES)
    integer_balances = [0]*14
    for (u, v), value in zip(EDGES, integer_flow):
        integer_balances[u] += value
        integer_balances[v] -= value
    need(not any(integer_balances) and all(0 < abs(v) < 3 for v in integer_flow),
         'explicit integral 3-flow')
    torus_rows = embeddings[0][1]
    first_seven_color_flow = tuple(sum(i*row[e] for i, row in enumerate(torus_rows)) % 7
                                  for e in range(21))
    need(all(first_seven_color_flow) and check_flow(first_seven_color_flow, Group(7)),
         'first possible torus face-color threshold has an explicit flow')
    print('Threshold controls: offset1 matching carries -2, others+1: integral3-flow; constant2 is Z6 flow;')
    print('  seven distinct torus face labels give a Z7 boundary flow.')

    for group in (Group(3), Group(4, True), Group(5), Group(8, True)):
        normalized = group.q == 8
        multiplier = 42 if normalized else 1
        bases = [rref(rows, 2 if group.xor else group.q) for genus, rows in embeddings]
        count, boundary_counts = 0, [0, 0]
        sectors = [Counter(), Counter()]
        first_flow = None
        for flow in nowhere_zero_flows(group, normalized):
            count += 1
            need(check_flow(flow, group) and all(flow), 'independent full-edge conservation')
            if first_flow is None:
                first_flow = flow
            for i, basis in enumerate(bases):
                residue = reduce_mod_boundaries(flow, basis, group)
                boundary_counts[i] += not any(residue)
                if not normalized:
                    sectors[i][residue] += 1
        full_count = multiplier*count
        full_boundaries = [multiplier*x for x in boundary_counts]
        expected = [prod(group.q-i for i in range(1, 7)), (group.q-1)*(group.q-2)]
        need(full_boundaries == expected, ('dual-color count', group.q, full_boundaries, expected))
        print('Group', 'F2^r' if group.xor else 'Z'+str(group.q), 'order', group.q,
              'nowhere-zero flows', full_count, 'boundary counts genus1/genus3', full_boundaries)
        if normalized:
            print('  Enumerated all', 7**6, 'normalized cotree assignments; GL(3,2) orbit multiplier42.')
        else:
            print('  Enumerated all', (group.q-1)**8, 'nonzero cotree assignments.')
            print('  Homology sector occupancy histograms:', [dict(sorted(Counter(s.values()).items())) for s in sectors])
        print('  Triangle lift count over each Heawood nowhere-zero flow:', (group.q-3)**7,
              '; number of descending nowhere-zero K7 flows:', full_count*(group.q-3)**7)
        choices = []
        for t in range(7):
            x, y, z = [(t+d) % 7 for d in D]
            bx, by = [first_flow[EINDEX[(p, 7+t)]] for p in (x, y)]
            choices.append([a for a in range(group.q) if all((a, group.add(a, by), group.sub(a, bx)))])
        # Complete lifts of one witness:1,128,78125 choices respectively.
        witness_lifts = 0
        for parameters in product(*choices):
            lift = lift_to_paley(first_flow, parameters, group)
            need(all(lift.values()) and project_from_paley(lift, group) == first_flow,
                 'nowhere-zero lift and retained projection')
            witness_lifts += 1
        need(witness_lifts == (group.q-3)**7, 'complete witness fibre')
        print('  Complete independently checked lifts of one witness:', witness_lifts)

    print('Conditional cubic Euler bound: all face lengths>=7 gives V<=28(g-1).')
    print('Negative curvature alone is not a no-flow certificate: both embeddings have identical flows.')
    print('PASS: exact local fibres, two embeddings, complete flow counts and boundary obstructions.')


if __name__ == '__main__':
    main()
