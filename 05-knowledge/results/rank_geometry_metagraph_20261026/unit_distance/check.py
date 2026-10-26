"""Exact boundary tables and planar metric gates; standard library only."""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations, product
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rank(rows):
    a = [[Q(x) for x in row] for row in rows]
    r = 0
    for j in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        divisor = a[r][j]
        a[r] = [v/divisor for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                multiplier = a[i][j]
                a[i] = [x-multiplier*y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def norm(a):
    # a0+a1*alpha, alpha=(1+i*sqrt(3))/2.
    x, y = a
    return x*x+x*y+y*y


def dist(a, b):
    return norm((a[0]-b[0], a[1]-b[1]))


def dot(a, b):
    return Q(norm((a[0]+b[0], a[1]+b[1]))-norm(a)-norm(b), 2)


def frame(perm, shift, x):
    y = shift
    for j in range(3):
        y ^= ((x >> j) & 1) << perm[j]
    return y


def palette(perm, shift):
    return tuple(2*perm[j] + (b ^ ((shift >> perm[j]) & 1))
                 for j in range(3) for b in range(2))


def isometric(points, p):
    return all(dist(points[i], points[j]) == dist(points[p[i]], points[p[j]])
               for i in range(len(points)) for j in range(i))


def gram(points):
    differences = [(points[2*j+1][0]-points[2*j][0],
                    points[2*j+1][1]-points[2*j][1]) for j in range(3)]
    return [[dot(a, b) for b in differences] for a in differences]


def gram_equations(perm, shift):
    variables = [(i, j) for i in range(3) for j in range(i, 3)]
    equations = []
    for i, j in variables:
        row = [0]*6
        row[variables.index((i, j))] += 1
        target = tuple(sorted((perm[i], perm[j])))
        sign = (-1)**(((shift >> perm[i]) & 1) + ((shift >> perm[j]) & 1))
        row[variables.index(target)] -= sign
        equations.append(row)
    return equations


def matrix_product(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for j in range(1,len(vertices)):
        for rest in matchings(vertices[1:j]+vertices[j+1:]):
            yield ((first,vertices[j]),)+rest


def pair_transport(matching, motion):
    perm, shift = [], 0
    for u,v in matching:
        target = next((j for j,pair in enumerate(matching)
                       if set(pair) == {motion[u],motion[v]}), None)
        if target is None:
            return None
        perm.append(target)
        if motion[u] == matching[target][1]:
            shift |= 1 << target
    return tuple(perm),shift


def boundary_counts(n, edges, boundary, q):
    """Enumerate all proper labelled colorings, retaining complete port tuples."""
    previous = [[] for _ in range(n)]
    for a, b in edges:
        a, b = sorted((a, b))
        previous[b].append(a)
    colors = [-1]*n
    counts = Counter()
    def visit(v):
        if v == n:
            counts[tuple(colors[j] for j in boundary)] += 1
            return
        forbidden = {colors[j] for j in previous[v]}
        for color in range(q):
            if color not in forbidden:
                colors[v] = color
                visit(v+1)
    visit(0)
    return counts


def diamond_strip(length):
    # (a,b) represents (sqrt(3)*a/2,b/2), with squared norm (3a^2+b^2)/4.
    points = [(2*j, 0) for j in range(length+1)]
    points += [(2*j+1, b) for j in range(length) for b in (-1, 1)]
    edges = set()
    for i in range(len(points)):
        for j in range(i):
            a, b = points[i][0]-points[j][0], points[i][1]-points[j][1]
            if 3*a*a+b*b == 4:
                edges.add((j, i))
    intended = set()
    for j in range(length):
        u, v = length+1+2*j, length+2+2*j
        intended.update((min(a,b), max(a,b)) for a,b in
                        ((u,v), (j,u), (j,v), (j+1,u), (j+1,v)))
    need(edges == intended and len(edges) == 5*length, 'complete induced strip edges')
    need(len(points) == len(set(points)), 'injective strip')
    return points, edges


def main():
    # The inherited two-loop frame maps.
    pi = (1, 2, 0)
    g1 = tuple(frame(pi, 0, x) for x in range(8))
    g2 = tuple(frame(pi, 5, x) for x in range(8))
    fixed1 = [x for x in range(8) if g1[x] == x]
    fixed2 = [x for x in range(8) if g2[x] == x]
    need(fixed1 == [0,7] and fixed2 == [3,4], 'literal F10 fixed sets')
    even = [0,3,5,6]
    maps = []
    for perm in permutations(range(3)):
        for shift in range(8):
            if shift.bit_count() % 2 == 0:
                maps.append(tuple(frame(perm, shift, x) for x in even))
    need(len(set(maps)) == 24, 'full residual S4')
    table_pairs = 0
    for f in maps:
        for g in maps:
            a = [[int(f[i] == even[j]) for j in range(4)] for i in range(4)]
            b = [[int(g[i] == even[j]) for j in range(4)] for i in range(4)]
            composed = matrix_product(a, b)
            need(composed == [[int(g[even.index(f[i])] == even[j]) for j in range(4)]
                              for i in range(4)], 'full relation composition')
            need(sum(a[i][i]*b[i][i] for i in range(4)) ==
                 sum(f[i] == even[i] == g[i] for i in range(4)), 'common-frame gluing')
            table_pairs += 1

    # Exact planar metrics for each loop separately, but not jointly.
    directions = [(Q(1),Q(0)), (Q(-1),Q(1)), (Q(0),Q(-1))]
    hexagon = [p for d in directions for p in (d,(-d[0],-d[1]))]
    p1, p2 = palette(pi,0), palette(pi,5)
    relabel = palette((0,1,2),4)
    second = [hexagon[relabel[i]] for i in range(6)]
    need(isometric(hexagon,p1) and not isometric(hexagon,p2), 'first local metric')
    need(isometric(second,p2) and not isometric(second,p1), 'second local metric')
    need(rank(gram(hexagon)) == rank(gram(second)) == 2, 'two planar local Grams')
    eq1, eq2 = gram_equations(pi,0), gram_equations(pi,5)
    need(rank(eq1) == rank(eq2) == 4 and rank(eq1+eq2) == 5, 'joint Gram equations')
    identity = [1,0,0,1,0,1]
    need(all(sum(x*y for x,y in zip(row,identity)) == 0 for row in eq1+eq2),
         'only scalar identity Gram remains')
    # Unit octahedron realizes the same two permutations in dimension three.
    octahedron = [tuple((1 if b == 0 else -1) if j == k else 0 for k in range(3))
                  for j in range(3) for b in range(2)]
    for p in (p1,p2):
        need(all(sum((a-b)**2 for a,b in zip(octahedron[i],octahedron[j])) ==
                 sum((a-b)**2 for a,b in zip(octahedron[p[i]],octahedron[p[j]]))
                 for i in range(6) for j in range(i)), 'sharp three-dimensional model')
    geometric = []
    for perm in permutations(range(3)):
        for shift in range(8):
            if isometric(hexagon,palette(perm,shift)):
                geometric.append((perm,shift))
                e0 = frame(perm,shift,0)
                need((e0 ^ (7 if e0.bit_count()%2 else 0)) == 0,
                     'hexagon residual has common frame')
    need(len(geometric) == 12, 'exact hexagon palette isometry group')

    # Re-pairing the very same metric creates a genuine fixed-free even transport.
    clock = [(Q(1),Q(0)),(Q(0),Q(1)),(Q(-1),Q(1)),
             (Q(-1),Q(0)),(Q(0),Q(-1)),(Q(1),Q(-1))]
    dihedral = [tuple((sign*i+offset)%6 for i in range(6))
                for sign in (1,-1) for offset in range(6)]
    need(all(isometric(clock,motion) for motion in dihedral), 'literal hexagon motions')
    pairing_census = Counter()
    for matching in matchings(tuple(range(6))):
        actions = [pair_transport(matching,motion) for motion in dihedral]
        actions = [a for a in actions if a is not None]
        parity_closed = [(p,a) for p,a in actions if a.bit_count()%2 == 0]
        fixed = [x for x in even if all(frame(p,a,x)==x for p,a in parity_closed)]
        length_type = tuple(sorted(int(dist(clock[a],clock[b])) for a,b in matching))
        pairing_census[(length_type,len(actions),len(parity_closed),len(fixed))] += 1
    need(sum(pairing_census.values()) == 15, 'all matchings retained')
    mirror_pairs = ((0,3),(1,5),(2,4))
    reflection = tuple((-i)%6 for i in range(6))
    mirror_action = pair_transport(mirror_pairs,reflection)
    need(mirror_action == ((0,1,2),6), 'actual even double-flip reflection')
    need(all(frame(*mirror_action,x) != x for x in range(8)), 'fixed-free geometric transport')
    need(sum(count for (_,_,_,fixed),count in pairing_census.items() if fixed==0)==3,
         'three geometric pairings have no residual common frame')

    # Complete unprescribed boundary relation: none of this removes other colorings.
    rim_edges = {(min(i,(i+1)%6),max(i,(i+1)%6)) for i in range(6)}
    wheel_points = clock+[(Q(0),Q(0))]
    wheel_edges = {(j,i) for i in range(7) for j in range(i)
                   if dist(wheel_points[i],wheel_points[j]) == 1}
    need(wheel_edges == rim_edges | {(i,6) for i in range(6)}, 'induced wheel metric')
    wheel_table = boundary_counts(7,wheel_edges,tuple(range(6)),6)
    for colors in product(range(6),repeat=6):
        proper = all(colors[a] != colors[b] for a,b in rim_edges)
        need(wheel_table[colors] == int(proper)*(6-len(set(colors))), 'full wheel table')
    need(len(wheel_table) == 14910 and sum(wheel_table.values()) == 24600,
         'wheel full relation census')
    for x in range(8):
        colors = [-1]*6
        for j,(u,v) in enumerate(mirror_pairs):
            colors[u],colors[v] = 2*j+((x>>j)&1),2*j+(1-((x>>j)&1))
        need(wheel_table[tuple(colors)] == 0, 'rainbow frame chart is empty on wheel')
    need(all(([0,1,0,1,0,1,2][a] != [0,1,0,1,0,1,2][b]) for a,b in wheel_edges),
         'actual three-coloring outside rainbow frame chart')

    # Every edge of the explicit infinite strip family is checked through length 20.
    strip_colorings = 0
    for length in range(1,21):
        points, edges = diamond_strip(length)
        if length <= 2:
            for q in range(3,7):
                counts = boundary_counts(len(points),edges,(0,length),q)
                d = [[(q-1)*(q-2) if i==j else (q-2)*(q-3)
                      for j in range(q)] for i in range(q)]
                expected = d if length == 1 else matrix_product(d,d)
                need(all(counts[i,j] == expected[i][j] for i in range(q) for j in range(q)),
                     'complete strip boundary counts')
                strip_colorings += sum(counts.values())
    six = [[20 if i==j else 12 for j in range(6)] for i in range(6)]
    power = [[int(i==j) for j in range(6)] for i in range(6)]
    for length in range(1,21):
        power = matrix_product(power,six)
        need(all(power[i][j] == 8**length*int(i==j)+(80**length-8**length)//6
                 for i in range(6) for j in range(6)), 'closed six-color transfer power')
        need(all(value > 0 for row in power for value in row), 'full six-color support')

    # An actual unit-star has all label permutations as graph symmetries,
    # but generally not as symmetries of its complete Euclidean metric.
    rational = [(Q(1-j*j,1+j*j),Q(2*j,1+j*j)) for j in range(1,7)]
    def euclid(a,b): return sum((x-y)**2 for x,y in zip(a,b))
    need(all(euclid(p,(0,0)) == 1 for p in rational), 'unit star spokes')
    need(all(euclid(rational[i],rational[j]) not in (0,1)
             for i in range(6) for j in range(i)), 'induced rational star')
    star_hist = Counter()
    total = 0
    for colors in product(range(6),repeat=6):
        extensions = sum(c not in colors for c in range(6))
        need(extensions == 6-len(set(colors)), 'full star relation')
        star_hist[extensions] += 1
        total += extensions
    need(star_hist[0] == 720 and total == 6*5**6, 'prescribed versus unprescribed star')
    need(all(6-len({2*j+((x>>j)&1) for j in range(3)}) == 3 for x in range(8)),
         'all frame choices extend across three-port star')

    # Two independently legal unit edges create an undeclared cross edge at 60 degrees.
    triangle = [(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))]
    need(all(dist(triangle[i],triangle[j]) == 1 for i in range(3) for j in range(i)),
         'unintended cross-unit contact')
    path_counts = boundary_counts(3,{(0,1),(0,2)},(0,),6)
    triangle_counts = boundary_counts(3,{(0,1),(0,2),(1,2)},(0,),6)
    need(all(path_counts[c,] == 25 and triangle_counts[c,] == 20 for c in range(6)),
         'missing cross factor changes exact count')
    print(json.dumps({
        'status':'PASS', 'abstract_relation_pairs':table_pairs,
        'two_loop_fixed_frames':[fixed1,fixed2],
        'individual_Gram_constraint_ranks':[rank(eq1),rank(eq2)],
        'joint_Gram_constraint_rank':rank(eq1+eq2),
        'hexagon_geometric_palette_maps':len(geometric),
        'all15_pairings_census':{str(key):value for key,value in sorted(pairing_census.items())},
        'mirror_pairing_even_action':mirror_action,
        'wheel_boundary_support':len(wheel_table),
        'wheel_total_six_colorings':sum(wheel_table.values()),
        'strip_exact_geometry_lengths':20, 'strip_enumerated_colorings':strip_colorings,
        'six_color_single_diamond':six,
        'six_color_two_diamonds_diagonal_offdiagonal':[1120,1056],
        'unit_star_weight_histogram':dict(sorted(star_hist.items())),
        'unit_star_total_colorings':total,
        'cross_unit_contact_fixed_center_counts':[25,20],
        'scope':'exact finite controls; all-order proofs and geometry limitations in RESULT.md'
    },sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
