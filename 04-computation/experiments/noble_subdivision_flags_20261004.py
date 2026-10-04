"""Exact controls for flag versus triangular-lattice refinement.

Python 3 + networkx. Two pinned Hill OFF models are downloaded only if absent
from the existing temporary model cache; their Git blob hashes are verified.
No floating-point geometry or assertions are used. See the companion note.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from hashlib import sha1
from itertools import combinations, permutations, product
from pathlib import Path
import tempfile
import urllib.parse
import urllib.request

import networkx as nx


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def edges_of(faces):
    return {frozenset((f[i], f[(i + 1) % len(f)]))
            for f in faces for i in range(len(f))}


def counts(faces):
    return len(set().union(*map(set, faces))), len(edges_of(faces)), len(faces)


def flags_of(faces):
    """Literal incidences (v,e,f), independent of the older corner-offset code."""
    incidences = [(v, e, j) for j, face in enumerate(faces)
                  for e in edges_of([face]) for v in e]
    flags = sorted(incidences, key=repr)
    buckets = [defaultdict(list) for _ in range(3)]
    for i, flag in enumerate(flags):
        for rank in range(3):
            buckets[rank][tuple(x for k, x in enumerate(flag) if k != rank)].append(i)
    transitions = []
    for rank in range(3):
        row = [None] * len(flags)
        for pair in buckets[rank].values():
            check(len(pair) == 2, 'closed diamond condition')
            a, b = pair
            row[a], row[b] = b, a
        transitions.append(row)
    return flags, transitions


def flag_automorphisms(faces):
    flags, transitions = flags_of(faces)
    automorphisms = []
    for target in range(len(flags)):
        image, work, valid = {0: target}, deque([0]), True
        while work and valid:
            source = work.popleft()
            for row in transitions:
                a, b = row[source], row[image[source]]
                if a in image:
                    if image[a] != b:
                        valid = False
                        break
                else:
                    image[a] = b
                    work.append(a)
        if valid and len(image) == len(flags) and len(set(image.values())) == len(flags):
            automorphisms.append(image)
    return flags, automorphisms


def barycentric(faces):
    return [((0, v), (1, edge), (2, j)) for j, face in enumerate(faces)
            for edge in edges_of([face]) for v in edge]


def graph_of(faces):
    graph = nx.Graph()
    graph.add_edges_from(tuple(e) for e in edges_of(faces))
    return graph


def graph_automorphisms(faces, typed=False):
    graph = graph_of(faces)
    if typed:
        nx.set_node_attributes(graph, {v: v[0] for v in graph}, 'rank')
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph, node_match=(lambda a, b: a['rank'] == b['rank']) if typed else None)
    return list(matcher.isomorphisms_iter())


def orbit_sizes(objects, action_maps):
    objects = set(objects)
    sizes = []
    while objects:
        representative = next(iter(objects))
        orbit = {mapping(representative) for mapping in action_maps}
        check(orbit <= objects, 'disjoint group orbits')
        sizes.append(len(orbit))
        objects -= orbit
    return sorted(sizes)


def triples(total):
    return [(a, b, total-a-b) for a in range(total+1) for b in range(total-a+1)]


def edgewise(r):
    check(isinstance(r, int) and r >= 1, 'positive integer subdivision degree')
    faces, kinds = [], []
    for b in triples(r-1):
        faces.append(tuple(tuple(b[j] + (i == j) for j in range(3)) for i in range(3)))
        kinds.append('up')
    if r >= 2:
        for b in triples(r-2):
            faces.append(tuple(tuple(b[j] + (i != j) for j in range(3)) for i in range(3)))
            kinds.append('down')
    return faces, kinds


def boundary_count(faces):
    multiplicities = Counter(e for f in faces for e in edges_of([f]))
    check(set(multiplicities.values()) <= {1, 2}, 'surface edge incidence')
    return sum(m == 1 for m in multiplicities.values())


def edgewise_surface(faces, r):
    refined = []
    local, _ = edgewise(r)
    for face in faces:
        check(len(face) == 3, 'edgewise surface must be triangular')
        def point(weights):
            return tuple(sorted((face[i], weights[i]) for i in range(3) if weights[i]))
        refined.extend(tuple(point(w) for w in cell) for cell in local)
    return refined


def compose_edgewise(r, s):
    outer, _ = edgewise(r)
    inner, _ = edgewise(s)
    return {frozenset(tuple(sum(weights[i]*cell[i][j] for i in range(3))
                            for j in range(3)) for weights in micro)
            for cell in outer for micro in inner}


def hill_faces(name, expected_hash):
    commit = 'a801da7582fa927ba07e0af74d7c9c39445b7264'
    path = Path(tempfile.gettempdir()) / ('noble-tools-revised-' + commit[:7]) / (name+'.off')
    if not path.exists():
        url = ('https://raw.githubusercontent.com/Plasmath/noble-tools-revised/' + commit + '/'
               + urllib.parse.quote('library/0 Degrees of Freedom/D/'+name+'.off'))
        raw = urllib.request.urlopen(url, timeout=30).read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    raw = path.read_bytes()
    digest = sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    check(digest == expected_hash, name+' pinned Git blob')
    lines = [line.strip() for line in raw.decode().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    check(lines.pop(0) == 'OFF', 'OFF header')
    nv, nf, _ = map(int, lines.pop(0).split())
    faces = []
    for line in lines[nv:nv+nf]:
        row = list(map(int, line.split()))
        check(row[0] >= 3 and len(row[1:]) == row[0], 'polygon face')
        faces.append(tuple(row[1:]))
    return faces


def main():
    tetra = list(combinations(range(4), 3))
    cube_vertices = list(product((-1, 1), repeat=3))
    cube = []
    for axis in range(3):
        for sign in (-1, 1):
            varying = [a for a in range(3) if a != axis]
            face = []
            for u, v in [(-1,-1),(-1,1),(1,1),(1,-1)]:
                point = [0]*3
                point[axis], point[varying[0]], point[varying[1]] = sign, u, v
                face.append(cube_vertices.index(tuple(point)))
            cube.append(tuple(face))
    octa = [(a,b,c) for a in (0,1) for b in (2,3) for c in (4,5)]
    maps = [('tetrahedron', tetra, 24), ('cube', cube, 48), ('octahedron', octa, 48)]
    print('EXACT FLAG / EDGEWISE REFINEMENT CONTROLS')
    for name, faces, expected_aut in maps:
        v,e,f = counts(faces)
        flags, automorphisms = flag_automorphisms(faces)
        check(len(automorphisms) == expected_aut, name+' flag automorphisms')
        refined = barycentric(faces)
        check(counts(refined) == (v+e+f, 6*e, 4*e), name+' barycentric count')
        degrees = Counter(dict(graph_of(refined).degree()).values())
        check(degrees[4] == e and len(degrees) >= 2, name+' no vertex transitivity')
        typed = graph_automorphisms(refined, True)
        untyped = graph_automorphisms(refined)
        check(len(typed) == expected_aut, name+' reconstruction')
        check(len(untyped) == expected_aut*(2 if name == 'tetrahedron' else 1), name+' duality')
        print(name, 'base', (v,e,f), 'barycentric', counts(refined),
              'degree census', dict(sorted(degrees.items())),
              'typed/untyped Aut', (len(typed),len(untyped)))
    # Exact regular tetrahedron metric filters out abstract dualities.
    regular = [(Q(1),Q(1),Q(1)), (Q(1),Q(-1),Q(-1)),
               (Q(-1),Q(1),Q(-1)), (Q(-1),Q(-1),Q(1))]
    refined = barycentric(tetra)
    coordinates = {}
    for rank, cell in graph_of(refined):
        ids = [cell] if rank == 0 else (list(cell) if rank == 1 else tetra[cell])
        coordinates[(rank,cell)] = tuple(sum(regular[i][j] for i in ids)/len(ids) for j in range(3))
    def distance(a,b):
        return sum((coordinates[a][i]-coordinates[b][i])**2 for i in range(3))
    pairs = list(combinations(coordinates,2))
    geometric = [g for g in graph_automorphisms(refined)
                 if all(distance(a,b) == distance(g[a],g[b]) for a,b in pairs)]
    check(len(geometric) == 24, 'metric excludes barycentric abstract dualities')
    print('regular tetrahedron barycentric geometric Aut', len(geometric))

    # Pinned data controls recover the nonregular noble examples.
    models = [('D-4','eac135ab499faf4437ca05e7f9d05be2d954e5dc',(20,120,60),240),
              ('D-5','5978dbe15987115b8264c871d4203cfa3e5c8cd1',(20,90,60),120)]
    for name, digest, expected_counts, expected_aut in models:
        faces = hill_faces(name,digest)
        flags, automorphisms = flag_automorphisms(faces)
        check(counts(faces) == expected_counts and len(automorphisms) == expected_aut, name+' model census')
        refined = barycentric(faces)
        v,e,f = counts(faces)
        check(counts(refined) == (v+e+f,6*e,4*e), name+' refinement census')
        rank_orbits = []
        for rank in range(3):
            induced = []
            for automorphism in automorphisms:
                mapping = {flags[i][rank]: flags[j][rank] for i,j in automorphism.items()}
                check(all(mapping[flags[i][rank]] == flags[j][rank]
                          for i,j in automorphism.items()), 'well-defined cell action')
                induced.append(lambda cell,mapping=mapping: mapping[cell])
            rank_orbits.append(orbit_sizes((flag[rank] for flag in flags),induced))
        check(len(rank_orbits[0]) == len(rank_orbits[2]) == 1, name+' abstract V/F transitivity')
        print(name, 'base', counts(faces), 'flags',len(flags),'abstract Aut',len(automorphisms),
              'inherited abstract face orbits',len(flags)//len(automorphisms),
              'barycentric',counts(refined))
        print(name, 'base V/E/F orbit sizes',rank_orbits,
              'barycentric degree census',dict(sorted(Counter(dict(graph_of(refined).degree()).values()).items())))

    permutations3 = list(permutations(range(3)))
    local_orbit_census = []
    for r in range(1,13):
        faces,kinds = edgewise(r)
        check(counts(faces) == ((r+1)*(r+2)//2,3*r*(r+1)//2,r*r), 'edgewise local counts')
        check(boundary_count(faces) == 3*r, 'edgewise boundary count')
        # Independently reconstruct the graph from root differences.
        vertices = triples(r)
        root_edges = {frozenset((a,b)) for a,b in combinations(vertices,2)
                      if sorted(x-y for x,y in zip(a,b)) == [-1,0,1]}
        check(edges_of(faces) == root_edges, 'independent A2 adjacency')
        actions = [(lambda cell,p=p: frozenset(tuple(w[j] for j in p) for w in cell))
                   for p in permutations3]
        orbits = orbit_sizes(map(frozenset,faces),actions)
        partition_count = sum(1 for total in (r-1,r-2) if total >= 0
                              for a,b,c in triples(total) if a <= b <= c)
        check(len(orbits) == partition_count, 'partition orbit formula')
        local_orbit_census.append(len(orbits))
        for name, coarse, _ in maps:
            if any(len(face) != 3 for face in coarse):
                continue
            refined = edgewise_surface(coarse,r)
            v,e,f = counts(coarse)
            check(counts(refined) == (v+(r-1)*e+(r-1)*(r-2)*f//2,r*r*e,r*r*f), 'global edgewise counts')
            degrees = Counter(dict(graph_of(refined).degree()).values())
            if r >= 2:
                old_degree = 3 if name == 'tetrahedron' else 4
                check(degrees == {old_degree:v,6:counts(refined)[0]-v}, 'edgewise old/new degree obstruction')
    print('edgewise r=1..12 local face-orbit counts',local_orbit_census)
    for r in range(1,7):
        for s in range(1,7):
            direct,_ = edgewise(r*s)
            check(compose_edgewise(r,s) == set(map(frozenset,direct)), 'composition r*s')
    print('edgewise composition: all 36 pairs r,s=1..6 passed')
    faces,_ = edgewise(2)
    actions = [(lambda cell,p=p: frozenset(tuple(w[j] for j in p) for w in cell)) for p in permutations3]
    check(orbit_sizes(map(frozenset,faces),actions) == [1,3], 'central/corner hostile')
    local = [(0,1,2)]
    for _ in range(2):
        local = barycentric(local)
    direct,_ = edgewise(6)
    check(counts(local) == (25,60,36) and boundary_count(local) == 12, 'sd squared control')
    check(counts(direct) == (28,63,36) and boundary_count(direct) == 18, 'degree6 control')
    print('equal face-count hostile: sd^2',counts(local),'boundary',boundary_count(local),
          'versus edgewise6',counts(direct),'boundary',boundary_count(direct))
    midpoint, centroid = (3,3,0),(2,2,2)
    difference = tuple(a-b for a,b in zip(centroid,midpoint))
    check(sorted(difference) != [-1,0,1] and all(difference), 'median is not root-parallel')
    print('r6 shared-vertex hostile: midpoint-to-centroid direction',difference)
    print('Universes: 3 regular maps; 2 pinned noble maps; local r=1..12; tetra/octa r=1..12; 36 compositions.')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
