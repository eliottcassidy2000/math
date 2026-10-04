"""Exact tetrahedral pair-midpoint charts and their minimal common refinement.

Fraction coordinates throughout; no floating geometry or inferred embedding.
Run normally and with Python -O. The note states all minimality hypotheses.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json


def need(ok,message):
    if not ok:
        raise ValueError(message)


def cell(points):
    return tuple(sorted(points))


E = tuple(tuple(F(int(i == j)) for i in range(4)) for j in range(4))
CENTER = (F(1,4),)*4
PAIRS = tuple(combinations(range(4),2))
MID = {pair:tuple((E[pair[0]][i]+E[pair[1]][i])/2 for i in range(4)) for pair in PAIRS}
MATCHINGS = (((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
CORNER_CELLS = {cell((E[i],)+(tuple(MID[p] for p in PAIRS if i in p))) for i in range(4)}
H = ((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
KLEIN_LABELS = ((0,0),(0,1),(1,0),(1,1))
WALSH_CHARACTERS = ((0,0),(1,0),(0,1),(1,1))


def hadamard(point):
    return tuple(sum(a*b for a,b in zip(row,point)) for row in H)


def inverse_hadamard(point):
    return tuple(F(a)/4 for a in hadamard(point))


def hadamard_integer_image(point):
    return len({a % 2 for a in point}) == 1 and sum(point) % 4 == 0


def average(points):
    return tuple(sum(p[i] for p in points)/len(points) for i in range(4))


def barycentric_refinement(mesh):
    parents = {}
    for tet in mesh:
        for p in permutations(tet):
            small = cell(average(p[:k]) for k in range(1,5))
            need(small not in parents, "unique barycentric tetrahedron interior")
            need(all(contains(tet,v) for v in small), "exact barycentric parent containment")
            parents[small] = tet
    return parents


def det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def relative_volume(tet):
    a,b,c,d = tet
    return abs(det3([[v[i]-a[i] for v in (b,c,d)] for i in range(3)]))


def barycentric_coordinates(tet,point):
    a,b,c,d = tet
    matrix = [[v[i]-a[i] for v in (b,c,d)] for i in range(3)]
    determinant = det3(matrix)
    need(determinant != 0, "nondegenerate containing tetrahedron")
    rhs = [point[i]-a[i] for i in range(3)]
    weights = []
    for j in range(3):
        changed = [row[:] for row in matrix]
        for i in range(3):
            changed[i][j] = rhs[i]
        weights.append(det3(changed)/determinant)
    return (1-sum(weights),)+tuple(weights)


def contains(tet,point):
    return min(barycentric_coordinates(tet,point)) >= 0


def signed_endpoint(axis,sign):
    return MID[MATCHINGS[axis][0 if sign == 1 else 1]]


def diagonal_chart(axis):
    others = [j for j in range(3) if j != axis]
    central = {signs:cell((MID[MATCHINGS[axis][0]],MID[MATCHINGS[axis][1]],
                          signed_endpoint(others[0],signs[0]),signed_endpoint(others[1],signs[1])))
               for signs in product((-1,1),repeat=2)}
    return CORNER_CELLS | set(central.values()),central


STAR_CELLS = {signs:cell((CENTER,)+tuple(signed_endpoint(i,signs[i]) for i in range(3)))
              for signs in product((-1,1),repeat=3)}
STAR_MESH = CORNER_CELLS | set(STAR_CELLS.values())


def closure(mesh):
    return {size:{sub for tet in mesh for sub in combinations(tet,size)} for size in range(1,5)}


def stats(mesh):
    faces = closure(mesh)
    facets = Counter(f for tet in mesh for f in combinations(tet,3))
    need(set(facets.values()) <= {1,2}, "manifold triangular face incidence")
    need(sum(relative_volume(tet) for tet in mesh) == 1, "total exact relative volume")
    f_vector = tuple(len(faces[i]) for i in range(1,5))
    need(sum((-1)**i*n for i,n in enumerate(f_vector)) == 1, "ball Euler characteristic")
    return dict(f_vector=f_vector,boundary_triangles=sum(n == 1 for n in facets.values()),
                relative_volume_histogram={str(v):n for v,n in sorted(Counter(relative_volume(tet) for tet in mesh).items())})


def graph_automorphisms(vertices,edges):
    vertices = tuple(sorted(vertices))
    lookup = {v:i for i,v in enumerate(vertices)}
    adj = [set() for _ in vertices]
    for u,v in edges:
        i,j = lookup[u],lookup[v]
        adj[i].add(j)
        adj[j].add(i)
    signatures = [(len(adj[i]),tuple(sorted(len(adj[j]) for j in adj[i]))) for i in range(len(vertices))]
    domains = [[j for j in range(len(vertices)) if signatures[j] == signatures[i]] for i in range(len(vertices))]
    answer,mapping,used = [],{},set()
    def visit():
        if len(mapping) == len(vertices):
            answer.append(tuple(mapping[i] for i in range(len(vertices))))
            return
        i = max((i for i in range(len(vertices)) if i not in mapping),
                key=lambda i:(sum(j in mapping for j in adj[i]),-len(domains[i])))
        for j in domains[i]:
            if j in used:
                continue
            if any((a in adj[i]) != (b in adj[j]) for a,b in mapping.items()):
                continue
            mapping[i] = j
            used.add(j)
            visit()
            used.remove(j)
            del mapping[i]
    visit()
    return vertices,answer


def mesh_automorphisms(mesh):
    faces = closure(mesh)
    vertices,maps = graph_automorphisms({v[0] for v in faces[1]},faces[2])
    for p in maps:
        need({cell(vertices[p[vertices.index(v)]] for v in tet) for tet in mesh} == mesh,
             "computed graph automorphism preserves all tetrahedral faces")
    return vertices,maps


def coordinate_action(point,p):
    return tuple(point[i] for i in p)


def coordinate_mesh_action(mesh,p):
    return {cell(coordinate_action(v,p) for v in tet) for tet in mesh}


def cumulative_edge(u,v):
    total = 0
    delta = []
    for a,b in zip(u,v):
        total += 2*(b-a)
        delta.append(total)
    return all(a in (0,1) for a in delta) or all(a in (0,-1) for a in delta)


def main():
    charts = [diagonal_chart(i) for i in range(3)]
    original_vertices = set(E)|set(MID.values())
    edgewise = {cell(tet) for tet in combinations(sorted(original_vertices),4)
                if all(cumulative_edge(u,v) for u,v in combinations(tet,2))}
    need(edgewise == charts[1][0], "independent cumulative rule selects the 02|13 matching")
    chart_keys = {frozenset(mesh) for mesh,_ in charts}
    coordinate_images = Counter(frozenset(coordinate_mesh_action(edgewise,p)) for p in permutations(range(4)))
    need(set(coordinate_images) == chart_keys and sorted(coordinate_images.values()) == [8,8,8],
         "the 24 coordinate orders give exactly three eight-cell charts")
    for matching in MATCHINGS:
        need(tuple((MID[matching[0]][i]+MID[matching[1]][i])/2 for i in range(4)) == CENTER,
             "all opposite-pair diagonals intersect at the center")
    need(tuple(tuple(sum(H[i][k]*H[k][j] for k in range(4)) for j in range(4)) for i in range(4)) ==
         tuple(tuple(4*int(i == j) for j in range(4)) for i in range(4)), "H squared is four times identity")
    need(H == tuple(tuple((-1)**(a*x+b*y) for x,y in KLEIN_LABELS) for a,b in WALSH_CHARACTERS),
         "H is the exact Klein-four character table in the declared labels")
    for i,j in product(range(4),repeat=2):
        need(tuple(H[i][k]*H[j][k] for k in range(4)) == H[i^j],
             "Walsh rows multiply by the Klein-four group law")
    hdet = sum((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
               *H[0][p[0]]*H[1][p[1]]*H[2][p[2]]*H[3][p[3]] for p in permutations(range(4)))
    need(hdet == -16, "Hadamard lattice index is sixteen, not unimodular")
    for point in product(range(-3,4),repeat=4):
        decoded = inverse_hadamard(point)
        need(hadamard_integer_image(point) == all(x.denominator == 1 for x in decoded),
             "parity plus mod-four condition exactly characterizes integer image")
        need(inverse_hadamard(hadamard(point)) == point, "cut coordinate round trip")
    need(not hadamard_integer_image((2,0,0,0)) and inverse_hadamard((2,0,0,0)) == (F(1,2),)*4,
         "equal parity alone misses degree-two center obstruction")
    need(inverse_hadamard((4,0,0,0)) == (1,1,1,1), "center is an integer degree-four address")
    for axis in range(3):
        for sign in (-1,1):
            need(hadamard(signed_endpoint(axis,sign)) == (1,)+tuple(sign*int(i == axis) for i in range(3)),
                 "three pair matchings are the three cut-coordinate axes")
    need({hadamard(v)[1:] for v in E} == {signs for signs in product((-1,1),repeat=3) if signs[0]*signs[1]*signs[2] == 1},
         "original corners occupy the positive-product sign class")
    chart_stats = [stats(mesh) for mesh,_ in charts]
    need(all(row['f_vector'] == (10,25,24,8) and row['relative_volume_histogram'] == {'1/8':8}
             for row in chart_stats), "edgewise chart counts and equal volumes")
    star_stats = stats(STAR_MESH)
    need(star_stats['f_vector'] == (11,30,32,12), "center repair f-vector")
    need(star_stats['relative_volume_histogram'] == {'1/16':8,'1/8':4}, "symmetry repair has unequal volumes")
    boundaries = [{f for f,n in Counter(f for tet in mesh for f in combinations(tet,3)).items() if n == 1}
                  for mesh in [*(mesh for mesh,_ in charts),STAR_MESH]]
    need(all(b == boundaries[0] for b in boundaries), "all charts and repair have exactly the same boundary")

    parents,central_codes = {},[]
    for fine in STAR_MESH:
        selected = []
        for mesh,_ in charts:
            containing = [tet for tet in mesh if all(contains(tet,p) for p in fine)]
            need(len(containing) == 1, "one full-dimensional parent per chart")
            selected.append(containing[0])
        parents[fine] = tuple(selected)
    for i,j in combinations(range(3),2):
        need(len({(p[i],p[j]) for p in parents.values()}) == 12,
             "any two charts already distinguish all twelve common cells")
    for signs,fine in STAR_CELLS.items():
        codes = []
        for axis,(_,central) in enumerate(charts):
            code = tuple(signs[i] for i in range(3) if i != axis)
            need(parents[fine][axis] == central[code], "geometric parent map is sign-coordinate deletion")
            codes.append(code)
        central_codes.append(dict(octant=signs,three_chart_codes=codes))
    for axis in range(3):
        need(set(Counter(row['three_chart_codes'][axis] for row in central_codes).values()) == {2},
             "each chart forgets exactly one binary sign")

    # Independent finite classification with the prescribed boundary edges.
    boundary_edges = {edge for face in boundaries[0] for edge in combinations(face,2)}
    diagonal_edges = {cell(MID[p] for p in matching) for matching in MATCHINGS}
    allowed_edges = boundary_edges|diagonal_edges
    candidates = {cell(tet) for tet in combinations(sorted(original_vertices),4)
                  if all(cell(edge) in allowed_edges for edge in combinations(tet,2)) and relative_volume(tet)>0}
    need(len(boundary_edges) == 24 and len(candidates) == 16, "all allowable nondegenerate tetrahedra")
    need(candidates == set().union(*(mesh for mesh,_ in charts)), "candidate inventory is four corners plus twelve diagonal cells")
    for i in range(4):
        incident = {tet for tet in candidates if E[i] in tet}
        need(len(incident) == 1 and incident <= CORNER_CELLS, "each old corner forces its unique corner tetrahedron")
    central_candidates = candidates-CORNER_CELLS
    for tet in central_candidates:
        need(sum(edge in set(combinations(tet,2)) for edge in diagonal_edges) == 1,
             "each central nondegenerate cell uses exactly one opposite-pair diagonal")

    octa_vertices = set(MID.values())
    octa_edges = {cell((MID[p],MID[q])) for p,q in combinations(PAIRS,2) if set(p)&set(q)}
    _,octa_maps = graph_automorphisms(octa_vertices,octa_edges)
    _,chosen_octa_maps = graph_automorphisms(octa_vertices,octa_edges|{next(iter(diagonal_edges))})
    _,star_octa_maps = graph_automorphisms(octa_vertices|{CENTER},octa_edges|{cell((CENTER,v)) for v in octa_vertices})
    _,full_chart_maps = mesh_automorphisms(charts[0][0])
    _,full_star_maps = mesh_automorphisms(STAR_MESH)
    group_counts = dict(bare_midpoint_octahedron=len(octa_maps),octahedron_with_one_diagonal=len(chosen_octa_maps),
        bare_center_star=len(star_octa_maps),full_edgewise_chart=len(full_chart_maps),full_center_repair=len(full_star_maps))
    need(group_counts == dict(bare_midpoint_octahedron=48,octahedron_with_one_diagonal=16,
        bare_center_star=48,full_edgewise_chart=8,full_center_repair=24), "independent automorphism census")
    coordinate_permutations = tuple(permutations(range(4)))
    need(all(coordinate_mesh_action(STAR_MESH,p) == STAR_MESH for p in coordinate_permutations),
         "all 24 original coordinate symmetries preserve the repair")
    orbits,remaining = [],set(STAR_MESH)
    while remaining:
        representative = min(remaining)
        orbit = {cell(coordinate_action(v,p) for v in representative) for p in coordinate_permutations}
        need(orbit <= remaining, "cell orbits form a partition")
        orbits.append(orbit)
        remaining -= orbit
    need(sorted(map(len,orbits)) == [4,4,4], "restored S4 has three tetrahedron orbits")
    cut_action_determinants = Counter()
    for p in coordinate_permutations:
        columns = [hadamard(coordinate_action(signed_endpoint(i,1),p))[1:] for i in range(3)]
        matrix = [[columns[j][i] for j in range(3)] for i in range(3)]
        signs = [next(x for x in col if x) for col in columns]
        need(signs[0]*signs[1]*signs[2] == 1, "tetrahedral action preserves product sign")
        cut_action_determinants[int(det3(matrix))] += 1
    need(cut_action_determinants == {-1:12,1:12}, "tetrahedral S4 is not merely the rotational octahedral subgroup")

    flag_parents = barycentric_refinement(STAR_MESH)
    flag_mesh = set(flag_parents)
    flag_stats = stats(flag_mesh)
    need(flag_stats['f_vector'] == (85,420,624,288), "center repair flag-completion f-vector")
    base_flag_mesh = {p:cell(average(tuple(E[i] for i in p[:k])) for k in range(1,5))
                      for p in coordinate_permutations}
    for fine in flag_mesh:
        centroid = average(fine)
        order = tuple(sorted(range(4),key=lambda i:(-centroid[i],i)))
        need(all(v[order[i]] >= v[order[i+1]] for v in fine for i in range(3)),
             "each flag cell lies within one original coordinate-order chamber")
        need(all(contains(base_flag_mesh[order],v) for v in fine),
             "independent determinant check of original barycentric parent")
        for chart_parent in parents[flag_parents[fine]]:
            need(all(contains(chart_parent,v) for v in fine), "flag completion retains all three edgewise parents")
    report = dict(status="PROVED with exact finite controls; minimality is for common refinement of the stated embedded charts",
        pair_midpoint_labels=[dict(pair=p,point=list(map(str,MID[p]))) for p in PAIRS],
        diagonal_matchings=MATCHINGS,coordinate_orders_per_chart=sorted(coordinate_images.values()),
        hadamard_cut_coordinates=dict(matrix=H,determinant=hdet,lattice_index=abs(hdet),
            inverse="H/4",integer_image="all four coordinates have the same parity and their sum is 0 mod 4",
            klein_four_vertex_labels=KLEIN_LABELS,walsh_character_labels=WALSH_CHARACTERS,
            degree_two_center_cut=(2,0,0,0),degree_four_center_cut=(4,0,0,0),
            original_corner_cuts=[list(map(int,hadamard(v))) for v in E],cut_action_determinants=dict(sorted(cut_action_determinants.items()))),
        chart_counts=chart_stats,center=list(map(str,CENTER)),center_repair=star_stats,
        original_vertex_only_candidate_tetrahedra=len(candidates),
        original_vertex_only_central_candidate_tetrahedra=len(central_candidates),
        common_refinement_cells=len(parents),central_joint_codes=central_codes,
        flag_completion=flag_stats,flag_completion_parent_checks=len(flag_mesh),
        automorphism_counts=group_counts,full_center_repair_tetrahedron_orbit_sizes=sorted(map(len,orbits)),
        scope="Standard midpoint boundary fixed for ten-vertex classification; any two distinct charts force the unique twelve-cell minimal common refinement.",
        loss="A single diagonal chart forgets one octant sign; the bare octahedron forgets which triangular facets came from old vertices versus old faces.")
    print(json.dumps(report,indent=2))
    print("PASS: all checks remain active under -O")


if __name__ == '__main__':
    main()
