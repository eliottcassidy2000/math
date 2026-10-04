"""Exact cumulative-coordinate edgewise subdivision and hostile controls.

Only standard-library exact arithmetic is used. Run normally and with -O.
The accompanying note separates classical subdivision facts from proofs
and finite checks below. No abstract graph embedding is assumed.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


@lru_cache(None)
def compositions(r,m):
    if m == 1:
        return ((r,),)
    return tuple((a,)+tail for a in range(r+1) for tail in compositions(r-a,m-1))


def cumulative(v):
    out,total = [],0
    for a in v:
        total += a
        out.append(total)
    return tuple(out)


def uncumulative(v):
    return (v[0],)+tuple(v[i]-v[i-1] for i in range(1,len(v)))


def adjacent(u,v,literal=False):
    if u == v:
        return False
    delta = tuple(b-a for a,b in zip(cumulative(u),cumulative(v)))
    if literal:
        return all(abs(a) <= 1 for a in delta)
    return all(a in (0,1) for a in delta) or all(a in (0,-1) for a in delta)


@lru_cache(None)
def complex_faces(d,r,literal=False):
    vertices = compositions(r,d+1)
    neighbors = [set() for _ in vertices]
    for i,j in combinations(range(len(vertices)),2):
        if adjacent(vertices[i],vertices[j],literal):
            neighbors[i].add(j)
            neighbors[j].add(i)
    faces = []
    def visit(prefix,candidates):
        for i in sorted(candidates):
            face = prefix+(i,)
            faces.append(tuple(vertices[j] for j in face))
            visit(face,{j for j in candidates if j > i and j in neighbors[i]})
    visit((),set(range(len(vertices))))
    return tuple(faces)


def face_vector(faces):
    counts = Counter(len(face) for face in faces)
    return tuple(counts[k] for k in range(1,max(counts)+1))


def face_formula(d,r,k):
    return sum((-1)**(k-j)*comb(k,j)*comb(r*(j+1)+d,d) for j in range(k+1))


def determinant(matrix):
    a = [list(map(F,row)) for row in matrix]
    out = F(1)
    for col in range(len(a)):
        pivot = next((i for i in range(col,len(a)) if a[i][col]),None)
        if pivot is None:
            return F(0)
        if pivot != col:
            a[col],a[pivot] = a[pivot],a[col]
            out = -out
        value = a[col][col]
        out *= value
        for i in range(col+1,len(a)):
            ratio = a[i][col]/value
            for j in range(col,len(a)):
                a[i][j] -= ratio*a[col][j]
    return out


def top_cells(d,r):
    return tuple(face for face in complex_faces(d,r) if len(face) == d+1)


def refine_cell(cell,s,reverse_middle=False):
    ordered = sorted(cell,key=lambda v:cumulative(v))
    if reverse_middle and len(ordered) >= 4:
        ordered[1],ordered[2] = ordered[2],ordered[1]
    d = len(ordered)-1
    out = set()
    for small in top_cells(d,s):
        mapped = tuple(sorted(tuple(sum(w[i]*ordered[i][j] for i in range(d+1))
                                  for j in range(len(ordered[0]))) for w in small))
        out.add(mapped)
    return out


def refined_top_cells(d,r,s,reverse_middle=False):
    out = set()
    for cell in top_cells(d,r):
        out.update(refine_cell(cell,s,reverse_middle))
    return out


def point_cell(barycentric,r):
    """Lossless exact point->minimal mesh face plus positive local weights.

    Vertices have integer coordinates summing to r; their geometric
    positions are obtained by dividing by r.
    """
    barycentric = tuple(map(F,barycentric))
    need(all(a >= 0 for a in barycentric) and sum(barycentric) == 1,
         "point lies in the closed simplex")
    q = cumulative(tuple(r*a for a in barycentric))[:-1]
    floors = tuple(a.numerator//a.denominator for a in q)
    fractions = tuple(a-b for a,b in zip(q,floors))
    cuts = sorted({F(0),F(1),*fractions})
    out = {}
    for left,right in zip(cuts,cuts[1:]):
        midpoint = (left+right)/2
        s = tuple(a+int(b > midpoint) for a,b in zip(floors,fractions))+(r,)
        vertex = uncumulative(s)
        need(min(vertex) >= 0, "decoder retains the ordered cumulative chamber")
        out[vertex] = out.get(vertex,F(0))+(right-left)
    need(all(weight > 0 for weight in out.values()), "minimal-face weights are positive")
    need(sum(out.values()) == 1, "local weights sum to one")
    need(all(sum(weight*v[j] for v,weight in out.items())/r == barycentric[j]
             for j in range(len(barycentric))), "point reconstructs exactly")
    need(all(adjacent(u,v) for u,v in combinations(out,2)), "decoder face is a repaired clique")
    return out


def twice_point_cell(barycentric,r,s):
    coarse = point_cell(barycentric,r)
    ordered = sorted(coarse,key=cumulative)
    local = point_cell(tuple(coarse[v] for v in ordered),s)
    out = {}
    for w,weight in local.items():
        vertex = tuple(sum(w[i]*ordered[i][j] for i in range(len(ordered)))
                       for j in range(len(barycentric)))
        out[vertex] = weight
    return out


def coordinate_symmetries(d,r,literal=False):
    vertices = compositions(r,d+1)
    edges = {tuple(sorted((u,v))) for u,v in combinations(vertices,2) if adjacent(u,v,literal)}
    good = []
    for p in permutations(range(d+1)):
        changed = {tuple(sorted((tuple(u[i] for i in p),tuple(v[i] for i in p)))) for u,v in edges}
        if changed == edges:
            good.append(p)
    return tuple(good)


def dihedral_permutations(m):
    return {tuple((sign*i+shift)%m for i in range(m))
            for sign in (-1,1) for shift in range(m)}


def main():
    triangles = []
    triangle_distance_checks = 0
    for r in range(1,13):
        faces = complex_faces(2,r)
        literal = complex_faces(2,r,True)
        expected = (comb(r+2,2),3*r*(r+1)//2,r*r)
        expected_bad = (comb(r+2,2),2*r*r+r,2*r*r-r,r*(r-1)//2)
        need(face_vector(faces) == expected, "triangle repaired f-vector")
        need(face_vector(literal) == tuple(x for x in expected_bad if x), "literal clique f-vector")
        for u,v in combinations(compositions(r,3),2):
            delta = tuple(b-a for a,b in zip(u,v))
            root_edge = sorted(delta) == [-1,0,1]
            need(adjacent(u,v) == root_edge, "triangle edge iff unit root transfer")
        vertices = compositions(r,3)
        neighbors = {u:{v for v in vertices if adjacent(u,v)} for u in vertices}
        corners = {tuple(r if i == j else 0 for i in range(3)) for j in range(3)}
        need({v for v in vertices if len(neighbors[v]) == 2} == corners,
             "triangle corners are intrinsically the degree-two vertices")
        for source in vertices:
            distances,frontier = {source:0},[source]
            while frontier:
                current = frontier.pop()
                for other in neighbors[current]:
                    if other not in distances or distances[other] > distances[current]+1:
                        distances[other] = distances[current]+1
                        frontier.append(other)
            for target in vertices:
                need(distances[target] == sum(abs(a-b) for a,b in zip(source,target))//2,
                     "independent graph distance equals half the barycentric L1 distance")
                triangle_distance_checks += 1
        triangles.append(dict(r=r,repaired_faces=expected,literal_clique_faces=face_vector(literal),
                              extra_edges=r*(r-1)//2))
    u,v,p,q = (0,2,0),(1,0,1),(0,1,1),(1,1,0)
    need(adjacent(u,v,True) and not adjacent(u,v), "smallest long-chord hostile")
    need(adjacent(p,q), "crossed legitimate mesh edge")
    need(tuple(u[i]+v[i] for i in range(3)) == tuple(p[i]+q[i] for i in range(3)),
         "two diagonals meet at their common midpoint")
    need(all(adjacent(a,b,True) for a,b in combinations((u,v,p,q),2)), "coplanar literal K4")

    dimensions = []
    for d in range(1,5):
        for r in range(1,5):
            faces = complex_faces(d,r)
            f = face_vector(faces)
            need(f == tuple(face_formula(d,r,k) for k in range(d+1)), "unimodular Ehrhart face formula")
            need(f[-1] == r**d and sum((-1)**k*x for k,x in enumerate(f)) == 1,
                 "volume count and Euler characteristic")
            for face in top_cells(d,r):
                matrix = [[face[j+1][i]-face[0][i] for j in range(d)] for i in range(d)]
                need(abs(determinant(matrix)) == 1, "each mesh simplex is unimodular")
            dimensions.append(dict(d=d,r=r,f_vector=f))

    refinements = []
    for d in (1,2,3):
        for r,s in product(range(1,4),repeat=2):
            refined = refined_top_cells(d,r,s)
            direct = set(top_cells(d,r*s))
            need(refined == direct, "compatible edgewise degrees multiply")
            refinements.append(dict(d=d,r=r,s=s,top_cells=len(direct)))
    bad_refinement = refined_top_cells(3,2,2,True)
    direct = set(top_cells(3,4))
    need(bad_refinement != direct, "arbitrary tetrahedral child ordering changes refinement")
    wrong_cell = min(bad_refinement-direct)

    symmetry_rows = []
    for d in range(1,6):
        need(len(coordinate_symmetries(d,1)) == factorial(d+1), "unrefined simplex has every coordinate symmetry")
        for r in (2,3):
            actual = set(coordinate_symmetries(d,r))
            need(actual == dihedral_permutations(d+1), "coordinate symmetries are exactly dihedral")
            symmetry_rows.append(dict(d=d,r=r,coordinate_symmetries=len(actual),unrefined_symmetries=factorial(d+1)))
    need(len(coordinate_symmetries(2,2,True)) == 2, "literal absolute rule breaks triangle coordinate symmetry")
    tetra_u,tetra_v = (1,0,1,0),(0,1,0,1)
    changed_u,changed_v = (1,1,0,0),(0,0,1,1)
    need(adjacent(tetra_u,tetra_v) and not adjacent(changed_u,changed_v), "tetrahedral chart-order hostile")

    decoder_checks = nested_checks = 0
    for d in (1,2,3):
        for denominator in range(1,9):
            for numerator in compositions(denominator,d+1):
                point = tuple(F(a,denominator) for a in numerator)
                for r in range(1,7):
                    point_cell(point,r)
                    decoder_checks += 1
                for r,s in ((2,2),(2,3),(3,2)):
                    need(twice_point_cell(point,r,s) == point_cell(point,r*s),
                         "hierarchical address flattens with exact weights and inherited order")
                    nested_checks += 1
    a,b = (F(1,10),F(2,10),F(7,10)),(F(2,10),F(1,10),F(7,10))
    code_a,code_b = point_cell(a,2),point_cell(b,2)
    need(set(code_a) == set(code_b) and code_a != code_b, "cell address alone loses local geometric position")
    for signed_vector in product(range(-2,3),repeat=4):
        need(uncumulative(cumulative(signed_vector)) == signed_vector, "cumulative integer shear is lossless")

    # Combinatorial barycentric counts for a triangle at recursion depth k.
    barycentric_rows = []
    V,E,T = 3,3,1
    for depth in range(6):
        need(T == 6**depth, "barycentric face count is exponential in depth")
        need((V,E,T) == ((6**depth+3*2**depth)//2+1,(3*6**depth+3*2**depth)//2,6**depth),
             "barycentric recursion closed forms")
        barycentric_rows.append(dict(depth=depth,f_vector=(V,E,T)))
        V,E,T = V+E+T,2*E+6*T,6*T
    report = dict(status="PROVED scoped statements with FINITE-EXACT controls; classical edgewise subdivision",
        smallest_literal_hostile=dict(r=2,u=u,v=v,cumulative_difference=(1,-1,0),
            crossed_edge=(p,q),squared_lengths=dict(extra=6,mesh=2)),
        triangular_census=triangles,general_simplex_census=dimensions,
        triangle_graph_distance_checks=triangle_distance_checks,
        compatible_refinement_checks=refinements,
        arbitrary_child_order_hostile=dict(d=3,r=2,s=2,cell_absent_from_direct=wrong_cell,
            spurious_top_cells=len(bad_refinement-direct),missing_top_cells=len(direct-bad_refinement)),
        coordinate_symmetry_census=symmetry_rows,
        tetrahedron_chart_hostile=dict(edge=(tetra_u,tetra_v),nonedge_after_swap=(changed_u,changed_v)),
        point_decoder_checks=decoder_checks,hierarchical_point_decoder_checks=nested_checks,
        cell_only_loss=dict(first_point=list(map(str,a)),second_point=list(map(str,b)),
            common_vertices=sorted(code_a),first_weights={str(k):str(v) for k,v in code_a.items()},
            second_weights={str(k):str(v) for k,v in code_b.items()}),
        barycentric_depth_census=barycentric_rows)
    print(json.dumps(report,indent=2))
    print("PASS: all checks remain active under -O")


if __name__ == "__main__":
    main()
