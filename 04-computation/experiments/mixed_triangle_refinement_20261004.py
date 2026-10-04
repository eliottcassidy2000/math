"""Exact triangle meshes: scale, flags, and common-refinement counterexamples.

Universe and all checks are printed by main; Fraction coordinates throughout.
Run: python -X utf8 -B 04-computation/experiments/mixed_triangle_refinement_20261004.py
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def face(points):
    return tuple(sorted(points))


BASE = {face(((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1))))}


def average(points):
    return tuple(sum(p[i] for p in points)/len(points) for i in range(3))


def barycentric(mesh):
    return {face((a, average((a,b)), average((a,b,c))))
            for tri in mesh for a,b,c in permutations(tri)}


def edgewise(mesh, r):
    need(isinstance(r,int) and r >= 1, "positive integral scale")
    result = set()
    for a,b,c in mesh:
        def point(i,j):
            return tuple(((r-i-j)*a[k]+i*b[k]+j*c[k])/r for k in range(3))
        for i in range(r):
            for j in range(r-i):
                result.add(face((point(i,j),point(i+1,j),point(i,j+1))))
                if i+j < r-1:
                    result.add(face((point(i+1,j),point(i,j+1),point(i+1,j+1))))
    return result


def edges(mesh):
    return Counter(edge for tri in mesh for edge in combinations(tri,2))


def vertices(mesh):
    return set(v for tri in mesh for v in tri)


def signed_area(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def contains(tri,p):
    a,b,c = tri
    signs = [signed_area(a,b,p), signed_area(b,c,p), signed_area(c,a,p)]
    return all(x>=0 for x in signs) or all(x<=0 for x in signs)


def outside_cell(fine, coarse):
    """First fine triangle not contained in any one coarse triangle."""
    coarse = sorted(coarse)
    for tri in sorted(fine):
        if not any(all(contains(cell,p) for p in tri) for cell in coarse):
            return tri
    return None


def stats(mesh):
    es = edges(mesh)
    vs = vertices(mesh)
    degree = Counter(v for e in es for v in e)
    return dict(V=len(vs), E=len(es), F=len(mesh), boundary=sum(n==1 for n in es.values()),
                degree_histogram=dict(sorted(Counter(degree.values()).items())))


def validate(mesh):
    es = edges(mesh)
    need(set(es.values()) <= {1,2}, "edge incidence")
    need(sum(abs(signed_area(*tri)) for tri in mesh)==1, "total exact area")
    need(all(abs(signed_area(*tri))>0 for tri in mesh), "nondegenerate cells")
    s = stats(mesh)
    need(s['V']-s['E']+s['F']==1, "disk Euler characteristic")


def serial(point):
    return [str(x) for x in point]


def oriented_polygon(points):
    points = list(dict.fromkeys(points))
    if len(points)<3:
        return ()
    if sum(signed_area((Q(0),Q(0),Q(1)),a,b)
           for a,b in zip(points,points[1:]+points[:1]))<0:
        points.reverse()
    changed = True
    while changed and len(points)>=3:
        changed = False
        for i in range(len(points)):
            if signed_area(points[i-1],points[i],points[(i+1)%len(points)])==0:
                points.pop(i)
                changed = True
                break
    if len(points)<3:
        return ()
    offset = min(range(len(points)),key=points.__getitem__)
    return tuple(points[offset:]+points[:offset])


def intersect_triangles(first,second):
    """Sutherland-Hodgman clipping, exact and retaining the barycentric plane."""
    polygon = list(oriented_polygon(first))
    clip = oriented_polygon(second)
    for a,b in zip(clip,clip[1:]+clip[:1]):
        output = []
        if not polygon:
            return ()
        for p,q in zip(polygon,polygon[1:]+polygon[:1]):
            sp,sq = signed_area(a,b,p),signed_area(a,b,q)
            if sp>=0:
                output.append(p)
            if (sp<0<sq) or (sq<0<sp):
                t = sp/(sp-sq)
                output.append(tuple(p[k]+t*(q[k]-p[k]) for k in range(3)))
        polygon = list(dict.fromkeys(output))
    return oriented_polygon(polygon)


def overlay(first,second):
    """Full-dimensional cells carry BOTH parent triangles, not just area."""
    result = {}
    for a in sorted(first):
        for b in sorted(second):
            poly = intersect_triangles(a,b)
            if poly:
                need(poly not in result, "unique interior parent pair")
                need(all(contains(a,p) and contains(b,p) for p in poly), "intersection carrier")
                result[poly] = (a,b)
    return result


def overlay_barycentric(cells):
    triangles = {}
    for poly,parents in cells.items():
        center = average(poly)
        for a,b in zip(poly,poly[1:]+poly[:1]):
            midpoint = average((a,b))
            for endpoint in (a,b):
                tri = face((endpoint,midpoint,center))
                need(tri not in triangles, "unique overlay child")
                need(all(contains(parent,p) for p in tri for parent in parents), "both parent decoders")
                triangles[tri] = parents
    return triangles


def main():
    print("FINITE-EXACT: Fraction meshes; no numerical tolerance; triangle coordinate unit vectors.")
    sd = barycentric(BASE)
    for r in range(1,13):
        er = edgewise(BASE,r)
        left,right = edgewise(sd,r),barycentric(er)
        for mesh in (er,left,right):
            validate(mesh)
        need((stats(er)['V'], stats(er)['E'], stats(er)['F']) ==
             ((r+1)*(r+2)//2, 3*r*(r+1)//2, r*r), "pure scale formula")
        for mixed in (left,right):
            s = stats(mixed)
            need((s['V'],s['E'],s['F']) == (3*r*r+3*r+1,9*r*r+3*r,6*r*r), "mixed f-vector")
        need(outside_cell(right,sd) is None, "S E_r refines S")
        need(outside_cell(right,er) is None, "S E_r refines E_r")
    print("r=1..12: pure and mixed counts, manifold edges, area, and S(E_r) common refinement verified.")
    for r in range(1,6):
        for s in range(1,6):
            need(edgewise(edgewise(BASE,r),s)==edgewise(BASE,r*s), "edgewise multiplication")
    print("r,s=1..5: E_s(E_r)=E_(rs), exact equality of triangle sets.")
    current = BASE
    for k in range(5):
        s = stats(current)
        need((s['V'],s['E'],s['F'],s['boundary']) ==
             ((6**k+3*2**k)//2+1,(3*6**k+3*2**k)//2,6**k,3*2**k), "barycentric recurrence")
        validate(current)
        current = barycentric(current)
    print("k=0..4: barycentric f-vector and exact geometric validation.")
    e2s,s_e2 = edgewise(sd,2),barycentric(edgewise(BASE,2))
    print("E_2(S):", json.dumps(stats(e2s),sort_keys=True))
    print("S(E_2):", json.dumps(stats(s_e2),sort_keys=True))
    need(stats(e2s)['degree_histogram'] != stats(s_e2)['degree_histogram'], "abstract nonisomorphism")
    witness1 = outside_cell(e2s,edgewise(BASE,2))
    need(witness1 is not None, "opposite mixed order not a common refinement")
    witness2 = outside_cell(barycentric(s_e2),barycentric(sd))
    need(witness2 is not None, "one-level common refinement is not blindly iterable")
    print("E_2(S) not refining E_2; witness:", [serial(p) for p in witness1])
    print("S^2(E_2) not refining S^2; witness:", [serial(p) for p in witness2])
    print("S^2:", json.dumps(stats(barycentric(sd)),sort_keys=True))
    print("E_6:", json.dumps(stats(edgewise(BASE,6)),sort_keys=True))
    for r in range(1,5):
        for k in range(3):
            recursive = BASE
            for _ in range(k):
                recursive = barycentric(recursive)
            cells = overlay(edgewise(BASE,r),recursive)
            common = overlay_barycentric(cells)
            validate(common)
            # Independent exact carrier check is built into overlay_barycentric.
            if (r,k)==(2,2):
                print("Overlay E_2 with S^2: polygon cells",len(cells),
                      "cell-size histogram",dict(sorted(Counter(map(len,cells)).items())),
                      "barycentric refinement",json.dumps(stats(common),sort_keys=True))
    print("r=1..4,k=0..2: intersection overlay and barycentric common refinement pass; both parent maps retained.")
    print("All checks passed. General claims require the accompanying proofs, not this finite range.")


if __name__ == '__main__':
    main()
