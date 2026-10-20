"""Independent all-score-fibre and icosahedral-face audit.
Uses only the separate audit's arc-set routines, not either producer checker.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from independent_audit import tournament, code_of, words, degree, complement_path, require


def matmul(A, B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]


def characteristic_coefficients(A):
    # Exact Faddeev-LeVerrier: coefficients in descending degree order.
    n = len(A)
    B = [[int(i == j) for j in range(n)] for i in range(n)]
    out = [1]
    for k in range(1, n+1):
        AB = matmul(A, B)
        tr = sum(AB[i][i] for i in range(n))
        require(tr % k == 0, 'characteristic coefficient integrality')
        c = -tr // k
        out.append(c)
        B = [[AB[i][j] + c*(i == j) for j in range(n)] for i in range(n)]
    require(not any(x for row in B for x in row), 'Cayley-Hamilton closure')
    return out


def poly_times(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def main():
    summaries = []
    for n in range(1, 6):
        bank = [tournament(n,c) for c in range(1 << (n*(n-1)//2))]
        hps = [words(t,n) for t in bank]
        ds = [degree(t,n) for t in bank]
        by_score = defaultdict(list)
        for code,d in enumerate(ds):
            by_score[d].append(code)
        incident = []
        step = []
        for code,t in enumerate(bank):
            columns = {}
            for p in hps[code]:
                s = code_of(complement_path(t,p),n)
                require(s not in columns, 'path-to-column injectivity')
                columns[s] = p
            incident.append(set(columns))
            step.append(columns)
        direct = []
        for t in range(len(bank)):
            weight = Counter()
            for s in incident[t]:
                for u,q in step[s].items():
                    if ds[t] == ds[u]:
                        require((step[t][s][0],step[t][s][-1]) == (q[0],q[-1]), 'common incidence endpoints')
                        weight[u] += 1
            require(weight[t] == len(hps[t]), 'raw-H diagonal')
            direct.append(weight)
        checked = 0
        for fibre in by_score.values():
            for t in fibre:
                for u in fibre:
                    require(direct[t][u] == len(incident[t] & incident[u]), 'Gram identity')
                    checked += 1
        summaries.append(dict(n=n, tournaments=len(bank), score_fibres=len(by_score),
                              Hamilton_path_incidences=sum(map(len,hps)), ordered_same_score_pairs=checked,
                              same_score_return_presentations=sum(sum(w.values()) for w in direct)))

    regular = by_score[(2,)*5]
    column = defaultdict(set)
    for t in regular:
        for s in incident[t]:
            column[s].add(t)
    require(Counter(map(len,column.values())) == {1:240,3:40}, 'column degree census')
    private_count = Counter()
    face_columns = {}
    for s, support in column.items():
        if len(support) == 1:
            private_count[next(iter(support))] += 1
            continue
        require(len(support) == 3, 'column support size')
        ds_s = ds[s]
        x = ds_s.index(3)
        y = ds_s.index(1)
        inside = set(range(5)) - {x,y}
        unanimous = set.intersection(*(set(bank[t]) for t in support))
        require(unanimous == {(x,y)}, 'unique unanimous arc')
        majority = {uv for t in support for uv in bank[t] if sum(uv in bank[u] for u in support) >= 2}
        interior = {uv for uv in majority if set(uv) <= inside}
        reconstructed = frozenset({(y,x)} | {(x,v) for v in inside} | {(v,y) for v in inside} | interior)
        require(reconstructed == bank[s], 'complete intermediate reconstruction')
        require(all(sum(a == v for a,b in interior) == 1 for v in inside), 'interior is C3')
        triple = tuple(sorted(support))
        require(triple not in face_columns, 'face-map injectivity')
        face_columns[triple] = s
    require(set(private_count.values()) == {10} and len(private_count) == 24, 'ten private per source')
    graph = {t:{u for u in regular if u != t and direct[t][u]} for t in regular}
    all_triangles = {tuple(sorted((t,u,v))) for t in regular for u,v in combinations(graph[t],2) if v in graph[u]}
    require(set(face_columns) == all_triangles and len(all_triangles) == 40, 'onto all icosahedral faces')

    # Verify incidence factorization, both kernels, and scalar orbit compression.
    entries = 0
    for t in regular:
        two_stage = defaultdict(Fraction)
        for s in incident[t]:
            for u in column[s]:
                two_stage[u] += Fraction(1,15*len(column[s]))
        require(sum(two_stage.values()) == 1, 'second stochastic row sum')
        for u in regular:
            face_product = sum(t in face and u in face for face in face_columns)
            require(direct[t][u] == 10*(t == u) + face_product, 'private-plus-face Gram')
            require(direct[t][u] == 15*(t == u) + 2*(u in graph[t]), 'icosahedral weighted adjacency')
            require(Fraction(direct[t][u],25) == Fraction(3,5)*(t == u) + Fraction(2,25)*(u in graph[t]), 'pair-uniform kernel')
            require(two_stage[u] == Fraction(7,9)*(t == u) + Fraction(2,45)*(u in graph[t]), 'path-then-continuation kernel')
            entries += 1
    compression = Fraction(sum(sum(direct[t][u] for u in regular) for t in regular), 24)
    require(compression == 25, 'isomorphism orbit compression')

    vertices = [947,76,613,217,710,341,668,313,682,355,410,806]
    A = [[int(u in graph[t]) for u in vertices] for t in vertices]
    A2 = matmul(A,A)
    require(all(A2[i][j] == 3*(i == j) + 2 - 2*(vertices[i]^vertices[j] == 1023)
                for i in range(12) for j in range(12)), 'A-square identity')
    antipodal_pairs = {tuple(sorted((v,v^1023))) for v in vertices}
    require(len(antipodal_pairs) == 6, 'six antipodal pairs')
    for t in vertices:
        for pair in antipodal_pairs:
            require(sum(u in graph[t] for u in pair) == int(t not in pair), 'even quotient K6')
    expected = [1,-5]
    for _ in range(5): expected = poly_times(expected,[1,1])
    for _ in range(3): expected = poly_times(expected,[1,0,-5])
    actual = characteristic_coefficients(A)
    require(actual == expected, 'exact characteristic polynomial')
    W = [[15*(i == j)+2*A[i][j] for j in range(12)] for i in range(12)]
    expected_W = [1,-25]
    for _ in range(5): expected_W = poly_times(expected_W,[1,-13])
    for _ in range(3): expected_W = poly_times(expected_W,[1,-30,205])
    require(characteristic_coefficients(W) == expected_W, 'exact W characteristic polynomial')
    require(face_columns[(217,613,947)] == 654, 'explicit face representative')
    result = dict(status='PASS', all_score_fibre_Gram_checks=summaries,
                  regular5_intermediate_columns=len(column), degree_counts={1:240,3:40},
                  full_intermediate_reconstructions=40, reconstructed_graph_faces=40,
                  incidence_factorization_and_stochastic_matrix_entries=entries,
                  orbit_compression=str(compression),
                  A_characteristic_coefficients=actual,
                  W_characteristic_coefficients=expected_W,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  independent_audit_dependency_sha256=sha256(Path(__file__).with_name('independent_audit.py').read_bytes()).hexdigest())
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
