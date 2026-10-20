from collections import Counter, defaultdict
from itertools import combinations
import hashlib
import json
from pathlib import Path
from check_regular_transport import rows, encode, paths, complement_path, need


def main():
    summary = []
    for n in range(1, 6):
        bank = [rows(n, c) for c in range(1 << (n*(n-1)//2))]
        scores = [tuple(row.bit_count() for row in a) for a in bank]
        fibres = defaultdict(list)
        for code, score in enumerate(scores):
            fibres[score].append(code)
        incidences, direct = {}, {}
        path_total = return_total = 0
        for code, a in enumerate(bank):
            incident, weights = set(), Counter()
            h = 0
            for p in paths(a):
                h += 1
                s = complement_path(a, p)
                sid = encode(s)
                need(sid not in incident, 'unique path per intermediate')
                incident.add(sid)
                for q in paths(s):
                    u = complement_path(s, q)
                    uid = encode(u)
                    if scores[uid] == scores[code]:
                        weights[uid] += 1
            need(weights[code] == h, 'Gram diagonal raw H')
            incidences[code], direct[code] = incident, weights
            path_total += h
            return_total += sum(weights.values())
        pair_checks = 0
        for fibre in fibres.values():
            for t in fibre:
                for u in fibre:
                    need(direct[t][u] == len(incidences[t] & incidences[u]), 'direct return vs Gram')
                    pair_checks += 1
        summary.append({'n': n, 'tournaments': len(bank), 'score_fibres': len(fibres),
            'Hamilton_path_incidences': path_total, 'ordered_same_score_pairs': pair_checks,
            'same_score_return_presentations': return_total})
    regular = [t for t, score in enumerate(scores) if score == (2,2,2,2,2)]
    column_support = defaultdict(set)
    for t in regular:
        for s in incidences[t]:
            column_support[s].add(t)
    degrees = Counter(map(len,column_support.values()))
    need(degrees == {1:240,3:40}, 'private/face column counts')
    faces = [tuple(sorted(v)) for v in column_support.values() if len(v)==3]
    need(len(faces)==len(set(faces))==40, 'distinct face triples')
    for s, triple in column_support.items():
        if len(triple)!=3:
            continue
        need(all((u^v).bit_count()==6 for u,v in combinations(triple,2)), 'face triangle adjacency')
        a=rows(5,s); high=[v for v,r in enumerate(a) if r.bit_count()==3]
        low=[v for v,r in enumerate(a) if r.bit_count()==1]
        need(len(high)==len(low)==1,'face exceptional endpoints')
        x,y=high[0],low[0]; inside=[v for v in range(5) if v not in (x,y)]
        need(a[y]>>x&1 and all(a[x]>>v&1 and a[v]>>y&1 for v in inside),'face blowup arcs')
        need(all(sum(a[u]>>v&1 for v in inside if v!=u)==1 for u in inside),'face interior C3')
        unanimous=[]
        for i,j in combinations(range(5),2):
            values=[bool(rows(5,t)[i]>>j&1) for t in triple]
            if len(set(values))==1:
                unanimous.append((i,j) if values[0] else (j,i))
        need(unanimous==[(x,y)],'unique unanimous edge reconstructs endpoints')
    vertices = [947,76,613,217,710,341,668,313,682,355,410,806]
    A = [[int((u^v).bit_count()==6) for v in vertices] for u in vertices]
    for i in range(12):
        for j in range(12):
            square = sum(A[i][k]*A[k][j] for k in range(12))
            expected = 3*(i==j)+2-2*(vertices[i]^vertices[j]==1023)
            need(square == expected, 'icosahedron exact A-square identity')
    result = {'status':'PASS', 'all_score_fibre_Gram_checks':summary,
        'icosahedron_integer_matrix_checks':144,
        'regular5_intermediate_columns':len(column_support),
        'regular5_column_degrees':dict(sorted(degrees.items())),
        'regular5_distinct_face_columns':len(set(faces)),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    Path(__file__).with_suffix('.json').write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
