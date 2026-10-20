"""Independent C5 source-context lift replay: no imports from the main proof code."""
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

C5 = (6, 12, 24, 17, 3)
D = (4, 29, 8, 1, 13)  # b1 -> b4 -> C3(b0,b2,b3)


def need(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def join(a, b):
    m, k = len(a), len(b)
    return tuple(x | (((1 << k)-1) << m) for x in a) + tuple(x << m for x in b)


def hp(t, p):
    return all(t[a] >> b & 1 for a,b in zip(p,p[1:]))


def blue(t, p):
    n = len(p)
    return all((t[p[i]] >> p[j] & 1) == (t[p[n-1-j]] >> p[n-1-i] & 1)
               for i in range(n) for j in range(i+1,n))


def reverse(t, p):
    target = list(t)
    for a,b in zip(p,p[1:]):
        target[a] ^= 1 << b
        target[b] ^= 1 << a
    return tuple(target)


def four_paths(a, m):
    b0,b1,b2,b3,b4 = range(m,m+5)
    return (a+(b0,b1,b2,b3,b4),
            (b0,)+a[::-1]+(b4,b3,b2,b1),
            (b3,b0,b2,b4)+a+(b1,),
            (b2,b0,b3,b4,b1)+a[::-1])


def source_tournaments(m):
    pairs = list(combinations(range(1,m),2))
    for code in range(1 << len(pairs)):
        t = [(1 << m)-2] + [0]*(m-1)
        for k,(i,j) in enumerate(pairs):
            if code >> k & 1:
                t[i] |= 1 << j
            else:
                t[j] |= 1 << i
        yield tuple(t)


def main():
    rows = []
    for m in range(1,7):
        tournaments = presentations = 0
        for a in source_tournaments(m):
            tournaments += 1
            for tail in permutations(range(1,m)):
                p = (0,)+tail
                if not hp(a,p):
                    continue
                presentations += 1
                t = join(a,C5)
                for j,path in enumerate(four_paths(p,m),1):
                    need(sorted(path) == list(range(m+5)), ('permutation',m,p,j))
                    need(hp(t,path), ('not HP',m,a,p,j))
                    need(not blue(t,path), ('blue',m,a,p,j))
                    t = reverse(t,path)
                need(t == join(a,D), ('target',m,a,p))
        rows.append({'m':m,'labelled_source_contexts':tournaments,'Hamiltonian_context_presentations':presentations})
    symbolic = []
    for m in range(5,101):
        paths = four_paths(tuple(range(m)),m)
        u,v,r,s = 0,m-3,m+4,m-4
        for step in (1,2):
            p = paths[step]
            reflection = dict(zip(p,p[::-1]))
            need((reflection[v],reflection[u]) == (r,s), ('reflected pair',m,step))
            # Original source arc u->v is true; original cross arc b4->s is false.
            # Neither unordered pair is toggled by either prior R move.
            for previous in paths[:step]:
                edges = {frozenset((a,b)) for a,b in zip(previous,previous[1:])}
                need(frozenset((u,v)) not in edges and frozenset((r,s)) not in edges,
                     ('known arc toggled',m,step))
        symbolic.append(m)
    payload = {'status':'FINITE-EXACT supporting the all-order proof in five_atom_audit.md',
               'universe':'All labelled source contexts on 0..m-1 with source0, m1..6, every directed Hamiltonian path; no isomorphism filter',
               'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
               'C5_out_masks':C5,'target_local_out_masks':D,'rows':rows,
               'total_contexts':sum(r['labelled_source_contexts'] for r in rows),
               'total_presentations':sum(r['Hamiltonian_context_presentations'] for r in rows),
               'symbolic_reflection_and_untouched_pair_checks_m':symbolic,
               'limitations':['Source context only, no right context', 'The sampled symbolic sizes supplement but do not replace the uniform index proof']}
    dest = Path(__file__).with_suffix('.json')
    dest.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'all_checks':'PASS','contexts':payload['total_contexts'],'presentations':payload['total_presentations'],
                      'receipt_sha256':sha256(dest.read_bytes()).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
