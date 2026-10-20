"""Complete small boundary for noncomplementary Hamilton-cycle edges.

Universe: every fixed-path blue tournament at even orders 4,6,8.
Filter: strong and exactly two outdegree values. Each retained matrix
gets an explicit Hamiltonian cycle with a noncomplementary-score edge.
"""
from pathlib import Path
import json
from moves import need, blue_tilings, is_path
from moves_five import strong
from moves_endpoint_probe import paths


def main():
    rows = []
    expected = {4:(4,3),6:(64,21),8:(4096,537)}
    for n in (4,6,8):
        universe = 0
        certificates = []
        for a in blue_tilings(n):
            universe += 1
            scores = [r.bit_count() for r in a]
            if len(set(scores)) != 2 or not strong(a):
                continue
            witness = None
            for p in paths(a,0):
                if not (a[p[-1]] >> p[0] & 1):
                    continue
                edge = next(((u,v) for u,v in zip(p,p[1:]+p[:1])
                             if scores[u]+scores[v] != n-1),None)
                if edge is not None:
                    witness = (p,edge)
                    break
            need(witness is not None, 'unresolved small alternating-cycle residual')
            p,(u,v) = witness
            need(is_path(a,p) and a[p[-1]] >> p[0] & 1, 'cycle certificate failed')
            need(a[u] >> v & 1 and scores[u]+scores[v] != n-1, 'edge certificate failed')
            certificates.append(dict(a=a,cycle=p,edge=(u,v),
                                     imbalance=scores[u]+scores[v]-(n-1)))
        need((universe,len(certificates)) == expected[n], 'universe count changed')
        rows.append(dict(n=n,blue_universe=universe,
                         strong_two_score_count=len(certificates),certificates=certificates))
        print(n,universe,len(certificates),'PASS')
    Path(__file__).with_suffix('.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__ == '__main__':
    main()
