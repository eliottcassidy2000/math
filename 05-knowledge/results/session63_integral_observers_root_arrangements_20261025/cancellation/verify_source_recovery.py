"""Exact comparison with the concurrently recovered source witness; stdlib only."""
from itertools import product
from math import comb
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    points = list(product(range(5), repeat=2))
    index = {v:i for i,v in enumerate(points)}
    here = {(0,1),(0,2),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(4,0)}
    source = {(0,3),(0,4),(1,0),(3,0),(3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(4,3),(4,4)}
    negate = lambda v: tuple(-a % 5 for a in v)
    need(source == {negate(v) for v in here}, 'source connection set is exactly negative of reconstructed set')
    phi = [index[negate(v)] for v in points]
    def arc(S, i, j):
        return tuple((points[j][k]-points[i][k]) % 5 for k in range(2)) in S
    need(all(arc(here,i,j) == arc(source,phi[i],phi[j]) for i in range(25) for j in range(25)),
         'negation preserves all 625 core adjacency entries')
    phi50 = [2*phi[i//2]+i%2 for i in range(50)]
    def arc50(S, i, j):
        return (i%2 < j%2) if i//2 == j//2 else arc(S,i//2,j//2)
    need(all(arc50(here,i,j) == arc50(source,phi50[i],phi50[j]) for i in range(50) for j in range(50)),
         'blockwise negation preserves all 2500 inflated adjacency entries')
    P = [(comb(25,j)+24*(comb(5,j//5) if j%5==0 else 0))//25 for j in range(26)]
    local = [sum(P[i]*comb(i,j)*(-1)**(i-j) for i in range(j,26)) for j in range(26)]
    Q = [(-1)**(25+j)*a for j,a in enumerate(local)]
    need(all(a>=0 for a in Q), 'positive source-jet convention')
    need([a%5 for a in Q] == [0]*21+[4,3,3,4,1], 'source positive-jet mod5 formula agrees exactly')
    need(sum(P)==1342208, 'same subset-orbit count')
    print(json.dumps({'status':'PASS',
                      'source_commit':'a971f9b4fbdf409c31ba61ad9e8b3b00648536cf',
                      'source_path':'05-knowledge/results/unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md',
                      'source_section':4,'source_connection_set':sorted(source),
                      'reconstructed_connection_set':sorted(here),
                      'core_isomorphism':'u -> -u over F5^2',
                      'inflated_isomorphism':'(u,b) -> (-u,b)',
                      'core_adjacency_entries_checked':625,
                      'inflated_adjacency_entries_checked':2500,
                      'positive_jet_convention':'Q(u) = -P(-1-u)',
                      'mod5_source_nonzero_coefficients':{'21':4,'22':3,'23':3,'24':4,'25':1},
                      'subset_orbits':sum(P)},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
