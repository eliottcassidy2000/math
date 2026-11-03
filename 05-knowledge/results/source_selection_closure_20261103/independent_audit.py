#!/usr/bin/env python3
"""Independent cofactors, literal native-component enumeration, and controller audit.
No producer imports; reads frozen data and independently recomputes assertions.
"""
from fractions import Fraction as Q
from itertools import permutations
from functools import lru_cache
from math import gcd, isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok: raise RuntimeError(why)


def det(A):
    if len(A) == 1: return A[0][0]
    return sum((-1)**j*A[0][j]*det([row[:j]+row[j+1:] for row in A[1:]])
               for j in range(len(A)))


def inverse(A):
    d = det(A); n = len(A)
    return [[Q((-1)**(i+j)*det([row[:i]+row[i+1:] for row in A[:j]+A[j+1:]]), d)
             for j in range(n)] for i in range(n)]


def mul(A, B):
    return [[sum(row[k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for row in A]


@lru_cache(None)
def paths(T):
    return tuple(p for p in permutations(range(len(T)))
                 if all(T[a] >> b & 1 for a, b in zip(p, p[1:])))


def native(T, path):
    need(path in paths(T), 'legal native Hamilton word')
    full = (1 << len(T))-1
    rows = [full ^ (1 << i) ^ row for i, row in enumerate(T)]
    for a, b in zip(path, path[1:]):
        rows[a] ^= 1 << b; rows[b] ^= 1 << a
    return tuple(rows)


def run(T, ports, word):
    for letter in word:
        mask, drop = divmod(letter, 3); n = len(T)
        incoming = {ports[i] for i in range(2) if mask >> i & 1}
        T = tuple(row | ((1 << n) if i in incoming else 0) for i, row in enumerate(T))
        T += (sum(1 << i for i in range(n) if i not in incoming),)
        if drop < 2:
            ports = tuple(n if i == drop else p for i, p in enumerate(ports))
    return T, ports


def ratio(pair, p=11):
    b, t = (x % p for x in pair)
    need(b or t, 'primitive projective pair')
    return None if b == 0 else t*pow(b, -1, p) % p


def renewal(pair, k):
    b, t = pair
    return ((3*2**k+4*3**(k-1)+1)*b-2*t, (2**k+1)*b)


def prime(n):
    return n >= 2 and (n == 2 or (n % 2 and all(n % d for d in range(3, isqrt(n)+1, 2))))


def main():
    old = json.loads((HERE.parent/'two_port_return_induction_20261103/rooted_resets.json').read_text())
    certificate = json.loads((HERE/'grammar_inverse.json').read_text())
    P = [[1,0,0,0],[1,1,0,0],[1,0,1,0],[0,0,0,1]]
    Pi = inverse(P)
    matrices = {''.join(map(str, item['word'])):mul(mul(P, item['matrix']), Pi)
                for item in old['grammar']['alphabet']}
    need(matrices == certificate['matrices'], 'independent coordinate conjugation')
    need([det(matrices[w]) for w in ('041','403')] == [4,-4], 'cofactor determinants')
    for row in certificate['separation_certificates']:
        value = mul(mul(matrices[row['other']], inverse(matrices[row['target']])),
                    [[x] for x in row['dual']])
        need(all(x >= 0 for x in row['dual']) and all(x[0] <= 0 for x in value)
             and value[2][0] < 0, 'dual exact separation')
        need([str(x[0]) for x in value] == row['contraction'], 'supplied contraction')
    for M in matrices.values():
        out = mul(M, [[1],[1],[2],[0]])
        need(all(x >= 0 for row in M for x in row), 'positive physical matrices')
        need(all(out[i][0] >= [2,3,3,1][i] for i in range(4)), 'first-birth inverse bound')
        need(all(M[i][3] <= 3*out[i][0] for i in range(4)), 'scalar defect bound')
    branch = json.loads((HERE/'branch_probe.json').read_text())
    nodes = branch['named_native_component']; rows = [tuple(node['rows']) for node in nodes]
    edges = []
    for i, T in enumerate(rows):
        need(len(paths(T)) == 31, 'node count by all permutations')
        found = []
        for p in paths(T):
            U = native(T, p)
            if len(paths(U)) == 31:
                need(U in rows, 'complete neutral component')
                found.append((rows.index(U), p))
        expected = [(e['to'], tuple(e['path'])) for e in nodes[i]['edges']]
        need(sorted(found) == sorted(expected) and len(found) == 2, 'all neutral edges')
        edges.append([i, sorted(j for j, _ in found)])
    repair = branch['native_fixed_pair_repair']
    T = tuple(repair['source']); ports = tuple(repair['ports'])
    U = native(T, tuple(repair['path']))
    need(U == tuple(repair['target']), 'fixed-name native source')
    for source, key, target_h in ((T,'branch95',95),(U,'branch97_after_native',97)):
        endpoint = repair[key]
        V, pair = run(source, ports, endpoint['word'])
        need(V == tuple(endpoint['rows']) and pair == tuple(endpoint['ports']), 'literal retained ports')
        need(len(paths(V)) == target_h, 'independent factorial endpoint count')
    # Projective closure computed from matrices, not a supplied transition table.
    todo = [(3,1)]; seen = {}; transitions = []
    while todo:
        pair = todo.pop(); r = ratio(pair)
        if r in seen: continue
        seen[r] = pair
        for k in (3,5):
            if k == 5 and r == 7: continue
            child = renewal(pair,k); target = ratio(child)
            transitions.append([r,k,target]); todo.append(child)
    need(set(seen) == {0,1,3,4,5,6,7,None}, 'eight projective controller states')
    need(len(transitions) == 15, 'one forbidden fifth-block edge')
    b,t = renewal(renewal((3,1),3),5)
    need((b,t) == (76147,5973) and gcd(b,t) == 1, 'two-return source')
    values = [(2*b-t,20903),(226*b-t,2457607)]
    need(all(h == 7*p and prime(p) for h,p in values), 'deterministic trial division controls')
    print(json.dumps(dict(status='PASS', grammar_matrices=28, cofactor_separators=54,
                         native_nodes=6, native_edges=edges, positive_branches=[95,97],
                         controller_states=sorted((str(x) for x in seen)), controller_edges=transitions,
                         seven_prime_controls=[{'H':h,'prime':p,'trial_division_bound':isqrt(p)}
                                               for h,p in values]),indent=2,sort_keys=True))


if __name__ == '__main__': main()
