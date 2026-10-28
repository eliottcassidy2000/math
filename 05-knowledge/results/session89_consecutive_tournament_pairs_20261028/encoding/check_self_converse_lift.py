#!/usr/bin/env python3
"""Lossless lift to odd self-converse tournaments and anti-involution decoder.
The independent generation/audit lanes check unmarked-class completeness.
This checker uses explicit active checks, including under python -O.
"""
from itertools import permutations
import json
from check_encoding import alpha, pairs, relabel_adj

CHECKS = 0

def need(condition, message="lift verification failed"):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def lift(m, n):
    """Vertices root=0, p_i=1+i, q_i=1+n+i; return adjacency rows, tau."""
    size = 2*n + 1
    rows = [0] * size
    def arc(a,b,bit):
        if bit:
            rows[a] |= 1 << b
        else:
            rows[b] |= 1 << a
    for i in range(n):
        p, q = 1+i, 1+n+i
        arc(p,q,1)
        d = (m >> (i*n+i)) & 1
        arc(0,p,d)
        arc(0,q,1-d)
    for i,j in pairs(n):
        x = (m >> (i*n+j)) & 1
        y = (m >> (j*n+i)) & 1
        r = x ^ y
        arc(1+i,1+j,x)
        arc(1+n+i,1+n+j,1-x)
        arc(1+i,1+n+j,r)
        arc(1+j,1+n+i,r)
    tau = (0,) + tuple(range(1+n,1+2*n)) + tuple(range(1,1+n))
    return tuple(rows), tau


def is_tournament(rows):
    size = len(rows)
    return all(0 <= row < (1 << size) for row in rows) and all(
        not ((rows[i] >> i) & 1) for i in range(size)) and all(
        ((rows[i] >> j) & 1) ^ ((rows[j] >> i) & 1) for i,j in pairs(size))


def is_anti(rows, p):
    return all(((rows[i] >> j) & 1) == ((rows[p[j]] >> p[i]) & 1)
               for i,j in pairs(len(rows)))


def decode(rows, tau):
    """Decode any odd tournament with a supplied anti-involution.
    Pair ordering uses minimum original vertex only as an arbitrary frame.
    Return a matrix bitmask and n; output is canonical only up to alpha.
    """
    size = len(rows)
    need(size % 2 == 1 and is_tournament(rows), "input must be an odd tournament")
    need(sorted(tau) == list(range(size)), "tau must be a permutation")
    need(all(tau[tau[i]] == i for i in range(size)), "tau must be an involution")
    need(is_anti(rows,tau), "tau must reverse every arc")
    fixed = [i for i in range(size) if tau[i] == i]
    need(len(fixed) == 1, "anti-involution must have exactly one fixed vertex")
    root = fixed[0]
    pq = []
    for i in range(size):
        j = tau[i]
        if i < j:
            pq.append((i,j) if (rows[i] >> j) & 1 else (j,i))
    n = len(pq)
    m = 0
    for i,(p,q) in enumerate(pq):
        m |= ((rows[root] >> p) & 1) << (i*n+i)
    for i,j in pairs(n):
        p_i,q_i = pq[i]
        p_j,q_j = pq[j]
        x = (rows[p_i] >> p_j) & 1
        r = (rows[p_i] >> q_j) & 1
        need(r == ((rows[p_j] >> q_i) & 1), "cross bits must be symmetric")
        m |= x << (i*n+j)
        m |= (x ^ r) << (j*n+i)
    return m,n


def canonical_pair(m,n):
    return min(alpha(m,p) for p in permutations(range(n)))


def run():
    report = {"schema":1}
    # Every input matrix through n=3; no filtering of the two core tournaments.
    examples = 0
    equivariances = 0
    for n in range(4):
        for m in range(1 << (n*n)):
            rows,tau = lift(m,n)
            need(is_tournament(rows), "lift is a tournament")
            need(is_anti(rows,tau), "declared tau is an anti-automorphism")
            need(all(tau[tau[i]] == i for i in range(2*n+1)), "tau squares to identity")
            need(decode(rows,tau) == (m,n), "labelled round trip")
            for p in permutations(range(n)):
                large = (0,) + tuple(1+p[i] for i in range(n)) + tuple(1+n+p[i] for i in range(n))
                need(lift(alpha(m,p),n)[0] == relabel_adj(rows,large), "equivariant lift")
                equivariances += 1
            examples += 1
    report["exhaustive_labelled_round_trips"] = {"max_n":3,"matrices":examples}
    report["exhaustive_lift_equivariance"] = {"max_n":3,"matrix_permutation_pairs":equivariances}
    # All alternative anti-involutions, plus all vertex relabelings, through order5.
    anti_checks = 0
    label_checks = 0
    for n in range(3):
        size = 2*n+1
        ps = list(permutations(range(size)))
        for m in range(1 << (n*n)):
            rows,tau = lift(m,n)
            canonical = canonical_pair(m,n)
            for p in ps:
                if all(p[p[i]] == i for i in range(size)) and is_anti(rows,p):
                    decoded,order = decode(rows,p)
                    need(order == n and canonical_pair(decoded,n) == canonical,
                         "different anti-involutions give same pair class")
                    anti_checks += 1
                inverse = [0]*size
                for i,j in enumerate(p): inverse[j]=i
                transported = tuple(p[tau[inverse[j]]] for j in range(size))
                decoded,order = decode(relabel_adj(rows,p),transported)
                need(order == n and canonical_pair(decoded,n) == canonical,
                     "arbitrary vertex relabeling preserves decoded pair class")
                label_checks += 1
    report["all_alternative_anti_involutions"] = {"max_n":2,"checks":anti_checks}
    report["all_lift_vertex_relabelings"] = {"max_n":2,"checks":label_checks}
    report["active_checks"] = CHECKS
    report["status"] = "PASS"
    return report

if __name__ == "__main__":
    print(json.dumps(run(),indent=2,sort_keys=True))
