#!/usr/bin/env python3
"""Exact guarded path count without using the height-transfer identities.

Enumerate the core's two ordered lanes; stream each ordinal chain through a
nine-state DP (empty / good / bad per lane).  The prefix and suffix need
not be small.  The order-11 constant block is an independent F_3 replay.
"""
from collections import Counter
from itertools import permutations, product


def path_count(word, tail, mask):
    good = 1
    # L records last before apex, R records first after apex; 0 = empty.
    def join(old, block):
        result = Counter()
        for (l, r), a in old.items():
            for (m, s), b in block.items():
                result[m or l, r or s] += a*b
        return result

    def chain(contacts):
        res = Counter({(0, 0): 1})
        for b in map(int, contacts):
            res = join(res, {(2 if b else good, 0): 1,
                             (0, good if b else 2): 1})
        return res

    arc = set()
    k = 0
    out = {0, 1, 4}
    for i in range(6):
        for j in range(i+1, 6):
            arc.add((j, i) if (mask >> k) & 1 else (i, j))
            k += 1
    core = Counter()
    for p in permutations(range(6)):
        for i in range(7):
            if all(j+1 == i or (p[j], p[j+1]) in arc for j in range(5)):
                l = 0 if not i else (2 if p[i-1] in out else good)
                r = 0 if i == 6 else (good if p[i] in out else 2)
                core[l, r] += 1
    state = join(join(chain(word), core), chain(tail))
    return sum(state[l,r] for l in (0,1) for r in (0,1))


def q(n):
    result = 0
    while n:
        p, d = divmod(n, 2)
        result += d*((p+1)//2)
        n = p
    return result


def height(v, u, h, a, b, c, k):
    return 2*(a*v*u+b*v+c*u-h*(q(v)+q(u)))+k


def order11():
    def mul(A,B):
        return [[sum(a*b for a,b in zip(row,col)) % 3 for col in zip(*B)] for row in A]

    for d, e in ((2,1),(3,-1)):
        A = [[2,e],[e,2]]
        P = [[1]*5,[-1]*d+[1]*(5-d)]
        PPt = mul(P,list(zip(*P)))
        PtP = mul(list(zip(*P)),P)
        assert mul(A,A) == PPt
        AP = mul(A,P)
        survivors = set()
        free = [(i,j) for i in range(1,4) for j in range(i,4)]
        for z in product(range(3),repeat=6):
            S = [[0]*5 for _ in range(5)]
            for (i,j),v in zip(free,z):
                S[i][j] = S[j][i] = v
            for j in (1,2,3,0,4):
                u = -AP[0][j]-sum(S[i][j] for i in (1,2,3))
                v = -AP[1][j]-sum(P[1][i]*S[i][j] for i in (1,2,3))
                S[0][j], S[4][j] = 2*(u-v)%3, 2*(u+v)%3
                if j in (1,2,3):
                    S[j][0], S[j][4] = S[0][j], S[4][j]
            if S[0][4] != S[4][0] or mul(S,S) != PtP:
                continue
            assert mul(P,S) == [[-x%3 for x in row] for row in AP]
            survivors.add(tuple(map(tuple, S)))
        form = Counter(min(''.join(str(S[i][j]) for i in p for j in p)
                           for p in permutations(range(5))
                           if all((i < d) == (p[i] < d) for i in range(5)))
                       for S in survivors)
        assert sorted(form.values()) == [2,2,2,6]
        print(d, sorted(form.items()))


if __name__ == '__main__':
    word = '101000000100011011010110010100101110'
    v, u = int(word,2)//2, 1
    assert (v,q(v),len(word)+9) == (21511975575,9802336406,45)
    E = 2**22*(2**23-1)
    pair = [path_count(word,'10',m) for m in (177,181)]
    assert pair == [E+1,E+3]
    assert pair[0] == height(v,u,29,640,191,158,111)
    print('45 sites:',pair)
    E = 2**18*(2**19-1)
    assert height(2422,44328,29,640,191,158,111) == E-39
    assert height(0,27887181,651,17616,6768,2691,2313) == E+401
    order11()
