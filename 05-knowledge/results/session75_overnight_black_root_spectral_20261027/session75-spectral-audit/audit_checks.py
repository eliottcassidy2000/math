#!/usr/bin/env python3
"""Independent Session75 checks; standard library only; assertions survive -O.

Finite exact tests certify only the printed finite universes. The all-order
claims are proved in REPORT.md, not inferred from these tests.
"""
from functools import lru_cache
from itertools import permutations, product
from math import log, sqrt
from random import Random
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def matrix(n, value):
    return [[value(i, j) for j in range(n)] for i in range(n)]


def matvec(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]


def gram(a):
    return [[sum(x*y for x, y in zip(r, s)) for s in a] for r in a]


def rank_mod(a, prime=1000003):
    a = [[x % prime for x in row] for row in a]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [x*inv % prime for x in a[rank]]
        for r in range(rank+1, len(a)):
            c = a[r][col]
            if c:
                a[r] = [(x-c*y) % prime for x, y in zip(a[r], a[rank])]
        rank += 1
    return rank


# A leaf is a singleton tuple (integer label,); a node is (children, quotient).
def leaf(label):
    return (label,)


def leaves(t):
    if len(t) == 1:
        return [t[0]]
    return [x for c in t[0] for x in leaves(c)]


def size(t):
    return len(leaves(t))


def suppress(t):
    if len(t) == 1:
        return t
    cs = [suppress(c) for c in t[0]]
    if len(cs) == 1:
        return cs[0]
    return (cs, t[1])


def evaluate(t):
    labels = leaves(t)
    arcs = {}
    def visit(u):
        if len(u) == 1:
            return
        cs, q = u
        for c in cs:
            visit(c)
        for i, c in enumerate(cs):
            for j, d in enumerate(cs):
                if i != j:
                    for x in leaves(c):
                        for y in leaves(d):
                            arcs[x, y] = q[i][j]
    visit(t)
    return labels, [[0 if x == y else arcs[x, y] for y in labels] for x in labels]


def restrict(t, keep):
    if len(t) == 1:
        return t if t[0] in keep else None
    cs, q = t
    indexed = [(i, restrict(c, keep)) for i, c in enumerate(cs)]
    indexed = [(i, c) for i, c in indexed if c is not None]
    if not indexed:
        return None
    if len(indexed) == 1:
        return indexed[0][1]
    return ([c for i, c in indexed], [[q[i][j] for j, d in indexed] for i, c in indexed])


def cut(t, b):
    t = suppress(t)
    n = size(t)
    require(n >= 2 and b >= 2, "invalid cut domain")
    u = t
    while True:
        require(len(u) == 2, "large node unexpectedly leaf")
        large = [c for c in u[0] if size(c)*(b+1) > b*n]
        if not large:
            m = max(u[0], key=size)
            require(n < size(m)*(b+1) <= b*n, "balanced-cut interval")
            return m
        u = large[0]


@lru_cache(None)
def shapes(n, b):
    if n == 1:
        return ((),)
    out = []
    def compositions(total, length):
        if length == 1:
            if total >= 1:
                yield (total,)
            return
        for z in range(1, total-length+2):
            for rest in compositions(total-z, length-1):
                yield (z,)+rest
    for arity in range(2, min(n, b)+1):
        for parts in compositions(n, arity):
            for cs in product(*(shapes(x, b) for x in parts)):
                out.append(cs)
    return tuple(out)


def realize(shape, rng, counter):
    if not shape:
        z = counter[0]
        counter[0] += 1
        return leaf(z)
    cs = [realize(s, rng, counter) for s in shape]
    q = [[0]*len(cs) for _ in cs]
    for i in range(len(cs)):
        for j in range(i):
            q[i][j] = rng.choice((-1, 1))
            q[j][i] = -q[i][j]
    return (cs, q)


def check_cut(t, b, check_all_deletions=False):
    ls, a = evaluate(t)
    n = len(ls)
    mset = set(leaves(cut(t, b)))
    mis = [i for i, x in enumerate(ls) if x in mset]
    ois = [i for i, x in enumerate(ls) if x not in mset]
    require(all(len({a[i][j] for i in mis}) == 1 for j in ois), "selected set not module")
    # Exact cross-block Gram spectrum: eigenvalues m(n-m), twice, and zero.
    x = matrix(n, lambda i, j: a[i][j] if ((ls[i] in mset) != (ls[j] in mset)) else 0)
    g = gram(x)
    lam = len(mis)*len(ois)
    require(sum(g[i][i] for i in range(n)) == 2*lam, "cross Gram trace")
    require(all(sum(g[i][h]*g[h][j] for h in range(n)) == lam*g[i][j]
                for i in range(n) for j in range(n)), "cross Gram idempotence")
    keeps = [set(ls)-mset, mset]
    if check_all_deletions:
        keeps += [{ls[i] for i in range(n) if mask >> i & 1} for mask in range(1, 1 << n)]
    for keep in keeps:
        r = restrict(t, keep)
        rls, ra = evaluate(r)
        require(ra == [[a[ls.index(i)][ls.index(j)] for j in rls] for i in rls], "quotient restriction changed induced arcs")


def field27():
    # F_3[x]/(x^3-x+1). No roots in F_3, so the cubic is irreducible.
    def digits(a):
        return (a % 3, a//3 % 3, a//9)
    def number(v):
        return sum((x % 3)*3**i for i, x in enumerate(v))
    def sub(a, b):
        return number([x-y for x, y in zip(digits(a), digits(b))])
    def mul(a, b):
        aa, bb = digits(a), digits(b)
        c = [sum(aa[i]*bb[j] for i in range(3) for j in range(3) if i+j == d) % 3 for d in range(5)]
        for d in (4, 3):
            c[d-2] = (c[d-2]+c[d]) % 3
            c[d-3] = (c[d-3]+2*c[d]) % 3
        return number(c[:3])
    require(all((x**3-x+1) % 3 for x in range(3)), "cubic not irreducible")
    require(all(sorted(mul(a, b) for b in range(27)) == list(range(27)) for a in range(1, 27)), "nonzero multiplication not bijective")
    squares = {mul(a, a) for a in range(1, 27)}
    require(len(squares) == 13, "field square set cardinality")
    def chi(a):
        return 0 if a == 0 else 1 if a in squares else -1
    require(all(chi(mul(a,b)) == chi(a)*chi(b) for a in range(27) for b in range(27)), "character multiplicativity")
    k = matrix(27, lambda i, j: chi(sub(j, i)))
    require(all(k[i][j] == -k[j][i] for i in range(27) for j in range(27)), "Paley antisymmetry")
    require(all(sum(r) == 0 for r in k), "Paley regularity")
    require(gram(k) == matrix(27, lambda i,j: 26 if i == j else -1), "Paley Gram")
    # Hostile: integers modulo 27 are not the field of order 27.
    bad_squares = {a*a % 27 for a in range(1,27)}
    bad = matrix(27, lambda i,j: 0 if i == j else 1 if (j-i) % 27 in bad_squares else -1)
    require(any(bad[i][j] != -bad[j][i] for i in range(27) for j in range(27)), "ring hostile unexpectedly tournament")
    return k


def ternary_digits(x, k):
    return tuple(x//3**i % 3 for i in range(k-1,-1,-1))


def tower(k):
    n = 3**k
    labels = [ternary_digits(x,k) for x in range(n)]
    def arc(i,j):
        for a,b in zip(labels[i],labels[j]):
            if a != b:
                return 1 if (b-a) % 3 == 1 else -1
        return 0
    return matrix(n, arc), labels


def tower_checks():
    rows = []
    for k in (1,2,3):
        a, labels = tower(k)
        n = len(a)
        require(all(sum(row) == 0 for row in a), "tower not regular")
        basis = [[1]*n]
        mults = {}
        for p in range(k):
            lam = 3**(2*(k-p)-1)
            mults[lam] = 2*3**p
            for prefix in product(range(3), repeat=p):
                for w in ((1,-1,0),(1,0,-1)):
                    v = [w[z[p]] if z[:p] == prefix else 0 for z in labels]
                    require(matvec(a,matvec(a,v)) == [-lam*x for x in v], "tower contrast eigenvector")
                    basis.append(v)
        require(rank_mod(basis) == n, "tower contrast basis deficient")
        require(sum(mults.values()) == n-1, "tower multiplicity")
        # sqrt(lam)=3^(k-p-1)*sqrt(3), exact coefficient of sqrt(3).
        root3_coefficient = sum(2*3**p*3**(k-p-1) for p in range(k))
        require(root3_coefficient == 2*k*3**(k-1), "tower trace coefficient")
        rows.append({"k":k,"n":n,"gram_eigenvalue_multiplicities":mults,"trace_sqrt3_coefficient":root3_coefficient})
    bounds=[]
    for k in (1,3,5,7):
        n=3**k
        numerator=n*((n-1)-2*k*3**((k-1)//2))
        require(numerator % 4 == 0, "tower Paley bound integral")
        m=max(0,numerator//4)
        t=(m+n-2)//(n-1)
        bounds.append({"k":k,"n":n,"arc_lower_bound":m,"native_word_lower_bound":t})
    require(bounds[1]["arc_lower_bound"] == 54 and bounds[1]["native_word_lower_bound"] == 3,"n27 bounds")
    require(bounds[2]["arc_lower_bound"] == 9234 and bounds[2]["native_word_lower_bound"] == 39,"n243 bounds")
    return rows,bounds


def legal_paths(a):
    return [p for p in permutations(range(len(a))) if all(a[p[i]][p[i+1]] == 1 for i in range(len(a)-1))]


def native(a,p):
    n=len(a)
    out=[[-x for x in row] for row in a]
    for i,j in zip(p,p[1:]):
        require(a[i][j] == 1,"native path not directed")
        out[i][j],out[j][i]=a[i][j],a[j][i]
    return out


def disagreement(a,b):
    return sum(a[i][j] != b[i][j] for i in range(len(a)) for j in range(i))


def native_checks(rng):
    moves=0
    for n in (2,3,4):
        pairs=[(i,j) for i in range(n) for j in range(i)]
        for mask in range(1 << len(pairs)):
            a=matrix(n,lambda i,j:0)
            for bit,(i,j) in enumerate(pairs):
                a[i][j]=1 if mask >> bit & 1 else -1
                a[j][i]=-a[i][j]
            for p in legal_paths(a):
                c=native(a,p)
                require(disagreement(c,[[-x for x in row] for row in a]) == n-1,"native support tariff")
                require(native(c,p) == a,"native involution")
                moves+=1
    a=matrix(4,lambda i,j:0 if i==j else 1 if i<j else -1)
    c=native(a,(0,1,2,3))
    correction=[[c[i][j]+a[i][j] for j in range(4)] for i in range(4)]
    require(rank_mod(correction) == 4,"native correction rank-four hostile")
    # Genuine legal words with interleaved independent free relabelling/converse.
    for trial in range(100):
        n=6
        a=matrix(n,lambda i,j:0 if i==j else 1 if i<j else -1)
        b=[row[:] for row in a]
        for t in range(1,7):
            p=rng.choice(legal_paths(a))
            a=native(a,p)
            b=[[-x for x in row] for row in b]
            perm=list(range(n));rng.shuffle(perm)
            sign=rng.choice((-1,1))
            a=[[sign*a[i][j] for j in perm] for i in perm]
            b=[[sign*b[i][j] for j in perm] for i in perm]
            require(disagreement(a,b) <= t*(n-1),"native word tariff with symmetries")
    return {"all_native_moves_n2_n3_n4":moves,"random_symmetry_words":100,"native_rank_two_hostile_rank":4}


def main():
    rng=Random(750103)
    shape_count=0
    for b in (2,3,4):
        for n in range(2,8):
            for shape in shapes(n,b):
                t=realize(shape,rng,[0])
                # Unary wrappers test their transparent removal.
                t=([t],[[0]])
                check_cut(t,b,check_all_deletions=(n<=4))
                shape_count+=1
    # Deliberately maximally unbalanced binary chain.
    t=leaf(0)
    for n in range(2,33):
        t=([t,leaf(n-1)],[[0,1],[-1,0]])
        check_cut(t,2)
    # n <= b, including a star with every arc pattern at n=4.
    for b in (4,5,9,100):
        t=realize(((),(),(),()),rng,[0])
        check_cut(t,b,True)
    ratio_checks=0
    previous=float("inf")
    for i in range(1,10001):
        x=i/20000
        h=-x*log(x)-(1-x)*log(1-x)
        ratio=2*sqrt(x*(1-x))*log(2)/h
        require(ratio <= previous+1e-12,"entropy ratio numerical monotonicity control")
        previous=ratio;ratio_checks+=1
    pk=field27()
    tw,_=tower(3)
    sampled=[]
    for _ in range(128):
        p=list(range(27));rng.shuffle(p)
        sign=rng.choice((-1,1))
        b=[[sign*tw[i][j] for j in p] for i in p]
        m=disagreement(pk,b)
        pairing=sum(pk[i][j]*b[i][j] for i in range(27) for j in range(27))
        require(pairing == 27*26-4*m,"arc pairing factor four")
        require(m>=54,"sampled Paley-tower arc bound")
        sampled.append(m)
    spectra,bounds=tower_checks()
    result={
        "status":"PASS",
        "scope":"Independent finite exact controls; only entropy grid is floating-point; all-order proof separate",
        "parse_plane_shapes_n2_through_n7_b2_b3_b4":shape_count,
        "maximally_unbalanced_binary_chain_max_n":32,
        "entropy_ratio_grid_samples_numeric_only":ratio_checks,
        "GF27_irreducible_polynomial":"x^3-x+1 over F3",
        "GF27_Gram":"27I-J (exact)",
        "integer_mod27_hostile":"not antisymmetric, rejected",
        "tower_spectrum_exact_controls":spectra,
        "tower_Paley_bounds_exact_arithmetic":bounds,
        "sampled_Paley27_tower27_relabelling_converse_count":len(sampled),
        "sampled_arc_disagreement_range":[min(sampled),max(sampled)],
        "native_controls":native_checks(rng)
    }
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
