"""Unequal cyclic modules: exact rational kernel and full named controls.

Standard library only. All checks survive python -O.
"""
import itertools
import json
from collections import defaultdict
from math import factorial, prod


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def opposite(a):
    return tuple(((1 << len(a))-1) ^ row ^ (1 << i) for i,row in enumerate(a))


def tournament(n, bits):
    out = [0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v = (j,i) if bits >> bit & 1 else (i,j)
        out[u] |= 1 << v
    return tuple(out)


def transitive(a):
    return sorted(row.bit_count() for row in a) == list(range(len(a)))


def hamilton_subsets(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            unused = a[v] & ~mask
            while unused:
                bit = unused & -unused
                unused -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return [1]+[sum(row) for row in dp[1:]]


def hamilton(a):
    return hamilton_subsets(a)[-1]


def covers(a):
    n = len(a)
    H = hamilton_subsets(a)
    dp = [[0]*(n+1) for _ in range(1 << n)]
    dp[0][0] = 1
    for mask in range(1,1 << n):
        sub = mask
        while sub:
            rest = mask ^ sub
            for k in range(1,n+1):
                dp[mask][k] += H[sub]*dp[rest][k-1]
            sub = (sub-1)&mask
    return tuple(dp[-1][1:])


def substitute(quotient, blocks):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1]+len(block))
    out = [0]*offsets[-1]
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            out[offsets[i]+u] |= row << offsets[i]
            for j in range(len(blocks)):
                if quotient[i] >> j & 1:
                    out[offsets[i]+u] |= ((1 << len(blocks[j]))-1) << offsets[j]
    return tuple(out)


def ambient(core):
    # Directed triangle root->outside1->outside2->root; core names come first.
    return substitute((2,4,1),(core,(0,),(0,)))


def interval_data(k, eta):
    pairs = list(itertools.combinations(range(k),2))
    selected = [pair for bit,pair in enumerate(pairs) if eta >> bit & 1]
    delta = sum((sum(a <= i and j <= b for a,b in selected)%2) << bit
                for bit,(i,j) in enumerate(pairs))
    eps = tuple(sum(a <= i <= b for a,b in selected)%2 for i in range(k))
    return delta, eps, selected


def inverse_intervals(k, delta):
    pairs = list(itertools.combinations(range(k),2))
    table = {pair: (delta >> bit & 1) for bit,pair in enumerate(pairs)}
    return sum((table[a,b] ^ table.get((a-1,b),0) ^ table.get((a,b+1),0)
                ^ table.get((a-1,b+1),0)) << bit for bit,(a,b) in enumerate(pairs))


def literal_xor(a, blocks, intervals):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1]+len(block))
    out = list(a)
    for lo,hi in intervals:
        for u,v in itertools.combinations(range(offsets[lo],offsets[hi+1]),2):
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    return tuple(out)


def extend(a, x, c):
    n = len(a)
    out = list(a)+[1 << (n+1),0]
    for v in range(n):
        if x >> v & 1:
            out[v] |= 1 << n
        else:
            out[n] |= 1 << v
        if c >> v & 1:
            out[v] |= 1 << (n+1)
        else:
            out[n+1] |= 1 << v
    return tuple(out)


def polynomial(a):
    n = len(a)
    result = defaultdict(int)
    for word in itertools.permutations(range(n+2)):
        terms = {0:1}
        for u,v in zip(word,word[1:]):
            if u < n and v < n:
                if not a[u] >> v & 1:
                    break
            elif u >= n and v >= n:
                if u != n:
                    break
            else:
                bit = 1 << ((v-n)*n+u if u < n else (u-n)*n+v)
                updated = defaultdict(int)
                for mask, value in terms.items():
                    if u < n:
                        updated[mask | bit] += value
                    else:
                        updated[mask] += value
                        updated[mask | bit] -= value
                terms = updated
        else:
            for mask,value in terms.items():
                result[mask] += value
    return {mask:value for mask,value in result.items() if value}


def difference(p, q):
    return {mask:p.get(mask,0)-q.get(mask,0) for mask in p.keys() | q.keys()
            if p.get(mask,0) != q.get(mask,0)}


# Three-variable polynomial arithmetic: no third-party symbolic algebra.
ONE = {(0,0,0):1}


def add(*polys):
    out = defaultdict(int)
    for poly in polys:
        for powers,value in poly.items():
            out[powers] += value
    return {p:v for p,v in out.items() if v}


def scale(poly, value):
    return {p:v*value for p,v in poly.items() if v*value}


def multiply(*polys):
    out = ONE
    for poly in polys:
        new = defaultdict(int)
        for p,v in out.items():
            for q,w in poly.items():
                new[tuple(a+b for a,b in zip(p,q))] += v*w
        out = {p:v for p,v in new.items() if v}
    return out


def bilinear(left, matrix, right):
    return add(*(scale(matrix[i][j],left[i]*right[j])
                 for i in range(3) for j in range(3)))


def port_numerator(G, d, alpha, gamma):
    e = (1,1,1)
    oa = tuple(1-x for x in alpha)
    oc = tuple(1-x for x in gamma)
    return add(
        multiply(add(d,bilinear(e,G,alpha)),
                 add(d,bilinear(oa,G,gamma)),
                 add(d,bilinear(oc,G,e))),
        multiply(add(d,bilinear(e,G,gamma)),
                 bilinear(oc,G,alpha),
                 add(d,bilinear(oa,G,e))))


def kernel_identity():
    a,b,c = ({(1,0,0):1},{(0,1,0):1},{(0,0,1):1})
    p = multiply(a,b,c)
    d = add(ONE,scale(p,-1))
    G = [[a,multiply(a,b),p],[p,b,multiply(b,c)],[multiply(a,c),p,c]]
    Gt = list(map(list,zip(*G)))
    coords = ((1,0,0),(0,1,0),(0,0,1))
    lhs = {}
    for i,j in ((0,1),(1,2),(2,0)):
        lhs = add(lhs,port_numerator(G,d,coords[i],coords[j]),
                  scale(port_numerator(Gt,d,coords[i],coords[j]),-1))
    variables = (a,b,c)
    S = add(*(multiply(variables[i],variables[i],variables[j],variables[j])
              for i,j in ((0,1),(1,2),(2,0))),
            scale(multiply(p,add(a,b,c)),-1),
            *(multiply(variables[i],variables[i],variables[j])
              for i,j in itertools.permutations(range(3),2)),scale(p,-6))
    need(lhs == multiply(d,S), "cyclic rational identity")
    return S


def q_at(q, k):
    return q[k] if k < len(q) else 0


def cyclic_formula(profiles):
    terms = []
    for k in range(max(map(len,profiles))):
        v = [[q_at(p,k+j) for j in range(3)] for p in profiles]
        E = -6*prod(x[1] for x in v)
        for i,j,l in ((0,1,2),(1,2,0),(2,0,1)):
            E += v[i][2]*v[j][2]*v[l][0]-v[i][2]*v[j][1]*v[l][1]
        for i,j,l in itertools.permutations(range(3)):
            E += v[i][2]*v[j][1]*v[l][0]
        terms.append((k+1)*E)
    return sum(terms),terms


def cross_condition(profiles):
    return all(q_at(profiles[i],k+2)*q_at(profiles[j],k)
               <= q_at(profiles[i],k+1)*q_at(profiles[j],k+1)
               for i,j in itertools.permutations(range(3),2)
               for k in range(1,max(map(len,profiles))))


def uniform_bits(blocks):
    bits,offset = [],0
    for block in blocks:
        bits.append(((1 << len(block))-1) << offset)
        offset += len(block)
    return bits


def literal_three_queries(blocks):
    T = substitute((2,4,1),blocks)
    U = substitute((4,1,2),blocks)
    masks = uniform_bits(blocks)
    records = []
    for i,j in ((0,1),(1,2),(2,0)):
        a,c = masks[i],masks[j]
        ft,fu = hamilton(extend(T,a,c)),hamilton(extend(U,a,c))
        records.append({"a":a,"c":c,"T":ft,"U":fu,"gap":ft-fu})
    return T,U,records


def main():
    S = kernel_identity()
    small = [(0,),tournament(2,0),tournament(3,0),(2,4,1)]
    condition_checks,formula_checks = 0,0
    for blocks in itertools.product(small,repeat=3):
        profiles = [(0,)+covers(b) for b in blocks]
        formula,terms = cyclic_formula(profiles)
        _,_,queries = literal_three_queries(blocks)
        need(formula == sum(r["gap"] for r in queries), "unequal cyclic identity")
        formula_checks += 1
        if cross_condition(profiles):
            need(formula < 0 and terms[0] < 0 and all(v <= 0 for v in terms[1:]),
                 "cross-profile sign criterion")
            condition_checks += 1
        if min(map(len,blocks)) >= 2:
            need(cross_condition(profiles), "orders-two-or-three theorem")

    blocks = (tournament(2,0),tournament(2,0),(2,4,1))
    T,U,queries = literal_three_queries(blocks)
    n = len(T)
    need(hamilton(T) == hamilton(U) == 159, "mixed7 scalar")
    need(queries[0]["T"] == 1959 and queries[0]["U"] == 2043, "mixed7 witness")
    need(sum(r["gap"] for r in queries) == -324, "mixed7 three-query sum")
    pt,pu = polynomial(T),polynomial(U)
    contrast = difference(pt,pu)
    need(max(m.bit_count() for m in pt.keys() | pu.keys()) <= 4, "port degree")
    checks = 0
    for degree in range(5):
        for variables in itertools.combinations(range(2*n),degree):
            mask = sum(1 << v for v in variables)
            a,c = mask & ((1 << n)-1),mask >> n
            value = sum(v for m,v in contrast.items() if m & mask == m)
            need(value == hamilton(extend(T,a,c))-hamilton(extend(U,a,c)),
                 "full coefficient determining query")
            checks += 1
    uniform = []
    block_masks = uniform_bits(blocks)
    for alpha in range(8):
        for gamma in range(8):
            a = sum(m for i,m in enumerate(block_masks) if alpha >> i & 1)
            c = sum(m for i,m in enumerate(block_masks) if gamma >> i & 1)
            ft,fu = hamilton(extend(T,a,c)),hamilton(extend(U,a,c))
            mask = a | c << n
            need(ft-fu == sum(v for m,v in contrast.items() if m & mask == m), "uniform polynomial")
            uniform.append({"alpha":alpha,"gamma":gamma,"a":a,"c":c,"T":ft,"U":fu})
    need(sum(r["T"] != r["U"] for r in uniform) == 44, "uniform disagreement count")
    need(all((r["T"]-r["U"])%4 == 0 for r in uniform), "uniform mod4 hostile")

    boundary_blocks = tuple(tournament(n,0) for n in (2,3,5))
    boundary_profiles = [(0,)+covers(b) for b in boundary_blocks]
    boundary_gap,boundary_terms = cyclic_formula(boundary_profiles)
    _,_,boundary_queries = literal_three_queries(boundary_blocks)
    need(all(all(q_at(q,k)**2 >= q_at(q,k-1)*q_at(q,k+1)
                 for k in range(1,len(q))) for q in boundary_profiles), "boundary separate LC")
    need(not cross_condition(boundary_profiles), "boundary cross failure")
    need(boundary_gap == sum(r["gap"] for r in boundary_queries) == 4540, "boundary positive gap")
    record = {
        "status":"PASS",
        "kernel": [{"powers":list(p),"coefficient":v} for p,v in sorted(S.items())],
        "small_formula_checks":formula_checks,
        "small_cross_condition_checks":condition_checks,
        "mixed7":{
            "T":T,"U":U,"profiles":[(0,)+covers(b) for b in blocks],
            "scalar":159,"three_queries":queries,"cyclic_gap":-324,
            "coefficient_determining_queries":checks,
            "degree_counts":{str(d):sum(m.bit_count()==d for m in contrast) for d in range(5)},
            "all_unlisted_coefficients_are_zero":True,
            "full_named_contrast":[{"mask":m,"coefficient":v} for m,v in sorted(contrast.items())],
            "all_64_uniform_queries":uniform},
        "separate_log_concavity_sign_boundary":{
            "block_sizes":[2,3,5],"profiles":boundary_profiles,
            "cross_condition":False,"cyclic_gap":boundary_gap,
            "terms":boundary_terms,"literal_queries":boundary_queries},
        "scope":"Exact three-query identity and conditional unequal-profile separation; no unrestricted collision theorem."
    }
    print(json.dumps(record,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()

