"""Exact simultaneous-interval controls and full named polynomial witnesses."""
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


def main():
    interval_counts = {}
    scalar_checks, weighted_checks, ambient_checks = 0,0,0
    for k in range(1,7):
        total = 1 << (k*(k-1)//2)
        seen, preserving = set(),0
        baseline = covers(tournament(k,0)) if k <= 5 else None
        for eta in range(total):
            delta, eps, selected = interval_data(k,eta)
            need(inverse_intervals(k,delta) == eta, "interval Mobius inverse")
            seen.add(delta)
            K = tournament(k,delta)
            preserving += transitive(K)
            if k <= 5:
                q = covers(K)
                need(all(a >= b for a,b in zip(q,baseline)), "cover domination")
                need((q == baseline) == transitive(K), "scalar classification")
                scalar_checks += 1
        need(len(seen) == total and preserving == factorial(k), "interval bijection count")
        interval_counts[str(k)] = {"patterns":total,"scalar_preserving":preserving}

    fixtures = [((2,4,1),(0,),(2,0)), ((2,0),(0,),(0,),(0,)),
                ((2,4,1),(2,4,1),(0,))]
    contexts = [(0,), (2,0), (0,1), (2,4,1), (2,12,9,1)]
    for blocks in fixtures:
        k = len(blocks)
        original = substitute(tournament(k,0),blocks)
        original_q = covers(original)
        block_product = prod(hamilton(a) for a in blocks)
        for eta in range(1 << (k*(k-1)//2)):
            delta,eps,selected = interval_data(k,eta)
            K = tournament(k,delta)
            new_blocks = tuple(opposite(a) if bit else a for a,bit in zip(blocks,eps))
            changed = substitute(K,new_blocks)
            need(changed == literal_xor(original,blocks,selected), "literal interval action")
            new_q = covers(changed)
            need(all(a >= b for a,b in zip(new_q,original_q)), "weighted q domination")
            need((new_q == original_q) == transitive(K), "weighted equality classification")
            for Q in contexts:
                old_ambient = substitute(Q,(original,)+((0,),)*(len(Q)-1))
                new_ambient = substitute(Q,(changed,)+((0,),)*(len(Q)-1))
                gap = hamilton(new_ambient)-hamilton(old_ambient)
                need(gap == 0 if transitive(K) else gap >= 2*hamilton(Q)*block_product,
                     "ambient strict lower bound")
                ambient_checks += 1
            weighted_checks += 1

    replacements = {"nontransitive_scalar":0,"transitive_named":0}
    for m in range(2,5):
        old = ambient(tournament(m,0))
        for code in range(1,1 << (m*(m-1)//2)):
            U = tournament(m,code)
            new = ambient(U)
            if not transitive(U):
                need(hamilton(new)-hamilton(old) >= 2*3, "nontransitive replacement")
                replacements["nontransitive_scalar"] += 1
            else:
                # Both ports mimic root incidences to exterior; x beats core,
                # and y loses only to the one selected physical core vertex.
                x = 1 << (m+1)
                values = [(hamilton(extend(old,x,x | (1 << v))),
                           hamilton(extend(new,x,x | (1 << v)))) for v in range(m)]
                need(any(a != b for a,b in values), "named transitive reorder")
                replacements["transitive_named"] += 1

    # Full arbitrary contrast and full whole-converse controls for two strong
    # ambient examples. All 2^(2n) queries, not only a chosen marginal, replay.
    original_core = tournament(3,0)
    interval_families = {"crossing":[(0,1),(1,2)], "nested":[(0,2),(0,1)]}
    old = ambient(original_core)
    old_poly = polynomial(old)
    whole = difference(old_poly,polynomial(opposite(old)))
    need(all(mask.bit_count() <= 3 for mask in whole), "whole converse cubic")
    full_records = {}
    query_checks = 0
    whole_query_checks = 0
    for name,family in interval_families.items():
        new_core = literal_xor(original_core,((0,),(0,),(0,)),family)
        new = ambient(new_core)
        contrast = difference(old_poly,polynomial(new))
        need(contrast and all(mask.bit_count() <= 4 for mask in contrast), "full arbitrary contrast")
        n = len(old)
        for x in range(1 << n):
            for c in range(1 << n):
                mask = x | (c << n)
                got = sum(value for term,value in contrast.items() if term & mask == term)
                want = hamilton(extend(old,x,c))-hamilton(extend(new,x,c))
                need(got == want, "full contrast truth table")
                query_checks += 1
                if name == "crossing":
                    whole_got = sum(value for term,value in whole.items() if term & mask == term)
                    whole_want = hamilton(extend(old,x,c))-hamilton(extend(opposite(old),x,c))
                    need(whole_got == whole_want, "full whole-converse truth table")
                    whole_query_checks += 1
        full_records[name] = {"intervals":family,"original":old,"changed":new,
                              "scalar_difference":hamilton(old)-hamilton(new),
                              "full_degree":max(mask.bit_count() for mask in contrast),
                              "coefficient_omission":"Every unlisted full Boolean coefficient is zero.",
                              "complete_contrast":[[mask & ((1 << n)-1),mask >> n,value]
                                                   for mask,value in sorted(contrast.items())]}
    return {"status":"PASS","interval_counts":interval_counts,
            "quotient_cover_checks":scalar_checks,"weighted_interval_checks":weighted_checks,
            "ambient_H_checks":ambient_checks,"replacement_checks":replacements,
            "full_named_query_checks":query_checks,"whole_converse_query_checks":whole_query_checks,
            "full_contrast_records":full_records,
            "whole_converse_control":{"original":old,"converse":opposite(old),
                                      "complete_contrast":[[mask & 31,mask >> 5,value]
                                                           for mask,value in sorted(whole.items())]},
            "scope":"All-order interval scalar classification and arbitrary replacement of transitive source modules; nontransitive nested source cores remain outside the full no-collision theorem."}


if __name__ == '__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
