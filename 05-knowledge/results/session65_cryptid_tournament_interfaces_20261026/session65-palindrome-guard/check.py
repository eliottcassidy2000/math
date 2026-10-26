#!/usr/bin/env python3
"""Independent finite controls; no producer/repository imports, no asserts."""
import hashlib
import itertools as it
import json
import math
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tape(word):
    n = 3 * len(word)
    a = [[int(i < j) for j in range(n)] for i in range(n)]
    for j, bit in enumerate(word):
        if bit == '1':
            x, y = 3*j, 3*j+2
            a[x][y], a[y][x] = 0, 1
    return a


def hamilton(a, q):
    return len(q) == len(a) and set(q) == set(range(len(a))) and all(a[x][y] for x,y in zip(q,q[1:]))


def product_words(word):
    factors = []
    for j, bit in enumerate(word):
        base = 3*j
        local = [(base,base+1,base+2)] if bit == '0' else [
            (base,base+1,base+2), (base+1,base+2,base), (base+2,base,base+1)]
        factors.append(local)
    for parts in it.product(*factors):
        yield sum(parts, ())


def dfs_words(a):
    n = len(a)
    def visit(q, unused):
        if not unused:
            yield q
        for x in sorted(unused):
            if not q or a[q[-1]][x]:
                yield from visit(q+(x,), unused-{x})
    yield from visit((), set(range(n)))


def reversal(q):
    rho = [None]*len(q)
    for i,x in enumerate(q):
        rho[x] = q[-1-i]
    return rho


def is_anti(a, rho):
    return all(a[x][y] == a[rho[y]][rho[x]]
               for x in range(len(a)) for y in range(len(a)) if x != y)


def blue(a, q):
    return is_anti(a, reversal(q))


def native_k(a, q):
    require(hamilton(a,q), 'K word must be directed in actual source')
    path = {frozenset((x,y)) for x,y in zip(q,q[1:])}
    n = len(a)
    return [[0 if x == y else (a[x][y] if frozenset((x,y)) in path else a[y][x])
             for y in range(n)] for x in range(n)]


def scc_sizes(a):
    n = len(a)
    reach = [[int(x == y or a[x][y]) for y in range(n)] for x in range(n)]
    for z in range(n):
        for x in range(n):
            if reach[x][z]:
                for y in range(n):
                    reach[x][y] |= reach[z][y]
    pending = set(range(n))
    components = []
    while pending:
        x = min(pending)
        part = {y for y in pending if reach[x][y] and reach[y][x]}
        components.append(part)
        pending -= part
    # Canonical order by the number of predecessor components, not input labels.
    unsorted_components = components[:]
    components.sort(key=lambda part: sum(a[next(iter(other))][next(iter(part))]
                                        for other in unsorted_components if other is not part))
    require(all(a[x][y] for i,c in enumerate(components) for d in components[i+1:]
                for x in c for y in d), 'SCC order check')
    return [len(c) for c in components]


def decode(sizes):
    out, run = [], 0
    for size in sizes+[3]: # sentinel flushes last singleton run, then removed
        if size == 1:
            run += 1
        else:
            require(size == 3 and run % 3 == 0, 'invalid tape SCC signature')
            out.extend('0'*(run//3))
            out.append('1')
            run = 0
    return ''.join(out[:-1])


def relabel(a, pi):
    n = len(a)
    b = [[0]*n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            b[pi[x]][pi[y]] = a[x][y]
    return b


def adj_from_bits(n,bits):
    a = [[0]*n for _ in range(n)]
    for k,(x,y) in enumerate(it.combinations(range(n),2)):
        a[x][y] = (bits >> k) & 1
        a[y][x] = 1-a[x][y]
    return a


def changed_pairs(a,b):
    return [[x,y] for x,y in it.combinations(range(len(a)),2) if a[x][y] != b[x][y]]


def main():
    rows, records = [], []
    guard_total = 0
    brute_sc, brute_perms, dfs_total = 0, 0, 0
    for m in range(7):
        row = dict(m=m, tapes=0, palindromes=0, hamiltonian_words=0, blue_words=0, black_words=0)
        for bits in it.product('01',repeat=m):
            w = ''.join(bits)
            a = tape(w)
            pal = w == w[::-1]
            qs = list(product_words(w))
            require(len(set(qs)) == len(qs) == 3**w.count('1'), 'product count')
            sizes = scc_sizes(a)
            require(decode(sizes) == w, 'SCC decode original')
            require(decode(list(reversed(sizes))) == w[::-1], 'SCC decode reverse')
            revlabels = list(reversed(range(len(a))))
            require(decode(scc_sizes(relabel(a,revlabels))) == w, 'SCC decode relabel')
            require((sizes == sizes[::-1]) == pal, 'SCC palindrome criterion')
            blues = []
            for q in qs:
                require(hamilton(a,q), 'product path legality')
                b = blue(a,q)
                require(b == pal, 'tape color theorem')
                if b:
                    blues.append(q)
                k = native_k(a,q)
                require(hamilton(k,q) and blue(k,q) == b and native_k(k,q) == a, 'tape K checks')
            # Independent DFS covers n <= 9; beyond this enumeration is product-based.
            if m <= 3:
                actual = list(dfs_words(a))
                require(set(actual) == set(qs), 'DFS/product equality')
                dfs_total += len(actual)
            # All candidate anti-isomorphisms up to order 6, without early exit.
            if m <= 2:
                anti_count = 0
                for pi in it.permutations(range(len(a))):
                    anti_count += int(is_anti(a,pi))
                    brute_perms += 1
                require(bool(anti_count) == pal, 'brute self-converse criterion')
                brute_sc += 1
            gw = '0'+w+'1'
            ga = tape(gw)
            gqs = list(product_words(gw))
            for q in gqs:
                require(hamilton(ga,q) and not blue(ga,q), 'nonpalindromic guard theorem')
            require(len(gqs) == 3**(w.count('1')+1), 'guard product count')
            guard_total += len(gqs)
            record = dict(word=w, scc_sizes=sizes, hamiltonian_words=len(qs), blue_words=len(blues),
                          black_words=len(qs)-len(blues), first_word=list(qs[0]),
                          first_word_reversal=reversal(qs[0]), first_word_blue=blue(a,qs[0]))
            records.append(record)
            row['tapes'] += 1
            row['palindromes'] += int(pal)
            row['hamiltonian_words'] += len(qs)
            row['blue_words'] += len(blues)
            row['black_words'] += len(qs)-len(blues)
        require(row['hamiltonian_words'] == 4**m, 'sum product paths')
        require(row['blue_words'] == (1 if m%2 == 0 else 4)*10**(m//2), 'weighted palindrome count')
        require(row['palindromes'] == 2**((m+1)//2), 'palindrome tape count')
        rows.append(row)

    krows, totals = [], dict(tournaments=0, candidate_words=0, legal_pairs=0, blue_pairs=0, black_pairs=0,
                            relabel_controls=0, pure_blue_tournaments=0)
    for n in range(1,6):
        row = dict(n=n, tournaments=2**math.comb(n,2), candidate_words=0, legal_pairs=0,
                   blue_pairs=0, black_pairs=0, relabel_controls=0, pure_blue_tournaments=0)
        pure_blue_ids, black_targets = set(), set()
        pi = list(range(n))
        if n >= 2:
            pi[0],pi[1] = pi[1],pi[0]
        for bits in range(row['tournaments']):
            a = adj_from_bits(n,bits)
            has_black = False
            for q in it.permutations(range(n)):
                row['candidate_words'] += 1
                if not hamilton(a,q):
                    continue
                b = blue(a,q)
                k = native_k(a,q)
                require(hamilton(k,q), 'K retained legal path')
                require(blue(k,q) == b, 'K own-word blue invariance')
                require(native_k(k,q) == a, 'K involutivity')
                require(len(changed_pairs(a,k)) == (n-1)*(n-2)//2, 'K changed pair count')
                aq = tuple(pi[x] for x in q)
                ra = relabel(a,pi)
                require(hamilton(ra,aq) and blue(ra,aq) == b, 'color relabel invariant')
                require(native_k(ra,aq) == relabel(k,pi), 'K relabel equivariance')
                row['relabel_controls'] += 1
                row['legal_pairs'] += 1
                row['blue_pairs' if b else 'black_pairs'] += 1
                has_black |= not b
                if not b:
                    kbits = sum(k[x][y] << i for i,(x,y) in enumerate(it.combinations(range(n),2)))
                    black_targets.add(kbits)
            if not has_black:
                pure_blue_ids.add(bits)
        row['pure_blue_tournaments'] = len(pure_blue_ids)
        require(not (black_targets & pure_blue_ids), 'pure-blue target has no black edge')
        require(row['legal_pairs'] == math.factorial(n)*2**((n-1)*(n-2)//2), 'legal pair count')
        for key in totals:
            totals[key] += row[key]
        krows.append(row)

    u,v = '0001','0000'
    a,b = tape(u),tape(v)
    target_qs = list(product_words(v))
    require(changed_pairs(a,b) == [[9,11]], 'order 12 one-bit write')
    require(all(not blue(a,q) for q in product_words(u)), 'source pure black')
    require(len(target_qs) == 1 and blue(b,target_qs[0]), 'target pure blue')
    n = len(a)
    require(n == 12 and (n-1)*(n-2)//2 == 55, 'odd-C example')
    identity = tuple(range(n))
    require(native_k(a,identity) != b, 'one legal K is not target')
    example = dict(source=u,target=v,n=n,C_n=55,changed_pairs=changed_pairs(a,b),
                   source_hamiltonian_words=3,source_black_words=3,target_hamiltonian_words=1,
                   target_blue_words=1,target_blue_word=list(target_qs[0]),
                   target_anti_reversal=reversal(target_qs[0]))

    # Distinct fixed guards can be padded on the left by up to three zero cells
    # to enter either order residue with odd C, while preserving nonpalindromicity.
    padding_rows = []
    for m in range(7):
        choices = [p for p in range(4) if (m+2+p)%4 in (0,1)]
        p = min(choices)
        for bits in it.product('01',repeat=m):
            w = ''.join(bits)
            padded = '0'*(1+p)+w+'1'
            require(padded != padded[::-1], 'padding preserves guard')
            nn = 3*len(padded)
            require(((nn-1)*(nn-2)//2)%2 == 1, 'padding odd C')
        padding_rows.append(dict(payload_m=m,additional_left_zeros=p,order=3*(m+2+p)))

    output = dict(status='PASS', enumeration='Complete product Hamiltonian words for all binary tapes m=0..6; not n! enumeration beyond small controls',
                  tape_rows=rows, tape_totals={key:sum(r[key] for r in rows) for key in rows[0] if key != 'm'},
                  tape_records=records, independent_dfs=dict(max_m=3, tapes=15, hamiltonian_words=dfs_total),
                  brute_self_converse=dict(max_m=2,tapes=brute_sc,candidate_permutations=brute_perms),
                  guarded_tapes=dict(payload_max_m=6,tapes=127,hamiltonian_words=guard_total,all_black=True),
                  native_k_rows=krows,native_k_totals=totals,order12_write=example,parity_padding=padding_rows,
                  limits=['All-length claims are proof consequences, not finite-search conclusions.',
                          'No black-K reachability or simulation sufficiency is claimed.',
                          'No repository edits, commits, or pushes.'])
    print(json.dumps(output,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
