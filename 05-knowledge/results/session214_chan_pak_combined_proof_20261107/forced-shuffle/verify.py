#!/usr/bin/env python3
"""Exact finite controls for the two-colour shuffle proof. No asserts."""
import hashlib
import itertools
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def edge(rows, u, v):
    return bool(rows[u] >> v & 1)


def merge(rows, a, b):
    """Head-comparison merge of two directed words, preserving both orders."""
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        if edge(rows, a[i], b[j]):
            out.append(a[i]); i += 1
        else:
            out.append(b[j]); j += 1
    return tuple(out) + tuple(a[i:]) + tuple(b[j:])


def nonordered_shuffle(rows, a, b):
    i = max(i for i, u in enumerate(a) if any(edge(rows, u, v) for v in b))
    j = next(j for j, v in enumerate(b) if edge(rows, a[i], v))
    return merge(rows, a[:i], b[:j]) + (a[i],) + tuple(b[j:]) + tuple(a[i+1:])


def valid(rows, word):
    return len(set(word)) == len(word) and all(edge(rows, u, v) for u, v in zip(word, word[1:]))


def tournament(n, bits):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> k & 1: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
            k += 1
    return rows


def endpoints(rows):
    n = len(rows); full = (1 << n) - 1
    incoming = [full ^ rows[v] ^ (1 << v) for v in range(n)]
    E = [[0] * n for _ in range(1 << n)]
    H = [0] * (1 << n); H[0] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1): continue
            rest = mask ^ (1 << v)
            if not rest: E[mask][v] = 1
            else:
                pred = rest & incoming[v]
                total = 0
                while pred:
                    bit = pred & -pred; pred -= bit
                    total += E[rest][bit.bit_length() - 1]
                E[mask][v] = total
        H[mask] = sum(E[mask])
    return E, H


def strong(rows):
    n = len(rows); full = (1 << n) - 1
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1: reach[i] |= reach[k]
    return all(r == full for r in reach)


def extension(rows, A):
    n = len(rows); full = (1 << n) - 1; B = full ^ A
    return [r | ((1 << n) if B >> i & 1 else 0) for i, r in enumerate(rows)] + [A]


def words(rows, mask):
    vv = [v for v in range(len(rows)) if mask >> v & 1]
    return [p for p in itertools.permutations(vv) if valid(rows, p)]


def matrix_controls():
    cases = 0; digest = hashlib.sha256()
    for s in range(1, 5):
        for r in range(1, 5):
            a = tuple(range(s)); b = tuple(range(s, s+r))
            for mask in range(1, 1 << (s*r)):
                rows = [(sum(1 << j for j in range(i+1, s)) if i < s else
                         sum(1 << j for j in range(i+1, s+r))) for i in range(s+r)]
                for i in range(s):
                    for j in range(r):
                        if mask >> (i*r+j) & 1: rows[i] |= 1 << (s+j)
                        else: rows[s+j] |= 1 << i
                w = nonordered_shuffle(rows, a, b)
                require(valid(rows, w), ('shuffle invalid', s, r, mask, w))
                require(tuple(v for v in w if v < s) == a, 'A order lost')
                require(tuple(v for v in w if v >= s) == b, 'B order lost')
                require(w != b+a, 'shuffle ordered')
                digest.update(bytes([s, r]) + mask.to_bytes(2, 'big') + bytes(w))
                cases += 1
    return {'cross_matrix_dimensions': '1<=|A|,|B|<=4', 'nonzero_matrices': cases,
            'output_sha256': digest.hexdigest()}


def exhaustive_controls():
    stats = []; cores = []; total_injection_pairs = 0
    for n in range(1, 7):
        full = (1 << n) - 1
        tournaments = cuts = positive = zero = ext = literal = 0
        for bits in range(1 << (n*(n-1)//2)):
            rows = tournament(n, bits); E, H = endpoints(rows); h = H[full]
            rev = [full ^ r ^ (1 << i) for i, r in enumerate(rows)]
            S, unused = endpoints(rev)
            require(h % 2 == 1, 'Redei parity control')
            if n <= 5:
                pp = words(rows, full)
                require(len(pp) == h, 'literal H mismatch')
                literal += len(pp)
            if n >= 3 and strong(rows):
                for v in range(n):
                    h0 = H[full ^ (1 << v)]; delta = h-h0
                    if 0 < delta <= 8:
                        cores.append({'n':n, 'rows':rows, 'root':v, 'h':h0, 'child':h,
                                      'k':delta//2})
            for A in range(1 << n):
                B = full ^ A
                active = any(rows[a] & B for a in range(n) if A >> a & 1)
                p = H[A]*H[B]
                q = 0
                if A and B:
                    q = sum(E[B][b]*S[A][a] for b in range(n) if B >> b & 1
                            for a in range(n) if A >> a & 1 and edge(rows,b,a))
                else: q = h
                if active:
                    positive += 1
                    require(h >= p+q, ('product surplus', n,bits,A,h,p,q))
                    require(2*q+1 <= h, ('half bound',n,bits,A,h,q))
                    if n <= 4:
                        images = []
                        for aa in words(rows,A):
                            for bb in words(rows,B):
                                w = nonordered_shuffle(rows,aa,bb)
                                require(valid(rows,w), 'literal image invalid')
                                require(w != bb+aa, 'literal image ordered')
                                images.append(w)
                        require(len(images) == p == len(set(images)), 'injection collision')
                        total_injection_pairs += p
                else:
                    zero += 1
                    require(q == h == p, 'dominance boundary')
                if n <= 5:
                    y = extension(rows,A); hy = endpoints(y)[1][-1]
                    require(hy >= 2*h-q, 'slot bound')
                    if active:
                        require(2*hy >= 3*h+1, 'growth bound')
                        require(hy >= max(h+p, 2*h-p), 'profile growth bound')
                        require(hy > h, 'active cut not positive')
                    else: require(hy == h, 'inactive cut changed H')
                    ext += 1
                cuts += 1
            tournaments += 1
        stats.append({'n':n,'labelled_tournaments':tournaments,'all_cuts':cuts,
                      'active_cuts':positive,'inactive_cuts':zero,
                      'direct_extension_counts':ext,'literal_H_words':literal})
    levels = {}
    for k in range(1,5):
        witnesses = {}
        for c in cores:
            d = c['k']
            if k % d: continue
            g = k//d
            if g not in (1,3): continue
            h = g*c['h']; child = g*c['child']
            require(h <= 4*k-1 and child <= 6*k-1, 'core scaling bound')
            if h in witnesses: continue
            rows = c['rows'][:]; root = c['root']
            if g == 3:
                n = len(rows); extra = 7 << n
                rows = [r|extra for r in rows] + [2 << n,4 << n,1 << n]
            Hfull = endpoints(rows)[1]
            require(Hfull[-1] == child, 'witness child mismatch')
            require(Hfull[((1 << len(rows))-1) ^ (1 << root)] == h, 'witness parent mismatch')
            witnesses[h] = {'rows':rows,'deleted_vertex':root,'parent_H':h,'child_H':child,
                            'strong_core_order':c['n'],'external_product':g}
        want = {1:[1,3],2:[1,5],3:[3,5,9],4:[1,5,9,11]}[k]
        require(sorted(witnesses) == want, ('level mismatch',k,sorted(witnesses)))
        levels[str(k)] = {'increment':2*k,'all_order_parent_H':want,
                         'witnesses':[witnesses[h] for h in want]}
    return {'universe':stats,'explicit_injection_pairs':total_injection_pairs,
            'increment_levels':levels}


def main():
    result = {'status':'PASS','scope':'Finite controls; all-order proof is in PROOF.md.',
              'matrix_controls':matrix_controls(),'exhaustive_controls':exhaustive_controls()}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__ == '__main__': main()
