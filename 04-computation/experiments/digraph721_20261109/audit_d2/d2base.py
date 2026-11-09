"""audit D2 (THM-5056), item 5: the base below 512, recomputed independently.

The census (dom_K511.out) is used ONLY as a list of candidate connected graphs: its d and parity columns are ignored.
For every odd k < 512 this script
  * looks for an odd graph on exactly L(k) = bit length of k vertices with d = k, of the form C u T, where C is a
    connected census graph with d(C) | k, T is a threshold graph realising k / d(C) on L(k / d(C)) vertices (empty if
    the cofactor is 1), d(C) and the oddness of C are recomputed here (oddness: an explicit automorphism with an odd
    number of reversed edges, found among transpositions or by a full backtracking enumeration of Aut(C));
  * writes each realiser in wcheck format (for d2check, which re-verifies d, order and oddness of the whole union);
  * for k in {1, 5, 37}: verifies that every connected census graph whose d divides k is EVEN (full Aut enumeration),
    so no union of components can be odd (d is multiplicative, 5 and 37 are prime, a component with d = 1 is K_1).
usage: python3 d2base.py ../dom_K511.out d2base_witnesses.txt
"""
import sys, random
random.seed(20261109)

def bitlen(k): return k.bit_length()

def parse(line):
    p = line.split(); n = int(p[0]); e = p[5] if len(p) > 5 else ''
    adj = [0] * n
    for i in range(0, len(e), 2):
        u = int(e[i], 16); v = int(e[i + 1], 16); adj[u] |= 1 << v; adj[v] |= 1 << u
    return adj

def dcount(adj):                       # inclusion-exclusion: d = sum_U (-1)^|U| 2^(n - |N[U]|)
    n = len(adj); cl = [adj[v] | (1 << v) for v in range(n)]; tot = 0
    NU = [0] * (1 << n)
    for U in range(1 << n):
        if U:
            b = (U & -U).bit_length() - 1; NU[U] = NU[U & (U - 1)] | cl[b]
        tot += (-1) ** bin(U).count('1') * 2 ** (n - bin(NU[U]).count('1'))
    return tot

def rev_edges(adj, p):
    n = len(adj)
    return sum(1 for a in range(n) for b in range(a + 1, n) if (adj[a] >> b) & 1 and p[a] > p[b])

def is_aut(adj, p):
    n = len(adj)
    return all((((adj[a] >> b) & 1) == ((adj[p[a]] >> p[b]) & 1)) for a in range(n) for b in range(n))

def odd_transposition(adj):
    n = len(adj)
    for u in range(n):
        for v in range(u + 1, n):
            p = list(range(n)); p[u], p[v] = v, u
            if is_aut(adj, p) and rev_edges(adj, p) % 2: return p
    return None

def all_auts(adj):
    n = len(adj); dg = [bin(a).count('1') for a in adj]; img = [-1] * n; used = [False] * n; out = []
    def bt(v):
        if v == n: out.append(img[:]); return
        for c in range(n):
            if used[c] or dg[c] != dg[v]: continue
            if any(((adj[v] >> j) & 1) != ((adj[c] >> img[j]) & 1) for j in range(v)): continue
            img[v] = c; used[c] = True; bt(v + 1); used[c] = False; img[v] = -1
    bt(0); return out

def odd_cert(adj):
    p = odd_transposition(adj)
    if p is not None: return p, 'T'
    for g in all_auts(adj):
        if rev_edges(adj, g) % 2:
            assert is_aut(adj, g); return g, 'A'
    return None, 'E'

def threshold(m, N):                   # graph on N vertices with d = m (m odd, 1 <= m < 2^N), top-down (Lemma B1)
    if N <= 1:                         # empty graph or K_1 (the cone formula needs at least one vertex below)
        assert m == 1; return [0] * N
    if m >= 2 ** (N - 1):
        g = threshold(m - 2 ** (N - 1), N - 1); k = len(g)
        return [a | (1 << k) for a in g] + [(1 << k) - 1]
    return threshold(m, N - 1) + [0]

def union(a, b):
    na = len(a); return a + [x << na for x in b]

def main():
    cen = {}
    for line in open(sys.argv[1]):
        p = line.split(); n = int(p[0]); d = int(p[2])
        cen.setdefault(d, {}).setdefault(n, []).append(line)
    out = open(sys.argv[2], 'w')
    found = {}; kinds = {'T': 0, 'A': 0}; tested = 0
    for k in range(1, 512, 2):
        L = bitlen(k); best = None
        for a in sorted(cen):
            if k % a: continue
            m = k // a; Lm = 0 if m == 1 else bitlen(m)
            need = L - Lm                       # an order below need is impossible (d <= 2^n - 1)
            cands = cen[a].get(need, [])
            random.shuffle(cands)
            for line in cands:
                C = parse(line); tested += 1
                assert dcount(C) == a, ('census d wrong', line)
                p, kind = odd_cert(C)
                if p is None: continue
                G = union(C, threshold(m, Lm)) if m > 1 else C
                best = (G, kind); break
            if best: break
        if best:
            G, kind = best; kinds[kind] += 1; found[k] = len(G)
            n = len(G)
            edges = ' '.join(f'{u}-{v}' for u in range(n) for v in range(u + 1, n) if (G[u] >> v) & 1)
            out.write(f'base{k} O {n} {k} | {edges}\n')
    out.close()
    missing = [k for k in range(1, 512, 2) if k not in found]
    print('odd k < 512 with an odd realiser on exactly L(k) vertices:', len(found), 'of', 256)
    print('no odd realiser found on L(k) vertices:', missing)
    print('oddness certificate of the connected factor: odd transposition', kinds['T'], '; other automorphism', kinds['A'])
    # exclusion: every connected census graph with d | k, k in {1, 5, 37}, is even
    for k in (1, 5, 37):
        divs = [a for a in cen if k % a == 0]
        allg = [line for a in divs for nn in cen[a] for line in cen[a][nn]]
        par = []
        for line in allg:
            C = parse(line); assert dcount(C) == int(line.split()[2])
            auts = all_auts(C); odd = any(rev_edges(C, g) % 2 for g in auts)
            par.append((len(C), int(line.split()[2]), len(auts), 'ODD' if odd else 'even'))
        print(f'k={k}: connected census graphs with d | k:', len(allg), par)

if __name__ == '__main__':
    main()
