#!/usr/bin/env python3
"""Exact controls: marked module covers, shared-slot divisors, CRT height."""
import itertools as it
import json
import math
import random
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'digit_memory_all_bases_20261103'))
from tournament_check import path_counts, profile, SEEDS


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def val(n, p):
    if not n:
        return math.inf
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def factors(n):
    ans = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            ans[p] = ans.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        ans[n] = 1
    return ans


def cyclic(p):
    steps = {i*i % p for i in range(1, p)} if p % 4 == 3 else set(range(1, (p+1)//2))
    return tuple(sum(1 << j for j in range(p) if (j-i) % p in steps) for i in range(p))


@lru_cache(None)
def marked_covers(rows):
    """C_j plus specified-run start/end-at-root covers, for j=0,1,2,3."""
    n = len(rows)
    h = path_counts(rows)
    covers = [[0]*4 for _ in h]
    covers[0][0] = 1
    start = [[0]*n for _ in h]
    start[1][0] = 1
    for mask in range(1, 1 << n):
        for end in range(n):
            if mask >> end & 1:
                avail = rows[end] & ~mask
                while avail:
                    b = avail & -avail
                    v = b.bit_length()-1
                    start[mask | b][v] += start[mask][end]
                    avail ^= b
        bit = mask & -mask
        sub = mask
        while sub:
            if sub & bit:
                for j in range(1, 4):
                    covers[mask][j] += h[sub]*covers[mask ^ sub][j-1]
            sub = (sub-1) & mask
    # Root-ending paths are root-starting paths in the converse.
    rev = tuple(sum(1 << j for j in range(n) if rows[j] >> i & 1) for i in range(n))
    end_dp = [[0]*n for _ in h]
    end_dp[1][0] = 1
    for mask in range(1, 1 << n):
        for end in range(n):
            avail = rev[end] & ~mask if mask >> end & 1 else 0
            while avail:
                b = avail & -avail
                end_dp[mask | b][b.bit_length()-1] += end_dp[mask][end]
                avail ^= b
    C = [math.factorial(j)*covers[-1][j] for j in range(4)]
    S, E = [0]*4, [0]*4
    full = (1 << n)-1
    for sub in range(1, 1 << n, 2):
        for j in range(1, 4):
            weight = math.factorial(j-1)*covers[full ^ sub][j-1]
            S[j] += sum(start[sub])*weight
            E[j] += sum(end_dp[sub])*weight
    return C, S, E


def coupled(blocks, split=1):
    """Literal ordinal blocks plus feedback c and root-selecting reader z."""
    N = sum(map(len, blocks))
    rows = [0]*(N+2)
    roots, offsets = [], []
    off = 0
    for b in blocks:
        roots.append(off)
        offsets.append(off)
        for v, row in enumerate(b):
            rows[off+v] |= row << off
            rows[off+v] |= ((1 << N)-1) ^ ((1 << (off+len(b)))-1)
        off += len(b)
    cut = sum(map(len, blocks[:split]))
    for v in range(N):
        if v < cut:
            rows[N] |= 1 << v
        else:
            rows[v] |= 1 << N
        if v in roots:
            rows[N+1] |= 1 << v
        else:
            rows[v] |= 1 << (N+1)
    rows[N] |= 1 << (N+1)
    return tuple(rows)


def compressed(blocks, split=1):
    """Bijective three-segment enumeration after deleting named c,z."""
    M = len(blocks)
    profiles = [marked_covers(tuple(b)) for b in blocks]
    total, accepted = 0, 0
    for order in [('c', 'z'), ('z', 'c')]:
        for masks in it.product(range(1, 8), repeat=M):
            segments = [[i for i in range(M) if masks[i] >> j & 1] for j in range(3)]
            word = segments[0]+[order[0]]+segments[1]+[order[1]]+segments[2]
            good = True
            for a, b in zip(word, word[1:]):
                if a == 'c' and isinstance(b, int):
                    good &= b < split
                elif b == 'c' and isinstance(a, int):
                    good &= a >= split
                elif a == 'z' and b == 'c':
                    good = False
            if not good:
                continue
            pos = word.index('z')
            before = word[pos-1] if pos else None
            after = word[pos+1] if pos+1 < len(word) else None
            weight = 1
            for i, mask in enumerate(masks):
                C, S, E = profiles[i]
                j = mask.bit_count()
                if after == i:
                    # If also before==i, these are distinct runs; the root
                    # cannot occur in the incoming run as well.
                    weight *= S[j]
                elif before == i:
                    weight *= C[j]-E[j]
                else:
                    weight *= C[j]
            total += weight
            accepted += bool(weight)
    return total, accepted


def divisor(contents, q):
    k = max(len(contents)-2*q, 0)
    return math.gcd(*(math.prod(s) for s in it.combinations(contents, k)))


def slot(P, r):
    return sum(math.comb(r, j)*math.factorial(j)*v for j, v in enumerate(P) if j <= r)


def ordinal(A, B):
    a, b = len(A), len(B)
    return tuple(row | (((1 << b)-1) << a) for row in A)+tuple(row << a for row in B)


def valid(rows):
    return all(not (row >> i & 1) and row < 1 << len(rows) for i,row in enumerate(rows)) and all(
        ((rows[i] >> j & 1)+(rows[j] >> i & 1)) == 1 for i,j in it.combinations(range(len(rows)),2))


def cut_polynomial(rows):
    """Exact H of every one-birth cut from two-path endpoint covers."""
    n=len(rows);full=(1<<n)-1
    def endpoints(arcs):
        dp=[[0]*n for _ in range(1<<n)]
        for v in range(n):dp[1<<v][v]=1
        for mask in range(1,1<<n):
            for u in range(n):
                if not dp[mask][u]:continue
                avail=arcs[u]&~mask
                while avail:
                    b=avail&-avail
                    dp[mask|b][b.bit_length()-1]+=dp[mask][u]
                    avail^=b
        return dp
    e=endpoints(rows)
    rev=tuple(sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n))
    s=endpoints(rev)
    K=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        for u in range(n):
            if not e[mask][u]:continue
            for v in range(n):K[u][v]+=e[mask][u]*s[full^mask][v]
    H=sum(e[full])
    L=[s[full][v]-e[full][v]+sum(K[u][v] for u in range(n)) for v in range(n)]
    Q={(u,v):-(K[u][v]+K[v][u]) for u,v in it.combinations(range(n),2)}
    g=math.gcd(H,*L,*Q.values())
    def evaluate(mask):
        return H+sum(L[v] for v in range(n) if mask>>v&1)+sum(a for (u,v),a in Q.items() if mask>>u&1 and mask>>v&1)
    return H,L,Q,g,evaluate


def add_cut(rows, mask):
    n=len(rows)
    return tuple(row | (0 if mask>>i&1 else 1<<n) for i,row in enumerate(rows))+(mask,)


def main():
    A, B = (8,1,3,6), (24,1,3,6,14)
    Pa, Pb = profile(A), profile(B)
    Pc = profile(ordinal(A, B))
    content = lambda P: math.gcd(*(math.factorial(j)*x for j, x in enumerate(P)))
    need((content(Pa), content(Pb), content(Pc)) == (1,1,5), 'emergent 5')
    need(all(slot(Pa, r) == r**4+4*r*r and slot(Pb, r) == r**5+8*r**3+2*r for r in range(11)), 'polynomials')
    census, profiles = 0, {}
    for n in range(1, 6):
        pairs = list(it.combinations(range(n), 2))
        bank = {}
        for bits in range(1 << len(pairs)):
            rows = [0]*n
            for k, (i,j) in enumerate(pairs):
                u,v = (i,j) if bits >> k & 1 else (j,i)
                rows[u] |= 1 << v
            P = profile(rows)
            bank[tuple(P)] = [P[1], slot(P, 2), content(P)]
            census += 1
        profiles[n] = sorted(bank.values())
    need(census == 1099, 'census universe')
    need(all(H % 7 for H,_,_ in profiles[5]), 'no H multiple7 through5')
    need(all(F2 % 5 for n in range(1,5) for _,F2,_ in profiles[n]), 'F2 through4')
    records = []
    small = [(2,0), cyclic(3), A, *map(tuple, SEEDS)]
    for blocks in it.product(small, repeat=2):
        graph = coupled(blocks)
        need(valid(graph), 'literal tournament')
        literal = path_counts(graph)[-1]
        compressed_H, terms = compressed(blocks)
        need(literal == compressed_H, ('compressed pair', blocks, literal, compressed_H))
        cs = [content(profile(b)) for b in blocks]
        need(literal % divisor(cs, 1) == 0, 'pair divisor')
        records.append([list(map(len, blocks)), literal, terms])
    for word in it.product((0,1), repeat=3):
        blocks = [tuple(SEEDS[i]) for i in word]
        h, terms = compressed(blocks)
        need(h == path_counts(coupled(blocks))[-1] and h % 3 == 0, ('triple', word))
        records.append([list(word), h, terms])
    modules = []
    for p in (3,5,7,11):
        rows = cyclic(p)
        P = profile(rows)
        modules.append({'p':p, 'rows':rows, 'P':P, 'H':P[1], 'content':content(P), 'marked':marked_covers(rows)})
    blocks = [tuple(m['rows']) for m in modules]
    c = [m['content'] for m in modules]
    h, terms = compressed(blocks)
    D = divisor(c,1)
    need(valid(coupled(blocks)) and c == [3,15,63,55] and D == 3 and h % D == 0, 'heterogeneous divisor')
    # Arbitrary cuts, including two births, rather than only root selectors.
    random_controls = []
    rng = random.Random(307115)
    for M,q,repeat in [(3,1,5),(4,1,3),(5,2,2)]:
        base = ()
        for _ in range(M):
            base = ordinal(base,cyclic(3))
        for _ in range(repeat):
            rows = list(base)
            masks=[]
            for __ in range(q):
                n=len(rows); mask=rng.randrange(1<<n);masks.append(mask)
                rows = [row | (0 if mask>>i&1 else 1<<n) for i,row in enumerate(rows)]+[mask]
            need(valid(rows), 'arbitrary-cut tournament')
            H=path_counts(rows)[-1]
            bound=divisor([3]*M,q)
            need(H % bound == 0, ('arbitrary port bound',M,q,masks,H,bound))
            random_controls.append([M,q,masks,H,bound])
    minor_values = [math.prod(pair) for pair in it.combinations(c,2)]
    smith_deltas = [math.gcd(*(math.prod(s) for s in it.combinations(c,k))) for k in range(5)]
    smith = [smith_deltas[k]//smith_deltas[k-1] for k in range(1,5)]
    all_prime_checks = 0
    for m in range(1,6):
        for cs in it.product((1,3,5,9,15),repeat=m):
            for q in range(4):
                k=max(m-2*q,0)
                expected=math.prod(p**sum(sorted(val(c,p) for c in cs)[:k]) for p in (3,5))
                need(divisor(cs,q)==expected, ('prime formula',cs,q))
                all_prime_checks += 1
    # Every finite CRT projection of a non-diagonal idempotent is legal.
    crt, previous, modulus = [], 0, 1
    for j in range(2,16):
        M = math.factorial(j)
        v = val(M,2)
        two, odd = 2**v, M//2**v
        r = (two*pow(two,-1,odd)) % M if odd > 1 else 0
        need(r % modulus == previous and r*r % M == r and r % two == 0 and (r-1)%odd == 0, 'CRT coherent idempotent')
        signed = r if r <= M//2 else r-M
        crt.append([j,M,r,signed])
        previous, modulus = r, M
    # Digit rules in every tested base and modulus; separate nonunit tail.
    digit_cases=0
    for base in range(2,17):
        for m in range(2,65):
            u=m
            while (g:=math.gcd(u,base))>1:
                u//=g
            s=m//u
            L=0
            while base**L % s:
                L+=1
            for n in (0,1,m-1,m,m+1,base**5+731):
                a,d=divmod(n,base)
                need((n % m == 0) == (n % base**L % s == 0 and n % u == 0), 'CRT digit split')
                if u > 1:
                    need((n % u == 0) == ((a+pow(base,-1,u)*d)%u == 0), 'trim')
                digit_cases+=1
    cut_records=[];cut_tests=0
    for rows in [cyclic(3),A,B,ordinal(A,B),ordinal(cyclic(3),cyclic(3))]:
        H,L,Q,g,ev=cut_polynomial(rows)
        n=len(rows)
        masks=range(1<<n) if n<=6 else [0,*[1<<i for i in range(n)],*[(1<<i)|(1<<j) for i,j in Q],(1<<n)-1,85,170]
        values=[]
        for mask in masks:
            actual=path_counts(add_cut(rows,mask))[-1]
            need(ev(mask)==actual, ('cut polynomial',rows,mask,ev(mask),actual))
            values.append(actual);cut_tests+=1
        need(math.gcd(*values)==g,'sample ideal equality')
        all_values=sorted({ev(mask) for mask in range(1<<n)})
        cut_records.append({'rows':rows,'H':H,'gcd_all_one_births':g,'linear':L,'quadratic':[[i,j,v] for (i,j),v in Q.items()],
                            'distinct_counts':len(all_values),'minimum':min(all_values),'maximum':max(all_values),
                            'target_membership':{v:v in all_values for v in (1,111,167,169)}})
        if H==55:
            need(len(all_values)==354 and min(all_values)==55 and max(all_values)==2513 and all(v not in all_values for v in (111,167,169)), 'H55 hostile cut bank')
    triangle_gcds=[];rows=()
    for M in range(1,6):
        rows=ordinal(rows,cyclic(3))
        H,L,Q,g,ev=cut_polynomial(rows)
        need(g==3**max(M-2,0),('triangle ideal',M,g))
        triangle_gcds.append([M,H,g])
    print(json.dumps({'status':'FINITE-EXACT controls; see PROOFS for universal scoped claims',
      'emergence':{'A':A,'B':B,'P_A':Pa,'P_B':Pb,'P_sum':Pc,'contents':[1,1,5], 'zeros_mod5':[[r for r in range(5) if slot(P,r)%5==0] for P in (Pa,Pb)]},
      'census_count':census,'small_profiles_H_F2_content':profiles,
      'literal_compressed_checks':len(records),'small_reader_records':records,'modules':modules,
      'heterogeneous':{'order':sum(map(len,blocks))+2,'H':h,'valuations_H':{p:val(h,p) for p in (2,3,5,7,11,13)},'positive_skeletons':terms,'contents':c,'K4_pair_products':minor_values,'Smith':smith,'D2':D},
      'prime_budget_checks':all_prime_checks,'arbitrary_cut_controls':random_controls,
      'CRT_idempotent':crt,'digit_rule_cases':digit_cases,
      'cut_polynomial_literal_tests':cut_tests,'cut_polynomials':cut_records,'triangle_gcd_controls':triangle_gcds},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
