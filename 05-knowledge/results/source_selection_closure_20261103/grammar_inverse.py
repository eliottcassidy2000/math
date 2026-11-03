#!/usr/bin/env python3
"""Exact physical-cone suffix decoder; bounded structural controls, no H census."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
from math import gcd
import hashlib
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'two_port_return_induction_20261103'/'rooted_resets.json'
FULL = ('041', '403')
DUAL = {
'001>041':[0,0,0,1], '001>403':[0,0,1,0], '003>041':[1,0,0,0], '003>403':[0,0,0,1],
'01>041':[0,0,0,1], '01>403':[0,0,1,0], '021>041':[0,0,1,1], '021>403':[1,0,0,0],
'03>041':[1,0,0,0], '03>403':[0,0,1,0], '041>403':[1,0,0,0], '043>041':[0,0,0,1],
'043>403':[1,0,0,0], '053>041':[0,1,0,0], '053>403':[1,0,0,0], '1>041':[0,4,3,1],
'1>403':[3,0,0,1], '201>041':[0,0,0,1], '201>403':[0,0,1,0], '203>041':[1,0,0,0],
'203>403':[0,0,0,1], '21>041':[0,0,1,0], '21>403':[11,1,0,0], '221>041':[0,0,1,0],
'221>403':[0,1,0,0], '241>041':[0,1,0,0], '241>403':[0,1,0,0], '3>041':[0,0,3,1],
'3>403':[1,0,0,0], '401>041':[0,0,0,1], '401>403':[0,0,1,0], '403>041':[1,0,0,0],
'41>041':[0,1,0,0], '41>403':[0,1,0,0], '421>041':[0,0,1,0], '421>403':[0,1,0,0],
'43>041':[0,0,0,1], '43>403':[0,0,1,0], '441>041':[1,0,0,0], '441>403':[0,1,0,0],
'443>041':[0,0,0,1], '443>403':[0,0,1,0], '453>041':[0,0,0,1], '453>403':[1,0,0,0],
'503>041':[1,0,0,0], '503>403':[0,0,1,0], '53>041':[0,0,3,1], '53>403':[1,0,0,0],
'541>041':[0,1,0,0], '541>403':[0,1,0,0], '543>041':[0,0,0,1], '543>403':[0,0,1,0],
'553>041':[0,0,1,0], '553>403':[1,0,0,0]}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def transpose(A):
    return [list(row) for row in zip(*A)]


def mm(A, B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def step(x, M):
    return tuple(mm([x], M)[0])


def rref(A, coefficient_columns=None):
    A = [list(map(Q,row)) for row in A]
    pivots = []
    for j in range(len(A[0]) if coefficient_columns is None else coefficient_columns):
        r = len(pivots)
        k = next((k for k in range(r,len(A)) if A[k][j]), None)
        if k is None:
            continue
        A[r],A[k] = A[k],A[r]
        z = A[r][j]
        A[r] = [x/z for x in A[r]]
        for k in range(len(A)):
            if k != r and A[k][j]:
                z = A[k][j]
                A[k] = [x-z*y for x,y in zip(A[k],A[r])]
        pivots.append(j)
        if len(pivots) == len(A):
            break
    return A,pivots


def inverse(A):
    n = len(A)
    R,p = rref([list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(A)],n)
    need(p == list(range(n)), 'invertible matrix')
    return [row[n:] for row in R]


def determinant(A):
    A = [list(map(Q,row)) for row in A]
    value = Q(1)
    for j in range(len(A)):
        k = next((k for k in range(j,len(A)) if A[k][j]),None)
        if k is None:
            return Q(0)
        if k != j:
            A[j],A[k] = A[k],A[j]
            value = -value
        z = A[j][j]
        value *= z
        for k in range(j+1,len(A)):
            t = A[k][j]/z
            A[k] = [x-t*y for x,y in zip(A[k],A[j])]
    return value


def h(x):
    return x[0]+x[1]+2*x[2]


def physical_cone(x):
    return all(z >= 0 for z in x) and x[2] > 0


def integral(x):
    return all(Q(z).denominator == 1 for z in x)


def birth(T,root):
    n = len(T)
    return tuple(row|((1<<n) if i==root else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<root),)


def compile_source(entry, macros):
    T,slots,controller = (2,4,1),[0,1],[0,1]
    raw,roots = [],[]
    for word in (entry,)+tuple(macros):
        for character in word:
            r,discard = divmod(int(character),3)
            actual_r = controller[r]
            actual_d = controller[discard] if discard<2 else 2
            before,root = T,slots[actual_r]
            T = birth(T,root)
            raw.append(3*actual_r+actual_d)
            roots.append(root)
            if actual_d<2:
                slots[actual_d] = len(before)
        if word[-1]=='3':
            controller.reverse()
    need(slots[controller[1]] == len(T)-1, 'last interface newborn')
    need(slots[controller[0]] == root, 'last interface old root')
    return before,root,tuple(raw),tuple(roots[:-1])


def literal_profile(T,root):
    """All-subset endpoint DP, separate from the source matrix construction."""
    n = len(T)
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for mask in range(1,1<<n):
        for v,count in enumerate(dp[mask]):
            choices = T[v]&~mask
            while count and choices:
                bit = choices&-choices
                choices -= bit
                dp[mask|bit][bit.bit_length()-1] += count
    H = [sum(row) for row in dp]
    H[0] = 1
    full = (1<<n)-1
    start = [[0]*n for _ in range(1<<n)]
    start[1<<root][root] = 1
    for mask in range(1,1<<n):
        for v,count in enumerate(start[mask]):
            choices = T[v]&~mask
            while count and choices:
                bit = choices&-choices
                choices -= bit
                start[mask|bit][bit.bit_length()-1] += count
    s = sum(start[full])
    q = H[full^(1<<root)]
    F = sum(dp[mask][root]*H[full^mask] for mask in range(1,1<<n))
    return H[full],s,q,F


def coordinates(z):
    H,s,q,F = z
    return H-s-q,q-s,s,F-H-q+s


def finite_fibre(y,M):
    """All cone/integer candidates obeying the necessary first-birth bound."""
    H = int(h(y))
    R,p = rref([row+[value] for row,value in zip(transpose(M),y)],4)
    if any(all(x==0 for x in row[:4]) and row[4] for row in R):
        return
    free = [i for i in range(4) if i not in p]
    for values in product(range(H+1),repeat=len(free)):
        x = [Q(0)]*4
        for i,value in zip(free,values):
            x[i] = Q(value)
        for i,j in enumerate(p):
            x[j] = R[i][4]-sum(R[i][k]*x[k] for k in free)
        if not integral(x) or not physical_cone(x):
            continue
        x = tuple(map(int,x))
        if 2*x[0]+3*x[1]+3*x[2]+x[3] <= H:
            need(step(x,M)==tuple(y), 'exact inverse-fibre candidate')
            yield x


def main():
    source = json.loads(SOURCE.read_text())
    P = [[1,0,0,0],[1,1,0,0],[1,0,1,0],[0,0,0,1]]
    Pi = inverse(P)
    matrices = {''.join(map(str,z['word'])):mm(mm(P,z['matrix']),Pi)
                for z in source['grammar']['alphabet']}
    entries = {tuple(map(int,step(z['profile'],Pi))):''.join(map(str,z['entry']))
               for z in source['grammar']['saturations'][1]['entry_states']}
    inverses = {w:inverse(matrices[w]) for w in FULL}
    need(len(entries)==21 and len(matrices)==28, 'inherited grammar size')
    need(all(x>=0 for M in matrices.values() for row in M for x in row), 'all matrices positive in physical cone')
    need([determinant(matrices[w]) for w in FULL]==[4,-4], 'full-rank determinants')
    certificate = []
    for target in FULL:
        for other,M in matrices.items():
            if other==target:
                continue
            lam = DUAL[other+'>'+target]
            value = tuple(mm(mm(M,inverses[target]),[[z] for z in lam])[i][0] for i in range(4))
            need(all(z>=0 for z in lam) and all(z<=0 for z in value) and value[2]<0,
                 'strict physical last-letter separation')
            certificate.append({'other':other,'target':target,'dual':lam,'contraction':list(map(str,value))})
    need(len(certificate)==54, 'complete ordered separation universe')
    hc = [1,1,2,0]
    for word,M in matrices.items():
        out = [row[0]+row[1]+2*row[2] for row in M]
        need(all(2*out[i]>=3*hc[i] for i in range(4)), 'all macros grow H by at least3/2')
        need(all(M[i][3]<=3*out[i] for i in range(4)), 'delta output at most3H')
    for word,factor in zip(FULL,(7,6)):
        need(all(matrices[word][i][0]+matrices[word][i][1]+2*matrices[word][i][2]>=factor*hc[i]
                 for i in range(4)), 'full-rank H growth')
    def peel(y):
        found = []
        for word,I in inverses.items():
            x = step(y,I)
            if physical_cone(x):
                found.append((word,x))
        need(len(found)<=1, 'at most one full-rank predecessor cone')
        return found[0] if found else None
    def binary_decode(y):
        reverse = []
        while y not in entries:
            result = peel(y)
            if result is None or not integral(result[1]):
                return None
            word,x = result
            need(h(x)<h(y), 'inverse descent')
            reverse.append(word)
            y = tuple(map(int,x))
        return entries[y],tuple(reversed(reverse))
    literal = []
    seen = set()
    for x,entry in entries.items():
        need(physical_cone(x), 'strong entry cone')
        need(x[3]<=3*h(x), 'entry defect bound')
        for length in range(3):
            for words in product(FULL,repeat=length):
                y = x
                for word in words:
                    y = step(y,matrices[word])
                need(y not in seen, 'unique entry/full-rank word')
                seen.add(y)
                need(binary_decode(y)==(entry,words), 'exact binary inverse')
                T,root,raw,roots = compile_source(entry,words)
                actual = coordinates(literal_profile(T,root))
                need(actual==y, 'literal graph rooted profile')
                need(gcd(*map(int,y))==gcd(*x), 'odd content preserved')
                literal.append({'entry':entry,'macros':words,'x':list(map(int,y)),
                                'rows':T,'root':root,'raw_with_final_birth':raw,'parent_roots':roots})
    @lru_cache(None)
    def membership(y):
        if not integral(y) or not physical_cone(y):
            return None
        y = tuple(map(int,y))
        if y in entries:
            return entries[y],()
        result = peel(y)
        if result is not None:
            word,x = result
            if not integral(x):
                return None
            previous = membership(tuple(map(int,x)))
            return None if previous is None else (previous[0],previous[1]+(word,))
        for word,M in matrices.items():
            if word in FULL:
                continue
            for x in finite_fibre(y,M):
                previous = membership(x)
                if previous is not None:
                    return previous[0],previous[1]+(word,)
        return None
    base = (1,0,1,0)
    fibres = []
    for word,M in matrices.items():
        y = tuple(map(int,step(base,M)))
        result = membership(y)
        need(result is not None, 'every one-macro actual target accepted')
        T,root,_,_ = compile_source(*result)
        need(coordinates(literal_profile(T,root))==y, 'inverse witness literal replay')
        if word not in FULL:
            candidates = list(finite_fibre(y,M))
            need(base in candidates, 'true parent in complete bounded singular fibre')
            fibres.append({'word':word,'target':y,'candidate_count':len(candidates)})
    boundary = []
    for x,entry in entries.items():
        for word,I in inverses.items():
            predecessor = step(x,I)
            if all(z>=0 for z in predecessor):
                boundary.append({'entry':entry,'word':word,'weak_predecessor':list(map(int,predecessor))})
                need(tuple(predecessor)==(0,1,0,0), 'sole nonnegative entry predecessor is weak TT2 sink')
    need(len(boundary)==2, 'both entry boundary controls')
    need(membership((0,1,0,0)) is None, 'weak source is not grammar entry')
    print(json.dumps({
        'scope':'structural cone and actual-program controls; no target-support saturation',
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'coordinate_order':['h-s-q','q-s','s','F-h-q+s'],
        'matrices':{w:[[int(x) for x in row] for row in M] for w,M in matrices.items()},
        'fullrank_inverse':{w:[[str(x) for x in row] for row in I] for w,I in inverses.items()},
        'separation_certificates':certificate,'entry_count':21,
        'primitive_entry_count':sum(gcd(*x)==1 for x in entries),
        'literal_binary_programs':literal,'maximum_literal_order':max(len(z['rows']) for z in literal),
        'singular_fibre_controls':fibres,'membership_one_macro_controls':28,
        'entry_weak_boundaries':boundary,'membership_cache':membership.cache_info()._asdict(),
    },indent=2,sort_keys=True))


if __name__=='__main__':
    main()
