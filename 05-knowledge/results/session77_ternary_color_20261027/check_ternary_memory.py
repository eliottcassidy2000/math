#!/usr/bin/env python3
"""Exact path-cover composition and scalar-memory hostile controls.

Standalone standard-library computation. No external files or downloads.
Run with ordinary Python or python -O; all checks use explicit exceptions.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def tournament(n, code):
    require(n >= 1 and 0 <= code < (1 << (n*(n-1)//2)), 'bad tournament code')
    out = [0]*n
    for k, (i,j) in enumerate(combinations(range(n), 2)):
        if (code >> k) & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


def path_cover_profile(out):
    """Subset Hamiltonian DP, then root-block set-partition recurrence."""
    n = len(out)
    require(all(not (out[v] >> v & 1) for v in range(n)), 'loop')
    require(all(((out[v] >> w)&1) + ((out[w] >> v)&1) == 1
                for v in range(n) for w in range(v)), 'not a tournament')
    dp = [[0]*n for _ in range(1 << n)]
    h = [0]*(1 << n)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        h[mask] = sum(dp[mask])
        for v in range(n):
            rem = out[v] & ~mask
            while rem:
                bit = rem & -rem
                rem -= bit
                dp[mask | bit][bit.bit_length()-1] += dp[mask][v]
    pc = [[0]*(n+1) for _ in range(1 << n)]
    pc[0][0] = 1
    for mask in range(1, 1 << n):
        first = mask & -mask
        sub = mask
        while sub:
            if sub & first:
                for d in range(1, n+1):
                    pc[mask][d] += h[sub]*pc[mask ^ sub][d-1]
            sub = (sub-1) & mask
    return pc[-1][1:]


def permutation_cut_profile(out):
    """Independent ordered-permutation-and-cut enumeration, divided by d!."""
    n = len(out)
    ordered = [0]*(n+1)
    for word in permutations(range(n)):
        compulsory = 0
        for j in range(n-1):
            if not (out[word[j]] >> word[j+1] & 1):
                compulsory |= 1 << j
        # Optional cuts can occur at every other position.
        free = ((1 << (n-1))-1) ^ compulsory
        subset = free
        while True:
            cuts = compulsory | subset
            ordered[cuts.bit_count()+1] += 1
            if subset == 0:
                break
            subset = (subset-1) & free
    answer = []
    for d in range(1,n+1):
        require(ordered[d] % factorial(d) == 0, 'decorated free-action failure')
        answer.append(ordered[d]//factorial(d))
    return answer


def hamiltonian_permutations(out):
    n = len(out)
    return sum(all(out[w[j]] >> w[j+1] & 1 for j in range(n-1))
               for w in permutations(range(n)))


def substitute(quotient, blocks):
    require(len(quotient) == len(blocks), 'quotient/block mismatch')
    offsets = [0]
    for block in blocks:
        require(len(block) > 0, 'empty block')
        offsets.append(offsets[-1]+len(block))
    out = [0]*offsets[-1]
    for i, block in enumerate(blocks):
        for u in range(len(block)):
            out[offsets[i]+u] = block[u] << offsets[i]
            for j, other in enumerate(blocks):
                if quotient[i] >> j & 1:
                    out[offsets[i]+u] |= ((1 << len(other))-1) << offsets[j]
    return out


def walk_series(quotient, bounds):
    """Bounded coefficient DP retaining the final quotient vertex."""
    q = len(quotient)
    states = {}
    for i in range(q):
        c = tuple(int(j == i) for j in range(q))
        if bounds[i]:
            states[(c,i)] = 1
    total = defaultdict(int)
    for length in range(1, sum(bounds)+1):
        next_states = defaultdict(int)
        for (c,last), count in states.items():
            total[c] += count
            for j in range(q):
                if quotient[last] >> j & 1 and c[j] < bounds[j]:
                    new = list(c)
                    new[j] += 1
                    next_states[(tuple(new),j)] += count
        states = next_states
    return dict(total)


def multiply(left, right, bounds):
    result = defaultdict(int)
    for a, va in left.items():
        for b, vb in right.items():
            c = tuple(x+y for x,y in zip(a,b))
            if all(c[i] <= bounds[i] for i in range(len(bounds))):
                result[c] += va*vb
    return dict(result)


def profile_from_kernel(quotient, profiles):
    bounds = tuple(map(len,profiles))
    W = walk_series(quotient,bounds)
    power = {(0,)*len(bounds): 1}
    result = []
    checked = 0
    for d in range(1,sum(bounds)+1):
        power = multiply(power,W,bounds)
        total = 0
        for c, coeff in power.items():
            if not all(c):
                continue
            numerator = coeff
            factor = 1
            for i,r in enumerate(c):
                numerator *= factorial(r)
                factor *= profiles[i][r-1]
            require(numerator % factorial(d) == 0, 'nonintegral decorated kernel')
            total += numerator//factorial(d)*factor
            checked += 1
        result.append(total)
    return result, checked


def is_strong(out):
    full=(1 << len(out))-1
    for start in range(len(out)):
        reached=1 << start
        while True:
            new=reached
            for v in range(len(out)):
                if reached >> v & 1:
                    new |= out[v]
            if new == reached:
                break
            reached=new
        if reached != full:
            return False
    return True


def automorphism_order(out):
    n=len(out)
    return sum(all(((out[i] >> j)&1) == ((out[p[i]] >> p[j])&1)
                   for i in range(n) for j in range(n))
               for p in permutations(range(n)))


def c3_coeff(c):
    if len(set(c)) == 1:
        return 3
    return int(max(c)-min(c) == 1)


def main():
    c3 = [2,4,1]
    A, B = tournament(5,69), tournament(5,76)
    pA, pB = path_cover_profile(A), path_cover_profile(B)
    require(pA == [15,39,33,10,1], 'A profile')
    require(pB == [15,45,35,10,1], 'B profile')
    require(pA == permutation_cut_profile(A), 'A independent profile')
    require(pB == permutation_cut_profile(B), 'B independent profile')
    require(pA[0] == pB[0] and pA != pB, 'scalar hostile vanished')
    require(is_strong(A) and is_strong(B), 'factor strongness')
    require(automorphism_order(A) == 3 and automorphism_order(B) == 5, 'factor automorphisms')
    nested = substitute(c3,[c3,[0],[0]])
    relabel = [0,3,2,1,4]
    require(all(((nested[i] >> j)&1) == ((A[relabel[i]] >> relabel[j])&1)
                for i in range(5) for j in range(5)), 'nested C3 identification')
    regular5 = [sum(1 << ((v+d)%5) for d in (1,2)) for v in range(5)]
    regular_maps = [p for p in permutations(range(5))
                    if all(((regular5[i] >> j)&1) == ((B[p[i]] >> p[j])&1)
                           for i in range(5) for j in range(5))]
    require(len(regular_maps) == 5, 'regular-five identification')
    contexts = [substitute(c3,[T,[0],[0]]) for T in (A,B)]
    require(all(is_strong(T) for T in contexts), 'context strongness')
    contextual_profiles = []
    for T, p in zip(contexts,(pA,pB)):
        direct = path_cover_profile(T)
        independent = permutation_cut_profile(T)
        transferred, _ = profile_from_kernel(c3,[p,[1],[1]])
        require(direct == independent == transferred, 'context profile mismatch')
        require(direct[0] == 3*p[0]+2*p[1], 'one-block context identity')
        contextual_profiles.append(direct)
    require([p[0] for p in contextual_profiles] == [123,135], 'context counts')

    # Exhaustive minimal-order search: every labelled tournament on orders 1..4.
    small_count = 0
    for n in range(1,5):
        by_h = {}
        for code in range(1 << (n*(n-1)//2)):
            p = path_cover_profile(tournament(n,code))
            require(p == permutation_cut_profile(tournament(n,code)), 'small profile')
            require(p[0] not in by_h or by_h[p[0]] == p, 'smaller witness found')
            by_h[p[0]] = p
            small_count += 1

    W = walk_series(c3,(6,6,6))
    for c in product(range(1,7),repeat=3):
        require(W.get(c,0) == c3_coeff(c), 'C3 balanced-shell kernel')
    blocks = [[0],tournament(2,1),tournament(3,7),c3]
    profiles = [path_cover_profile(T) for T in blocks]
    complete_cases = cells = 0
    for ids in product(range(4),repeat=3):
        bs = [blocks[i] for i in ids]
        transferred, checked = profile_from_kernel(c3,[profiles[i] for i in ids])
        require(transferred == path_cover_profile(substitute(c3,bs)), 'C3 full-profile transfer')
        complete_cases += 1
        cells += checked
    # General quotient controls include all labelled q=3 tournaments and
    # all singleton/edge block choices: 64 complete profiles, not just H.
    general_cases = 0
    for code in range(8):
        Q = tournament(3,code)
        for ids in product(range(2),repeat=3):
            bs = [blocks[i] for i in ids]
            transferred, checked = profile_from_kernel(Q,[profiles[i] for i in ids])
            require(transferred == path_cover_profile(substitute(Q,bs)), 'general quotient transfer')
            general_cases += 1
            cells += checked
    tower = substitute(c3,[c3,c3,c3])
    tower_profile = path_cover_profile(tower)
    require(tower_profile[0] == hamiltonian_permutations(tower) == 3159, 'C3 tower')

    # Bare quotient words need not be distinct. In W_C3^2, x^2 has coefficient1;
    # dividing by2! before assigning two actual singleton components gives1/2.
    Wsmall = walk_series(c3,(2,1,1))
    square = multiply(Wsmall,Wsmall,(2,1,1))
    require(Fraction(square[(2,0,0)],factorial(2)) == Fraction(1,2), 'bare-word hostile')
    require(Fraction(square[(2,0,0)]*factorial(2),factorial(2)) == 1, 'decorated-word repair')
    fourth = multiply(multiply(square,Wsmall,(2,1,1)),Wsmall,(2,1,1))
    require(fourth[(2,1,1)] == 12, 'all-positive singleton walk tuples')
    require(Fraction(fourth[(2,1,1)],factorial(4)) == Fraction(1,2), 'all-positive bare-word hostile')
    require(Fraction(fourth[(2,1,1)]*factorial(2),factorial(4)) == 1, 'all-positive decoration repair')
    # Deliberately false scalar rule and shell multiplicity both fail.
    rejected = 0
    for cond in (contextual_profiles[0][0] == contextual_profiles[1][0],
                 W[(1,1,1)] == 1,
                 Fraction(square[(2,0,0)],2).denominator == 1):
        try:
            require(cond,'intentional hostile')
        except ValueError:
            rejected += 1
    require(rejected == 3,'hostile rejection count')
    report = {
        'status':'PASS',
        'convention':'bit k=1 means i->j for lexicographic pairs i<j',
        'witness_order':5,
        'A_code':69,'B_code':76,
        'A_C3_nested_relabelling':relabel,
        'B_regular5_relabelling':regular_maps[0],
        'A_outdegrees':[row.bit_count() for row in A],
        'B_outdegrees':[row.bit_count() for row in B],
        'A_profile':pA,'B_profile':pB,
        'A_arcs':[(i,j) for i in range(5) for j in range(5) if A[i] >> j & 1],
        'B_arcs':[(i,j) for i in range(5) for j in range(5) if B[i] >> j & 1],
        'A_strong':is_strong(A),'B_strong':is_strong(B),
        'A_automorphism_order':automorphism_order(A),
        'B_automorphism_order':automorphism_order(B),
        'contexts_strong':[is_strong(T) for T in contexts],
        'context':'C3[T,K1,K1]',
        'context_profiles':contextual_profiles,
        'small_tournaments_independently_checked':small_count,
        'C3_kernel_cells':216,
        'C3_complete_profile_cases':complete_cases,
        'general_complete_profile_cases':general_cases,
        'integral_kernel_cells':cells,
        'C3_tower_profile':tower_profile,
        'C3_tower_H_independent_permutations':3159,
        'bare_word_division':'1/2 before decoration; 1 after factorial assignment',
        'all_positive_bare_word_coefficient':12,
        'explicit_hostile_rejections':rejected,
        'scope':'Exact scalar-H failure and inherited path-cover substitution; no QCD identification or complexity collapse.'
    }
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
