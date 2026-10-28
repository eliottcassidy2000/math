#!/usr/bin/env python3
"""Exact tournament counts, not tournament-class generation (stdlib only).

The mathematics is inherited from Davis/Burnside and repository THM-514/516.
This module tests integer-kernel implementations and batching, not formula novelty.
"""
from collections import defaultdict
from math import comb, factorial, gcd


def _check_n(n):
    if type(n) is not int or n < 0:
        raise ValueError('n must be a nonnegative integer')


def odd_partition_numbers(n, minimum=1):
    _check_n(n)
    a = [1] + [0] * n
    for k in range(minimum, n + 1, 2):
        for m in range(k, n + 1):
            a[m] += a[m-k]
    return a


def core_walk(n, minimum, visit):
    """Visit all odd cores of mass <= n, with parts >= odd minimum.

    Visit signature: (mass, length, three_divisible, energy, centralizer).
    Grouping equal parts and gcd interaction are incremental. No tuple leaves
    or rational numbers are built. Each core is visited exactly once.
    """
    _check_n(n)
    if minimum < 3 or minimum % 2 == 0:
        raise ValueError('minimum must be odd and at least 3')
    parts = list(range(minimum, n+1, 2))
    table = {k: {j: gcd(k, j) for j in parts} for k in parts}
    stack = []
    nodes = 0
    gcd_terms = 0

    def rec(max_part, mass, length, c3, energy, z):
        nonlocal nodes, gcd_terms
        nodes += 1
        visit(mass, length, c3, energy, z)
        start = min(max_part, n-mass)
        start -= 1 - start % 2
        for k in range(start, minimum-1, -2):
            row = table[k]
            cross = sum(m * row[j] for j, m in stack)
            gcd_terms += len(stack)
            new_z, new_e = z, energy
            delta = cross + (k-1)//2
            for mult in range(1, (n-mass)//k+1):
                new_z *= k * mult
                new_e += delta
                delta += k
                stack.append((k, mult))
                rec(k-2, mass+k*mult, length+mult,
                    c3+(mult if k % 3 == 0 else 0), new_e, new_z)
                stack.pop()

    rec(n, 0, 0, 0, 0, 1)
    return {'core_nodes': nodes, 'gcd_group_terms': gcd_terms,
            'gcd_table_entries': len(parts)**2}


def count_direct(n, buckets=True):
    """One n: analytically fill 1-tail; optionally bucket before shifting."""
    _check_n(n)
    facts = [factorial(k) for k in range(n+1)]
    total = 0
    bins = defaultdict(int)

    def visit(m, t, c3, energy, z):
        nonlocal total
        r = n-m
        exponent = energy + r*t + r*(r-1)//2
        weight = facts[n] // (z*facts[r])
        if buckets:
            bins[exponent] += weight
        else:
            total += weight << exponent

    stats = core_walk(n, 3, visit)
    if buckets:
        total = sum(v << e for e, v in bins.items())
    answer, rem = divmod(total, facts[n])
    if rem:
        raise ArithmeticError('Burnside numerator is not divisible by n!')
    stats['exponent_buckets'] = len(bins) if buckets else None
    return answer, stats


def build_kernel(n, peel_three=False):
    """K[m,(t,c3)]=sum m!/z * 2^e over cores (c3 unused if not peeling3)."""
    _check_n(n)
    facts = [factorial(k) for k in range(n+1)]
    rows = [defaultdict(int) for _ in range(n+1)]

    def visit(m, t, c3, energy, z):
        rows[m][(t, c3 if peel_three else 0)] += (facts[m] // z) << energy

    stats = core_walk(n, 5 if peel_three else 3, visit)
    stats['kernel_cells'] = sum(map(len, rows))
    return rows, stats


def restore_unit_kernel(rows):
    """Fold the 3-tail once, retaining total mass and total cycle count."""
    N = len(rows)-1
    facts = [factorial(k) for k in range(N+1)]
    K = [defaultdict(int) for _ in range(N+1)]
    three_factors = [3**s*facts[s] for s in range(N//3+1)]
    terms = 0
    for m,row in enumerate(rows):
        for (t,c3),value in row.items():
            for s in range((N-m)//3+1):
                M = m+3*s
                coefficient = facts[M]//(facts[m]*three_factors[s])
                exponent = s*(t+2*c3+1)+3*s*(s-1)//2
                K[M][(t+s,0)] += (coefficient*value) << exponent
                terms += 1
    return K,terms


def counts_from_kernel(rows, peel_three=False, colors=1):
    """Recover the entire prefix 0..N from an integer kernel."""
    N = len(rows)-1
    if type(colors) is not int or colors < 1:
        raise ValueError('colors must be a positive literal integer')
    facts = [factorial(k) for k in range(N+1)]
    powers = [colors**k for k in range(N+1)]
    K = rows
    transfer_terms = 0
    if peel_three:
        # Merge the 3-tail into the ordinary unit-tail kernel first.  This
        # avoids a 3-tail sum for each n, and remains integer term by term.
        K,transfer_terms = restore_unit_kernel(rows)
    result = []
    for n in range(N+1):
        total = 0
        for m, row in enumerate(K[:n+1]):
            r = n-m
            binomial = comb(n,m)
            for (t, _), value in row.items():
                total += (binomial*value*powers[t+r]) << (r*t+r*(r-1)//2)
                transfer_terms += 1
        value, rem = divmod(total, facts[n])
        if rem:
            raise ArithmeticError('kernel numerator is not divisible by n!')
        result.append(value)
    return result, {'transfer_terms': transfer_terms,
                    'unit_kernel_cells': sum(map(len, K))}


def count_prefix(n, peel_three=True):
    rows, stats = build_kernel(n, peel_three)
    result, second = counts_from_kernel(rows, peel_three)
    return result, dict(stats, **second)


def palette_prefixes(n, colors):
    """Amortize an integer K3 build over many vertex-color evaluations."""
    colors = tuple(colors)
    if not colors or any(type(c) is not int or c < 1 for c in colors):
        raise ValueError('a nonempty collection of positive literal colors is required')
    rows,stats = build_kernel(n,True)
    kernel,transfer = restore_unit_kernel(rows)
    values = {c:counts_from_kernel(kernel,colors=c)[0] for c in colors}
    stats.update({'three_tail_transfer_terms':transfer,
                  'unit_kernel_cells':sum(map(len,kernel)),
                  'palette_count':len(values)})
    return values,stats


def rooted_prefix(n):
    """One distinguished vertex; A_root(n)=sum m1*2^e/z.

    Remove the marked fixed point: permutation types on n-1 remain odd, and
    every cycle chooses one direction relative to the root, giving 2^length.
    """
    _check_n(n)
    facts = [factorial(k) for k in range(n+1)]
    answer = [0]
    for order in range(1,n+1):
        total = 0
        def visit(m,t,c3,e,z):
            nonlocal total
            r = order-1-m
            total += (facts[order-1]//(z*facts[r])) << (
                e+r*t+r*(r-1)//2+t+r)
        core_walk(order-1, 3, visit)
        value, rem = divmod(total, facts[order-1])
        if rem:
            raise ArithmeticError('rooted numerator is not divisible')
        answer.append(value)
    return answer


def packed_integer_dp(bound):
    """Degree-local ordinary DP with one integer key for degree+response.

    The response quotient and odd-factorial normalization are inherited from
    Session21/22. This version only repacks their state, folds no-op positive
    filtering, and precomputes transfer tuples. No new asymptotic claim.
    """
    _check_n(bound)
    width = max(1, bound.bit_length())
    digit = (1 << width)-1
    masks = [((1 << (width*c))-1) << width
             for c in range(bound//2+2)]
    facts = [factorial(k) for k in range(bound+1)]
    twos = [(f & -f).bit_length()-1 for f in facts]
    odds = [f >> v for f, v in zip(facts, twos)]
    states = {0: 1}
    rows = []
    precomputed = 0
    for ell in range(bound-1+bound%2, 0, -2):
        shift = width*((ell-1)//2+1)
        row = sum(gcd(ell,j) << (width*(i+1))
                  for i,j in enumerate(range(1,ell,2)))
        tables = [None]*(bound+1)
        for key in states:
            d = key & digit
            if tables[d] is not None:
                continue
            entries = []
            for k in range((bound-d)//ell+1):
                nd = d+ell*k
                h = ell**k*odds[k]
                factor, rem = divmod(odds[nd],odds[d]*h)
                if rem:
                    raise ArithmeticError('nonintegral degree-local transfer')
                base = k*(ell-1)//2+ell*k*(k-1)//2-twos[k]
                count = max(0,(min(ell-2,bound-nd)+1)//2)
                entries.append((k,nd,base,factor,masks[count],k*row))
                precomputed += 1
            tables[d] = entries
        new = defaultdict(int)
        transitions = 0
        for key, weight in states.items():
            charge = (key >> shift) & digit
            for k,nd,base,factor,mask,offset in tables[key & digit]:
                new[((key+offset)&mask)|nd] += (weight*factor) << (base+k*charge)
                transitions += 1
        # Ordinary weights are positive, so no zero filtering is necessary.
        states = new
        rows.append({'cycle_length':ell,'states':len(states),'transitions':transitions})
    result = [0]*(bound+1)
    for key, weight in states.items():
        if key > digit:
            raise ArithmeticError('unfinished packed response')
        result[key], rem = divmod(weight,odds[key])
        if rem:
            raise ArithmeticError('nonintegral tournament count')
    return result, {'rows':rows,'precomputed_transfers':precomputed,
                    'transitions':sum(r['transitions'] for r in rows),
                    'peak_states':max([1]+[r['states'] for r in rows])}


def degree_bucket_dp(bound, colors=1, terminal_fold=True):
    """Session22's exact DP regrouped by degree, with positive vertex colors.

    C_n(colors) counts tournaments whose vertices carry one of `colors` colors;
    isomorphisms preserve colors. The state quotient is UNCHANGED. Storing
    degree in the containing array removes it from every dictionary key, and
    loop inversion shares one transfer across all responses of that degree.
    """
    _check_n(bound)
    if type(colors) is not int or colors < 1:
        raise ValueError('colors must be a positive literal integer')
    width = max(1,bound.bit_length())
    digit = (1 << width)-1
    masks = [(1 << (width*c))-1 for c in range(bound//2+2)]
    facts = [factorial(k) for k in range(bound+1)]
    twos = [(f & -f).bit_length()-1 for f in facts]
    odds = [f >> v for f,v in zip(facts,twos)]
    powers = [colors**k for k in range(bound+1)]
    states = [{} for _ in range(bound+1)]
    states[0][0] = 1
    rows = []
    precomputed = 0
    for ell in range(bound-1+bound%2,0,-2):
        if ell == 1 and terminal_fold:
            # Once only 1-cycles remain, response=t is the core cycle count.
            # Evaluate the unit-tail polynomial before multiplying the shared
            # degree transfer. No response dictionary is needed afterward.
            totals = [0]*(bound+1)
            transitions = 0
            for d,bucket in enumerate(states):
                if not bucket:
                    continue
                limit = bound-d
                transitions += len(bucket)*(limit+1)
                precomputed += limit+1
                for k in range(limit+1):
                    nd = d+k
                    factor,rem = divmod(odds[nd],odds[d]*odds[k])
                    if rem:
                        raise ArithmeticError('nonintegral terminal transfer')
                    polynomial = sum(weight << (k*response)
                                     for response,weight in bucket.items())
                    base = k*(k-1)//2-twos[k]
                    totals[nd] += (polynomial*factor*powers[k]) << base
            result = []
            for d,total in enumerate(totals):
                value,rem = divmod(total,odds[d])
                if rem:
                    raise ArithmeticError('nonintegral terminal count')
                result.append(value)
            rows.append({'cycle_length':1,'states':bound+1,'transitions':transitions})
            return result, {'rows':rows,'precomputed_transfers':precomputed,
                            'transitions':sum(r['transitions'] for r in rows),
                            'peak_states':max([1]+[r['states'] for r in rows])}
        shift = width*((ell-1)//2)
        row = sum(gcd(ell,j) << (width*i)
                  for i,j in enumerate(range(1,ell,2)))
        new = [defaultdict(int) for _ in range(bound+1)]
        transitions = 0
        for d,bucket in enumerate(states):
            if not bucket:
                continue
            limit = (bound-d)//ell
            transitions += len(bucket)*(limit+1)
            precomputed += limit+1
            # Zero multiplicity has identity arithmetic; only the response
            # expiration can merge states, so it must still be applied.
            mask = masks[max(0,(min(ell-2,bound-d)+1)//2)]
            target = new[d]
            for response,weight in bucket.items():
                target[response & mask] += weight
            for k in range(1,limit+1):
                nd = d+ell*k
                factor,rem = divmod(odds[nd],odds[d]*ell**k*odds[k])
                if rem:
                    raise ArithmeticError('nonintegral degree-local transfer')
                factor *= powers[k]
                base = k*(ell-1)//2+ell*k*(k-1)//2-twos[k]
                mask = masks[max(0,(min(ell-2,bound-nd)+1)//2)]
                offset = k*row
                target = new[nd]
                for response,weight in bucket.items():
                    exponent = base+k*((response >> shift)&digit)
                    target[(response+offset)&mask] += (weight*factor) << exponent
        states = new
        rows.append({'cycle_length':ell,'states':sum(map(len,states)),
                     'transitions':transitions})
    result = [0]*(bound+1)
    for d,bucket in enumerate(states):
        if bucket:
            if set(bucket) != {0}:
                raise ArithmeticError('unfinished response')
            result[d],rem = divmod(bucket[0],odds[d])
            if rem:
                raise ArithmeticError('nonintegral colored tournament count')
    return result, {'rows':rows,'precomputed_transfers':precomputed,
                    'transitions':sum(r['transitions'] for r in rows),
                    'peak_states':max([1]+[r['states'] for r in rows])}


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, nargs='?', default=50)
    parser.add_argument('--single', action='store_true')
    parser.add_argument('--one-tail', action='store_true')
    args = parser.parse_args()
    values, stats = (count_direct(args.n) if args.single else
                     count_prefix(args.n, not args.one_tail))
    print(json.dumps({'n':args.n,'values':values,'stats':stats},indent=2))
