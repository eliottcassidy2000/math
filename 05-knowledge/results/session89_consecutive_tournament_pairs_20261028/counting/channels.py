#!/usr/bin/env python3
"""Exact positive tournament-layer channels on one shared response graph.

F_n(r,s) = sum_(odd partitions lambda of n) 2^(r E(lambda)+s length(lambda))/z_lambda.
Requires r >= 1, s >= 0: these are counts, not tournament representatives.
The quotient/arithmetic are inherited from Sessions 21/22, loop layout and
terminal folding from Session 87. New here: independent channel weights share
one traversal and one transfer calculation. Standard library only.
"""
from collections import defaultdict
from math import factorial, gcd


def _validate(bound, channels):
    if type(bound) is not int or bound < 0:
        raise ValueError('nonnegative literal integer bound required')
    channels = tuple(channels)
    if not channels or any(not isinstance(c, (tuple, list)) or len(c) != 2 or
            type(c[0]) is not int or c[0] < 1 or
            type(c[1]) is not int or c[1] < 0 for c in channels):
        raise ValueError('nonempty channels (r,s), literal integers r>=1, s>=0 required')
    return tuple(tuple(c) for c in channels)


def _setup(bound):
    width = max(1, bound.bit_length())
    facts = [factorial(k) for k in range(bound+1)]
    twos = [(f & -f).bit_length()-1 for f in facts]
    return (width, (1 << width)-1,
            [(1 << (width*c))-1 for c in range(bound//2+2)],
            twos, [f >> v for f,v in zip(facts,twos)])


def _finish(totals, odds):
    result = []
    for d,total in enumerate(totals):
        value,rem = divmod(total,odds[d])
        if rem:
            raise ArithmeticError('nonintegral final channel count')
        result.append(value)
    return result


def _stats(rows, precomputed, q):
    return {'rows':rows, 'precomputed_transfers':precomputed,
            'transitions':sum(r['transitions'] for r in rows),
            'peak_states':max([1]+[r['states'] for r in rows]),
            'channels':q}


def scalar_prefix(bound, r=1, s=0, terminal_fold=True):
    """Strong S87 scalar baseline, minimally generalized from colors to r,s.

    r=1 agrees with S87 degree_bucket_dp(bound, colors=2**s), stage for stage.
    """
    ((r,s),) = _validate(bound, [(r,s)])
    width,digit,masks,twos,odds = _setup(bound)
    states = [{} for _ in range(bound+1)]
    states[0][0] = 1
    rows = []
    precomputed = 0
    for ell in range(bound-1+bound%2,0,-2):
        if ell == 1 and terminal_fold:
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
                    polynomial = sum(weight << (r*k*response)
                                     for response,weight in bucket.items())
                    base = r*k*(k-1)//2-twos[k]+s*k
                    totals[nd] += (polynomial*factor) << base
            rows.append({'cycle_length':1,'states':bound+1,'transitions':transitions})
            return _finish(totals, odds), _stats(rows,precomputed,1)
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
            mask = masks[max(0,(min(ell-2,bound-d)+1)//2)]
            target = new[d]
            for response,weight in bucket.items():
                target[response & mask] += weight
            for k in range(1,limit+1):
                nd = d+ell*k
                factor,rem = divmod(odds[nd],odds[d]*ell**k*odds[k])
                if rem:
                    raise ArithmeticError('nonintegral degree-local transfer')
                base = r*(k*(ell-1)//2+ell*k*(k-1)//2)-twos[k]+s*k
                mask = masks[max(0,(min(ell-2,bound-nd)+1)//2)]
                offset = k*row
                target = new[nd]
                for response,weight in bucket.items():
                    exponent = base+r*k*((response >> shift)&digit)
                    target[(response+offset)&mask] += (weight*factor) << exponent
        states = new
        rows.append({'cycle_length':ell,'states':sum(map(len,states)),
                     'transitions':transitions})
    totals = []
    for bucket in states:
        if bucket and set(bucket) != {0}:
            raise ArithmeticError('unfinished response')
        totals.append(bucket.get(0,0))
    return _finish(totals,odds),_stats(rows,precomputed,1)


def fused_prefix(bound, channels=((1,0),(2,1))):
    """Return one full prefix per channel, sharing response/transfer metadata.

    Channel values are independent exact integers, never squared after merging.
    Duplicate channels are allowed and preserve caller order. Positive weights
    make the response graph identical for every supported channel.
    """
    channels = _validate(bound,channels)
    if len(channels) == 2:
        return _fused_pair(bound,channels)
    q = len(channels)
    width,digit,masks,twos,odds = _setup(bound)
    states = [{} for _ in range(bound+1)]
    states[0][0] = [1]*q
    rows = []
    precomputed = 0
    for ell in range(bound-1+bound%2,0,-2):
        if ell == 1:
            totals = [[0]*(bound+1) for _ in channels]
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
                    for i,(r,s) in enumerate(channels):
                        polynomial = sum(weights[i] << (r*k*response)
                                         for response,weights in bucket.items())
                        base = r*k*(k-1)//2-twos[k]+s*k
                        totals[i][nd] += (polynomial*factor) << base
            rows.append({'cycle_length':1,'states':bound+1,'transitions':transitions})
            return [_finish(t,odds) for t in totals],_stats(rows,precomputed,q)
        shift = width*((ell-1)//2)
        row = sum(gcd(ell,j) << (width*i)
                  for i,j in enumerate(range(1,ell,2)))
        new = [{} for _ in range(bound+1)]
        transitions = 0
        for d,bucket in enumerate(states):
            if not bucket:
                continue
            limit = (bound-d)//ell
            transitions += len(bucket)*(limit+1)
            precomputed += limit+1
            mask = masks[max(0,(min(ell-2,bound-d)+1)//2)]
            target = new[d]
            for response,weights in bucket.items():
                key = response & mask
                values = target.get(key)
                if values is None:
                    target[key] = weights.copy()
                else:
                    for i,w in enumerate(weights):
                        values[i] += w
            for k in range(1,limit+1):
                nd = d+ell*k
                factor,rem = divmod(odds[nd],odds[d]*ell**k*odds[k])
                if rem:
                    raise ArithmeticError('nonintegral degree-local transfer')
                inside = k*(ell-1)//2+ell*k*(k-1)//2
                bases = [(r*inside-twos[k]+s*k,r*k) for r,s in channels]
                mask = masks[max(0,(min(ell-2,bound-nd)+1)//2)]
                offset = k*row
                target = new[nd]
                for response,weights in bucket.items():
                    key = (response+offset)&mask
                    values = target.get(key)
                    charge = (response >> shift)&digit
                    if values is None:
                        target[key] = [(w*factor) << (base+rk*charge)
                                      for w,(base,rk) in zip(weights,bases)]
                    else:
                        for i,(base,rk) in enumerate(bases):
                            values[i] += (weights[i]*factor) << (base+rk*charge)
        states = new
        rows.append({'cycle_length':ell,'states':sum(map(len,states)),
                     'transitions':transitions})
    return [[1] for _ in channels],_stats(rows,precomputed,q)


def _fused_pair(bound, channels):
    """Two-channel specialization avoids a Python channel-loop per edge."""
    (r0,s0),(r1,s1) = channels
    width,digit,masks,twos,odds = _setup(bound)
    states = [{} for _ in range(bound+1)]
    states[0][0] = (1,1)
    rows = []
    precomputed = 0
    for ell in range(bound-1+bound%2,0,-2):
        if ell == 1:
            t0,t1 = [0]*(bound+1),[0]*(bound+1)
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
                    p0=p1=0
                    kr0,kr1=k*r0,k*r1
                    for response,(w0,w1) in bucket.items():
                        p0 += w0 << (kr0*response)
                        p1 += w1 << (kr1*response)
                    inside = k*(k-1)//2
                    t0[nd] += (p0*factor) << (r0*inside-twos[k]+s0*k)
                    t1[nd] += (p1*factor) << (r1*inside-twos[k]+s1*k)
            rows.append({'cycle_length':1,'states':bound+1,'transitions':transitions})
            return [_finish(t0,odds),_finish(t1,odds)],_stats(rows,precomputed,2)
        shift = width*((ell-1)//2)
        row = sum(gcd(ell,j) << (width*i)
                  for i,j in enumerate(range(1,ell,2)))
        new = [{} for _ in range(bound+1)]
        transitions = 0
        for d,bucket in enumerate(states):
            if not bucket:
                continue
            limit = (bound-d)//ell
            transitions += len(bucket)*(limit+1)
            precomputed += limit+1
            mask = masks[max(0,(min(ell-2,bound-d)+1)//2)]
            target = new[d]
            for response,(w0,w1) in bucket.items():
                key = response & mask
                v0,v1 = target.get(key,(0,0))
                target[key] = (v0+w0,v1+w1)
            for k in range(1,limit+1):
                nd = d+ell*k
                factor,rem = divmod(odds[nd],odds[d]*ell**k*odds[k])
                if rem:
                    raise ArithmeticError('nonintegral degree-local transfer')
                inside = k*(ell-1)//2+ell*k*(k-1)//2
                b0,b1 = r0*inside-twos[k]+s0*k,r1*inside-twos[k]+s1*k
                kr0,kr1 = k*r0,k*r1
                mask = masks[max(0,(min(ell-2,bound-nd)+1)//2)]
                offset = k*row
                target = new[nd]
                for response,(w0,w1) in bucket.items():
                    key = (response+offset)&mask
                    v0,v1 = target.get(key,(0,0))
                    charge = (response >> shift)&digit
                    target[key] = (v0+((w0*factor) << (b0+kr0*charge)),
                                   v1+((w1*factor) << (b1+kr1*charge)))
        states = new
        rows.append({'cycle_length':ell,'states':sum(map(len,states)),
                     'transitions':transitions})
    return [[1],[1]],_stats(rows,precomputed,2)
