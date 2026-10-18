#!/usr/bin/env python3
"""Exact whole-support audit of actual versus ideal ROOT coefficients.

The old matching scripts inspect only nonzero actual coefficients. This one
enumerates every primitive odd seed in the union of the two finite envelopes.
No convergence assumption or floating-point rounding enters the comparison.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from root_transport_energy_profile_20261017 import initial, step

K0 = np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2],
               [-3, 0, 3], [-1, -1, 2], [0, 0, 0]], dtype=np.int64)


def require(p, message):
    if not p:
        raise RuntimeError(message)


def ideal_step(k):
    bands, q = k.shape
    out = np.zeros((bands, 3*q), dtype=np.int64)
    r = np.arange(q, dtype=np.int64)
    for a in range(1, bands+1):
        short = 1 if a % 2 else 2
        ca = (3 + 2**short - 2**a)//3
        t = ((3*r + ca % (3*q))*pow(2**a, -1, 3*q)) % (3*q)
        out[a-1:, t] += (-1 if short == 1 else 1)*k[:bands-a+1, :]
    return out


def scalar_actual(n, depth):
    """Independent original odd-source accumulation, with its sign clock."""
    out = {}
    for m in range(1, 4*n+1, 2):
        chi = (0, 1, -1)[m % 3]
        c = -3*chi if m < 2*n else (2 if m % 3 == 2 else -1)
        s = m
        while s % 8 == 5:
            s = (s-1)//4
        for _ in range(depth):
            a = 2 if s % 8 == 1 else 1
            s = (3*s+1)//(2**a)
            c *= (-1)**a
            while s % 8 == 5:
                s = (s-1)//4
        out[s] = out.get(s, 0)+c
    return {s:c for s,c in out.items() if c}


def parent(s):
    a = 2 if s % 8 == 1 else 1
    s = (3*s+1)//(2**a)
    while s % 8 == 5:
        s = (s-1)//4
    return s


def boundary_seeds(n, depth):
    """Odd lattice points in (theta-1/3,theta+1/2], theta=3^d oddpart(N)/2^ell."""
    odd = n//(n & -n)
    top = 3**depth*odd
    out = set()
    den = 1
    while den <= 2*top:
        t = top//den
        for s in (t,t+1):
            if s>0 and s%2 and s%8!=5 and 3*(den*s-top)>-den and 2*(den*s-top)<=den:
                out.add(s)
        den *= 2
    return out


def scalar_ideal(n, d, k):
    out = {}
    top = 4*3**d*n
    for s in range(1,top//2**d+1,2):
        if s%8==5:
            continue
        b = (top//(2**d*s)).bit_length()-1
        value = int(k[b,s%k.shape[1]])
        if value:
            out[s] = value
    return np.array(list(out),dtype=np.int64),np.array(list(out.values()),dtype=np.int64)


def compare(n, d, s, c, k):
    upper = max(4*3**d*n//2**d, int(s[-1]) if len(s) else 0)
    universe = np.arange(1, upper+1, 2, dtype=np.int64)
    universe = universe[universe % 8 != 5]
    actual = np.zeros(len(universe), dtype=np.int64)
    positions = np.searchsorted(universe, s)
    require(np.array_equal(universe[positions], s), 'actual support escaped universe')
    actual[positions] = c
    # Integer comparisons avoid logarithm rounding on exact band boundaries.
    scaled = 2**d*universe
    ceiling = 4*3**d*n
    b = np.full(len(universe), -1, dtype=np.int64)
    for band in range(k.shape[0]):
        b[scaled <= ceiling] = band
        ceiling //= 2
        if ceiling == 0:
            break
    require(ceiling == 0, 'band truncation')
    ideal = np.zeros_like(actual)
    valid = b >= 0
    ideal[valid] = k[b[valid], universe[valid] % k.shape[1]]
    error = actual-ideal
    old = actual != 0
    off = universe != 1
    nz = error != 0
    omitted = (~old) & (ideal != 0)
    sq = lambda x: int(x @ x)
    maxc = max(int(np.abs(actual).max(initial=0)), int(np.abs(ideal).max(initial=0)))
    require(len(universe)*(2*maxc)**2 < 2**63, 'int64 sum safety')
    return dict(N=n, d=d, universe=len(universe), actual_support=int(old.sum()),
                full_error=sq(error), support_only_error=sq(error[old]),
                off_root_error=sq(error[off]), omitted_energy=sq(ideal[omitted]),
                omitted_seeds=universe[omitted].tolist()[:40],
                error_count=int(nz.sum()), off_root_error_count=int((nz & off).sum()),
                error_seeds=universe[nz].tolist()[:40],
                error_vector=list(zip(universe[nz].tolist()[:40],error[nz].tolist()[:40])),
                actual_energy=sq(actual), ideal_lattice_energy=sq(ideal),
                root_error=int(error[0]), full_l1=int(np.abs(error).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--depth', type=int, default=8)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    ns = [4095,4096,4097,7910,16383,16384,16385,65535,65536,65537]
    bands = (4*max(ns)*3**args.depth).bit_length()+1
    require(3*math.comb(bands+args.depth, args.depth) < 2**63, 'coefficient safety')
    k = K0[np.arange(bands) % 6]
    states = {n:initial(n) for n in ns}
    rows = []
    supports = {n:boundary_seeds(n,0) for n in ns}
    for d in range(args.depth+1):
        for n in ns:
            s,c = states[n]
            row = compare(n,d,s,c,k)
            require(set(row['error_seeds']) <= supports[n], 'sparse boundary-module support')
            require(row['error_count']<=40, 'all error seeds must be retained for this test')
            row['predicted_support'] = sorted(supports[n])
            row['predicted_support_size'] = len(supports[n])
            require(row['full_error'] == row['support_only_error']+row['omitted_energy'], 'support partition')
            rows.append(row)
            print(n,d,'full',row['full_error'],'old',row['support_only_error'],
                  'omitted',row['omitted_energy'],'seeds',row['error_count'],flush=True)
            states[n] = step(s,c)
            supports[n] = {parent(t) for t in supports[n]} | boundary_seeds(n,d+1)
        if d < args.depth:
            k = ideal_step(k)
    controls = local_controls = support_controls = 0
    small_ks = [K0[np.arange(bands)%6]]
    for d in range(1,6):
        small_ks.append(ideal_step(small_ks[-1]))
    for n in range(1,129):
        s,c = initial(n)
        ds = boundary_seeds(n,0)
        for d in range(6):
            require(dict(zip(s.tolist(),c.tolist())) == scalar_actual(n,d), 'independent source/clock control')
            row = compare(n,d,s,c,small_ks[d])
            require(row['error_count']<=40 and set(row['error_seeds'])<=ds,'small-source localization')
            support_controls += 1
            if d<5:
                si,ci = scalar_ideal(n,d,small_ks[d])
                si,ci = step(si,ci)
                local = compare(n,d+1,si,ci,small_ks[d+1])
                require(local['error_count']<=40 and set(local['error_seeds'])<=boundary_seeds(n,d+1),
                        'local forcing support')
                local_controls += 1
            controls += 1
            s,c = step(s,c)
            ds = {parent(t) for t in ds} | boundary_seeds(n,d+1)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(dict(status='FINITE-EXACT',depth=args.depth,bands=bands,
        independent_controls=controls,small_source_support_controls=support_controls,
        local_forcing_controls=local_controls,rows=rows),indent=2)+'\n')


if __name__ == '__main__':
    main()
