#!/usr/bin/env python3
"""Positive, clock-marked transport of the ideal ROOT band row.

Exact integer coefficient tables below a band cutoff; rigorous positive
majorant for the omitted tail. No eventual-root information is used.
"""
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def initial(r, b):
    # Interior point 3/2^b; multiply by 2^b to avoid floating endpoints.
    x = 3
    ans = 0
    while x <= 4 * 2**b:
        if x <= 2 * 2**b:
            ans += (0, -3, 3)[r % 3]
        else:
            ans += (-1, -1, 2)[r % 3]
        r += 1
        x *= 4
    return ans


def predecessor_residues(oldmod, A):
    mod = oldmod * 3
    a = 1 if A % 2 else 2
    j = (A-a)//2
    inverse2 = pow(2, -a, mod)
    inverse4 = pow(4, -j, mod)
    offset = ((4**j-1)//3) % mod
    old = np.arange(oldmod, dtype=np.int64)
    out = (((3*old+1)*inverse2 % mod - offset)*inverse4) % mod
    need(len(np.unique(out)) == oldmod, 'each branch is injective on old residues')
    # Independent reverse map modulo the new modulus.
    rho = (pow(4,j,mod)*out + offset) % mod
    need(np.all(rho % 3 == (2 if a == 1 else 1)), 'correct inverse odd branch')
    need(np.all(((2**a*rho-1)//3) % oldmod == old), 'literal inverse-fibre check')
    return out


def run(depth, bands):
    initial_period = [sum(abs(initial(r,b)) for r in range(3)) for b in range(6)]
    # M0(t) = integral t^b d|nu0|, evaluated at t=3/2.
    t = F(3,2)
    mgf0 = sum(F(v,4)*(t/2)**b for b,v in enumerate(initial_period))/(1-(t/2)**6)
    old = np.array([[(-1)**b*initial(r,b) for r in range(3)] for b in range(bands)], dtype=np.int64)
    rows = []
    for d in range(depth+1):
        if d:
            oldmod = old.shape[1]
            need(int(np.max(np.abs(old)))*bands < 2**63,
                 'all intermediate coefficient additions fit in int64')
            new = np.zeros((bands,oldmod*3), dtype=np.int64)
            for A in range(1,bands+1):
                residue = predecessor_residues(oldmod,A)
                new[A-1:,residue] += old[:bands-A+1,:]
            old = new
        max_abs = int(np.max(np.abs(old)))
        need(max_abs*old.shape[1] < 2**63, 'row reductions do not overflow int64')
        need(np.all(old.sum(axis=1)==0), 'every clock row has zero signed total')
        q = old.shape[1]//3
        need(np.all(old[:,:q]+old[:,q:2*q]+old[:,2*q:]==0), 'last ternary digit balance')
        mass = [int(x) for x in np.abs(old).sum(axis=1)]
        lower = sum(F(v,4*2**(d+b)) for b,v in enumerate(mass))
        tail = mgf0*(2-t)**(-d)*t**(-bands)
        upper = lower+tail
        # Compare with the quotient that discards the band clock. Integer
        # Horner evaluation keeps this diagnostic exact as well.
        marginal = old[0].astype(object)
        for band in old[1:]:marginal = 2*marginal+band
        collapsed = F(sum(abs(int(v)) for v in marginal),4*2**(d+bands-1))
        need(collapsed<=lower,'forgetting the clock can create additional cancellation')
        if d:
            need(lower <= F(rows[-1]['upper']), 'positive-kernel contraction control')
        rows.append(dict(depth=d,modulus=3**(d+1),bands=bands,
                         lower=str(lower),upper=str(upper),tail_bound=str(tail),
                         lower_decimal=float(lower),upper_decimal=float(upper),
                         clock_forgotten_lower=str(collapsed),clock_forgotten_decimal=float(collapsed),
                         max_integer_coefficient=max_abs,
                         coefficient_sha256=hashlib.sha256(old.tobytes()).hexdigest(),
                         absolute_band_sums=mass))
        print(json.dumps({k:v for k,v in rows[-1].items() if k not in ('lower','upper','tail_bound','absolute_band_sums','coefficient_sha256')},sort_keys=True),flush=True)
    return dict(status='FINITE-EXACT rational tariff brackets; positive-transport and tail proof in companion note',
                initial_absolute_period=initial_period,mgf_at_three_halves=str(mgf0),rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--depth',type=int,default=8);p.add_argument('--bands',type=int,default=104)
    p.add_argument('--output',required=True);a=p.parse_args()
    result=run(a.depth,a.bands)
    Path(a.output).write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
