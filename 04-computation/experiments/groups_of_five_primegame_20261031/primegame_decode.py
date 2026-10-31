#!/usr/bin/env python3
"""mac-mini-2026-10-31-fivegroups: exact decoding of PRIMEGAME over the first ten primes.

Two published lists (both run in procgen_atlas_20260924_primegame.py; opus-2026-10-31-S70 HYP-9448 / THM-4852
treats their equivalence and the phase law):
  CONWAY (Conway 1987; Guy 1983):     ..., 13/11, 15/14, 15/2, 55/1
  OWNER  (the owner's list of 2026-10-31, = the Wikipedia form): ..., 13/11, 15/2, 1/7, 55/1

Prints, for both lists:
  (1) every fraction factored over 2..29, its exponent change vector and its support;
  (2) the register/state split, with an inductive check that at most one state prime is ever present;
  (3) the machine table (state, tests, register moves, next state) = a 4-register Minsky machine;
  (4) the run from 2: the configurations (exponent vectors) through the emission of 2^2, 2^3, 2^5, and the step
      counts at which 2^p appears for the first primes (compared with OEIS A007546 and A007547);
  (5) the residue-class (generalized Collatz / RCWA) form: modulus lcm(denominators), and the exact density with
      which each fraction is selected on residues;
  (6) list- and label-dependence of the owner's {1/7, 1/17, 55/1}.

usage: python3 primegame_decode.py
"""
import itertools
import math
from fractions import Fraction

PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
REG = [2, 3, 5, 7]
STA = [11, 13, 17, 19, 23, 29]
HEAD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11)]
CONWAY = HEAD + [(15, 14), (15, 2), (55, 1)]
OWNER = HEAD + [(15, 2), (1, 7), (55, 1)]
# OEIS (JSON API, retrieved 2026-10-31): steps to produce 2^(p_k), first ten terms.
# A007546 ("fast version"; its comment names the ..., 15/14, 15/2, 55/1 program); A007547 ("slow version").
A007546_HEAD = [19, 69, 280, 707, 2363, 3876, 8068, 11319, 19201, 36866]
A007547_HEAD = [19, 69, 281, 710, 2375, 3893, 8102, 11361, 19268, 36981]


def fac(x):
    v = [0] * 10
    for i, p in enumerate(PR):
        while x % p == 0:
            x //= p
            v[i] += 1
    assert x == 1, 'not 29-smooth'
    return v


def fstr(v):
    parts = []
    for p, e in zip(PR, v):
        if e == 1:
            parts.append(str(p))
        elif e > 1:
            parts.append(f'{p}^{e}')
    return '*'.join(parts) if parts else '1'


def cfg(v):
    return fstr(v)


def step(prog, v):
    for j, (a, b) in enumerate(prog):
        den = fac(b)
        if all(v[i] >= den[i] for i in range(10)):
            num = fac(a)
            return j, [v[i] - den[i] + num[i] for i in range(10)]
    return None, v


def section(t):
    print()
    print('=' * 100)
    print(t)
    print('=' * 100)


def factor_table(name, prog):
    section(f'(1) {name}: the fractions over the first ten primes 2..29')
    print('  #   fraction    numerator        denominator      change vector (2,3,5,7,11,13,17,19,23,29)   support')
    used = set()
    for j, (a, b) in enumerate(prog):
        na, nb = fac(a), fac(b)
        ch = [na[i] - nb[i] for i in range(10)]
        sup = [PR[i] for i in range(10) if na[i] or nb[i]]
        used.update(sup)
        print(f'  {j + 1:<3} {a:>3}/{b:<5}  {fstr(na):16} {fstr(nb):16} {str(ch):45} {sup}')
    print(f'  primes used: {sorted(used)}  = the first ten primes: {sorted(used) == PR}')
    lcm = 1
    for a, b in prog:
        lcm = lcm * b // math.gcd(lcm, b)
    print(f'  lcm of denominators = {lcm} = 29# (primorial) : {lcm == math.prod(PR)}')
    triv = [f'{a}/{b}' for a, b in prog if a == 1 or b == 1]
    print(f'  fractions with a trivial side (numerator or denominator 1): {triv}  (count {len(triv)})')
    return triv


def state_split(name, prog):
    section(f'(2) {name}: register/state split and the one-hot state invariant')
    # every numerator and denominator carries at most one state prime, to exponent at most 1
    ok = True
    for a, b in prog:
        for x in (a, b):
            v = fac(x)
            if sum(v[4:]) > 1:
                ok = False
    print(f'  every numerator and denominator has at most one state prime (11..29), exponent <= 1: {ok}')
    # each state prime s has an exit fraction with denominator exactly s, before every no-state fraction
    nostate_idx = [j for j, (a, b) in enumerate(prog) if sum(fac(b)[4:]) == 0]
    exits = {}
    for s in STA:
        for j, (a, b) in enumerate(prog):
            if b == s:
                exits[s] = j
                break
    first_nostate = min(nostate_idx)
    ok2 = all(s in exits and exits[s] < first_nostate for s in STA)
    print(f'  every state s has an unconditional exit s -> (numerator) listed before the first no-state fraction '
          f'(#{first_nostate + 1}): {ok2}  exits: {[(s, prog[exits[s]]) for s in STA]}')
    adds = [(a, b) for a, b in prog if sum(fac(a)[4:]) == 1 and sum(fac(b)[4:]) == 0]
    print(f'  fractions that add a state without removing one: {adds} (no-state fractions, so they act only when no '
          f'state is present)')
    print('  => by induction from 2 (no state), the state part of every configuration is zero or one prime to the '
          'first power (PROVED by the two checks above).')
    return ok and ok2


def machine_table(name, prog):
    section(f'(3) {name}: the machine table (a 4-register Minsky machine, states = state primes, 1 = no state)')
    rows = {}
    for j, (a, b) in enumerate(prog):
        na, nb = fac(a), fac(b)
        sin = [STA[i - 4] for i in range(4, 10) if nb[i]]
        sout = [STA[i - 4] for i in range(4, 10) if na[i]]
        s_in = sin[0] if sin else 1
        s_out = sout[0] if sout else 1
        dec = [f'r{PR[i]}' for i in range(4) if nb[i]]
        inc = [f'r{PR[i]}' for i in range(4) if na[i]]
        rows.setdefault(s_in, []).append((j + 1, dec, inc, s_out))
    for s in [1] + STA:
        print(f'  state {s}:')
        for (j, dec, inc, so) in rows.get(s, []):
            test = (' and '.join(f'{d}>0' for d in dec)) if dec else 'always'
            mv = ', '.join([f'{d}-=1' for d in dec] + [f'{u}+=1' for u in inc]) or 'no register change'
            print(f'      fraction #{j:<2} if {test:22} then {mv:28} -> state {so}')
    # state transition digraph
    arcs = sorted({(s, so) for s in rows for (_, _, _, so) in rows[s]})
    print(f'  state digraph arcs: {arcs}')
    return arcs


def run(name, prog, nprimes):
    section(f'(4) {name}: the run from 2 (exponent vectors), emissions 2^p')
    v = fac(2)
    steps = 0
    found = []
    trace = [(0, None, v)]
    while len(found) < nprimes:
        j, v = step(prog, v)
        steps += 1
        if len(found) < 2:
            trace.append((steps, j, v))
        if v[0] > 0 and all(x == 0 for x in v[1:]):
            found.append((v[0], steps))
    print('  step  fraction  configuration (prime factorization)          exponents (2,3,5,7 | 11,13,17,19,23,29)')
    for (t, j, w) in trace:
        fr = '-' if j is None else f'{prog[j][0]}/{prog[j][1]}'
        print(f'  {t:>4}  {fr:>7}   {cfg(w):42} {w[:4]} | {w[4:]}')
    print(f'  ... (trace shown through the emission of 2^3)')
    print(f'  emissions (p, step): {found}')
    return found


def rcwa(name, prog):
    section(f'(5) {name}: residue-class form; exact selection densities on residues mod 29#')
    # n mod 29# determines which fraction applies; divisibility by distinct primes is independent on residues.
    dens = [Fraction(0)] * len(prog)
    for pattern in itertools.product((0, 1), repeat=10):  # pattern[i] = 1 iff p_i | n
        w = Fraction(1)
        for i, p in enumerate(PR):
            w *= Fraction(1, p) if pattern[i] else Fraction(p - 1, p)
        for j, (a, b) in enumerate(prog):
            nb = fac(b)
            if all(pattern[i] for i in range(10) if nb[i]):
                dens[j] += w
                break
    for j, (a, b) in enumerate(prog):
        print(f'  #{j + 1:<2} {a:>3}/{b:<4} selected on a fraction {str(dens[j]):>28} = {float(dens[j]):.6f} of residues')
    print(f'  total = {sum(dens)} (the last fraction has denominator 1, so some fraction always applies)')
    return dens


def main():
    triv = {}
    for name, prog in (('CONWAY list (Conway 1987, Guy 1983)', CONWAY), ('OWNER list (= Wikipedia form)', OWNER)):
        triv[name] = factor_table(name, prog)
        assert state_split(name, prog)
        machine_table(name, prog)
    fc = run('CONWAY list', CONWAY, 10)
    fo = run('OWNER list', OWNER, 10)
    primes10 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    print()
    print(f'  CONWAY emissions = first ten primes: {[p for p, _ in fc] == primes10}; steps = {[s for _, s in fc]}')
    print(f'  OEIS A007546 (fast version, first ten terms):      {A007546_HEAD}  match: '
          f'{[s for _, s in fc] == A007546_HEAD}')
    print(f'  OWNER emissions = first ten primes: {[p for p, _ in fo] == primes10}; steps = {[s for _, s in fo]}')
    print(f'  OEIS A007547 (slow version, first ten terms):      {A007547_HEAD}  match: '
          f'{[s for _, s in fo] == A007547_HEAD}')
    assert [s for _, s in fc] == A007546_HEAD and [s for _, s in fo] == A007547_HEAD
    assert [p for p, _ in fc] == primes10 and [p for p, _ in fo] == primes10
    rcwa('CONWAY list', CONWAY)
    rcwa('OWNER list', OWNER)
    section("(6) the owner's {1/7, 1/17, 55/1}: list dependence and label dependence")
    for k, v in triv.items():
        print(f'  {k}: trivial-sided fractions {v}')
    print('  The two lists emit the same primes (above; all n in opus-2026-10-31-S70 HYP-9448), yet the set of')
    print('  trivial-sided fractions has 2 elements in Conway\'s list and 3 in the owner\'s: it is not an invariant of')
    print('  the computation. Renaming the six state primes (any bijection of {11,...,29}) or the four registers gives')
    print('  an equivalent program; e.g. the transposition 13 <-> 17 turns 1/17 into 1/13. The numbers 7, 17, 55 are')
    print('  labels, not invariants.')
    print()
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
