#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: two more halts-iff-counterexample FRACTRAN programs on PRIMEGAME's divisor loop.

  FERMATGAME   halts iff some Fermat number F_m = 2^(2^m) + 1 with m >= 5 is prime (a sixth Fermat prime; THM-4847).
  LEGENDREGAME halts iff some n >= 1 has no prime p with n^2 < p < (n+1)^2 (Legendre's conjecture fails).

Both reuse the eleven PRIMEGAME fractions of GOLDBACHGAME's kernel (states 11..29, registers 2, 3, 5, 7; exit
1/17 -> X/17; contract: from state 11 with (r2, r3, r5, r7) = (0, M-1, M, 0), M >= 2, reach X with
(M, 0, 0, D(M) - 1)). Loops use a helper state s' with the single fraction s/s' (FRACTRAN has no self-loops).

FERMATGAME (registers K = 2^m, CK, K2, E, T): from CPK with E = 1: copy K into CK and 2K into K2; double E CK times
(DA/DB); move E into r5, r3 with r5 += 1 (M = 2^K + 1); kernel; at X: r7 > 0 (composite) -> clear r7, clear r2,
K := K2, E := 1, next m; r7 = 0 (prime) -> H (halt). Start for m = s: CPK * K^(2^s) * E.
LEGENDREGAME (registers C1 = C - 1, W = candidates left after C, L = 2n, T): test C; prime -> next window
(C1 += W + 2, L += 2, W := L - 1); composite -> next candidate, or H if the window is exhausted. Start (n = 1):
TEST * C1 * W * L^2  (window {2, 3}).

PRE-REGISTERED (before the first run of this file):
  F1  FERMATGAME has 36 fractions over 32 primes (9 registers, 23 states).
  F2  Started at m = s for s = 0, 1, 2, 3, 4 it halts after testing exactly F_s (prime).
  F3  The inverted control (halt at the first composite) started at m = 0 passes F_0 .. F_4 and reaches the start of
      the test of F_5 = 4294967297 without halting.
  L1  LEGENDREGAME has 33 fractions over 29 primes (8 registers, 21 states).
  L2  It does not halt for any window n = 1 .. 120, and the prime it accepts in each window is the least prime
      greater than n^2 (direct check).
  L3  Negative control (only the first three candidates of each window are tried): halts while processing n = 5.
  K1  For both programs, literal FRACTRAN execution and the state-dispatch simulator agree configuration by
      configuration on an initial segment (F: from m = 2 until the halt; L: windows n = 1 .. 3).

usage: python3 primekernel_programs.py
"""
import sys
from fractions import Fraction

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


KERNEL = [('S13', {'r7': 1}, {}, 'S17'), ('S17', {'r5': 1}, {'r2': 1, 'r3': 1}, 'S13'), ('S17', {'r3': 1}, {}, 'S19'),
          ('S19', {'r2': 1}, {}, 'S23'), ('S11', {'r3': 1}, {}, 'S29'), ('S29', {}, {'r7': 1}, 'S11'),
          ('S23', {}, {'r5': 1}, 'S19'), ('S19', {}, {'r7': 1}, 'S11'), ('S17', {}, {}, 'X'),
          ('S13', {}, {}, 'S11'), ('S11', {}, {}, 'S13')]
KPRIMES = {'r2': 2, 'r3': 3, 'r5': 5, 'r7': 7, 'S11': 11, 'S13': 13, 'S17': 17, 'S19': 19, 'S23': 23, 'S29': 29}
POOL = [31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139]


class Prog:
    def __init__(self, regs, states, ctrl):
        self.regs = ['r2', 'r3', 'r5', 'r7'] + regs
        self.prime = dict(KPRIMES)
        pool = iter(POOL)
        for name in regs + states:
            self.prime[name] = next(pool)
        self.I = KERNEL + ctrl
        for s, g, inc, nxt in self.I:
            assert s != nxt
        self.by = {}
        for s, g, inc, nxt in self.I:
            self.by.setdefault(s, []).append((g, inc, nxt))

    def mono(self, d):
        v = 1
        for k, e in d.items():
            v *= self.prime[k] ** e
        return v

    def fractions(self):
        out = []
        for s, g, inc, nxt in self.I:
            f = Fraction(self.prime[nxt] * self.mono(inc), self.prime[s] * self.mono(g))
            assert f.denominator % self.prime[s] == 0
            out.append(f)
        return out

    def nprimes(self):
        used = set()
        for f in self.fractions():
            for x in (f.numerator, f.denominator):
                for p in self.prime.values():
                    if x % p == 0:
                        used.add(p)
        return len(used)

    def pick(self, s, r):
        for idx, (g, inc, nxt) in enumerate(self.by.get(s, [])):
            if all(r[k] >= e for k, e in g.items()):
                return idx, g, inc, nxt
        return None

    def run(self, s, r, hook, maxsteps=10 ** 15, trace_limit=0):
        r = {k: r.get(k, 0) for k in self.regs}
        steps, tr = 0, []
        while steps < maxsteps:
            got = self.pick(s, r)
            if got is None:
                return ('HALT', s, steps, r), tr
            idx, g, inc, nxt = got
            if not trace_limit and nxt in self.by:
                mid = {k: inc.get(k, 0) - g.get(k, 0) for k in self.regs}
                r1 = {k: r[k] + mid[k] for k in self.regs}
                got2 = self.pick(nxt, r1)
                if got2 is not None and got2[3] == s:
                    idx2, g2, inc2, _ = got2
                    d = {k: mid[k] + inc2.get(k, 0) - g2.get(k, 0) for k in self.regs}
                    ok = all(all(d[k] <= 0 for k in self.by[s][j][0]) for j in range(idx)) and \
                        all(all(d[k] <= 0 for k in self.by[nxt][j][0]) for j in range(idx2))
                    dec = [k for k in self.regs if d[k] < 0 and (k in g or k in g2)]
                    if ok and dec:
                        kmax = min((r[k] - max(g.get(k, 0), g2.get(k, 0) - mid[k])) // (-d[k]) + 1 for k in dec)
                        if kmax > 1:
                            for k in self.regs:
                                r[k] += kmax * d[k]
                            steps += 2 * kmax
                            continue
            for k, e in g.items():
                r[k] -= e
            for k, e in inc.items():
                r[k] += e
            prev, s = s, nxt
            steps += 1
            if trace_limit and steps <= trace_limit:
                v = {k: x for k, x in r.items() if x}
                v[s] = 1
                tr.append(v)
            stop = hook(prev, s, r, steps)
            if stop:
                return ('STOP', s, steps, r, stop), tr
        return ('BUDGET', s, steps, r), tr

    def literal(self, start_state, regs, nsteps):
        n = self.prime[start_state] * self.mono(regs)
        fr = self.fractions()
        out = []
        for _ in range(nsteps):
            for f in fr:
                if n % f.denominator == 0:
                    n = n // f.denominator * f.numerator
                    break
            else:
                break
            out.append(self.decode(n))
        return out

    def decode(self, n):
        v = {}
        for name, p in self.prime.items():
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            if e:
                v[name] = e
        assert n == 1
        return v


def fermat(inverted=False):
    regs = ['K', 'CK', 'K2', 'E', 'T']
    states = ['CPK', 'CPKh', 'OUT', 'DA', 'DAh', 'DB', 'DBh', 'TM', 'TMh', 'X', 'CL', 'CLh', 'C2', 'C2h', 'MK',
              'MKh', 'H']
    c = [('CPK', {'K': 1}, {'CK': 1, 'K2': 2}, 'CPKh'), ('CPK', {}, {}, 'OUT'), ('CPKh', {}, {}, 'CPK'),
         ('OUT', {'CK': 1}, {}, 'DA'), ('OUT', {}, {}, 'TM'),
         ('DA', {'E': 1}, {'T': 2}, 'DAh'), ('DA', {}, {}, 'DB'), ('DAh', {}, {}, 'DA'),
         ('DB', {'T': 1}, {'E': 1}, 'DBh'), ('DB', {}, {}, 'OUT'), ('DBh', {}, {}, 'DB'),
         ('TM', {'E': 1}, {'r5': 1, 'r3': 1}, 'TMh'), ('TM', {}, {'r5': 1}, 'S11'), ('TMh', {}, {}, 'TM')]
    if not inverted:
        c += [('X', {'r7': 1}, {}, 'CL'), ('X', {}, {}, 'H')]
    else:
        c += [('X', {'r7': 1}, {}, 'H'), ('X', {}, {}, 'C2')]
    c += [('CL', {'r7': 1}, {}, 'CLh'), ('CL', {}, {}, 'C2'), ('CLh', {}, {}, 'CL'),
          ('C2', {'r2': 1}, {}, 'C2h'), ('C2', {}, {}, 'MK'), ('C2h', {}, {}, 'C2'),
          ('MK', {'K2': 1}, {'K': 1}, 'MKh'), ('MK', {}, {'E': 1}, 'CPK'), ('MKh', {}, {}, 'MK')]
    return Prog(regs, states, c)


def legendre(control=False):
    regs = ['C1', 'W', 'L', 'T']
    states = ['TEST', 'TESTh', 'X', 'Xh', 'X3', 'CL', 'CLh', 'NEXTC', 'NW', 'NWh', 'CP', 'CPh', 'RT', 'RTh', 'H']
    c = [('TEST', {'C1': 1}, {'r5': 1, 'r3': 1}, 'TESTh'), ('TEST', {}, {'r5': 1}, 'S11'), ('TESTh', {}, {}, 'TEST'),
         ('X', {'r2': 1}, {'C1': 1}, 'Xh'), ('X', {'C1': 1}, {}, 'X3'), ('Xh', {}, {}, 'X'),
         ('X3', {'r7': 1}, {}, 'CL'), ('X3', {}, {}, 'NW'),
         ('CL', {'r7': 1}, {}, 'CLh'), ('CL', {}, {}, 'NEXTC'), ('CLh', {}, {}, 'CL'),
         ('NEXTC', {'W': 1}, {'C1': 1}, 'TEST'), ('NEXTC', {}, {}, 'H'),
         ('NW', {'W': 1}, {'C1': 1}, 'NWh')]
    if not control:
        c += [('NW', {}, {'C1': 2, 'L': 2}, 'CP')]
    else:
        c += [('NW', {}, {'C1': 2, 'L': 2, 'W': 2}, 'TEST')]
    c += [('NWh', {}, {}, 'NW'),
          ('CP', {'L': 1}, {'W': 1, 'T': 1}, 'CPh'), ('CP', {}, {}, 'RT'), ('CPh', {}, {}, 'CP'),
          ('RT', {'T': 1}, {'L': 1}, 'RTh'), ('RT', {'W': 1}, {}, 'TEST'), ('RTh', {}, {}, 'RT')]
    return Prog(regs, states, c)


def is_prime(m):
    return m > 1 and all(m % p for p in range(2, int(m ** 0.5) + 1))


def main():
    # ---------------- FERMATGAME
    P = fermat()
    fr = P.fractions()
    print('FERMATGAME fractions:', ', '.join(str(f) for f in fr))
    check('F1: 36 fractions over 32 primes (9 registers, 23 states)', len(fr) == 36 and P.nprimes() == 32,
          f'{len(fr)} fractions, {P.nprimes()} primes')
    tested = []

    def hook_f(prev, s, r, steps):
        if s == 'S11' and prev == 'TM':
            tested.append(r['r5'])
        return None
    ok = True
    for sidx in range(5):
        tested.clear()
        st, _ = P.run('CPK', {'K': 2 ** sidx, 'E': 1}, hook_f)
        Fm = 2 ** (2 ** sidx) + 1
        good = st[0] == 'HALT' and st[1] == 'H' and tested == [Fm]
        ok &= good
        print(f'    start m={sidx}: status {st[0]} in state {st[1]} after {st[2]} steps; tested {tested}')
    check('F2: started at m = 0..4, halts after testing exactly F_m', ok)
    Pi = fermat(inverted=True)
    tested.clear()

    def hook_fi(prev, s, r, steps):
        if s == 'S11' and prev == 'TM':
            tested.append(r['r5'])
            if r['r5'] == 2 ** 32 + 1:
                return 'reached F5 test'
        return None
    st, _ = Pi.run('CPK', {'K': 1, 'E': 1}, hook_fi)
    check('F3: inverted control from m = 0 passes F_0..F_4 and reaches the test of F_5',
          st[0] == 'STOP' and tested == [3, 5, 17, 257, 65537, 4294967297], f'tested {tested}; status {st[0]}')
    # K1 for FERMATGAME: literal vs dispatch from m = 2 until halt
    st, tr = P.run('CPK', {'K': 4, 'E': 1}, lambda *a: None, trace_limit=10 ** 7)
    lit = P.literal('CPK', {'K': 4, 'E': 1}, len(tr) + 5)
    check('K1 (F): literal FRACTRAN = dispatch, every configuration, m = 2 until the halt',
          lit == tr and st[0] == 'HALT', f'{len(tr)} configurations')
    # ---------------- LEGENDREGAME
    L = legendre()
    fr = L.fractions()
    print('LEGENDREGAME fractions:', ', '.join(str(f) for f in fr))
    check('L1: 33 fractions over 29 primes (8 registers, 21 states)', len(fr) == 33 and L.nprimes() == 29,
          f'{len(fr)} fractions, {L.nprimes()} primes')
    accepted = []

    def hook_l(prev, s, r, steps, nmax=120):
        if prev == 'X3' and s == 'NW':
            accepted.append((r['L'] // 2, r['C1'] + 1))
            if len(accepted) >= nmax:
                return 'enough windows'
        return None
    st, _ = L.run('TEST', {'C1': 1, 'W': 1, 'L': 2}, hook_l)
    ok_least = all(p == min(x for x in range(n * n + 1, (n + 1) ** 2) if is_prime(x)) for n, p in accepted)
    check('L2: no halt for windows n = 1..120; accepted prime = least prime > n^2',
          st[0] == 'STOP' and [n for n, _ in accepted] == list(range(1, 121)) and ok_least,
          f'first accepted {accepted[:6]}; status {st[0]}, steps {st[2]}')
    Lc = legendre(control=True)
    accepted.clear()
    st, _ = Lc.run('TEST', {'C1': 1, 'W': 1, 'L': 2}, hook_l)
    halt_n = st[3]['L'] // 2 if st[0] == 'HALT' else None
    check('L3: negative control (first three candidates per window) halts while processing n = 5',
          st[0] == 'HALT' and st[1] == 'H' and [n for n, _ in accepted] == [1, 2, 3, 4],
          f'accepted windows {accepted}; halted in {st[1]} with L = {st[3]["L"]} (window n = {halt_n})')
    accepted.clear()
    st, tr = L.run('TEST', {'C1': 1, 'W': 1, 'L': 2}, lambda prev, s, r, steps: 'stop' if (prev == 'X3' and s == 'NW'
                   and r['L'] == 6) else None, trace_limit=10 ** 7)
    lit = L.literal('TEST', {'C1': 1, 'W': 1, 'L': 2}, len(tr))
    check('K1 (L): literal FRACTRAN = dispatch, every configuration, windows n = 1..3', lit == tr,
          f'{len(tr)} configurations')
    print()
    if FAIL:
        print('FAILED:', FAIL)
        return 1
    print('ALL PRE-REGISTERED CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
