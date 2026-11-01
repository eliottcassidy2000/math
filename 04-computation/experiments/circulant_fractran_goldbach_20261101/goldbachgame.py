#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: GOLDBACHGAME - a FRACTRAN program, built around PRIMEGAME's own divisor loop,
that halts iff the strong Goldbach conjecture fails.

Design (a Minsky machine compiled to FRACTRAN; one state prime is present at every step):
  * Primality subroutine = PRIMEGAME's first eleven fractions verbatim (registers 2, 3, 5, 7; states 11..29), with
    Conway's exit 1/17 replaced by X/17. Contract (Conway 1987; THM-4853 proof; session126 section 4): entering
    state 11 with (r2, r3, r5, r7) = (0, M-1, M, 0), M >= 2, it reaches X with (r2, r3, r5, r7) = (M, 0, 0, D(M)-1),
    D = largest proper divisor; so M is prime iff r7 = 0.
  * Controller registers: P' = P - 2 (prime 31), Q' = Q - 2 (prime 37), flag F (prime 41); the current target is
    N = P + Q. Start: the single state prime TP (N = 4, P = Q = 2).
      TP   : move P' into r5, r3; then add r5 += 2, r3 += 1 and enter 11            (test P)
      X    : F > 0 ? (F--, XQ) : XP
      XP   : move r2 back into P'; then P' -= 2; -> XP3
      XP3  : r7 > 0 ? (r7--, CL) : TQ                                                 (P composite / prime)
      TQ   : move Q' into r5, r3; then r5 += 2, r3 += 1, F += 1, enter 11              (test Q)
      XQ   : move r2 back into Q'; then Q' -= 2; -> XQ3
      XQ3  : r7 > 0 ? (r7--, CL) : SUCC                                               (Q composite / prime)
      CL   : clear r7; -> NEXT
      NEXT : Q' > 0 ? (Q'--, P'++, TP) : H                                            (next split, or HALT)
      SUCC : move P' into Q'; then Q' += 2; -> TP                                     (next target N + 2)
      H    : no fractions, so FRACTRAN halts.
    FRACTRAN cannot express a self-loop (s.a/(s.b) reduces to a/b and loses the state guard), so every loop
    s -> s uses a helper state s' with the single fraction s/s'.
  Every fraction carries its state prime in the denominator, so in state H nothing applies.

PRE-REGISTERED (before the first run of this file):
  G1  The compiled program has 37 fractions over 30 distinct primes (7 registers, 23 states).
  G2  Literal FRACTRAN execution (exact integers) and the state-dispatch simulator agree configuration by
      configuration through the first 6 targets (N = 4 .. 14).
  G3  The program never halts while processing every even N <= 1000, the targets are processed in increasing
      order, and for each N the accepted split (P, N - P) is the one with the least prime P (direct check).
  G4  Negative control: adding the single fraction H/(NEXT * 31) first in NEXT's group (only P in {2, 3} may be
      tried) makes the program halt, and it halts exactly while processing N = 12 (12 - 2 and 12 - 3 composite).
  G5  The per-target step count is reported; no prediction.

usage: python3 goldbachgame.py [NMAX]
"""
import sys
from fractions import Fraction

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


# ------------------------------------------------------------------ primes / names
REG = {'r2': 2, 'r3': 3, 'r5': 5, 'r7': 7, 'Pp': 31, 'Qp': 37, 'F': 41}
SUB = {'S11': 11, 'S13': 13, 'S17': 17, 'S19': 19, 'S23': 23, 'S29': 29}
CTRL_NAMES = ['TP', 'TPh', 'X', 'XP', 'XPh', 'XP3', 'TQ', 'TQh', 'XQ', 'XQh', 'XQ3', 'CL', 'CLh', 'NEXT',
              'SUCC', 'SUCCh', 'H']
CTRL_PRIMES = [43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
STATE = dict(SUB)
STATE.update(dict(zip(CTRL_NAMES, CTRL_PRIMES)))


def mono(d):
    v = 1
    for k, e in d.items():
        v *= (REG[k] if k in REG else STATE[k]) ** e
    return v


def program(control_variant=False):
    """list of (state, guard dict, inc dict, next state); FRACTRAN fraction = next*inc / (state*guard)"""
    I = []

    def ins(s, guard, inc, nxt):
        assert s != nxt, 'no self-loops in FRACTRAN'
        I.append((s, guard, inc, nxt))
    # PRIMEGAME's first eleven fractions, written as instructions (exit 1/17 -> X/17)
    ins('S13', {'r7': 1}, {}, 'S17')                     # 17/91
    ins('S17', {'r5': 1}, {'r2': 1, 'r3': 1}, 'S13')     # 78/85
    ins('S17', {'r3': 1}, {}, 'S19')                     # 19/51
    ins('S19', {'r2': 1}, {}, 'S23')                     # 23/38
    ins('S11', {'r3': 1}, {}, 'S29')                     # 29/33
    ins('S29', {}, {'r7': 1}, 'S11')                     # 77/29
    ins('S23', {}, {'r5': 1}, 'S19')                     # 95/23
    ins('S19', {}, {'r7': 1}, 'S11')                     # 77/19
    ins('S17', {}, {}, 'X')                              # X/17   (Conway: 1/17)
    ins('S13', {}, {}, 'S11')                            # 11/13
    ins('S11', {}, {}, 'S13')                            # 13/11
    # controller
    ins('TP', {'Pp': 1}, {'r5': 1, 'r3': 1}, 'TPh')
    ins('TP', {}, {'r5': 2, 'r3': 1}, 'S11')
    ins('TPh', {}, {}, 'TP')
    ins('X', {'F': 1}, {}, 'XQ')
    ins('X', {}, {}, 'XP')
    ins('XP', {'r2': 1}, {'Pp': 1}, 'XPh')
    ins('XP', {'Pp': 2}, {}, 'XP3')
    ins('XPh', {}, {}, 'XP')
    ins('XP3', {'r7': 1}, {}, 'CL')
    ins('XP3', {}, {}, 'TQ')
    ins('TQ', {'Qp': 1}, {'r5': 1, 'r3': 1}, 'TQh')
    ins('TQ', {}, {'r5': 2, 'r3': 1, 'F': 1}, 'S11')
    ins('TQh', {}, {}, 'TQ')
    ins('XQ', {'r2': 1}, {'Qp': 1}, 'XQh')
    ins('XQ', {'Qp': 2}, {}, 'XQ3')
    ins('XQh', {}, {}, 'XQ')
    ins('XQ3', {'r7': 1}, {}, 'CL')
    ins('XQ3', {}, {}, 'SUCC')
    ins('CL', {'r7': 1}, {}, 'CLh')
    ins('CL', {}, {}, 'NEXT')
    ins('CLh', {}, {}, 'CL')
    if control_variant:
        ins('NEXT', {'Pp': 1}, {}, 'H')                   # negative control: only P in {2, 3}
    ins('NEXT', {'Qp': 1}, {'Pp': 1}, 'TP')
    ins('NEXT', {}, {}, 'H')
    ins('SUCC', {'Pp': 1}, {'Qp': 1}, 'SUCCh')
    ins('SUCC', {}, {'Qp': 2}, 'TP')
    ins('SUCCh', {}, {}, 'SUCC')
    return I


def fractions_of(I):
    fr = []
    for s, g, inc, nxt in I:
        num = STATE[nxt] * mono(inc)
        den = STATE[s] * mono(g)
        f = Fraction(num, den)
        assert f.denominator % STATE[s] == 0, 'state guard lost by reduction'
        fr.append(f)
    return fr


# ------------------------------------------------------------------ simulators
def literal(fr, start, maxsteps, watch):
    n = start
    trace = []
    for t in range(maxsteps):
        for f in fr:
            if n % f.denominator == 0:
                n = n // f.denominator * f.numerator
                break
        else:
            return trace, ('HALT', t)
        trace.append(n)
        if watch(n, t):
            break
    return trace, ('STOP', len(trace))


def decode(n):
    v = {}
    for name, p in list(REG.items()) + list(STATE.items()):
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            v[name] = e
    assert n == 1
    return v


class Machine:
    def __init__(self, I):
        self.by = {}
        for s, g, inc, nxt in I:
            self.by.setdefault(s, []).append((g, inc, nxt))
        self.regs = list(REG)

    def pick(self, s, r):
        for idx, (g, inc, nxt) in enumerate(self.by[s]):
            if all(r[k] >= e for k, e in g.items()):
                return idx, g, inc, nxt
        return None

    def run(self, nmax, trace_limit=0):
        """exact step counts; with trace_limit = 0, runs of a 2-cycle s -> t -> s whose two chosen rules stay the
        first applicable ones are executed in one jump (the iteration count is computed exactly)."""
        r = {k: 0 for k in self.regs}
        s = 'TP'
        steps = 0
        events = []          # (N, P, steps) when XQ3 -> SUCC (a split accepted)
        tr = []
        while True:
            if s == 'H':
                return events, ('HALT', steps, r['Pp'] + r['Qp'] + 4), tr
            got = self.pick(s, r)
            if got is None:
                return events, ('STUCK', steps, s), tr
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
            if s == 'XQ3' and nxt == 'SUCC':
                N = r['Pp'] + r['Qp'] + 4
                events.append((N, r['Pp'] + 2, steps + 1))
                if N >= nmax:
                    return events, ('STOP', steps + 1, N), tr
            s = nxt
            steps += 1
            if trace_limit and steps <= trace_limit:
                v = {k: x for k, x in r.items() if x}
                v[s] = 1
                tr.append(v)


def is_prime(m):
    return m > 1 and all(m % p for p in range(2, int(m ** 0.5) + 1))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    I = program()
    fr = fractions_of(I)
    primes_used = set()
    for f in fr:
        for x in (f.numerator, f.denominator):
            for p in list(REG.values()) + list(STATE.values()):
                if x % p == 0:
                    primes_used.add(p)
    print('GOLDBACHGAME fractions (in order):')
    print('  ' + ', '.join(str(f) for f in fr))
    print(f'  start integer: {STATE["TP"]} (state TP, N = 4)')
    check('G1: 37 fractions over 30 primes (7 registers, 23 states)',
          len(fr) == 37 and len(primes_used) == 30, f'{len(fr)} fractions, {len(primes_used)} primes')
    check('first eleven fractions are PRIMEGAME\'s with 1/17 -> X/17',
          [str(f) for f in fr[:11]] == ['17/91', '78/85', '19/51', '23/38', '29/33', '77/29', '95/23', '77/19',
                                       f'{STATE["X"]}/17', '11/13', '13/11'])
    # G2: literal vs dispatch, configuration by configuration
    M = Machine(I)
    ev, st, tr = M.run(14, trace_limit=10 ** 6)
    lit, lst = literal(fr, STATE['TP'], len(tr), lambda n, t: False)
    same = all(decode(lit[i]) == tr[i] for i in range(len(tr)))
    check('G2: literal FRACTRAN = state-dispatch simulator, every configuration through N = 14',
          same and len(lit) == len(tr), f'{len(tr)} configurations, targets {[e[0] for e in ev]}')
    # G3
    ev, st, _ = M.run(nmax)
    Ns = [e[0] for e in ev]
    ok_order = Ns == list(range(4, nmax + 1, 2))
    ok_pairs = all(P == min(p for p in range(2, N - 1) if is_prime(p) and is_prime(N - p)) for N, P, _ in ev)
    check(f'G3: no halt for every even N <= {nmax}; targets in order; accepted split has the least prime P',
          st[0] == 'STOP' and ok_order and ok_pairs, f'status {st}')
    # G4
    Ic = program(control_variant=True)
    Mc = Machine(Ic)
    evc, stc, _ = Mc.run(nmax)
    check('G4: control variant (P in {2,3}) halts while processing N = 12',
          stc[0] == 'HALT' and stc[2] == 12 and [e[0] for e in evc] == [4, 6, 8, 10],
          f'status {stc}, accepted {[e[0] for e in evc]}')
    # G5
    cum = [(N, s) for N, P, s in ev]
    per = [(cum[i][0], cum[i][1] - (cum[i - 1][1] if i else 0)) for i in range(len(cum))]
    print(f'  G5: steps to accept N: {cum[:8]} ... {cum[-1]}')
    print(f'      per-target steps (first): {per[:10]}')
    big = [p for p in per if p[0] in (100, 200, 400, 800, 1000) or p[0] == nmax]
    print(f'      per-target steps at selected N: {big}')
    print()
    if FAIL:
        print('FAILED:', FAIL)
        return 1
    print('ALL PRE-REGISTERED CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
