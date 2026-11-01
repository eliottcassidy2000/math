#!/usr/bin/env python3
"""mac-mini-2026-11-01-lenses: exact translations between a Turing machine, a cellular automaton and a FRACTRAN
program, made executable and checked in lockstep on Leng's 25-state Goldbach machine (and on BB(2) as a halting
control).

(TM -> CA)  Standard radius-1 encoding (folklore; e.g. Smith 1971): cell alphabet = tape symbols {0,1}, head cells
            (q, g) for every state q and symbol g, and a halt marker HALT. Rule F(a,b,c):
              b = (q,g): delta(q,g) undefined -> HALT; else write w (the head leaves the cell);
              a = (q,g) with move R -> (q', b);   c = (q,g) with move L -> (q', b);   HALT persists; else b.
            One CA step is exactly one TM step.
(TM -> FRACTRAN)  Two-stack (Minsky) encoding with one-hot control: R = sum_{i>=0} tape[h+i] 2^i (head = bit 0),
            L = sum_{i>=1} tape[h-i] 2^(i-1); registers 2 (L), 3 (R), 5 (A), 7 (C); per TM state four fractions
            extract the head bit (A := R div 2), per instruction nine (move R: R := A, L := 2L + w) or ten (move L:
            L := L div 2 with its low bit pushed onto R, R := 4A + 2w + lowbit) fractions apply it. Loops use
            helper states (FRACTRAN has no self-loops). An undefined slot is a control state with no fractions.

PRE-REGISTERED (before the first run of this file):
  T1  The CA of Leng's machine has 53 symbols (2 + 25*2 + HALT) and agrees with the TM at every step for 200,000
      steps from the blank tape (full-window updates for the first 20,000 steps, then head-local updates, which
      equal full updates because F is the identity on neighbourhoods without a head).
  T2  The compiled FRACTRAN program LENGGAME has 25*4 + 22*9 + 27*10 = 568 fractions (Leng's table has 22 right
      moves and 27 left moves among its 49 instructions), and its register contents at the k-th entry into a
      TM-state prime equal the TM's two stacks after k steps, for k <= 3,000 (accelerated dispatch), and literal
      integer FRACTRAN agrees with the dispatch simulator configuration by configuration for the first 300 TM
      steps.
  T3  Control: BB(2) (A0:1RB A1:1LB B0:1LA B1:1RH) enters its halt state after 6 steps with four 1s; its CA writes
      HALT at step 7 (the attempted seventh step, where the TM finds no transition) and its compiled FRACTRAN halts
      with the same tape (four 1s).
usage: python3 machine_translations.py
"""
import collections
import os
import re
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
LEAN = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'session65_cryptid_tournament_interfaces_20261026',
                    'session65-twentyfive', 'inputs', 'goldbach_tm', 'GoldbachTm', 'Tm25', 'TuringMachine25.lean')
FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def leng_table():
    src = open(LEAN).read()
    pat = re.compile(r"\|\s*⟨\s*(\d+),\s*_⟩,\s*Γ\.(zero|one)\s*=>\s*(none|some\s*⟨⟨\s*(\d+),\s*by omega⟩,\s*⟨Turing\.Dir\.(left|right),\s*Γ\.(zero|one)⟩⟩)")
    tab = {}
    for m in pat.finditer(src):
        q, g = int(m.group(1)), 0 if m.group(2) == 'zero' else 1
        tab[(q, g)] = None if m.group(3) == 'none' else (int(m.group(4)), 'R' if m.group(5) == 'right' else 'L',
                                                          0 if m.group(6) == 'zero' else 1)
    assert len(tab) == 50
    return tab, list(range(25)), 0


def bb2_table():
    # states A=0, B=1, H=2 (no transitions); entries (next, move, write)
    tab = {(0, 0): (1, 'R', 1), (0, 1): (1, 'L', 1), (1, 0): (0, 'L', 1), (1, 1): (2, 'R', 1),
           (2, 0): None, (2, 1): None}
    return tab, [0, 1, 2], 0


# ---------------------------------------------------------------- TM
class TM:
    def __init__(self, tab, q0):
        self.tab, self.q, self.h, self.tape, self.halted = tab, q0, 0, collections.defaultdict(int), False

    def step(self):
        a = self.tab[(self.q, self.tape[self.h])]
        if a is None:
            self.halted = True
            return False
        nq, mv, w = a
        self.tape[self.h] = w
        self.h += 1 if mv == 'R' else -1
        self.q = nq
        return True

    def stacks(self):
        R = L = 0
        pos = [p for p, v in self.tape.items() if v]
        if pos:
            for p in pos:
                if p >= self.h:
                    R += 1 << (p - self.h)
                else:
                    L += 1 << (self.h - p - 1)
        return L, R


# ---------------------------------------------------------------- CA
def ca_rule(tab):
    def F(a, b, c):
        if b == 'HALT':
            return 'HALT'
        if isinstance(b, tuple):
            act = tab[b]
            return 'HALT' if act is None else act[2]
        if isinstance(a, tuple):
            act = tab[a]
            if act is not None and act[1] == 'R':
                return (act[0], b)
        if isinstance(c, tuple):
            act = tab[c]
            if act is not None and act[1] == 'L':
                return (act[0], b)
        return b
    return F


def ca_alphabet(tab, states):
    return [0, 1] + [(q, g) for q in states for g in (0, 1)] + ['HALT']


def ca_lockstep(tab, states, q0, nsteps, full_until, width):
    F = ca_rule(tab)
    cells = [0] * (2 * width + 1)       # index i <-> position i - width
    cells[width] = (q0, 0)
    tm = TM(tab, q0)
    for t in range(1, nsteps + 1):
        hpos = next(i for i, v in enumerate(cells) if isinstance(v, tuple) or v == 'HALT') if t <= full_until else hidx
        if t <= full_until:
            new = [cells[0]] + [F(cells[i - 1], cells[i], cells[i + 1]) for i in range(1, len(cells) - 1)] + [cells[-1]]
        else:
            new = list(cells) if t == full_until + 1 else cells
            lo, hi = hpos - 1, hpos + 1
            vals = [F(cells[i - 1], cells[i], cells[i + 1]) for i in range(lo, hi + 1)]
            for k, i in enumerate(range(lo, hi + 1)):
                new[i] = vals[k]
        cells = new
        ok_tm = tm.step()
        # locate head / halt
        hidx = next((i for i in range(max(0, (hpos if t > 1 else width) - 2), min(len(cells), (hpos if t > 1 else width) + 3))
                     if isinstance(cells[i], tuple) or cells[i] == 'HALT'), None)
        if hidx is None:
            return False, t, 'head lost'
        if not ok_tm:
            return cells[hidx] == 'HALT', t, 'halted'
        if hidx - width != tm.h or cells[hidx] != (tm.q, tm.tape[tm.h]):
            return False, t, 'head mismatch'
        if t <= full_until or t % 5000 == 0:
            for i, v in enumerate(cells):
                sym = v[1] if isinstance(v, tuple) else v
                if sym != tm.tape[i - width]:
                    return False, t, f'tape mismatch at {i - width}'
        if hidx < 3 or hidx > len(cells) - 4:
            return False, t, 'window too small'
    return True, nsteps, 'ran'


# ---------------------------------------------------------------- TM -> FRACTRAN
class Compiled:
    def __init__(self, tab, states):
        self.regs = {'L': 2, 'R': 3, 'A': 5, 'C': 7}
        primes = self._primes(5000)
        it = iter(p for p in primes if p > 7)
        self.sp = {}

        def S(name):
            if name not in self.sp:
                self.sp[name] = next(it)
            return name
        self.I = []

        def ins(s, guard, inc, nxt):
            assert s != nxt
            self.I.append((S(s), guard, inc, S(nxt)))
        self.entry = {}
        for q in states:
            P, Ph, B0, B1 = f'P{q}', f'Ph{q}', f'B0_{q}', f'B1_{q}'
            S(P); S(Ph); S(B0); S(B1)
            self.entry[q] = P
            ins(P, {'R': 2}, {'A': 1}, Ph)
            ins(P, {'R': 1}, {}, B1)
            ins(P, {}, {}, B0)
            ins(Ph, {}, {}, P)
        for (q, g), act in tab.items():
            B = f'B{g}_{q}'
            if act is None:
                continue
            nq, mv, w = act
            tgt = f'P{nq}'
            if mv == 'R':
                M1h, M2, M2h, M3, M3h = (f'{x}_{q}_{g}' for x in ('M1h', 'M2', 'M2h', 'M3', 'M3h'))
                ins(B, {'A': 1}, {'R': 1}, M1h); ins(B, {}, {}, M2); ins(M1h, {}, {}, B)
                ins(M2, {'L': 1}, {'A': 2}, M2h); ins(M2, {}, {}, M3); ins(M2h, {}, {}, M2)
                ins(M3, {'A': 1}, {'L': 1}, M3h); ins(M3, {}, {'L': w} if w else {}, tgt); ins(M3h, {}, {}, M3)
            else:
                N1h, N2, N2h, N3, N3h = (f'{x}_{q}_{g}' for x in ('N1h', 'N2', 'N2h', 'N3', 'N3h'))
                ins(B, {'L': 2}, {'C': 1}, N1h); ins(B, {'L': 1}, {'R': 1}, N2); ins(B, {}, {}, N2)
                ins(N1h, {}, {}, B)
                ins(N2, {'C': 1}, {'L': 1}, N2h); ins(N2, {}, {}, N3); ins(N2h, {}, {}, N2)
                ins(N3, {'A': 1}, {'R': 4}, N3h); ins(N3, {}, {'R': 2 * w} if w else {}, tgt); ins(N3h, {}, {}, N3)
        self.by = {}
        for s, g, inc, nxt in self.I:
            self.by.setdefault(s, []).append((g, inc, nxt))

    @staticmethod
    def _primes(n):
        return [p for p in range(2, n) if all(p % d for d in range(2, int(p ** 0.5) + 1))]

    def fractions(self):
        out = []
        for s, g, inc, nxt in self.I:
            num = self.sp[nxt]
            den = self.sp[s]
            for k, e in inc.items():
                num *= self.regs[k] ** e
            for k, e in g.items():
                den *= self.regs[k] ** e
            f = Fraction(num, den)
            assert f.denominator % self.sp[s] == 0 and f.numerator % self.sp[nxt] == 0
            out.append(f)
        return out

    def pick(self, s, r):
        for idx, (g, inc, nxt) in enumerate(self.by.get(s, [])):
            if all(r[k] >= e for k, e in g.items()):
                return idx, g, inc, nxt
        return None

    def run(self, start, nentries, accel=True, trace=False):
        """returns the list of (state, L, R) at each entry into a P-state, the number of FRACTRAN steps, halt flag"""
        r = {'L': 0, 'R': 0, 'A': 0, 'C': 0}
        s = start
        steps = 0
        entries = [(s, r['L'], r['R'], 0)]
        tr = []
        Pstates = set(self.entry.values())
        while len(entries) <= nentries:
            got = self.pick(s, r)
            if got is None:
                return entries, steps, True, tr
            idx, g, inc, nxt = got
            if accel and nxt in self.by:
                mid = {k: inc.get(k, 0) - g.get(k, 0) for k in r}
                r1 = {k: r[k] + mid[k] for k in r}
                got2 = self.pick(nxt, r1)
                if got2 is not None and got2[3] == s:
                    idx2, g2, inc2, _ = got2
                    d = {k: mid[k] + inc2.get(k, 0) - g2.get(k, 0) for k in r}
                    ok = all(all(d[k] <= 0 for k in self.by[s][j][0]) for j in range(idx)) and \
                        all(all(d[k] <= 0 for k in self.by[nxt][j][0]) for j in range(idx2))
                    dec = [k for k in r if d[k] < 0 and (k in g or k in g2)]
                    if ok and dec:
                        kmax = min((r[k] - max(g.get(k, 0), g2.get(k, 0) - mid[k])) // (-d[k]) + 1 for k in dec)
                        if kmax > 1:
                            for k in r:
                                r[k] += kmax * d[k]
                            steps += 2 * kmax
                            continue
            for k, e in g.items():
                r[k] -= e
            for k, e in inc.items():
                r[k] += e
            s = nxt
            steps += 1
            if trace:
                tr.append((s, dict(r)))
            if s in Pstates:
                entries.append((s, r['L'], r['R'], steps))
        return entries, steps, False, tr

    def literal(self, start, nsteps):
        n = self.sp[start]
        fr = self.fractions()
        inv = {p: k for k, p in self.regs.items()}
        sinv = {p: k for k, p in self.sp.items()}
        out = []
        for _ in range(nsteps):
            for f in fr:
                if n % f.denominator == 0:
                    n = n // f.denominator * f.numerator
                    break
            else:
                break
            m, regs, st = n, {'L': 0, 'R': 0, 'A': 0, 'C': 0}, None
            for p, k in inv.items():
                while m % p == 0:
                    m //= p
                    regs[k] += 1
            st = sinv[m]
            out.append((st, regs))
        return out


def main():
    tab, states, q0 = leng_table()
    print('(T1) Leng TM -> radius-1 CA')
    alpha = ca_alphabet(tab, states)
    check('T1a: CA alphabet has 53 symbols', len(alpha) == 53, f'{len(alpha)}')
    ok, t, why = ca_lockstep(tab, states, q0, 200_000, 20_000, width=320)
    check('T1b: CA = TM at every step for 200,000 steps', ok and why == 'ran', f'{why} at t={t}')
    print('(T2) Leng TM -> FRACTRAN (LENGGAME)')
    C = Compiled(tab, states)
    fr = C.fractions()
    nR = sum(1 for a in tab.values() if a and a[1] == 'R')
    nL = sum(1 for a in tab.values() if a and a[1] == 'L')
    check('T2a: 568 fractions (25*4 + 22*9 + 27*10)', len(fr) == 568 and (nR, nL) == (22, 27),
          f'{len(fr)} fractions, moves R={nR} L={nL}, {len(C.sp) + 4} primes, largest {max(C.sp.values())}')
    entries, steps, halted, _ = C.run(C.entry[q0], 3000)
    tm = TM(tab, q0)
    agree = True
    for k, (s, L, R, st) in enumerate(entries):
        if k > 0:
            tm.step()
        if s != C.entry[tm.q] or (L, R) != tm.stacks():
            agree = False
            print(f'    mismatch at TM step {k}: fractran {s, L, R} vs tm {tm.q, tm.stacks()}')
            break
    check('T2b: register stacks = TM stacks at each of the first 3,000 TM steps', agree and len(entries) == 3001,
          f'FRACTRAN steps for 3,000 TM steps: {steps}')
    # literal vs dispatch (no acceleration), first 300 TM steps
    entries_na, steps_na, _, trace = C.run(C.entry[q0], 300, accel=False, trace=True)
    lit = C.literal(C.entry[q0], len(trace))
    same = len(lit) == len(trace) and all(lit[i][0] == trace[i][0] and lit[i][1] == trace[i][1]
                                          for i in range(len(trace)))
    check('T2c: literal integer FRACTRAN = dispatch, every configuration, first 300 TM steps', same,
          f'{len(trace)} configurations')
    print('    first fractions:', ', '.join(str(f) for f in fr[:12]), '...')
    # step overhead growth
    sizes = [(k, entries[k][3]) for k in (100, 1000, 2000, 3000)]
    print(f'    cumulative FRACTRAN steps at TM steps 100/1000/2000/3000: {sizes}')
    print('(T3) control: BB(2)')
    btab, bstates, bq0 = bb2_table()
    tm = TM(btab, bq0)
    nsteps = 0
    while tm.step():
        nsteps += 1
    ones = sum(tm.tape.values())
    okca, tca, whyca = ca_lockstep(btab, bstates, bq0, 50, 50, width=20)
    Cb = Compiled(btab, bstates)
    eb, sb, hb, _ = Cb.run(Cb.entry[bq0], 100)
    lastL, lastR = eb[-1][1], eb[-1][2]
    check('T3: BB(2) halts after 6 steps with four 1s; CA gives HALT at step 6; FRACTRAN halts with four 1s',
          nsteps == 6 and ones == 4 and okca and whyca == 'halted' and tca == 7 and hb
          and bin(lastL).count('1') + bin(lastR).count('1') == 4,
          f'TM steps {nsteps}, ones {ones}; CA {whyca} at t={tca}; FRACTRAN halted={hb} after {sb} steps, '
          f'final stacks L={lastL} R={lastR}')
    print()
    if FAIL:
        print('FAILED:', FAIL)
        return 1
    print('ALL PRE-REGISTERED CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
