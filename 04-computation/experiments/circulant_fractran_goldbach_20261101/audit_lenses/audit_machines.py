#!/usr/bin/env python3
"""Independent audit of THM-4858 (Leng's 25-state Goldbach TM -> radius-1 CA, -> FRACTRAN "LENGGAME").

Own code only (stdlib).  Reads Leng's Lean table as data (write-then-move, Stmt = (move, write)).
  M1  table: 50 slots, 49 instructions, 22 right / 27 left moves, (24, 0) the only undefined slot.
  M2  TM -> CA: 53-symbol alphabet, exhaustive check of the local rule's quiescence (F(a,b,c) = b when no head or HALT
      is in the neighbourhood); FULL-window CA updates (every cell, every step) for 200,000 steps in lockstep with an
      own TM simulator; BB(2) control (HALT written at step 7).
  M3  TM -> FRACTRAN following the construction in THM-4858 (b): size 568 = 25*4 + 22*9 + 27*10, primes, the first
      12 fractions printed in machine_translations.out, one-hot property after reduction, unique halting state.
  M4  own exponent-vector FRACTRAN interpreter (literal first-applicable semantics inside the current control
      prime's fractions, justified by M3) with own loop acceleration: stacks at every instruction exit = TM stacks for
      the first 5,000 TM steps; FRACTRAN step count at 3,000 TM steps vs the claimed 24,405,662, and vs a closed-form
      per-step cost computed from the TM alone.
  M5  LITERAL big-integer FRACTRAN (all 568 fractions scanned in order) = vector interpreter at every configuration
      for every FRACTRAN step of the first 300 TM steps; how many TM steps the session's T2c really covered.
  M6  reproduction of the session's failed T2b / T3 numbers from the misaligned "entry" convention.
  M7  BB(2) compiled: halts with four 1s.
usage (from the repo root): python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit_lenses/audit_machines.py
"""
import math
import os
import re
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
LEAN = os.path.join(HERE, '..', '..', '..', '..', '05-knowledge', 'results',
                    'session65_cryptid_tournament_interfaces_20261026', 'session65-twentyfive', 'inputs',
                    'goldbach_tm', 'GoldbachTm', 'Tm25', 'TuringMachine25.lean')
FAILS = []


def report(name, ok, detail=''):
    print(f'  [{"OK" if ok else "MISMATCH"}] {name}' + (f' -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ the table
def parse_lean():
    tab = {}
    for line in open(LEAN, encoding='utf-8'):
        line = line.strip()
        if not line.startswith('| ⟨') or '=>' not in line or '_+25' in line:
            continue
        lhs, rhs = line.split('=>', 1)
        q = int(re.search(r'⟨\s*(\d+)\s*,', lhs).group(1))
        g = 1 if 'Γ.one' in lhs else 0
        if rhs.strip().startswith('none'):
            tab[(q, g)] = None
            continue
        nq = int(re.search(r'⟨⟨\s*(\d+)\s*,', rhs).group(1))
        mv = 'R' if 'Dir.right' in rhs else 'L'
        w = 1 if 'Γ.one' in rhs else 0
        tab[(q, g)] = (nq, mv, w)          # Stmt = (move, write): write w, then move mv, then state nq
    return tab


def bb2():
    # A = 0, B = 1, H = 2 (no transitions)
    return {(0, 0): (1, 'R', 1), (0, 1): (1, 'L', 1), (1, 0): (0, 'L', 1), (1, 1): (2, 'R', 1),
            (2, 0): None, (2, 1): None}


class TM:
    def __init__(self, tab):
        self.tab, self.q, self.h, self.tape = tab, 0, 0, defaultdict(int)

    def step(self):
        act = self.tab[(self.q, self.tape[self.h])]
        if act is None:
            return False
        nq, mv, w = act
        self.tape[self.h] = w
        self.h += 1 if mv == 'R' else -1
        self.q = nq
        return True

    def stacks(self):
        L = R = 0
        for p, v in self.tape.items():
            if v:
                if p >= self.h:
                    R |= 1 << (p - self.h)
                else:
                    L |= 1 << (self.h - 1 - p)
        return L, R


# ------------------------------------------------------------------ TM -> CA
def build_ca(tab, nstates):
    # symbols: 0, 1 tape; 2 + 2q + g head (q, g); HALT last
    HALT = 2 + 2 * nstates
    S = HALT + 1

    def head(s):
        return 2 <= s < HALT

    def qg(s):
        return (s - 2) // 2, (s - 2) % 2

    def F(a, b, c):
        if b == HALT:
            return HALT
        if head(b):
            act = tab[qg(b)]
            return HALT if act is None else act[2]
        if head(a):
            act = tab[qg(a)]
            if act is not None and act[1] == 'R':
                return 2 + 2 * act[0] + b
        if head(c):
            act = tab[qg(c)]
            if act is not None and act[1] == 'L':
                return 2 + 2 * act[0] + b
        return b
    table = [F(a, b, c) for a in range(S) for b in range(S) for c in range(S)]
    return table, S, HALT, head


def ca_run(tab, nstates, nsteps, check_every=100):
    table, S, HALT, head = build_ca(tab, nstates)
    S2 = S * S
    # quiescence lemma (exhaustive): no head/HALT among a, b, c  =>  F = b
    quiet = all(table[a * S2 + b * S + c] == b for a in (0, 1) for b in (0, 1) for c in (0, 1))
    off = 64
    cells = [0] * (2 * off + 1)
    cells[off] = 2 + 0            # head (q0 = 0, g = 0)
    tm = TM(tab)
    for t in range(1, nsteps + 1):
        # grow the window if anything non-blank comes within 2 cells of an end (background 0 is quiescent)
        if any(cells[:3]) or any(cells[-3:]):
            cells = [0] * 64 + cells + [0] * 64
            off += 64
        ext = [0] + cells + [0]
        cells = [table[ext[i] * S2 + ext[i + 1] * S + ext[i + 2]] for i in range(len(cells))]   # FULL update
        alive = tm.step()
        hs = [i for i in range(max(0, tm.h + off - 2), min(len(cells), tm.h + off + 3)) if cells[i] >= 2]
        if not alive:
            return quiet, S, t, 'halt' if any(cells[i] == HALT for i in hs) else 'halt-mismatch'
        if len(hs) != 1 or hs[0] - off != tm.h or cells[hs[0]] != 2 + 2 * tm.q + tm.tape[tm.h]:
            return quiet, S, t, 'head mismatch'
        if t % check_every == 0 or t == nsteps:
            nheads = sum(1 for x in cells if x >= 2)
            if nheads != 1:
                return quiet, S, t, 'extra head'
            for i, x in enumerate(cells):
                sym = x if x < 2 else (x - 2) % 2
                if sym != tm.tape.get(i - off, 0):
                    return quiet, S, t, 'tape mismatch'
    return quiet, S, nsteps, 'ran'


# ------------------------------------------------------------------ TM -> FRACTRAN (THM-4858 (b))
def primes_after7(count):
    out, p = [], 11
    while len(out) < count:
        if all(p % d for d in range(2, int(p ** 0.5) + 1)):
            out.append(p)
        p += 2
    return out


REG = {'L': 2, 'R': 3, 'A': 5, 'C': 7}


def compile_tm(tab, nstates):
    """returns (control prime map, ordered rule list [(state, guard, inc, next)])"""
    names = []

    def st(x):
        if x not in names:
            names.append(x)
        return x
    rules = []

    def rule(s, guard, inc, nxt):
        assert s != nxt
        rules.append((st(s), guard, inc, st(nxt)))
    for q in range(nstates):
        P, Ph, B0, B1 = f'P{q}', f'Ph{q}', f'B0_{q}', f'B1_{q}'
        for x in (P, Ph, B0, B1):
            st(x)
        rule(P, {'R': 2}, {'A': 1}, Ph)     # R >= 2: R -= 2, A += 1
        rule(P, {'R': 1}, {}, B1)           # R = 1: read bit 1
        rule(P, {}, {}, B0)                 # R = 0: read bit 0
        rule(Ph, {}, {}, P)
    for (q, g), act in tab.items():         # table order of the Lean file
        if act is None:
            continue
        nq, mv, w = act
        B, T = f'B{g}_{q}', f'P{nq}'
        k = f'{q}_{g}'
        if mv == 'R':                       # R := A ; L := 2L + w
            rule(B, {'A': 1}, {'R': 1}, 'X1h' + k); rule(B, {}, {}, 'X2' + k); rule('X1h' + k, {}, {}, B)
            rule('X2' + k, {'L': 1}, {'A': 2}, 'X2h' + k); rule('X2' + k, {}, {}, 'X3' + k); rule('X2h' + k, {}, {}, 'X2' + k)
            rule('X3' + k, {'A': 1}, {'L': 1}, 'X3h' + k); rule('X3' + k, {}, {'L': w} if w else {}, T)
            rule('X3h' + k, {}, {}, 'X3' + k)
        else:                               # L := L div 2, R := 4A + 2w + (L mod 2)
            rule(B, {'L': 2}, {'C': 1}, 'Y1h' + k); rule(B, {'L': 1}, {'R': 1}, 'Y2' + k); rule(B, {}, {}, 'Y2' + k)
            rule('Y1h' + k, {}, {}, B)
            rule('Y2' + k, {'C': 1}, {'L': 1}, 'Y2h' + k); rule('Y2' + k, {}, {}, 'Y3' + k); rule('Y2h' + k, {}, {}, 'Y2' + k)
            rule('Y3' + k, {'A': 1}, {'R': 4}, 'Y3h' + k); rule('Y3' + k, {}, {'R': 2 * w} if w else {}, T)
            rule('Y3h' + k, {}, {}, 'Y3' + k)
    pr = dict(zip(names, primes_after7(len(names))))
    return pr, rules


def to_fraction(pr, r):
    s, guard, inc, nxt = r
    num, den = pr[nxt], pr[s]
    for k, e in inc.items():
        num *= REG[k] ** e
    for k, e in guard.items():
        den *= REG[k] ** e
    g = math.gcd(num, den)
    return num // g, den // g


def factor_over(m, plist):
    out = {}
    for p in plist:
        while m % p == 0:
            m //= p
            out[p] = out.get(p, 0) + 1
    assert m == 1
    return out


class Vec:
    """exponent-vector FRACTRAN on a one-hot program: from control prime s only fractions with s in the denominator can
    apply (M3 checks every denominator has exactly one control prime), taken in global list order."""

    def __init__(self, fracs, ctrl, plist):
        self.f = []
        self.by = defaultdict(list)
        for idx, (num, den) in enumerate(fracs):
            fn, fd = factor_over(num, plist), factor_over(den, plist)
            sd = [p for p in fd if p in ctrl][0]
            sn = [p for p in fn if p in ctrl][0]
            guard = {p: e for p, e in fd.items() if p not in ctrl}
            inc = {p: e for p, e in fn.items() if p not in ctrl}
            self.by[sd].append((idx, guard, inc, sn))

    def first(self, s, reg):
        for k, (idx, guard, inc, sn) in enumerate(self.by.get(s, [])):
            if all(reg.get(p, 0) >= e for p, e in guard.items()):
                return k, idx, guard, inc, sn
        return None

    def apply(self, reg, guard, inc, times=1):
        for p, e in guard.items():
            reg[p] -= e * times
        for p, e in inc.items():
            reg[p] = reg.get(p, 0) + e * times

    def run(self, s, reg, accelerate, on_step=None, stop=None):
        """single or accelerated steps; on_step(prev, s, reg, steps) after each real step; stop() -> bool"""
        steps = 0
        while True:
            got = self.first(s, reg)
            if got is None:
                return s, reg, steps, True
            k, idx, guard, inc, sn = got
            if accelerate and k == 0 and len(self.by.get(sn, [])) == 1:
                _, g2, i2, back = self.by[sn][0]
                if back == s and not g2:
                    d = defaultdict(int)
                    for p, e in guard.items():
                        d[p] -= e
                    for p, e in inc.items():
                        d[p] += e
                    for p, e in i2.items():
                        d[p] += e
                    lim = [(reg.get(p, 0) - e) // (-d[p]) + 1 for p, e in guard.items() if d[p] < 0]
                    if lim:
                        kk = min(lim)
                        if kk > 1:
                            for p, e in d.items():
                                reg[p] = reg.get(p, 0) + e * kk
                            steps += 2 * kk
                            continue
            self.apply(reg, guard, inc)
            prev, s = s, sn
            steps += 1
            if on_step is not None:
                on_step(prev, s, reg, steps)
            if stop is not None and stop():
                return s, reg, steps, False


def per_step_cost(L, R, mv):
    A = R // 2
    extract = 2 * A + 1
    if mv == 'R':
        return extract + (2 * A + 1) + (2 * L + 1) + (4 * L + 1)
    return extract + (2 * (L // 2) + 1) + (2 * (L // 2) + 1) + (2 * A + 1)


def main():
    t0 = time.time()
    print('M1  Leng table')
    tab = parse_lean()
    nR = sum(1 for a in tab.values() if a and a[1] == 'R')
    nL = sum(1 for a in tab.values() if a and a[1] == 'L')
    undefined = [k for k, a in tab.items() if a is None]
    report('50 slots, 49 instructions, 22 right and 27 left moves, (24,0) the only undefined slot',
           len(tab) == 50 and (nR, nL) == (22, 27) and undefined == [(24, 0)], f'R={nR} L={nL} undefined={undefined}')

    print('M2  TM -> radius-1 CA (full-window updates)')
    quiet, S, t, why = ca_run(tab, 25, 200_000)
    report('alphabet size 53 = 2 + 25*2 + 1; rule quiescent away from heads (exhaustive)', S == 53 and quiet, f'{S}')
    report('CA = TM at every step for 200,000 steps (full updates; full tape compared every 100 steps)',
           why == 'ran' and t == 200_000, f'{why} at t={t}')
    q2, S2, t2, why2 = ca_run(bb2(), 3, 50, check_every=1)
    report('BB(2) CA (9 symbols) writes HALT at step 7', why2 == 'halt' and t2 == 7 and S2 == 9, f'{why2} at t={t2}')

    print('M3  TM -> FRACTRAN (THM-4858 (b) construction)')
    pr, rules = compile_tm(tab, 25)
    fr = [to_fraction(pr, r) for r in rules]
    ctrl = set(pr.values())
    plist = [2, 3, 5, 7] + sorted(ctrl)
    used = set()
    for num, den in fr:
        used |= set(factor_over(num, plist)) | set(factor_over(den, plist))
    report('568 fractions = 25*4 + 22*9 + 27*10; 349 primes (4 registers + 345 control); largest 2351',
           len(fr) == 568 == 25 * 4 + 22 * 9 + 27 * 10 and len(used) == 349 and len(ctrl) == 345 and max(used) == 2351,
           f'{len(fr)} fractions, {len(used)} primes, {len(ctrl)} control, largest {max(used)}')
    printed = '65/99, 19/33, 17/11, 11/13, 145/207, 37/69, 31/23, 23/29, 215/369, 53/123, 47/41, 41/43'
    mine = ', '.join(f'{a}/{b}' for a, b in fr[:12])
    report('first 12 fractions equal those printed in machine_translations.out', mine == printed, mine)
    onehot = True
    for (num, den), r in zip(fr, rules):
        fn, fd = factor_over(num, plist), factor_over(den, plist)
        cn = [p for p in fn if p in ctrl]
        cd = [p for p in fd if p in ctrl]
        regs_n = {p for p in fn if p not in ctrl}
        regs_d = {p for p in fd if p not in ctrl}
        if len(cn) != 1 or len(cd) != 1 or fn[cn[0]] != 1 or fd[cd[0]] != 1 or cn == cd or regs_n & regs_d:
            onehot = False
        # reduction must not have cancelled anything
        num0, den0 = pr[r[3]], pr[r[0]]
        for k, e in r[2].items():
            num0 *= REG[k] ** e
        for k, e in r[1].items():
            den0 *= REG[k] ** e
        if (num0, den0) != (num, den):
            onehot = False
    report('one-hot after reduction: exactly one control prime (exponent 1) in each numerator and denominator, '
           'different ones; no register on both sides (nothing cancels)', onehot)
    by = defaultdict(list)
    for r in rules:
        by[r[0]].append(r)
    halting = [s for s in pr if s not in by]
    last_uncond = all(not lst[-1][1] for lst in by.values())
    report('every state with fractions ends with an unconditional fraction; the only fraction-less state is B0_24',
           last_uncond and halting == ['B0_24'], f'halting states {halting}')

    print('M4  own vector interpreter, accelerated, vs own TM')
    V = Vec(fr, ctrl, plist)
    Pinv = {pr[f'P{q}']: q for q in range(25)}
    helpers = {pr[f'Ph{q}'] for q in range(25)}
    rec = []

    def on_step(prev, s, reg, steps):
        if s in Pinv and prev not in helpers:
            rec.append((Pinv[s], reg.get(2, 0), reg.get(3, 0), steps))
    N = 5000
    s, reg, steps, halted = V.run(pr['P0'], {}, True, on_step, lambda: len(rec) >= N)
    tm = TM(tab)
    ok = True
    cum = 0
    cost3000 = None
    for k in range(N):
        L0, R0 = tm.stacks()
        mv = tab[(tm.q, tm.tape[tm.h])][1]
        cum += per_step_cost(L0, R0, mv)
        tm.step()
        if (tm.q, *tm.stacks()) != rec[k][:3] or rec[k][3] != cum:
            ok = False
            print('    first disagreement at TM step', k + 1, rec[k], (tm.q, tm.stacks(), cum))
            break
        if k + 1 == 3000:
            cost3000 = (rec[k][3], cum)
    report(f'stacks and control state at every instruction exit = TM after k steps, k = 1..{N}', ok and len(rec) == N)
    report('FRACTRAN steps for the first 3,000 TM steps = 24,405,662 (interpreter count = closed-form count)',
           cost3000 == (24405662, 24405662), f'{cost3000}')
    print(f'    FRACTRAN steps for the first {N} TM steps: {rec[-1][3]}; largest register then: '
          f'{max(rec[-1][1], rec[-1][2])}')

    print('M5  literal big-integer FRACTRAN vs the vector interpreter (single steps)')
    tm = TM(tab)
    target = sum(0 for _ in range(0))
    cost300 = 0
    for k in range(300):
        L0, R0 = tm.stacks()
        cost300 += per_step_cost(L0, R0, tab[(tm.q, tm.tape[tm.h])][1])
        tm.step()
    trace = []
    V.run(pr['P0'], {}, False, lambda prev, s, reg, st: trace.append((s, dict(reg))), lambda: len(trace) >= cost300)
    n = pr['P0']
    lit_ok = True
    dens = [d for _, d in fr]
    nums = [a for a, _ in fr]
    for i in range(cost300):
        for j in range(568):
            if n % dens[j] == 0:
                n = n // dens[j] * nums[j]
                break
        else:
            lit_ok = False
            break
        s_i, reg_i = trace[i]
        want = s_i
        for p, e in reg_i.items():
            want *= p ** e
        if n != want:
            lit_ok = False
            print('    literal/vector disagreement at FRACTRAN step', i + 1)
            break
    report(f'literal integer FRACTRAN = vector interpreter at every one of the {cost300} configurations of the first '
           f'300 TM steps', lit_ok)
    # what the session's T2c covered: run until 300 P-arrivals INCLUDING helper returns, unaccelerated
    arrivals = [0]
    exits = [0]
    count_steps = [0]

    def t2c_step(prev, s, reg, st):
        if s in Pinv:
            arrivals[0] += 1
            if prev not in helpers:
                exits[0] += 1
        count_steps[0] = st
    V.run(pr['P0'], {}, False, t2c_step, lambda: arrivals[0] >= 300)
    report('the session\'s T2c ("first 300 TM steps") stopped after 1,948 FRACTRAN steps', count_steps[0] == 1948,
           f'{count_steps[0]} FRACTRAN steps = only {exits[0]} TM steps')

    print('M6  the failed T2b / T3 numbers from the misaligned "entry" convention (accelerated; helper returns '
          'recorded when the extraction loop runs once)')
    arr = []

    def mis_step(prev, s, reg, st):
        if s in Pinv:
            arr.append((Pinv[s], reg.get(2, 0), reg.get(3, 0), st, prev in helpers))
    V.run(pr['P0'], {}, True, mis_step, lambda: len(arr) >= 3000)
    true_steps = sum(1 for a in arr if not a[4])
    report('3,000 misaligned entries end at FRACTRAN step 20,789,496 (the failed T2b figure)', arr[-1][3] == 20789496,
           f'{arr[-1][3]}; they cover only {true_steps} TM steps')
    report('misaligned entry 3 is (P23, L=0, R=1), a helper return, vs TM step 3 = (23, L=0, R=6)',
           arr[2][:3] == (23, 0, 1) and arr[2][4], f'{arr[2]}')

    print('M7  BB(2) compiled')
    btab = bb2()
    tm = TM(btab)
    k = 0
    while tm.step():
        k += 1
    ones_tm = sum(tm.tape.values())
    pb, rb = compile_tm(btab, 3)
    fb = [to_fraction(pb, r) for r in rb]
    cb = set(pb.values())
    Vb = Vec(fb, cb, [2, 3, 5, 7] + sorted(cb))
    s, reg, steps, halted = Vb.run(pb['P0'], {}, False)
    inv = {v: k2 for k2, v in pb.items()}
    g = 1 if inv[s].startswith('B1') else 0
    ones = bin(reg.get(2, 0)).count('1') + g + bin(reg.get(5, 0)).count('1') + bin(reg.get(3, 0)).count('1')
    report('BB(2): TM halts after 6 steps with four 1s; compiled FRACTRAN halts in B1_2 with four 1s (L, head bit, A)',
           k == 6 and ones_tm == 4 and halted and inv[s] == 'B1_2' and ones == 4,
           f'TM {k} steps, {ones_tm} ones; FRACTRAN {steps} steps, state {inv[s]}, regs {reg}, ones {ones}')
    print()
    print(f'elapsed {time.time() - t0:.1f} s')
    print('ALL AUDIT CHECKS AGREE' if not FAILS else f'MISMATCHES: {FAILS}')
    return 0 if not FAILS else 1


if __name__ == '__main__':
    sys.exit(main())
