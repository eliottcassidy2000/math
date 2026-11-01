#!/usr/bin/env python3
"""Independent audit of THM-4856 (GOLDBACHGAME, FERMATGAME, LEGENDREGAME).

Written from scratch; imports nothing from the session. The fraction lists are PARSED from the session's
published text (THM-4856 for GOLDBACHGAME, primekernel_programs.out for the other two), never rebuilt from
the session's generator code. The state/register split is DERIVED from the fractions by constraint
propagation (each numerator and denominator must contain exactly one state prime), then checked.

Three interpreters:
  lit   - literal FRACTRAN on Python integers (n -> n*f for the first f with n*f integral);
  vec   - the same rule on exact prime-exponent vectors (cross-checked step by step against lit);
  fast  - vec plus exact acceleration of 2-cycles s -> t -> s whose rule choices provably stay fixed.
Semantic check: every kernel call (entry into state 11 from the controller, exit 17 -> X) is logged, the
kernel contract (pre- and post-conditions) is checked on every call, and the sequence of tested numbers and
verdicts is compared with the sequence predicted directly from number theory.

usage (from the repo root; the defaults reproduce audit_fractran.out in about 3 minutes):
  python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit/audit_fractran.py [GOLD_NMAX=1000] [LEG_NMAX=120] [FERMAT_M4=1]
Literal integer runs: GOLDBACHGAME to the acceptance of N = 100, FERMATGAME from m = 0..3 to the halt,
LEGENDREGAME windows n = 1..12. Accelerated runs: to the CLI bounds, with exact step counts.
"""
import re
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
EXP = ROOT / '04-computation/experiments/circulant_fractran_goldbach_20261101'
THM = ROOT / '01-canon/theorems/THM-4856-goldbachgame-fermatgame-legendregame-halts-iff-counterexample-fractran-programs-built-on-primegames-divisor-loop.md'
NOTE = ROOT / '05-knowledge/results/fractran_turing_circulant_goldbach_20261101.md'

FAILS, MISMATCH = [], []


def report(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


def claim(name, ok, detail=''):
    print(f'  [{"CLAIM OK" if ok else "CLAIM WRONG"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        MISMATCH.append(name)


# ----------------------------------------------------------------------------- number theory helpers
def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def lpd(n):
    """largest proper divisor D(n) = n / spf(n), n >= 2."""
    d = 2
    while d * d <= n:
        if n % d == 0:
            return n // d
        d += 1
    return 1


# ----------------------------------------------------------------------------- parsing
def parse_list(text):
    return [Fraction(int(a), int(b)) for a, b in re.findall(r'(\d+)/(\d+)', text)]


def programs_from_files():
    thm = THM.read_text()
    blocks = re.findall(r'```text\n(.*?)```', thm, re.S)
    gold_thm = [b for b in blocks if '705/1333' in b]
    note = NOTE.read_text()
    gold_note = [b for b in re.findall(r'```text\n(.*?)```', note, re.S) if '705/1333' in b]
    gout = (EXP / 'goldbachgame.out').read_text()
    gold_out = re.search(r'GOLDBACHGAME fractions \(in order\):\n(.*?)\n', gout).group(1)
    pk = (EXP / 'primekernel_programs.out').read_text()
    ferm = re.search(r'FERMATGAME fractions: (.*?)\n', pk).group(1)
    leg = re.search(r'LEGENDREGAME fractions: (.*?)\n', pk).group(1)
    G = parse_list(gold_thm[0])
    same = parse_list(gold_note[0]) == G and parse_list(gold_out) == G
    return G, parse_list(ferm), parse_list(leg), same, (gold_thm[0], ferm, leg)


# ----------------------------------------------------------------------------- program structure
class Program:
    def __init__(self, fracs, name):
        self.name = name
        self.fracs = fracs
        self.nf = [factor(f.numerator) for f in fracs]
        self.df = [factor(f.denominator) for f in fracs]
        self.primes = sorted(set(p for d in self.nf + self.df for p in d))

    def derive_states(self, seed_states=()):
        """constraint propagation: every side has exactly one state prime, with exponent 1."""
        state, reg = set(seed_states), set()
        sides = self.nf + self.df
        for d in sides:
            for p, e in d.items():
                if e >= 2:
                    reg.add(p)
        changed = True
        while changed:
            changed = False
            for d in sides:
                cand = [p for p, e in d.items() if e == 1 and p not in reg]
                known = [p for p in d if p in state]
                if len(known) == 0 and len(cand) == 1 and cand[0] not in state:
                    state.add(cand[0])
                    changed = True
                if len(known) == 1:
                    for p in d:
                        if p not in state and p not in reg:
                            reg.add(p)
                            changed = True
        undecided = [p for p in self.primes if p not in state and p not in reg]
        return state, reg, undecided

    def build(self, states, regs):
        self.states = sorted(states)
        self.regs = sorted(regs)
        self.sidx = {p: i for i, p in enumerate(self.states)}
        self.ridx = {p: i for i, p in enumerate(self.regs)}
        self.rules = [[] for _ in self.states]      # per state: (guard tuple, delta tuple, next, frac index)
        self.onehot_ok = True
        self.problems = []
        for i, (n, d) in enumerate(zip(self.nf, self.df)):
            sd = [p for p in d if p in states]
            sn = [p for p in n if p in states]
            if len(sd) != 1 or len(sn) != 1 or d[sd[0]] != 1 or n[sn[0]] != 1 or sd[0] == sn[0]:
                self.onehot_ok = False
                self.problems.append((str(self.fracs[i]), sd, sn))
                continue
            guard = tuple((self.ridx[p], e) for p, e in d.items() if p in regs)
            delta = [0] * len(self.regs)
            for p, e in d.items():
                if p in regs:
                    delta[self.ridx[p]] -= e
            for p, e in n.items():
                if p in regs:
                    delta[self.ridx[p]] += e
            self.rules[self.sidx[sd[0]]].append((guard, tuple(delta), self.sidx[sn[0]], i))

    def encode(self, state_prime, regvals):
        n = state_prime
        for p, e in regvals.items():
            n *= p ** e
        return n

    def encode_state(self, s, R):
        n = self.states[s]
        for p, e in zip(self.regs, R):
            if e:
                n *= p ** e
        return n

    def decode(self, n):
        st, R = None, [0] * len(self.regs)
        for p in self.states:
            if n % p == 0:
                if st is not None or (n // p) % p == 0:
                    raise ValueError('not one-hot')
                st = p
                n //= p
        for p in self.regs:
            while n % p == 0:
                R[self.ridx[p]] += 1
                n //= p
        if n != 1:
            raise ValueError('foreign prime')
        return st, R


# ----------------------------------------------------------------------------- interpreters
def lit_run(prog, n, nsteps):
    """literal FRACTRAN; yields successive integers."""
    fr = [(f.numerator, f.denominator) for f in prog.fracs]
    for _ in range(nsteps):
        for a, b in fr:
            if n % b == 0:
                n = n // b * a
                break
        else:
            return
        yield n


def vec_step(prog, s, R):
    for guard, delta, t, fi in prog.rules[s]:
        if all(R[r] >= a for r, a in guard):
            return t, [x + y for x, y in zip(R, delta)], fi
    return None


def lit_lockstep(prog, n, s, R, nsteps, hook=None, check_every=2000):
    """literal integer FRACTRAN in lockstep with the exponent-vector machine: at every step the index of the
    first applicable fraction (integer divisibility test on the full list) must equal the vector machine's
    chosen rule; every `check_every` steps and at the end the integer must equal the encoding of (s, R).
    Returns (status, s, R, steps, n, agree)."""
    fr = [(f.numerator, f.denominator) for f in prog.fracs]
    agree = True
    steps = 0
    status = 'BUDGET'
    while steps < nsteps:
        fi_lit = -1
        for i, (a, b) in enumerate(fr):
            if n % b == 0:
                fi_lit = i
                break
        got = vec_step(prog, s, R)
        if fi_lit < 0 or got is None:
            if not (fi_lit < 0 and got is None):
                agree = False
            status = 'HALT'
            break
        t, R2, fi_vec = got
        if fi_vec != fi_lit:
            agree = False
            break
        a, b = fr[fi_lit]
        n = n // b * a
        prev, s, R = s, t, R2
        steps += 1
        if steps % check_every == 0 and n != prog.encode_state(s, R):
            agree = False
            break
        if hook is not None and hook(prev, s, R, steps):
            status = 'STOP'
            break
    if n != prog.encode_state(s, R):
        agree = False
    return status, s, R, steps, n, agree


class Fast:
    """exact simulation with 2-cycle acceleration; calls hook(prev, new, R, steps) on every single
    (non-accelerated) transition into a watched state."""

    def __init__(self, prog, watched, hook, kernel):
        self.p = prog
        self.watched = watched
        self.hook = hook
        self.kernel = kernel          # set of state indices; 2-cycles never cross the kernel boundary

    def first(self, s, R):
        for j, (guard, delta, t, fi) in enumerate(self.p.rules[s]):
            ok = True
            for r, a in guard:
                if R[r] < a:
                    ok = False
                    break
            if ok:
                return j
        return -1

    def caps(self, s, i, X, D):
        """iterations l = 0..k-1 (state s, registers X + l D) keep rule i the first applicable one."""
        k = None
        rules = self.p.rules[s]
        for r, a in rules[i][0]:
            if D[r] < 0:
                c = (X[r] - a) // (-D[r]) + 1
                k = c if k is None else min(k, c)
        for j in range(i):
            guard = rules[j][0]
            safe = False
            need = 0
            for r, a in guard:
                if X[r] < a:
                    if D[r] <= 0:
                        safe = True
                        break
                    need = max(need, -((X[r] - a) // D[r]))     # ceil((a - X)/D)
            if not safe:
                k = need if k is None else min(k, need)
        return k

    def run(self, s, R, maxsteps, stop=None):
        p = self.p
        steps = 0
        while steps < maxsteps:
            i = self.first(s, R)
            if i < 0:
                return 'HALT', s, R, steps
            guard, delta, t, fi = p.rules[s][i]
            R1 = [x + y for x, y in zip(R, delta)]
            if (s in self.kernel) == (t in self.kernel):
                j = self.first(t, R1)
                if j >= 0 and p.rules[t][j][2] == s:
                    D = [x + y for x, y in zip(delta, p.rules[t][j][1])]
                    k1 = self.caps(s, i, R, D)
                    k2 = self.caps(t, j, R1, D)
                    ks = [k for k in (k1, k2) if k is not None]
                    if not ks:
                        return 'LOOP', s, R, steps          # a provably infinite 2-cycle
                    k = min(ks)
                    if k >= 2:
                        R = [x + k * d for x, d in zip(R, D)]
                        steps += 2 * k
                        continue
            prev, s, R = s, t, R1
            steps += 1
            if s in self.watched or (prev in self.kernel and s not in self.kernel):
                if self.hook(prev, s, R, steps):
                    return 'STOP', s, R, steps
        return 'BUDGET', s, R, steps


# ----------------------------------------------------------------------------- kernel-call logger
class KernelLog:
    def __init__(self, prog, kernel_states, exit_state, r):
        self.prog = prog
        self.K = kernel_states
        self.X = exit_state
        self.r = r                    # dict name -> register index for 2,3,5,7
        self.calls = []               # (M, verdict_prime, entry_step, exit_step)
        self.contract_bad = []
        self.cur = None

    def on(self, prev, s, R, steps):
        r2, r3, r5, r7 = (R[self.r[k]] for k in (2, 3, 5, 7))
        if s in self.K and prev not in self.K:                     # kernel entry
            M = r5
            if not (r2 == 0 and r7 == 0 and r3 == M - 1 and M >= 2) or self.cur is not None:
                self.contract_bad.append(('entry', M, r2, r3, r5, r7, steps))
            self.cur = (M, steps)
        elif prev in self.K and s not in self.K:                   # kernel exit
            M, t0 = self.cur
            self.cur = None
            ok = s == self.X and r2 == M and r3 == 0 and r5 == 0 and r7 == lpd(M) - 1
            if not ok:
                self.contract_bad.append(('exit', M, r2, r3, r5, r7, steps))
            self.calls.append((M, r7 == 0, t0, steps))
        return False


# ----------------------------------------------------------------------------- predictions
def goldbach_tests(Nmax):
    seq = []
    for N in range(4, Nmax + 1, 2):
        found = False
        for P in range(2, N - 1):
            seq.append(P)
            if is_prime(P):
                seq.append(N - P)
                if is_prime(N - P):
                    found = True
                    break
        if not found:
            return seq, N
    return seq, None


def legendre_tests(nmax):
    seq = []
    for n in range(1, nmax + 1):
        for C in range(n * n + 1, (n + 1) ** 2):
            seq.append(C)
            if is_prime(C):
                break
    return seq


def common_checks(name, prog, start_int, want_nf, want_np, gold_prefix=True):
    print(f'--- {name}: {len(prog.fracs)} fractions, {len(prog.primes)} primes')
    claim(f'{name}: {want_nf} fractions over {want_np} primes', len(prog.fracs) == want_nf and len(prog.primes) == want_np,
          f'{len(prog.fracs)} fractions, {len(prog.primes)} primes')
    reduced = all(f.numerator == int(re.sub(r'/.*', '', str(f))) for f in prog.fracs)
    kernel = [Fraction(17, 91), Fraction(78, 85), Fraction(19, 51), Fraction(23, 38), Fraction(29, 33),
              Fraction(77, 29), Fraction(95, 23), Fraction(77, 19), None, Fraction(11, 13), Fraction(13, 11)]
    ok = all(k is None or prog.fracs[i] == k for i, k in enumerate(kernel)) and prog.fracs[8].denominator == 17
    report(f'{name}: first eleven = PRIMEGAME\'s with 1/17 -> X/17 (X = {prog.fracs[8].numerator})', ok and reduced)
    states, regs, und = prog.derive_states()
    report(f'{name}: state/register split derived by propagation, nothing undecided', not und,
           f'states {sorted(states)}; registers {sorted(regs)}')
    prog.build(states, regs)
    report(f'{name}: one-hot: every reduced numerator and denominator has exactly one state prime (exp 1), '
           f'no self-loop', prog.onehot_ok, f'problems {prog.problems}' if prog.problems else
           f'{len(states)} states, {len(regs)} registers')
    st, R = prog.decode(start_int)
    report(f'{name}: start integer {start_int if start_int < 10**12 else "(large)"} is one-hot', st is not None)
    kern_idx = {prog.sidx[p] for p in (11, 13, 17, 19, 23, 29)}
    helpers = [s for s in range(len(prog.states)) if len(prog.rules[s]) == 1 and s not in kern_idx]
    hok = all(not prog.rules[s][0][0] and not any(prog.rules[s][0][1]) and
              any(r[2] == s for r in prog.rules[prog.rules[s][0][2]]) for s in helpers)
    report(f'{name}: the single-fraction controller states {[prog.states[s] for s in helpers]} are pure returns '
           f'(s/s\') to a state that calls them', hok)
    # every controller 2-cycle s -> helper -> s strictly decreases a guarded register (terminates)
    kern = {prog.sidx[p] for p in (11, 13, 17, 19, 23, 29)}
    loops_ok = True
    for s in range(len(prog.states)):
        for guard, delta, t, fi in prog.rules[s]:
            if len(prog.rules[t]) == 1 and prog.rules[t][0][2] == s:
                if not any(delta[r] < 0 for r, a in guard):
                    loops_ok = False
    report(f'{name}: every helper loop decrements one of its own guard registers (each loop terminates)', loops_ok)
    return kern


# ----------------------------------------------------------------------------- GOLDBACHGAME
def goldbach(G, nmax, nlit):
    prog = Program(G, 'GOLDBACHGAME')
    kern = common_checks('GOLDBACHGAME', prog, 43, 37, 30)
    S = prog.sidx
    X = S[prog.fracs[8].numerator]
    r = {p: prog.ridx[p] for p in (2, 3, 5, 7)}
    Pp, Qp = prog.ridx[31], prog.ridx[37]
    claim('GOLDBACHGAME registers are 2,3,5,7,31,37,41 and 23 states', prog.regs == [2, 3, 5, 7, 31, 37, 41]
          and len(prog.states) == 23, f'registers {prog.regs}')
    st, R = prog.decode(43)
    report('start 43 = state TP with P\' = Q\' = 0 (N = 4)', st == 43 and not any(R))
    # SUCC = the state whose first rule moves 31 into 37 (THM: Q' := P' + Q' + 2)
    succ = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((Pp, 1),)
            and prog.rules[s][0][1][Qp] == 1]
    succ = succ[0]
    succh = prog.rules[succ][0][2]
    # ---- literal integer run in lockstep with the exponent-vector machine, decoded semantically
    t0 = time.time()
    log = KernelLog(prog, kern, X, r)
    accepted_lit = []

    def hook_lit(prev, s, R, steps):
        log.on(prev, s, R, steps)
        if s == succ and prev != succh:
            accepted_lit.append((R[Pp] + R[Qp] + 4, R[Pp] + 2, steps))
            return R[Pp] + R[Qp] + 4 >= nlit
        return False
    st_, s_, R_, steps, n_end, agree = lit_lockstep(prog, 43, S[43], [0] * len(prog.regs), 10 ** 9, hook_lit)
    report(f'literal integer FRACTRAN = exponent-vector machine at every step (first applicable fraction identical; '
           f'integer re-encoded every 2000 steps)', agree, f'{steps} literal steps to the acceptance of N = {nlit} '
                                                          f'({time.time() - t0:.0f}s)')
    pred, _ = goldbach_tests(1000)
    obs = [m for m, v, a, b in log.calls]
    report(f'literal run: kernel contract holds on every call; tests = predicted sequence',
           not log.contract_bad and obs == pred[:len(obs)],
           f'{len(obs)} kernel calls, accepted N up to {accepted_lit[-1][0]} in {steps} literal steps')
    # ---- fast run to nmax
    log = KernelLog(prog, kern, X, r)
    acc = []

    def hook(prev, s, R, steps):
        log.on(prev, s, R, steps)
        if s == succ and prev != succh:
            N = R[Pp] + R[Qp] + 4                 # at SUCC entry: N = P + Q with P' = P - 2, Q' = Q - 2
            acc.append((N, R[Pp] + 2, steps))
            return N >= nmax
        return False
    watched = set(kern) | {succ}
    t0 = time.time()
    res = Fast(prog, watched, hook, kern).run(S[43], [0] * len(prog.regs), 10 ** 15)
    pred, halt_at = goldbach_tests(nmax)
    obs = [m for m, v, a, b in log.calls]
    verd = all(v == is_prime(m) for m, v, a, b in log.calls)
    least = all(P == min(p for p in range(2, N - 1) if is_prime(p) and is_prime(N - p)) for N, P, _ in acc)
    report(f'fast run to N = {nmax}: no halt, targets 4..{nmax} in order, least-prime splits, kernel contract on '
           f'every call, tests = predicted', res[0] == 'STOP' and [a[0] for a in acc] == list(range(4, nmax + 1, 2))
           and least and verd and not log.contract_bad and obs == pred[:len(obs)],
           f'{len(log.calls)} kernel calls, {res[3]} steps, {time.time() - t0:.0f}s')
    lit_prefix = [(N, st) for N, P, st in accepted_lit]
    fast_prefix = [(N, st) for N, P, st in acc[:len(lit_prefix)]]
    report('literal and fast step counts agree at every accepted N of the literal prefix', lit_prefix == fast_prefix)
    cum = {N: st for N, P, st in acc}
    sess = [(4, 50), (6, 269), (8, 670), (10, 1325), (12, 2735), (14, 4120), (16, 5975), (18, 9279)]
    claim('session G5 cumulative steps to accept N = 4..18', all(cum.get(N) == v for N, v in sess),
          str([(N, cum.get(N)) for N, _ in sess]))
    per = {N: cum[N] - cum[N - 2] for N in cum if N - 2 in cum}
    for N, v in ((100, 86987), (200, 354081), (400, 1428647), (800, 5738237), (1000, 8973223)):
        if N in per:
            claim(f'session per-target steps at N = {N}: {v}', per[N] == v, f'audit {per[N]}  (= {per[N] / N**2:.2f} N^2)')
    if 1000 in cum:
        claim('session G3 total steps to accept N = 1000: 2,991,801,422', cum[1000] == 2991801422, str(cum[1000]))
    # ---- negative control G4: H/(NEXT * 31) inserted first in NEXT's group
    nxt = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((Qp, 1),)
           and prog.rules[s][0][1][Pp] == 1]
    nxt = nxt[0]
    H = [s for s in range(len(prog.states)) if not prog.rules[s]][0]
    fi_first = prog.rules[nxt][0][3]
    Gc = list(G)
    Gc.insert(fi_first, Fraction(prog.states[H], prog.states[nxt] * 31))
    progc = Program(Gc, 'G4')
    progc.build(set(prog.states), set(prog.regs))
    logc = KernelLog(progc, kern, X, r)
    accc = []

    def hookc(prev, s, R, steps):
        logc.on(prev, s, R, steps)
        if s == succ and prev != succh:
            accc.append(R[Pp] + R[Qp] + 4)
        return False
    resc = Fast(progc, set(kern) | {succ}, hookc, kern).run(S[43], [0] * len(prog.regs), 10 ** 9)
    Rh = resc[2]
    last_pair = (logc.calls[-2][0], logc.calls[-1][0]) if len(logc.calls) >= 2 else None
    print(f'     G4 control: {resc[0]} in state {progc.states[resc[1]]} after {resc[3]} steps; accepted {accc}; '
          f'P\' = {Rh[Pp]}, Q\' = {Rh[Qp]} at the halt; last two tests {last_pair}')
    report('G4 control halts while processing N = 12 (accepted 4,6,8,10; last split (3,9); P\'+Q\'+4 = 11 because the '
           'guard consumed one unit of P\')', resc[0] == 'HALT' and progc.states[resc[1]] == prog.states[H]
           and accc == [4, 6, 8, 10] and last_pair == (3, 9) and Rh[Pp] + Rh[Qp] + 4 == 11)
    claim('session G4 halt after 2174 steps', resc[3] == 2174, str(resc[3]))


# ----------------------------------------------------------------------------- FERMATGAME
def fermat(F, do_m4):
    prog = Program(F, 'FERMATGAME')
    start5 = 53 * 31 ** 32 * 43
    kern = common_checks('FERMATGAME', prog, start5, 36, 32)
    S = prog.sidx
    X = S[prog.fracs[8].numerator]
    r = {p: prog.ridx[p] for p in (2, 3, 5, 7)}
    claim('FERMATGAME: 9 registers, 23 states', len(prog.regs) == 9 and len(prog.states) == 23,
          f'registers {prog.regs}')
    st, R = prog.decode(start5)
    report('start 53*31^32*43 = state 53 with register 31 = 32 = 2^5 and register 43 = 1',
           st == 53 and R[prog.ridx[31]] == 32 and R[prog.ridx[43]] == 1 and sum(R) == 33)
    sess = {0: 66, 1: 185, 2: 1911, 3: 400283, 4: 25772326379}
    for m in range(0, 5 if do_m4 else 4):
        start = 53 * 31 ** (2 ** m) * 43
        log = KernelLog(prog, kern, X, r)
        if m <= 3:
            t0 = time.time()
            st_, s_, R_, res_steps, n_end, agree = lit_lockstep(prog, start, S[53], prog.decode(start)[1], 10 ** 8,
                                                                log.on)
            halted_lit = st_ == 'HALT' and agree
            res_state = prog.states[s_]
            mode = f'literal integer run in lockstep, agree={agree}, {time.time() - t0:.0f}s'
        else:
            t0 = time.time()
            res = Fast(prog, set(kern), log.on, kern).run(S[53], prog.decode(start)[1], 10 ** 15)
            halted_lit = res[0] == 'HALT'
            res_steps, res_state = res[3], prog.states[res[1]]
            mode = f'fast, {time.time() - t0:.0f}s'
        Fm = 2 ** (2 ** m) + 1
        report(f'FERMATGAME from m = {m}: halts after testing exactly F_{m} = {Fm} (prime), kernel contract held',
               halted_lit and [c[0] for c in log.calls] == [Fm] and log.calls[0][1] and not log.contract_bad,
               f'{res_steps} steps, halted in state {res_state} ({mode})')
        claim(f'session step count from m = {m}: {sess[m]}', res_steps == sess[m], str(res_steps))
    # inverted control: swap X's two rules' targets (r7>0 -> H, else -> continue)
    Xs = prog.rules[X]
    H = [s for s in range(len(prog.states)) if not prog.rules[s]][0]
    cont = Xs[0][2]          # X: r7 > 0 -> CL (composite path) ; else -> H
    # the composite path is CL -> C2 -> MK; inverted: r7 > 0 -> H, else -> the state after CL's loop (C2)
    CL = cont
    C2 = prog.rules[CL][-1][2]
    Fi = list(F)
    Fi[Xs[0][3]] = Fraction(prog.states[H] * 1, prog.states[X] * 7)
    Fi[Xs[1][3]] = Fraction(prog.states[C2], prog.states[X])
    progi = Program(Fi, 'FERMAT-inverted')
    progi.build(set(prog.states), set(prog.regs))
    logi = KernelLog(progi, kern, X, r)

    def hooki(prev, s, R, steps):
        logi.on(prev, s, R, steps)
        return s == S[11] and prev not in kern and R[r[5]] == 2 ** 32 + 1
    res = Fast(progi, set(kern), hooki, kern).run(S[53], prog.decode(53 * 31 * 43)[1], 10 ** 15)
    tested = [c[0] for c in logi.calls]
    report('inverted control from m = 0 passes F_0..F_4 and reaches the kernel entry with M = F_5 = 4294967297',
           res[0] == 'STOP' and tested == [3, 5, 17, 257, 65537] and not logi.contract_bad,
           f'completed tests {tested}, then entry with r5 = {res[2][r[5]]}, r3 = {res[2][r[3]]}')


# ----------------------------------------------------------------------------- LEGENDREGAME
def legendre(L, nmax, nlit_windows):
    prog = Program(L, 'LEGENDREGAME')
    start = 47 * 31 * 37 * 41 ** 2
    kern = common_checks('LEGENDREGAME', prog, start, 33, 29)
    S = prog.sidx
    X = S[prog.fracs[8].numerator]
    r = {p: prog.ridx[p] for p in (2, 3, 5, 7)}
    claim('LEGENDREGAME: 8 registers, 21 states', len(prog.regs) == 8 and len(prog.states) == 21,
          f'registers {prog.regs}')
    st, R = prog.decode(start)
    report('start 47*31*37*41^2 = state 47 with C1 = 1, W = 1, L = 2 (n = 1, window {2,3})',
           st == 47 and R[prog.ridx[31]] == 1 and R[prog.ridx[37]] == 1 and R[prog.ridx[41]] == 2)
    pred = legendre_tests(max(nmax, 130))
    # literal integer run (lockstep) for the first windows
    log = KernelLog(prog, kern, X, r)
    target = len(legendre_tests(nlit_windows))
    t0 = time.time()

    def hook_l(prev, s, R, steps):
        log.on(prev, s, R, steps)
        return len(log.calls) >= target and prev in kern
    st_, s_, R_, steps, n_end, agree = lit_lockstep(prog, start, S[47], prog.decode(start)[1], 10 ** 9, hook_l)
    obs = [c[0] for c in log.calls]
    report(f'literal run (lockstep with vectors), windows n = 1..{nlit_windows}: tests = predicted, kernel contract '
           f'held', agree and obs[:target] == pred[:target] and not log.contract_bad,
           f'{len(obs)} calls, {steps} literal steps ({time.time() - t0:.0f}s)')
    # structure: NW (first rule moves W into C1; else C1 += 2, L += 2), X3 (r7 > 0 -> CL; else -> NW)
    C1, W, Lr = prog.ridx[31], prog.ridx[37], prog.ridx[41]
    NW = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((W, 1),)
          and prog.rules[s][0][1][C1] == 1 and len(prog.rules[s]) == 2 and prog.rules[s][1][1][Lr] == 2][0]
    X3 = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((r[7], 1),)
          and len(prog.rules[s]) == 2 and prog.rules[s][1][2] == NW][0]
    log = KernelLog(prog, kern, X, r)
    acc = []

    def hook(prev, s, R, steps):
        log.on(prev, s, R, steps)
        if prev == X3 and s == NW:
            acc.append((R[Lr] // 2, R[C1] + 1, steps))
            return len(acc) >= nmax
        return False
    t0 = time.time()
    res = Fast(prog, set(kern) | {NW}, hook, kern).run(S[47], prog.decode(start)[1], 10 ** 15)
    obs = [c[0] for c in log.calls]
    want = len(legendre_tests(nmax))
    least = all(p == min(c for c in range(n * n + 1, (n + 1) ** 2) if is_prime(c)) for n, p, _ in acc)
    report(f'fast run, windows n = 1..{nmax}: no halt, windows in order, least prime above n^2 accepted, tests = '
           f'predicted, kernel contract on every call', res[0] == 'STOP' and [a[0] for a in acc] ==
           list(range(1, nmax + 1)) and least and obs == pred[:want] and not log.contract_bad,
           f'{len(obs)} calls, {res[3]} steps, {time.time() - t0:.0f}s; last accepted {acc[-1][:2]}')
    if nmax == 120:
        claim('session L2 total 138,922,424,635 steps (at the 120th acceptance)', acc[-1][2] == 138922424635,
              str(acc[-1][2]))
    # ---- L3 as coded: NW's else-fraction replaced by C1 += 2, L += 2, W += 2 -> TEST
    C1, W, Lr = prog.ridx[31], prog.ridx[37], prog.ridx[41]
    NW = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((W, 1),)
          and prog.rules[s][0][1][C1] == 1 and len(prog.rules[s]) == 2 and prog.rules[s][1][1][Lr] == 2]
    NW = NW[0]
    Lc = list(L)
    Lc[prog.rules[NW][1][3]] = Fraction(47 * 31 ** 2 * 41 ** 2 * 37 ** 2, prog.states[NW])
    progc = Program(Lc, 'L3')
    progc.build(set(prog.states), set(prog.regs))
    logc = KernelLog(progc, kern, X, r)
    resc = Fast(progc, set(kern), logc.on, kern).run(S[47], prog.decode(start)[1], 10 ** 10)
    acc_c = [m for m, v, a, b in logc.calls if v]
    print(f'     L3 as coded: {resc[0]} in state {progc.states[resc[1]]}; accepted primes {acc_c}; L = {resc[2][Lr]} '
          f'at the halt; tests {[m for m, v, a, b in logc.calls]}')
    report('L3 as coded does NOT halt at n = 5 (accepts 2,5,11,13,17,23; halts with L = 14): the session\'s FAILED '
           'report is accurate', resc[0] == 'HALT' and acc_c == [2, 5, 11, 13, 17, 23] and resc[2][Lr] == 14)
    # ---- corrected control (own construction): budget register B (prime 127), 2 extra tests per window
    #   composite: CLR -> NEXTB: B>0 ? (B--, NEXTC) : H ;  prime: X3 -> ZB (clear B) -> NW ; RT's exit adds B += 2
    B = 127
    NEXTB, ZB, ZBh = 131, 137, 139
    Lb = []
    CL = None
    for i, f in enumerate(L):
        Lb.append(f)
    # locate the original states by structure
    X3 = [s for s in range(len(prog.states)) if prog.rules[s] and prog.rules[s][0][0] == ((r[7], 1),)
          and len(prog.rules[s]) == 2 and prog.rules[s][1][2] == NW][0]
    CLs = prog.rules[X3][0][2]
    NEXTC = prog.rules[CLs][-1][2]
    RT = [s for s in range(len(prog.states)) if len(prog.rules[s]) == 2 and prog.rules[s][1][0] == ((W, 1),)
          and prog.rules[s][1][2] == S[47]][0]
    # CL's exit now goes to NEXTB instead of NEXTC
    Lb[prog.rules[CLs][-1][3]] = Fraction(NEXTB, prog.states[CLs])
    # X3's prime branch goes to ZB instead of NW
    Lb[prog.rules[X3][1][3]] = Fraction(ZB, prog.states[X3])
    # RT's exit to TEST also sets B := 2 (B is 0 there)
    Lb[prog.rules[RT][1][3]] = Fraction(47 * B ** 2, prog.states[RT] * 37)
    Lb += [Fraction(prog.states[NEXTC], NEXTB * B), Fraction(prog.states[H_of(prog)], NEXTB),
           Fraction(ZBh, ZB * B), Fraction(prog.states[NW], ZB), Fraction(ZB, ZBh)]
    progb = Program(Lb, 'L3-corrected')
    stb, regb, und = progb.derive_states()
    progb.build(stb | set(prog.states), regb)
    logb = KernelLog(progb, {progb.sidx[p] for p in (11, 13, 17, 19, 23, 29)}, progb.sidx[prog.states[X]],
                     {p: progb.ridx[p] for p in (2, 3, 5, 7)})
    st0, R0 = progb.decode(start * B ** 2)
    resb = Fast(progb, set(logb.K), logb.on, logb.K).run(progb.sidx[st0], R0, 10 ** 10)
    acc_b = [m for m, v, a, b in logb.calls if v]
    tests_b = [m for m, v, a, b in logb.calls]
    report('corrected budget control (own construction, B = 2 extra tests per window): halts while processing n = 5 '
           'after accepting 2, 5, 11, 17 (tests 26, 27, 28)', resb[0] == 'HALT' and acc_b == [2, 5, 11, 17]
           and tests_b[-3:] == [26, 27, 28] and resb[2][progb.ridx[41]] == 10,
           f'{resb[0]}; tests {tests_b}; L = {resb[2][progb.ridx[41]]}')


def H_of(prog):
    return [s for s in range(len(prog.states)) if not prog.rules[s]][0]


def main():
    gold_nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    leg_nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    do_m4 = (sys.argv[3] != '0') if len(sys.argv) > 3 else True
    G, F, L, same, raw = programs_from_files()
    print('GOLDBACHGAME parsed from THM-4856; identical in the results note and goldbachgame.out:', same)
    report('the three published GOLDBACHGAME lists agree', same)
    goldbach(G, gold_nmax, 100)
    fermat(F, do_m4)
    legendre(L, leg_nmax, 12)
    print()
    print('SESSION CLAIMS FOUND WRONG: ' + ('; '.join(MISMATCH) if MISMATCH else 'none'))
    print('AUDIT CHECK FAILURES: ' + ('; '.join(FAILS) if FAILS else 'none'))
    return 1 if FAILS else 0


if __name__ == '__main__':
    sys.exit(main())
