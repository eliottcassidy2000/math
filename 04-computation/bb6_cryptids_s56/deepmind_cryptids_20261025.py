#!/usr/bin/env python3
"""The eleven new BB(6) Cryptid accelerations of Dos Santos et al. (arXiv:2610.22603v2, Appendix A), S56.

For each machine: the Turing machine (standard text format), the configuration family of the paper's simplified
system, and an implementation of that rule system (Python big integers, with the paper's self-loop rules jumped in
one step where they only shift parameters linearly).

  validate  : run the TM directly and check that the rule system's predicted configurations occur on the TM's tape,
              in order (catches transcription errors; the rules themselves are Lean-verified in the paper).
  run       : iterate the rule system for a budget of macro steps; record halting, sizes, and the halting-margin
              statistic of each machine.
Usage: python deepmind_cryptids_20261025.py validate [maxTMsteps] | run [macro_steps]
Session opus-2026-10-25-S56.
"""
import sys, json, time, os
HERE = os.path.dirname(os.path.abspath(__file__))

MACHINES = {
    'A1': '1RB1LD_1RC0RA_1LA1RE_0LA0LB_0RB1RF_---1RD',
    'A2': '1RB1RD_0LC1RA_0RA1LC_1RA1LE_1RF0LD_---0RE',
    'A3': '1RB1RB_1LC0RE_---1LD_0RA1LB_1RA1LF_0LE1LE',
    'A4': '1RB0LF_1LC1RA_0RE0RD_---1LE_1LF1RC_1LC1LA',
    'A5': '1RB1RD_0LC1RA_1RA1LC_0RB1LE_1RF0LD_---0RE',
    'A6': '1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA',
    'A7': '1RB0LC_0LC1RF_1LE1RD_0RC1RC_0LA1LA_---1RA',
    'A8': '1RB1LF_0RC1RD_1LA1LD_1LE0RA_---1LC_0LA1LA',
    'A9': '1RB1RC_1LC1LC_1RE0LD_1LB1RF_---1RA_0RD1RD',
    'A10': '1RB1LE_1RC1RA_0LD1RB_1LE1LD_1RF0LA_---0RE',
    'A11': '1RB1LA_1LC0RE_1LF1LD_1RE0LC_1LB1RA_---0LA',
}

HALT = 'HALT'


# ---------------------------------------------------------------- configuration families: (state, left, right)
# ">" patterns: head reads the first symbol of `right`.  ("<" for A3 is converted to ">" form.)
def fam(name, p):
    s = lambda x, n: x * n
    if name == 'A1':
        a, b, c, d = p
        return 'A', s('1', a) + s('10', b), s('01', c) + s('1', d)
    if name == 'A2':
        a, b, w = p
        return 'D', s('10', a) + s('1', b), s('1', w + 1)
    if name == 'A3':
        a, c, f = p
        left = s('1', a) + '0' + s('1', c) + '0' + '1'
        return 'D', left[:-1], left[-1] + s('1', f)
    if name == 'A4':
        a, b = p
        return 'E', '', s('10', a) + s('1', b)
    if name == 'A5':
        a, z, w = p
        return 'D', s('10', a) + s('1', z), '1' + s('1', w)
    if name == 'A6':
        u, v, w = p
        return 'A', s('1', 5) + s('10', u) + '1', '11' + s('01', v) + s('1', w) + '01'
    if name == 'A7':
        a, k, c = p
        return 'A', s('1', a), '0' + s('11', k) + '0' + s('1', c)
    if name == 'A8':
        a, b = p
        return 'A', '110', s('11', a) + '0' + s('1', b)
    if name == 'A9':
        a, b, c = p
        return 'E', s('1', a + 2), s('1', 2 * b + 2) + '0' + s('1', 2 * c + 2)
    if name == 'A10':
        a, b, c = p
        return 'E', s('10', b) + s('1', 2 * a + 2), s('1', c + 1)
    if name == 'A11':
        a, b, c, d = p
        return 'C', s('1', a + 1), '00' + s('1', b + 4) + s('01', c) + s('1', d) + '011'
    raise KeyError(name)


START = {'A1': (0, 0, 0, 0), 'A2': (0, 1, 1), 'A3': (6, 2, 2), 'A4': (0, 3), 'A5': (0, 2, 2), 'A6': (0, 0, 0),
         'A7': (3, 0, 2), 'A8': (0, 9), 'A9': (5, 2, 0), 'A10': (1, 0, 2), 'A11': (4, 0, 0, 3)}


# ---------------------------------------------------------------- rule systems (paper's simplified versions)
def step(name, p, jump=True):
    """One rule application (jump=True also collapses the linear self-loop rules)."""
    if name == 'A1':
        a, b, c, d = p
        if c > 0:
            return (a, b + c, 0, d) if jump else (a, b + 1, c - 1, d)
        if d == 0:
            return (a + 1, 0, 0, 2 * b + 2)
        if a >= 3:
            t = min(a // 3, d) if jump else 1
            return (a - 3 * t, b + 2 * t, 0, d - t)
        if a == 2:
            return (0, 0, b + 1, d - 1)
        if a == 1:
            return (1, 1, b + 1, d - 1) if b > 0 else (0, 0, 0, d + 2)
        if d == 1:
            return (0, 0, 2, 5) if b == 0 else (0, 0, 2, 2 * b + 1)
        if b >= d - 1:
            return (3 * d - 4, 1, 2, 2 * b - 2 * d + 3)
        if d == b + 2:
            return HALT
        if d == b + 3:
            return (3 * b + 2, 0, 0, 6)
        if d == b + 4:
            return (3 * b + 4, 0, 0, 4)
        return (3 * b + 3, 2, 0, d - b - 5)
    if name == 'A2':
        a, b, w = p
        if a == 0:
            if b % 2 == 1:
                k = (b - 1) // 2
                if w <= 1:
                    return (0, 1, 2 * k + 2 * w + 2)
                return (0, 2 * k + 4, w - 2)
            k = (b - 2) // 2
            if w == 0:
                return HALT
            if w <= 3:
                return (k + 1, 1, 2 * w - 1)
            return (k + 1, 5, w - 4)
        assert b % 2 == 1, ('A2 invariant', p)
        if w >= 2:
            t = min(a, w // 2) if jump else 1
            return (a - t, b + 4 * t, w - 2 * t)
        k = (b - 1) // 2
        return (a - 1, 1, 2 * k + 2 * w + 3)
    if name == 'A3':
        a, c, f = p
        K, M = a // 2, c // 2
        if a == 2 and c in (2, 4):
            return HALT
        if K >= 2 and M >= 2:
            if K < M:
                return (2, 2 * M - 2 * K + 2, f + 6 * K - 6)
            if M < K:
                return (2 * K - 2 * M + 2, 2, f + 6 * M - 6)
            return (2, 2, f + 6 * K - 6)
        if M == 1:
            if K >= 3:
                return (2 * K - 4, f + 8, 2) if f % 2 == 0 else (2 * K - 4, f - 1 + 8, 5)
            if K == 2:
                return (f + 4, 6, 5) if f % 2 == 0 else (f - 1 + 8, 2, 2)
            return HALT
        # K == 1, M >= 3
        if M >= 4:
            return (2 * M - 6, f + 10, 5) if f % 2 == 0 else (2 * M - 6, f - 1 + 12, 2)
        return (f + 10, 2, 2) if f % 2 == 0 else (f - 1 + 8, 6, 5)
    if name == 'A4':
        a, b = p
        if b > 2 * a + 2:
            if not jump:
                return (2 * a + 2, b - 2 * a - 2)
            # cascade a -> 2a+2 while b > 2a+2 (closed form): a_i + 2 = (a+2) 2^i, cost 2a_i + 2 each step
            while b > 2 * a + 2:
                b -= 2 * a + 2
                a = 2 * a + 2
            return (a, b)
        if b == 1:
            return (2, 2 * a + 2)
        if b == 2:
            return HALT
        if b % 2 == 0:
            k = b // 2
            x = 6 * a + 2 * k + 2
        else:
            k = (b - 1) // 2
            x = 6 * a + 2 * k + 6
        return (x // 4, x % 4 + 1)
    if name == 'A5':
        a, z, w = p
        if a == 0:
            if z % 2 == 1:
                k = (z - 1) // 2
                if w >= 2:
                    return (0, 2 * k + 4, w - 2)
                return (0, 2, 2 * k + 2) if w == 0 else (0, 2, 2 * k + 4)
            k = (z - 2) // 2
            if w >= 4:
                return (k + 1, 5, w - 4)
            if w == 0:
                return HALT
            return (k, 3, 2 * w - 1)
        assert z % 2 == 1, ('A5 invariant', p)
        k = (z - 1) // 2
        if w >= 2:
            t = min(a, w // 2) if jump else 1
            return (a - t, z + 4 * t, w - 2 * t)
        if a == 1:
            return (0, 2, 2 * k + 3) if w == 0 else (0, 2, 2 * k + 5)
        return (a - 2, 3, 2 * k + 3) if w == 0 else (a - 2, 3, 2 * k + 5)
    if name == 'A6':
        u, v, w = p
        if u == 0:
            if w == 0:
                return (0, 0, 2 * v + 12) if v % 2 == 0 else (2, 4, 2 * v + 7)
            if w == 1:
                return (2, 4, 2 * v + 9) if v % 2 == 0 else (0, 0, 2 * v + 14)
            if w == 2:
                return (2, v + 4, 1) if v % 2 == 0 else (0, v + 4, 1)
            return (v + 3, 0, w - 3)
        if w >= 2:
            t = min(u, w // 2) if jump else 1
            return (u - t, v + 2 * t, w - 2 * t)
        if w == 0:
            if u == 1:
                return (0, 0, 2 * v + 12)
            return (2, 4, 2 * u + 2 * v + 5) if u % 2 == 0 else (0, 0, 2 * u + 2 * v + 10)
        if u == 1:
            return HALT
        if u == 2:
            return (0, 0, 2 * v + 10)
        return (2, u, 2 * v + 5) if u % 2 == 0 else (0, u, 2 * v + 5)
    if name == 'A7':
        a, k, c = p
        if k == 0:
            if c == 0:
                return (3, 0, a + 2)
            m = a // 2
            if c % 2 == 1:
                return (3, m, c) if a % 2 == 0 else (4, m + 1, c)
            if c == 2:
                return (0, 1, 2 * m + 2) if a % 2 == 0 else (1, 1, 2 * m + 4)
            return (3, m + 1, c - 2) if a % 2 == 0 else (4, m + 2, c - 2)
        if k == 1:
            if c % 2 == 0:
                return (a + 3, c // 2, 0)
            if c == 1:
                return HALT
            if c == 3:
                return (3, 0, a + 7)
            m = a // 2
            if c == 5:
                return (1, 1, 2 * m + 8) if a % 2 == 0 else (0, 1, 2 * m + 8)
            return (4, m + 4, c - 5) if a % 2 == 0 else (3, m + 4, c - 5)
        K = k - 2
        if c == 0:
            m = a // 2
            return (3, m + 1, 2 * K + 1) if a % 2 == 0 else (4, m + 2, 2 * K + 1)
        if c % 2 == 1:
            if jump:
                t = k // 2            # E(a,k+2,odd) -> E(a+6,k,odd), repeat while k >= 2
                return (a + 6 * t, k - 2 * t, c)
            return (a + 6, K, c)
        if c == 2:
            return (a + 3, 1, 2 * K + 2)
        if jump:                      # E(a,k+2,2j+4) -> E(a+6,k+1,2j+2): k -> k-1, c -> c-2 while k >= 2, c >= 4
            t = min(k - 1, (c - 2) // 2)
            return (a + 6 * t, k - t, c - 2 * t)
        return (a + 6, K + 1, c - 2)
    if name == 'A8':
        a, b = p
        if b % 2 == 0:
            k = (b - 6) // 2
            if a == 0:
                return (k + 1, 9)
            ap = a - 1
            if k <= ap + 3:
                return (ap - k + 4, 6 * k + 12)
            if k in (ap + 4, ap + 5):
                return HALT
            return (k - ap - 5, 6 * ap + 39)
        k = (b - 1) // 2
        if k > a:
            return (k - a, 6 * a + 6)
        if a in (k, k + 1):
            return HALT
        return (a - k - 1, 6 * k + 9)
    if name == 'A9':
        a, b, c = p
        if b == 0:
            if c == 0:
                return HALT
            if c == 1:
                return (5, 2, a // 2 + 1) if a % 2 == 0 else (2, 0, (a - 1) // 2 + 3)
            return (2, a // 2 + 3, c - 2) if a % 2 == 0 else (5, (a - 1) // 2 + 3, c - 2)
        if c == 0:
            return (a + 3, 0, b - 1)
        t = min(b, c) if jump else 1
        return (a + 6 * t, b - t, c - t)
    if name == 'A10':
        a, b, c = p
        if c == 0:
            if a >= 1:
                return (a - 1, b, 1)
            if b == 0:
                return (1, 0, 1)
            if b == 1:
                return (1, 0, 2)
            if b == 2:
                return (1, 0, 6)
            return (1, b - 3, 4)
        if b == 0:
            if c <= 2:
                return (1, 0, 2 * a + 2 * c + 1)
            if c == 3:
                return HALT
            if c <= 6:
                return (1, a + 1, 2 * (c - 4) + 2)
            return (1, a + 3, c - 6)
        if c >= 3:
            t = min(b, (c - 1) // 2) if jump else 1
            return (a + 2 * t, b - t, c - 2 * t)
        if b == 1:
            return (1, 0, 2 * a + 2 * c + 2)
        if b == 2:
            return (1, 2, 2 * a + 2 * c)
        return (1, b - 3, 2 * a + 2 * c + 4)
    if name == 'A11':
        a, b, c, d = p
        if c == 0 and d == 6:
            b, d = b + 3, 3
        if a >= 2:
            t = a // 2 if jump else 1
            return (a - 2 * t, b + 3 * t, c, d)
        if c >= 2:
            return (a + b + 4, 3, c - 2, d)
        if c == 1 and d == 6:
            return (a + b + 4, 0, 3, 3)
        if c == 0:
            if b % 2 == 0:
                return (a, 0, b // 2 + 3, 3)
            return HALT if a == 0 else (6, 3, (b - 1) // 2, 3)
        # c == 1, d == 3
        if b % 2 == 1:
            return (a, 0, (b - 1) // 2 + 3, 6)
        if a == 0:
            return HALT
        return (6, 0, 3, 3) if b == 0 else (6, 3, (b - 2) // 2, 6)
    raise KeyError(name)


# ---------------------------------------------------------------- Turing machine
def parse_tm(code):
    rows = code.split('_')
    tab = {}
    for si, row in enumerate(rows):
        st = 'ABCDEF'[si]
        for sym in (0, 1):
            t = row[3 * sym:3 * sym + 3]
            if t == '---':
                tab[(st, sym)] = None
            else:
                tab[(st, sym)] = (int(t[0]), 1 if t[1] == 'R' else -1, t[2])
    return tab


def validate(name, max_steps=3_000_000, max_preds=60):
    tab = parse_tm(MACHINES[name])
    size = 1 << 21
    tape = bytearray(size)
    head = size // 2
    lo = hi = head
    state = 'A'
    preds = []
    p = START[name]
    for _ in range(max_preds):
        if p == HALT:
            preds.append(HALT)
            break
        preds.append(p)
        p = step(name, p, jump=False)
    targets = []
    for q in preds:
        if q == HALT:
            targets.append(HALT)
        else:
            st, L, R = fam(name, q)
            targets.append((st, L.encode().replace(b'1', b'\x01').replace(b'0', b'\x00'),
                            R.encode().replace(b'1', b'\x01').replace(b'0', b'\x00'), q))
    ti = 0
    matched = []
    for t in range(max_steps + 1):
        if ti < len(targets):
            tg = targets[ti]
            if tg != HALT and state == tg[0]:
                L, R = tg[1], tg[2]
                if head - len(L) >= 0 and tape[head - len(L):head] == L and tape[head:head + len(R)] == R:
                    if tape.find(b'\x01', lo, head - len(L)) == -1 and tape.find(b'\x01', head + len(R), hi + 1) == -1:
                        matched.append((t, tg[3]))
                        ti += 1
        sym = tape[head]
        tr = tab[(state, sym)]
        if tr is None:
            if ti < len(targets) and targets[ti] == HALT:
                matched.append((t, HALT))
                ti += 1
            break
        w, mv, ns = tr
        tape[head] = w
        head += mv
        if head < lo:
            lo = head
        if head > hi:
            hi = head
        state = ns
        if head <= 2 or head >= size - 3:
            break
    return {'machine': name, 'predicted': len(targets), 'matched_in_order': ti, 'tm_steps_used': t,
            'last_match': matched[-1] if matched else None}


def run(name, budget=200000, size_cap_bits=200000):
    p = START[name]
    n = 0
    maxbits = 0
    t0 = time.time()
    while n < budget:
        if p == HALT:
            return {'machine': name, 'halted': True, 'macro_steps': n}
        p = step(name, p, jump=True)
        n += 1
        if p != HALT:
            mb = max(x.bit_length() for x in p)
            maxbits = max(maxbits, mb)
            if mb > size_cap_bits:
                break
        if time.time() - t0 > 300:
            break
    return {'machine': name, 'halted': p == HALT, 'macro_steps': n, 'max_bits': maxbits,
            'final_bits': [x.bit_length() for x in p] if p != HALT else None}


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'validate'
    if mode == 'validate':
        ms = int(sys.argv[2]) if len(sys.argv) > 2 else 3_000_000
        res = []
        for nm in MACHINES:
            r = validate(nm, ms)
            res.append(r)
            print(f"{nm}: predicted {r['predicted']}, matched in order {r['matched_in_order']} "
                  f"(TM steps used {r['tm_steps_used']}; last match {r['last_match']})", flush=True)
        json.dump(res, open(os.path.join(HERE, 'validate.json'), 'w'), indent=1, default=str)
    elif mode == 'run':
        bud = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
        res = []
        for nm in MACHINES:
            r = run(nm, bud)
            res.append(r)
            print(r, flush=True)
        json.dump(res, open(os.path.join(HERE, 'run.json'), 'w'), indent=1, default=str)
