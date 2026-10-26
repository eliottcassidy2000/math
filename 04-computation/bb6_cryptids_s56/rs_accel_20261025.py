#!/usr/bin/env python3
"""Generic certified accelerator for the DeepMind BB(6) Cryptid rule systems (S56).

Labelled rule steps: every label is one branch of the paper's simplified rule system (Appendix A of
arXiv:2610.22603v2), with the linear self-loops collapsed and split by which bound stops them, so that each label
is an affine map on (polyhedron) cap (lattice coset).  Three-point lemma (see a6_accel_20261025.py): if b, b+D,
b+K*D all follow a fixed label sequence L with images b+D, b+2D, b+(K+1)D, then so does every b+k*D, 0<=k<=K; if
K = 2^E with E above three times the bit size involved, the chain is infinite (certified translated cycler: the
machine can never halt from there, since HALT is not among the labels).
Level 1: cycles of labelled steps inside a segment (segments start at an anchor configuration; history is fresh
for every segment, so the tokenisation of a segment is a function of its start state).  Level 2: cycles of
segments.  States are tuples of nonnegative integers.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))


class Halt(Exception):
    pass


class Capture(Exception):
    pass


# ---------------------------------------------------------------- labelled steps
def ls_A2(p):
    a, b, w = p
    if a == 0:
        if b % 2 == 1:
            if w <= 1:
                return (0, 1, b + 1 + 2 * w), 'o_small'
            return (0, b + 3, w - 2), 'o_big'
        if w == 0:
            return 'HALT', 'e_halt'
        if w <= 3:
            return (b // 2, 1, 2 * w - 1), 'e_small'
        return (b // 2, 5, w - 4), 'e_big'
    if w >= 2:
        if a <= w // 2:
            return (0, b + 4 * a, w - 2 * a), 'loop_a'
        h = w // 2
        return (a - h, b + 4 * h, w - 2 * h), ('loop_w', w % 2)
    return (a - 1, 1, b + 2 + 2 * w), 'a_small'


def ls_A5(p):
    a, z, w = p
    if a == 0:
        if z % 2 == 1:
            if w >= 2:
                return (0, z + 3, w - 2), 'o_big'
            return (0, 2, z + 1 + 2 * w), 'o_small'
        if w >= 4:
            return (z // 2, 5, w - 4), 'e_big'
        if w == 0:
            return 'HALT', 'e_halt'
        return (z // 2 - 1, 3, 2 * w - 1), 'e_small'
    if w >= 2:
        if a <= w // 2:
            return (0, z + 4 * a, w - 2 * a), 'loop_a'
        h = w // 2
        return (a - h, z + 4 * h, w - 2 * h), ('loop_w', w % 2)
    if a == 1:
        return (0, 2, z + 2 + 2 * w), 'a1_small'
    return (a - 2, 3, z + 2 + 2 * w), 'a_small'


def ls_A10(p):
    a, b, c = p
    if c == 0:
        if a >= 1:
            return (a - 1, b, 1), 'c0_a'
        if b <= 2:
            return ((1, 0, 1), (1, 0, 2), (1, 0, 6))[b], ('c0_b', b)
        return (1, b - 3, 4), 'c0_b3'
    if b == 0:
        if c <= 2:
            return (1, 0, 2 * a + 2 * c + 1), 'b0_small'
        if c == 3:
            return 'HALT', 'b0_halt'
        if c <= 6:
            return (1, a + 1, 2 * c - 6), 'b0_mid'
        return (1, a + 3, c - 6), 'b0_big'
    if c >= 3:
        if 2 * b + 1 <= c:
            return (a + 2 * b, 0, c - 2 * b), 'loop_b'
        h = (c - 1) // 2
        return (a + 2 * h, b - h, c - 2 * h), ('loop_c', c % 2)
    if b == 1:
        return (1, 0, 2 * a + 2 * c + 2), 'b1'
    if b == 2:
        return (1, 2, 2 * a + 2 * c), 'b2'
    return (1, b - 3, 2 * a + 2 * c + 4), 'b3'


def ls_A1(p):
    a, b, c, d = p
    if c > 0:
        return (a, b + c, 0, d), 'c'
    if d == 0:
        return (a + 1, 0, 0, 2 * b + 2), 'd0'
    if a >= 3:
        if 3 * d <= a:
            return (a - 3 * d, b + 2 * d, 0, 0), 'loop_d'
        q, r = divmod(a, 3)
        return (r, b + 2 * q, 0, d - q), ('loop_a', r)
    if a == 2:
        return (0, 0, b + 1, d - 1), 'a2'
    if a == 1:
        if b > 0:
            return (1, 1, b + 1, d - 1), 'a1b'
        return (0, 0, 0, d + 2), 'a1b0'
    if d == 1:
        if b == 0:
            return (0, 0, 2, 5), 'd1b0'
        return (0, 0, 2, 2 * b + 1), 'd1'
    if b >= d - 1:
        return (3 * d - 4, 1, 2, 2 * b - 2 * d + 3), 'bge'
    if d == b + 2:
        return 'HALT', 'halt'
    if d == b + 3:
        return (3 * b + 2, 0, 0, 6), 'db3'
    if d == b + 4:
        return (3 * b + 4, 0, 0, 4), 'db4'
    return (3 * b + 3, 2, 0, d - b - 5), 'db5'


LSTEP = {'A1': ls_A1, 'A2': ls_A2, 'A5': ls_A5, 'A10': ls_A10}
ANCHOR = {'A2': lambda p: p[1] == 1, 'A5': lambda p: p[1] <= 3, 'A10': lambda p: p[0] == 1,
          'A1': lambda p: p[2] == 0 and p[3] == 0}
START = {'A1': (0, 0, 0, 0), 'A2': (0, 1, 1), 'A5': (0, 2, 2), 'A10': (1, 0, 2)}


def vadd(x, D, k=1):
    return tuple(xi + k * di for xi, di in zip(x, D))


def vsub(x, y):
    return tuple(xi - yi for xi, yi in zip(x, y))


def maxbits(states):
    return max(abs(v).bit_length() for s in states for v in s)


def max_valid(valid, E):
    if valid(1 << E):
        return None
    lo, hi = 1, 2
    while valid(hi):
        lo, hi = hi, 2 * hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if valid(mid):
            lo = mid
        else:
            hi = mid
    return lo


class Engine:
    def __init__(self, name, gram=3, pmax=4000, seg_cap=5 * 10 ** 6):
        self.name = name
        self.ls = LSTEP[name]
        self.anchor = ANCHOR[name]
        self.gram, self.pmax, self.seg_cap = gram, pmax, seg_cap
        self.stats = {'steps': 0, 'l1_jumps': 0, 'l2_jumps': 0}

    def run_labels(self, s, labels):
        for L in labels:
            if min(s) < 0:
                return None
            q, lab = self.ls(s)
            if lab != L or q == 'HALT':
                return None
            s = q
        return s

    def segment(self, s0):
        """From state s0, run until the next anchor state (after >= 1 step). Returns (state, tokens)."""
        xs = [s0]
        labs = []
        toks = []
        occ = {}
        n_steps = 0
        while True:
            q, lab = self.ls(xs[-1])
            n_steps += 1
            self.stats['steps'] += 1
            if q == 'HALT':
                toks.append(lab)
                raise Halt((xs[-1], tuple(toks[-50:])))
            labs.append(lab)
            toks.append(lab)
            xs.append(q)
            if self.anchor(q):
                return q, tuple(toks)
            if n_steps > self.seg_cap:
                raise RuntimeError(('segment cap', q))
            n = len(labs)
            if n < 2 * self.gram:
                continue
            key = tuple(labs[-self.gram:])
            lst = occ.setdefault(key, [])
            jumped = False
            for i in reversed(lst):
                p = n - i
                if p > self.pmax or 2 * p > n:
                    break
                D = vsub(xs[n], xs[n - p])
                if D != vsub(xs[n - p], xs[n - 2 * p]) or not any(D):
                    continue
                if labs[n - p:n] != labs[n - 2 * p:n - p]:
                    continue
                if any(isinstance(l, tuple) and l and l[0] == 'cyc' for l in labs[n - 2 * p:n]):
                    continue
                cyc = tuple(labs[n - p:n])
                b = xs[n - 2 * p]
                E = max(1000, 3 * maxbits(xs[n - 2 * p:n + 1]) + 300)

                def valid(K):
                    return self.run_labels(vadd(b, D, K), cyc) == vadd(b, D, K + 1)
                Kmax = max_valid(valid, E)
                if Kmax is None:
                    raise Capture({'base': b, 'D': D, 'labels': cyc, 'E': E})
                # the chain from b may pass an anchor state: stop at the first anchor inside the jump range
                m = Kmax + 1
                new = vadd(b, D, m)
                self.stats['l1_jumps'] += 1
                del xs[n - 2 * p + 1:]
                del labs[n - 2 * p:]
                del toks[len(toks) - 2 * p:]
                toks.append(('cyc', cyc, m))
                labs.append(('cyc', cyc, m))
                xs.append(new)
                occ = {}
                jumped = True
                break
            if not jumped:
                lst.append(n)
                if len(lst) > 8:
                    del lst[0]
            elif self.anchor(xs[-1]):
                return xs[-1], tuple(toks)

    def run_segments(self, s, toks_list):
        for T in toks_list:
            if min(s) < 0:
                return None
            try:
                s2, t = self.segment(s)
            except (Halt, Capture, RuntimeError):
                return None
            if t != T:
                return None
            s = s2
        return s

    def orbit(self, s0, max_segments=10 ** 18, stop=None, gram=2, pmax=3000, on_seg=None):
        Ss = [s0]
        toks = []
        occ = {}
        nseg = 0
        while nseg < max_segments:
            s = Ss[-1]
            if stop is not None and stop(s, nseg):
                return 'stopped', s, nseg, None
            try:
                s2, t = self.segment(s)
            except Halt as h:
                return 'HALT', s, nseg, h.args[0]
            except Capture as c:
                return 'CAPTURE', s, nseg, c.args[0]
            nseg += 1
            toks.append(t)
            Ss.append(s2)
            if on_seg is not None:
                on_seg(s, s2, t, nseg)
            n = len(toks)
            if len(Ss) > 20000:
                del Ss[:10000]
                del toks[:10000]
                occ = {}
                continue
            if n < 2 * gram:
                continue
            key = tuple(toks[-gram:])
            lst = occ.setdefault(key, [])
            found = False
            for i in reversed(lst):
                p = n - i
                if p > pmax or 2 * p > n:
                    break
                D = vsub(Ss[n], Ss[n - p])
                if D != vsub(Ss[n - p], Ss[n - 2 * p]) or not any(D):
                    continue
                if toks[n - p:n] != toks[n - 2 * p:n - p]:
                    continue
                cyc = tuple(toks[n - p:n])
                b = Ss[n - 2 * p]
                E = max(1000, 3 * maxbits(Ss[n - 2 * p:n + 1]) + 300)

                def valid(K):
                    return self.run_segments(vadd(b, D, K), cyc) == vadd(b, D, K + 1)
                Kmax = max_valid(valid, E)
                if Kmax is None:
                    return 'INFINITE-SEGMENT-CYCLE', b, nseg, {'base': b, 'D': D, 'p': p, 'E': E}
                m = min(Kmax + 1, (max_segments - nseg + 2 * p) // p)
                new = vadd(b, D, m)
                self.stats['l2_jumps'] += 1
                nseg += m * p - 2 * p
                del Ss[n - 2 * p + 1:]
                del toks[n - 2 * p:]
                toks.append(('cyc2', len(cyc), m))
                Ss.append(new)
                occ = {}
                found = True
                break
            if not found:
                lst.append(n)
                if len(lst) > 8:
                    del lst[0]
        return 'max', Ss[-1], nseg, None


if __name__ == '__main__':
    # self-test: labelled steps reproduce deepmind_cryptids.step (jump=True) on the first 200000 macro steps
    sys.path.insert(0, HERE)
    from deepmind_cryptids_20261025 import step, START as DSTART, HALT
    for nm in ('A1', 'A2', 'A5', 'A10'):
        p = DSTART[nm]
        ok = True
        for i in range(200000):
            q1 = step(nm, p)
            q2, lab = LSTEP[nm](p)
            if q2 == 'HALT':
                q2 = HALT
            if q1 != q2:
                ok = False
                print(nm, 'MISMATCH at', p, q1, q2, lab)
                break
            if q1 == HALT:
                break
            p = q1
        print(nm, 'labelled steps agree on 200000 macro steps:', ok, flush=True)
