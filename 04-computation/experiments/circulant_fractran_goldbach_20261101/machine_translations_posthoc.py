#!/usr/bin/env python3
"""mac-mini-2026-11-01-lenses: POST-HOC (not pre-registered) anatomy of the failed checks T2b and T3 of
machine_translations.py. The verdicts there stand (FAILED as stated).

Anatomy: Compiled.run recorded an "entry into a TM-state prime P_q" every time control reached P_q, including the
returns P_q <- Ph_q of the head-bit extraction loop (made whenever that loop runs a single iteration and the
accelerator does not jump it). Those extra records shifted the step alignment. Here an entry is recorded only when
P_q is reached from a non-helper state (an instruction exit or the start), which is the intended "k-th TM step".
The FRACTRAN program itself is unchanged.

usage: python3 machine_translations_posthoc.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mt', os.path.join(HERE, 'machine_translations.py'))
mt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mt)


def run_entries(C, start, nentries):
    r = {'L': 0, 'R': 0, 'A': 0, 'C': 0}
    s, steps = start, 0
    entries = [(s, 0, 0, 0)]
    Pstates = set(C.entry.values())
    helpers = {f'Ph{q}' for q in C.entry}
    while len(entries) <= nentries:
        got = C.pick(s, r)
        if got is None:
            return entries, steps, True
        idx, g, inc, nxt = got
        if nxt in C.by:
            mid = {k: inc.get(k, 0) - g.get(k, 0) for k in r}
            r1 = {k: r[k] + mid[k] for k in r}
            got2 = C.pick(nxt, r1)
            if got2 is not None and got2[3] == s:
                idx2, g2, inc2, _ = got2
                d = {k: mid[k] + inc2.get(k, 0) - g2.get(k, 0) for k in r}
                ok = all(all(d[k] <= 0 for k in C.by[s][j][0]) for j in range(idx)) and \
                    all(all(d[k] <= 0 for k in C.by[nxt][j][0]) for j in range(idx2))
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
        prev, s = s, nxt
        steps += 1
        if s in Pstates and prev not in helpers:
            entries.append((s, r['L'], r['R'], steps))
    return entries, steps, False


def main():
    ok_all = True
    tab, states, q0 = mt.leng_table()
    C = mt.Compiled(tab, states)
    N = 3000
    entries, steps, halted = run_entries(C, C.entry[q0], N)
    tm = mt.TM(tab, q0)
    agree = True
    for k, (s, L, R, st) in enumerate(entries):
        if k > 0:
            tm.step()
        if s != C.entry[tm.q] or (L, R) != tm.stacks():
            agree = False
            print(f'  mismatch at TM step {k}: {s, L, R} vs {tm.q, tm.stacks()}')
            break
    print(f'T2b (post-hoc): stacks agree at each of the first {N} TM steps: {agree and len(entries) == N + 1}; '
          f'FRACTRAN steps {steps}')
    ok_all &= agree and len(entries) == N + 1
    marks = [(k, entries[k][3]) for k in (10, 100, 1000, 2000, 3000)]
    print(f'  cumulative FRACTRAN steps at TM steps 10/100/1000/2000/3000: {marks}')
    print(f'  largest register value at TM step 3000: {max(entries[-1][1], entries[-1][2])} '
          f'(= 2^{max(entries[-1][1], entries[-1][2]).bit_length() - 1}..)')
    btab, bstates, bq0 = mt.bb2_table()
    Cb = mt.Compiled(btab, bstates)
    eb, sb, hb = run_entries(Cb, Cb.entry[bq0], 100)
    tm = mt.TM(btab, bq0)
    seq = []
    for k, (s, L, R, st) in enumerate(eb):
        if k > 0:
            tm.step()
        seq.append((s, L, R) == (Cb.entry[tm.q], *tm.stacks()))
    L, R = eb[-1][1], eb[-1][2]
    ones = bin(L).count('1') + bin(R).count('1')
    ok_b = hb and len(eb) == 7 and all(seq) and ones == 4
    print(f'T3 (post-hoc): BB(2) compiled FRACTRAN: {len(eb) - 1} TM steps recorded, all agree {all(seq)}, '
          f'halted {hb} after {sb} FRACTRAN steps, final stacks L={L} R={R} ({ones} ones): {ok_b}')
    ok_all &= ok_b
    print()
    print('POST-HOC CHECKS PASSED' if ok_all else 'POST-HOC CHECK FAILED')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
