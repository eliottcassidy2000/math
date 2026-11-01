#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: POST-HOC (not pre-registered) anatomy of the failed controls G4 and L3.
Nothing here changes the pre-registered verdicts (G4 and L3 FAILED as stated).

G4 (goldbachgame.py): the control fraction H/(NEXT * 31) tests P' > 0 by consuming one unit of P' (FRACTRAN can only
    test by dividing). The check read the target as P' + Q' + 4 = 11 at the halt. Here the halt configuration is
    decoded: P' = 0 after the consumed unit, Q' = 7, so the target being processed was (P' + 1) + Q' + 4 = 12 and
    the last pair tested was (3, 9).
L3 (primekernel_programs.py): the control set W = 2 at each window, but W also drives the jump to the next window
    (C1 += W + 2), so windows drifted (it accepted 13 for n = 4). A correctly specified control keeps the full W
    and adds a separate budget register B (two further tests after the first) per window:
      NEXTC: B > 0 ? (B--, NEXTC2) : H;   NEXTC2: W > 0 ? (W--, C1++, TEST) : H;
      at a prime: clear B, then the usual next-window code, whose last fraction also sets B := 2.
    Prediction for the corrected control: halts while processing n = 5 (26, 27, 28 composite), after accepting
    the windows n = 1..4 with primes 2, 5, 11, 17.

usage: python3 posthoc_controls.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    gg = load('goldbachgame')
    Ic = gg.program(control_variant=True)
    M = gg.Machine(Ic)
    ev, st, _ = M.run(1000)
    # rerun to the halt and read the registers there
    r = {k: 0 for k in M.regs}
    s, steps = 'TP', 0
    while s != 'H':
        idx, g, inc, nxt = M.pick(s, r)
        for k, e in g.items():
            r[k] -= e
        for k, e in inc.items():
            r[k] += e
        s = nxt
        steps += 1
    print('G4 anatomy:')
    print(f'  accepted targets {[e[0] for e in ev]}; halt after {steps} steps with P\' = {r["Pp"]}, Q\' = {r["Qp"]}')
    print(f'  target being processed = (P\'+1) + Q\' + 4 = {r["Pp"] + 1 + r["Qp"] + 4}; last pair tested = '
          f'({r["Pp"] + 1 + 2}, {r["Qp"] + 2})')
    ok_g4 = [e[0] for e in ev] == [4, 6, 8, 10] and r['Pp'] + 1 + r['Qp'] + 4 == 12
    print(f'  consistent with halting while processing N = 12: {ok_g4}')
    print()
    pk = load('primekernel_programs')
    regs = ['C1', 'W', 'L', 'T', 'B']
    states = ['TEST', 'TESTh', 'X', 'Xh', 'X3', 'CL', 'CLh', 'NEXTC', 'NEXTC2', 'CB', 'CBh', 'NW', 'NWh', 'CP', 'CPh',
              'RT', 'RTh', 'H']
    c = [('TEST', {'C1': 1}, {'r5': 1, 'r3': 1}, 'TESTh'), ('TEST', {}, {'r5': 1}, 'S11'), ('TESTh', {}, {}, 'TEST'),
         ('X', {'r2': 1}, {'C1': 1}, 'Xh'), ('X', {'C1': 1}, {}, 'X3'), ('Xh', {}, {}, 'X'),
         ('X3', {'r7': 1}, {}, 'CL'), ('X3', {}, {}, 'CB'),
         ('CL', {'r7': 1}, {}, 'CLh'), ('CL', {}, {}, 'NEXTC'), ('CLh', {}, {}, 'CL'),
         ('NEXTC', {'B': 1}, {}, 'NEXTC2'), ('NEXTC', {}, {}, 'H'),
         ('NEXTC2', {'W': 1}, {'C1': 1}, 'TEST'), ('NEXTC2', {}, {}, 'H'),
         ('CB', {'B': 1}, {}, 'CBh'), ('CB', {}, {}, 'NW'), ('CBh', {}, {}, 'CB'),
         ('NW', {'W': 1}, {'C1': 1}, 'NWh'), ('NW', {}, {'C1': 2, 'L': 2}, 'CP'), ('NWh', {}, {}, 'NW'),
         ('CP', {'L': 1}, {'W': 1, 'T': 1}, 'CPh'), ('CP', {}, {}, 'RT'), ('CPh', {}, {}, 'CP'),
         ('RT', {'T': 1}, {'L': 1}, 'RTh'), ('RT', {'W': 1}, {'B': 2}, 'TEST'), ('RTh', {}, {}, 'RT')]
    P = pk.Prog(regs, states, c)
    accepted = []

    def hook(prev, s, r, steps):
        if prev == 'X3' and s == 'CB':
            accepted.append((r['L'] // 2, r['C1'] + 1))
        return None
    st, _ = P.run('TEST', {'C1': 1, 'W': 1, 'L': 2, 'B': 2}, hook)
    print('L3 corrected control (budget register B):')
    print(f'  accepted {accepted}; status {st[0]} in state {st[1]} with L = {st[3]["L"]} '
          f'(window n = {st[3]["L"] // 2})')
    ok_l3 = st[0] == 'HALT' and accepted == [(1, 2), (2, 5), (3, 11), (4, 17)] and st[3]['L'] == 10
    print(f'  halts while processing n = 5 after accepting 2, 5, 11, 17: {ok_l3}')
    print()
    print('POST-HOC ANATOMY DONE' if ok_g4 and ok_l3 else 'POST-HOC CHECK FAILED')
    return 0 if ok_g4 and ok_l3 else 1


if __name__ == '__main__':
    sys.exit(main())
