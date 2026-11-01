#!/usr/bin/env python3
"""Independent audit of the 25-state Goldbach Turing machine claims (results note section 2.3,
goldbach_tm25_structure.out). The table is PARSED from the Lean source
(session65 .../Tm25/TuringMachine25.lean); semantics as in that file's `step`:
(M q T.1).map fun (q', a) => (q', (T.write a.write).move a.move)   -- write, then move.
Stage n (Content.lean `never_halt_step`): state 21, head on a 0, exactly n ones immediately left of the head,
nothing else on the tape.

usage (from the repo root):
  python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit/audit_tm25.py [STEPS]
"""
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
LEAN = ROOT / ('05-knowledge/results/session65_cryptid_tournament_interfaces_20261026/session65-twentyfive/inputs/'
               'goldbach_tm/GoldbachTm/Tm25/TuringMachine25.lean')
FAILS, MISMATCH = [], []


def report(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


def claim(name, ok, detail=''):
    print(f'  [{"CLAIM OK" if ok else "CLAIM WRONG"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        MISMATCH.append(name)


def parse():
    src = LEAN.read_text()
    T = {}
    pat = re.compile(r'\|\s*⟨\s*(\d+),\s*_⟩,\s*Γ\.(zero|one)\s*=>\s*(none|some\s*⟨⟨\s*(\d+),\s*by omega⟩,\s*'
                     r'⟨Turing\.Dir\.(left|right),\s*Γ\.(zero|one)⟩⟩)')
    for m in pat.finditer(src):
        q, sym = int(m.group(1)), 1 if m.group(2) == 'one' else 0
        if m.group(3) == 'none':
            T[(q, sym)] = None
        else:
            T[(q, sym)] = (int(m.group(4)), -1 if m.group(5) == 'left' else 1, 1 if m.group(6) == 'one' else 0)
    return T


def scc(nodes, edges):
    """Tarjan, iterative."""
    index, low, on, stack, comps, counter = {}, {}, set(), [], [], [0]
    adj = {v: sorted(set(w for (a, w) in edges if a == v)) for v in nodes}
    for root in nodes:
        if root in index:
            continue
        work = [(root, 0)]
        while work:
            v, i = work.pop()
            if i == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on.add(v)
            recurse = False
            for j in range(i, len(adj[v])):
                w = adj[v][j]
                if w not in index:
                    work.append((v, j + 1))
                    work.append((w, 0))
                    recurse = True
                    break
                if w in on:
                    low[v] = min(low[v], index[w])
            if recurse:
                continue
            if low[v] == index[v]:
                comp = []
                while True:
                    w = stack.pop()
                    on.discard(w)
                    comp.append(w)
                    if w == v:
                        break
                comps.append(sorted(comp))
            if work:
                u = work[-1][0]
                low[u] = min(low[u], low[v])
    return comps


def main():
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    T = parse()
    report('Lean table parsed: 25 states x 2 symbols = 50 slots', len(T) == 50 and all((q, s) in T for q in range(25)
                                                                                         for s in (0, 1)))
    halts = [k for k, v in T.items() if v is None]
    claim('one halt slot (24, 0); 49 instructions', halts == [(24, 0)] and sum(v is not None for v in T.values()) == 49,
          f'halt slots {halts}')
    loops = sorted(k for k, v in T.items() if v is not None and v[0] == k[0])
    claim('14 self-loops', len(loops) == 14, str(loops))
    edges = set((q, v[0]) for (q, s), v in T.items() if v is not None and v[0] != q)
    comps = scc(list(range(25)), edges)
    claim('control digraph (self-loops dropped) is one strongly connected component', len(comps) == 1,
          f'{len(comps)} components')
    # ---- simulate
    tape = {}
    head, q, t = 0, 0, 0
    ones = 0
    visits = [0] * 25
    stages = []
    halted = False
    while t < nsteps:
        # stage test: state 21, head on 0, exactly n ones immediately left of the head, nothing else
        sym = tape.get(head, 0)
        if q == 21 and sym == 0 and ones > 0 and all(tape.get(head - i, 0) == 1 for i in range(1, ones + 1)):
            if not stages or stages[-1][1] != t:
                stages.append((ones, t))
        tr = T[(q, sym)]
        if tr is None:
            halted = True
            break
        visits[q] += 1
        nq, mv, wr = tr
        if wr != sym:
            ones += 1 if wr == 1 else -1
            tape[head] = wr
        head += mv
        q = nq
        t += 1
    report(f'{nsteps} steps from the blank tape, no halt', not halted, f'{t} steps')
    ns = [n for n, _ in stages]
    report('stages are exactly n = 4, 6, 8, ... in order (no other stage-shaped configuration)',
           ns == list(range(4, 4 + 2 * len(ns), 2)), f'last stage n = {ns[-1]}')
    st = dict(stages)
    claim('stage times 45, 132, 301, 596 and n = 74 at 1,793,334',
          [st.get(n) for n in (4, 6, 8, 10)] == [45, 132, 301, 596] and st.get(74) == 1793334,
          f'{[(n, st.get(n)) for n in (4, 6, 8, 10, 74)]}; last stage within the run: {stages[-1]}')
    tot = sum(visits)
    share = sorted(((visits[i] / tot, i) for i in range(25)), reverse=True)
    top4 = [i for _, i in share[:4]]
    claim('busiest states 11, 13, 10, 7 with ~88% of steps', top4 == [11, 13, 10, 7] and
          abs(sum(s for s, _ in share[:4]) - 0.88) < 0.01,
          f'top four {[(i, round(s, 3)) for s, i in share[:4]]}, sum {sum(s for s, _ in share[:4]):.3f}; '
          f'states used: {sum(v > 0 for v in visits)}')
    dts = [(stages[i][0], stages[i][1] - (stages[i - 1][1] if i else 0)) for i in range(len(stages))]
    sel = [(n, d) for n, d in dts if 14 <= n <= 74]
    xs = [math.log(n) for n, _ in sel]
    ys = [math.log(d) for _, d in sel]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    end = math.log(sel[-1][1] / sel[0][1]) / math.log(sel[-1][0] / sel[0][0])
    claim('least-squares log-log slope of the stage cost dt(n), n = 14..74, is ~3.26', abs(b - 3.26) < 0.01,
          f'LSQ slope {b:.3f}; endpoint slope {end:.3f}; local slopes vary (n=14..30: '
          f'{math.log(dict(dts)[30] / dict(dts)[14]) / math.log(30 / 14):.2f}, n=50..74: '
          f'{math.log(dict(dts)[74] / dict(dts)[50]) / math.log(74 / 50):.2f})')
    print()
    print('SESSION CLAIMS FOUND WRONG: ' + ('; '.join(MISMATCH) if MISMATCH else 'none'))
    print('AUDIT CHECK FAILURES: ' + ('; '.join(FAILS) if FAILS else 'none'))
    return 1 if FAILS else 0


if __name__ == '__main__':
    sys.exit(main())
