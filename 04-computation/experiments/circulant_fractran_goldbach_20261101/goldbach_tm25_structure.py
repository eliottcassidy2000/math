#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: the 25-state Goldbach Turing machine (Leng 2024, Lean-verified; table pinned in
05-knowledge/results/session65_cryptid_tournament_interfaces_20261026/session65-twentyfive/inputs/goldbach_tm/
GoldbachTm/Tm25/TuringMachine25.lean) - an independent simulator and a structural reading of its control graph.

The table is parsed directly from the Lean source (write, then move). Prints:
  (1) the table, the halt slot, the control digraph (state -> next state per symbol), self-loops ("scanners"),
      strongly connected components;
  (2) a run from the blank tape: the step at which each even target n first appears as the stage configuration
      (state 21 with exactly n ones to the left and an empty right side, as in Lean's cfg45 for n = 4), compared
      with session65's recorded first times (4:45, 6:132, 8:301, 10:596);
  (3) which control states are used between consecutive stages, and the per-stage step counts;
  (4) a coarse block decomposition of the control graph by the visit pattern.

usage: python3 goldbach_tm25_structure.py [STEPS]
"""
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEAN = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'session65_cryptid_tournament_interfaces_20261026',
                    'session65-twentyfive', 'inputs', 'goldbach_tm', 'GoldbachTm', 'Tm25', 'TuringMachine25.lean')


def parse():
    src = open(LEAN).read()
    pat = re.compile(r"\|\s*⟨\s*(\d+),\s*_⟩,\s*Γ\.(zero|one)\s*=>\s*(none|some\s*⟨⟨\s*(\d+),\s*by omega⟩,\s*⟨Turing\.Dir\.(left|right),\s*Γ\.(zero|one)⟩⟩)")
    tab = {}
    for m in pat.finditer(src):
        q, sym = int(m.group(1)), 0 if m.group(2) == 'zero' else 1
        if m.group(3) == 'none':
            tab[(q, sym)] = None
        else:
            tab[(q, sym)] = (int(m.group(4)), -1 if m.group(5) == 'left' else 1, 0 if m.group(6) == 'zero' else 1)
    assert len(tab) == 50, len(tab)
    return tab


def scc(nodes, edges):
    index, low, onstack, stack, out = {}, {}, set(), [], []
    counter = [0]
    sys.setrecursionlimit(10000)

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        onstack.add(v)
        for w in edges[v]:
            if w not in index:
                strong(w)
                low[v] = min(low[v], low[w])
            elif w in onstack:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack.discard(w)
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp))
    for v in nodes:
        if v not in index:
            strong(v)
    return out


def run(tab, steps):
    tape = collections.defaultdict(int)
    head, q = 0, 0
    stages = {}
    stage_states = []
    cur_states = collections.Counter()
    last_stage_t = 0
    for t in range(1, steps + 1):
        a = tab[(q, tape[head])]
        if a is None:
            return stages, stage_states, ('HALT', t)
        nq, mv, wr = a
        tape[head] = wr
        head += mv
        q = nq
        cur_states[q] += 1
        if q == 21:
            # stage test: exactly n ones immediately left of the head (contiguous), nothing else on the tape
            nz = [i for i, v in tape.items() if v]
            if nz and tape[head] == 0 and max(nz) == head - 1 and min(nz) == head - len(nz) \
                    and len(nz) % 2 == 0 and len(nz) >= 4:
                n = len(nz)
                if n not in stages:
                    stages[n] = t
                    stage_states.append((n, t, t - last_stage_t, dict(cur_states)))
                    cur_states = collections.Counter()
                    last_stage_t = t
    return stages, stage_states, ('RUNNING', steps)


def main():
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    tab = parse()
    print('(1) TABLE (state, read) -> (next, move, write); write then move')
    for q in range(25):
        row = []
        for s in (0, 1):
            a = tab[(q, s)]
            row.append(f'{s}: HALT' if a is None else f'{s}: ->{a[0]:>2} {"L" if a[1] < 0 else "R"} w{a[2]}')
        print(f'  {q:>2}  ' + '   '.join(row))
    halts = [k for k, v in tab.items() if v is None]
    print(f'  halt slots: {halts}; instructions: {sum(v is not None for v in tab.values())}')
    scanners = [(q, s) for (q, s), a in tab.items() if a and a[0] == q]
    print(f'  self-loops (scanners): {len(scanners)}: {sorted(scanners)}')
    print('    each rewrites its read symbol? ',
          {f'{q}/{s}': (tab[(q, s)][2] == s) for q, s in sorted(scanners)})
    edges = {q: sorted({a[0] for s in (0, 1) for a in [tab[(q, s)]] if a and a[0] != q}) for q in range(25)}
    comps = scc(range(25), edges)
    print(f'  strongly connected components of the control digraph (self-loops dropped): {comps}')
    indeg = collections.Counter(w for v in edges for w in edges[v])
    print(f'  out-degree / in-degree (distinct non-self targets): '
          f'{[(q, len(edges[q]), indeg[q]) for q in range(25)]}')
    print()
    print(f'(2) RUN from the blank tape for {steps} steps')
    stages, stage_states, status = run(tab, steps)
    print(f'  status: {status}')
    first = sorted(stages.items())
    print(f'  first stage times: {first[:12]} ... last {first[-1] if first else None}')
    s65 = {4: 45, 6: 132, 8: 301, 10: 596}
    print(f'  matches session65 first times {s65}: {all(stages.get(n) == t for n, t in s65.items())}')
    print()
    print('(3) control states used between consecutive stages')
    allused = set()
    for n, t, dt, cs in stage_states[:8]:
        allused |= set(cs)
        print(f'  stage n={n:>3} at t={t:>8} (+{dt:>7} steps): states visited {sorted(cs)}')
    tot = collections.Counter()
    for n, t, dt, cs in stage_states:
        tot.update(cs)
    print(f'  all {len(tot)} states used over the run; visit share of the 6 busiest: '
          f'{[(q, round(c / sum(tot.values()), 3)) for q, c in tot.most_common(6)]}')
    print()
    print('(4) stage step counts dt(n) and their growth')
    dts = [(n, dt) for n, t, dt, cs in stage_states]
    print(f'  {dts[:15]}')
    if len(dts) > 6:
        import math
        xs = [math.log(n) for n, _ in dts[5:]]
        ys = [math.log(dt) for _, dt in dts[5:]]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        print(f'  log-log slope of dt(n) over n = {dts[5][0]}..{dts[-1][0]}: {b:.3f}')
    print()
    print('ALL CHECKS DONE')


if __name__ == '__main__':
    main()
