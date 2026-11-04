#!/usr/bin/env python3
"""Optional depth-bounded BFS; each completed successor query is exhaustive.

The root six-atom receipt and source are frozen inputs. A nonempty frontier
at the depth cap is explicitly incomplete, never a complete return theorem.
"""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PINS = {
    'atom_six_probe.py': 'ed295d2f7ce1a7815a7c377ede9544f3118efe3736458d1fc921d460c21900ba',
    'atom_six_probe.json': '1a57397d99deacf460d6e7c01f80da705b0304240d0fec363108cedb5bc4d951',
}
for name, expected in PINS.items():
    if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != expected:
        raise RuntimeError('six-atom input pin mismatch')

from atom_six_probe import BoundedInverse
from component_inverse import Atom, direct
from shifted_returns import ShiftInverse, odd_part


def query(solver, grouped, c, R):
    searches, witnesses = [], []
    for (pre_h, pre_odd), prefix in sorted(solver.prefixes(c, R).items()):
        target_c, target_R = c // pre_h, R // pre_odd
        limit = solver.f_bound(target_c)
        B = max(solver.M ** 2, abs(target_c))
        f = target_R
        while f <= limit:
            candidates = 0
            for h in solver.inv.h_values(min(f, B)):
                for word in grouped.left(h, f):
                    if not solver.bank[word[0]].ports:
                        continue
                    candidates += 1
                    sigma = solver.inv.sigma(word)
                    if sigma[0] - sigma[4] != target_c:
                        continue
                    full_word = prefix + word
                    actual = solver.inv.sigma(full_word)
                    if actual[0] - actual[4] != c or odd_part(actual[2]) != R:
                        raise RuntimeError('bad reconstructed target')
                    rows, ports = solver.inv.graph(full_word)
                    checked = len(rows) <= 14
                    if checked and direct(rows, ports) != actual:
                        raise RuntimeError('literal witness mismatch')
                    witnesses.append({'word': [solver.bank[k].name for k in full_word],
                                      'rows': rows, 'ports': ports, 'sigma': actual,
                                      'literal_checked': checked})
            searches.append({'prefix': [pre_h, pre_odd], 'c': target_c, 'f': f,
                             'h_bound': B, 'f_bound': limit, 'candidates': candidates})
            f *= 2
    return {'c': c, 'R': R, 'searches': searches, 'witnesses': witnesses}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=8)
    parser.add_argument('--output', type=Path, default=HERE / 'six_tree_probe.json')
    args = parser.parse_args()
    if args.depth < 1:
        raise RuntimeError('positive depth required')
    receipt = json.loads((HERE / 'atom_six_probe.json').read_text())
    bank = tuple(Atom(r['name'], tuple(r['rows']), r['ports'], tuple(r['sigma']))
                 for r in receipt['bank'])
    solver, grouped = ShiftInverse(bank), BoundedInverse(bank)
    frontier = [('', -5, 7)]
    states = [{'address': '', 'c': -5, 'R': 7}]
    tests = []
    for depth in range(1, args.depth + 1):
        next_frontier = []
        for address, c, R in frontier:
            for digit in (2, 4):
                child = query(solver, grouped, 3 * c + digit, 3 * R)
                child['address'] = address + str(digit)
                tests.append(child)
                if child['witnesses']:
                    next_frontier.append((child['address'], child['c'], child['R']))
                    states.append({k: child[k] for k in ('address', 'c', 'R')})
                print(json.dumps({'target': child['address'], 'accepted': bool(child['witnesses']),
                                  'witnesses': len(child['witnesses'])}), flush=True)
        frontier = next_frontier
        output = {'scope': 'bounded BFS; each individual successor query is exhaustive at all word lengths and shifts',
                  'input_pins': PINS, 'requested_depth_cap': args.depth,
                  'completed_depth': depth, 'complete_tree': not frontier,
                  'frontier': frontier, 'states': states, 'tests': tests}
        text = json.dumps(output, indent=2, sort_keys=True) + '\n'
        args.output.write_text(text)
        print(json.dumps({'completed_depth': depth, 'frontier': [r[0] for r in frontier],
                          'complete_tree': not frontier}), flush=True)
        if not frontier:
            break
    print(json.dumps({'bytes': len(text.encode()), 'sha256': hashlib.sha256(text.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
