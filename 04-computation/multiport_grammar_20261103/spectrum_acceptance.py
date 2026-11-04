#!/usr/bin/env python3
"""Complete bounded scalar spectrum for the declared ordered-component grammar.

For target H and clock t>=1, j<=f-e implies H>=1+f and H>=2**t.
Consequently f<H is a proved exhaustive limit, not a search-order cutoff.
The exceptional seed1 is represented separately by the singleton tournament.
"""
import argparse
import hashlib
import json
from pathlib import Path

from component_inverse import Inverse, attach, direct, make_bank_five, require, strong


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bound', type=int, default=2000)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('spectrum_acceptance.json'))
    args = parser.parse_args()
    require(args.bound >= 1, 'positive bound required')
    bank = make_bank_five()
    inv = Inverse(bank)
    witnesses = {order: {1: None} for order in (1, 3, 4, 5)}
    checks = {'left_right_equalities': 0, 'candidates': 0, 'literal_attachments': 0}
    for f in range(1, args.bound):
        words = inv.words(f)
        require(words == inv.words_right(f), f'inverse mismatch at f={f}')
        checks['left_right_equalities'] += 1
        for word in words:
            sigma = inv.sigma(word)
            h, p, ff, g, j, e, s = sigma
            require(j <= f - e and ff == f, 'fixed-target inequality failed')
            checks['candidates'] += 1
            largest = max(len(bank[k].rows) for k in word)
            size = sum(len(bank[k].rows) for k in word)
            clock = 1
            while (target := h - j + (1 << clock) * f) <= args.bound:
                require(target % 2 == 1 and target > f, 'target parity or bound failed')
                for order in witnesses:
                    if largest > order:
                        continue
                    previous = witnesses[order].get(target)
                    if previous is None or size + clock + 1 < previous['order']:
                        witnesses[order][target] = {
                            'word': [bank[k].name for k in word], 'clock': clock,
                            'sigma': sigma, 'order': size + clock + 1}
                clock += 1
    # Independent actual-graph verification of every selected order-five witness.
    names = {a.name: k for k, a in enumerate(bank)}
    for target, rec in sorted(witnesses[5].items()):
        if target == 1:
            continue
        word = tuple(names[name] for name in rec['word'])
        rows, ports = inv.graph(word)
        full = attach(rows, ports, rec['clock'])
        require(strong(full), f'non-strong witness for {target}')
        observed = direct(full, 1 << (len(full) - 1))
        require(observed[0] == target, f'literal count mismatch for {target}')
        require(observed[2] - observed[0] == rec['sigma'][4], f'root defect mismatch for {target}')
        checks['literal_attachments'] += 1
    missing = {order: [target for target in range(1, args.bound + 1, 2)
                       if target not in values] for order, values in witnesses.items()}
    receipt = {
        'scope': 'FINITE-EXACT; complete bounded scalar acceptance, no source-preserving closure or unbounded coverage claimed',
        'bound': args.bound, 'C': str(inv.C), 'checks': checks,
        'atoms': [{'name': a.name, 'rows': a.rows, 'ports': a.ports, 'sigma': a.sigma} for a in bank],
        'missing_by_max_component_order': missing,
        'witnesses': witnesses[5],
        'boundary': 'H=1 is a separate singleton seed; all other witnesses use a nonempty transitive prefix and terminal-port attachment',
    }
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'bound': args.bound,
                      'missing_counts': {k: len(v) for k, v in missing.items()},
                      'missing_max_order_five': missing[5], 'checks': checks,
                      'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
