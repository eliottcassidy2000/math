#!/usr/bin/env python3
"""All-clock-shift successor decision via a proved fixed-intercept bound.

This is scalar replacement in an actual finite component language.  It does
not assert physical transport from the input word or a native return program.
"""
import argparse
from collections import deque
import hashlib
import json
from pathlib import Path

from component_inverse import Inverse, direct, make_bank_five, require, targets


def odd_part(x):
    require(x > 0, 'positive odd-part input required')
    return x // (x & -x)


class ShiftInverse:
    def __init__(self, bank):
        self.bank = bank
        self.inv = Inverse(bank)
        self.M = max(a.sigma[0] for a in bank)

    def f_bound(self, c):
        require(c != 0, 'zero intercept is outside this normalization theorem')
        B = max(self.M ** 2, abs(c))
        return max(max(a.sigma[2] for a in self.bank),
                   max(((a.sigma[2] + a.sigma[1]) * B - a.sigma[1] * c) // a.sigma[0]
                       for a in self.bank if a.ports))

    def prefixes(self, c, R):
        """One actual representative for every attainable bounded monoid pair."""
        found = {(1, 1): ()}
        queue = deque([(1, 1)])
        while queue:
            h, o = queue.popleft()
            for k, a in enumerate(self.bank):
                if a.ports:
                    continue
                H, P = a.sigma[:2]
                nxt = h * H, o * odd_part(P)
                if abs(c) % nxt[0] or R % nxt[1] or nxt in found:
                    continue
                found[nxt] = found[(h, o)] + (k,)
                queue.append(nxt)
        return found

    def successors(self, c, R):
        require(c != 0 and R > 0 and R % 2, 'bad normalized target')
        out, searches = [], []
        for (prefix_h, prefix_odd_p), prefix in sorted(self.prefixes(c, R).items()):
            target_c, target_R = c // prefix_h, R // prefix_odd_p
            limit = self.f_bound(target_c)
            f = target_R
            while f <= limit:
                words = self.inv.words(f)
                require(words == self.inv.words_right(f), 'shift inverse directions disagree')
                normalized = [w for w in words if self.bank[w[0]].ports]
                searches.append({'prefix_h': prefix_h, 'prefix_odd_p': prefix_odd_p,
                                 'normalized_c': target_c, 'f': f,
                                 'proved_f_bound': limit, 'candidates': len(normalized)})
                for word in normalized:
                    sigma = self.inv.sigma(word)
                    if sigma[0] - sigma[4] != target_c:
                        continue
                    full_word = prefix + word
                    full_sigma = self.inv.sigma(full_word)
                    require(full_sigma[0] - full_sigma[4] == c, 'prefix intercept error')
                    require(odd_part(full_sigma[2]) == R, 'prefix odd marginal error')
                    rows, ports = self.inv.graph(full_word)
                    literal = len(rows) <= 14
                    if literal:
                        require(direct(rows, ports) == full_sigma, 'shift witness literal mismatch')
                    out.append({'word': [self.bank[k].name for k in full_word],
                                'rows': rows, 'ports': ports, 'sigma': full_sigma,
                                'literal_checked': literal})
                f *= 2
        return {'c': c, 'R': R, 'searches': searches, 'witnesses': out}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-depth', type=int, default=8)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('shifted_returns.json'))
    args = parser.parse_args()
    solver = ShiftInverse(make_bank_five())
    pending = [('', -5, 7)]
    states = [{'address': '', 'c': -5, 'R': 7}]
    tests = []
    for depth in range(args.max_depth):
        nxt = []
        for address, c, R in pending:
            for d in (2, 4):
                child = solver.successors(3 * c + d, 3 * R)
                child['address'] = address + str(d)
                tests.append(child)
                if child['witnesses']:
                    nxt.append((child['address'], child['c'], child['R']))
                    states.append({k: child[k] for k in ('address', 'c', 'R')})
        print(f'depth={depth+1} accepted={[a for a,c,R in nxt]}', flush=True)
        pending = nxt
        if not pending:
            break
    receipt = {'scope': 'FINITE-EXACT all-clock-shift scalar returns, each successor decision exhaustive at all orders',
               'complete_tree': not pending, 'unexpanded_states': pending,
               'states': states, 'tests': tests,
               'normalization': 'bounded normalized intercept plus finite unmarked-prefix monoid',
               'max_atom_order': 5}
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'complete_tree': not pending, 'states': len(states),
                      'tests': len(tests), 'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
