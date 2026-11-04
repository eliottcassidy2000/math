#!/usr/bin/env python3
"""Complete fixed-marginal inverse for actual ordered component words.

The order bound is proved, not chosen: h <= f and p <= C f.  Two independent
inverse directions enumerate the same words.  Literal subset-path counting
checks signatures and selected attachments independently of the transfer.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def direct(rows, ports):
    """All seven observables from literal paths on complementary subsets."""
    n = len(rows)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    sp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
        sp[1 << v][v] = (ports >> v) & 1
    for mask in range(1, size):
        for v in range(n):
            if not dp[mask][v]:
                continue
            choices = rows[v] & (size - 1 - mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                w = bit.bit_length() - 1
                dp[mask | bit][w] += dp[mask][v]
                sp[mask | bit][w] += sp[mask][v]
    paths = [sum(row) for row in dp]
    paths[0] = 1
    ends = [sum(row[v] for v in range(n) if ports >> v & 1) for row in dp]
    starts = [sum(row) for row in sp]
    full = size - 1
    return (paths[full], sum(paths[m] * paths[full ^ m] for m in range(size)),
            sum(ends[m] * paths[full ^ m] for m in range(size)),
            sum(paths[m] * starts[full ^ m] for m in range(size)),
            sum(ends[m] * starts[full ^ m] for m in range(size)),
            ends[full], starts[full])


def product(x, y):
    h, p, f, g, j, e, s = x
    H, P, F, G, J, E, S = y
    return (h * H, p * P, f * H + p * F, g * P + h * G,
            j * H + h * J + e * S + g * F, h * E, s * H)


@dataclass(frozen=True)
class Atom:
    name: str
    rows: tuple[int, ...]
    ports: int
    sigma: tuple[int, ...]


def make_bank(include_c4=True):
    shapes = [('T', (0,), (0, 1)), ('C', (2, 4, 1), (0, 1, 3, 7))]
    if include_c4:
        shapes.append(('Q', (8, 1, 3, 6), range(16)))
    return tuple(Atom(f'{name}{ports:x}', rows, ports, direct(rows, ports))
                 for name, rows, masks in shapes for ports in masks)


def make_bank_five():
    """All labelled strong order-five atoms/all ports, modulo full signature.

    The quotient is used only for ordered-assembly scalar acceptance. Each
    representative retains actual rows and ports; it does not identify names.
    """
    bank = list(make_bank())
    seen = {a.sigma for a in bank}
    edges = tuple(itertools.combinations(range(5), 2))
    strong_count = 0
    for mask in range(1 << len(edges)):
        rows = [0] * 5
        for k, (i, j) in enumerate(edges):
            u, v = (i, j) if mask >> k & 1 else (j, i)
            rows[u] |= 1 << v
        if not strong(rows):
            continue
        strong_count += 1
        for ports in range(32):
            sigma = direct(rows, ports)
            if sigma not in seen:
                seen.add(sigma)
                bank.append(Atom(f'V5_{mask:x}_{ports:x}', tuple(rows), ports, sigma))
    require(strong_count == 544, 'strong order-five census changed')
    require(len(bank) == 170, 'signature quotient changed')
    return tuple(bank)


class Inverse:
    def __init__(self, bank):
        self.bank = bank
        self.C = max(Fraction(a.sigma[1], a.sigma[2]) for a in bank if a.ports)

    @lru_cache(None)
    def left(self, h, f):
        """Peel the first atom; no empty recursive tails are admitted."""
        out = []
        if h > f or h < 1 or f < 1:
            return ()
        for k, a in enumerate(self.bank):
            H, P, F = a.sigma[:3]
            if a.ports and (h, f) == (H, F):
                out.append((k,))
            if h % H:
                continue
            tail_h = h // H
            numerator = f - F * tail_h
            if numerator <= 0 or numerator % P:
                continue
            tail_f = numerator // P
            out.extend((k,) + word for word in self.left(tail_h, tail_f))
        return tuple(out)

    def h_values(self, f):
        seen, pending = {1}, [1]
        factors = set(a.sigma[0] for a in self.bank)
        while pending:
            h = pending.pop()
            for factor in factors:
                nxt = h * factor
                if nxt <= f and nxt not in seen:
                    seen.add(nxt)
                    pending.append(nxt)
        return sorted(seen)

    def words(self, f):
        return tuple(sorted(word for h in self.h_values(f) for word in self.left(h, f)))

    def hp_values(self, f):
        """Finite monoid in (h,p), independent of the left recursion."""
        seen, pending = {(1, 1)}, [(1, 1)]
        factors = set(a.sigma[:2] for a in self.bank)
        while pending:
            h, p = pending.pop()
            for H, P in factors:
                nxt = h * H, p * P
                if nxt[0] <= f and nxt[1] <= self.C * f and nxt not in seen:
                    seen.add(nxt)
                    pending.append(nxt)
        return sorted(seen)

    @lru_cache(None)
    def right_prefix(self, h, p, f):
        """Right inverse of arbitrary prefixes; their ports may all be empty."""
        if (h, p, f) == (1, 1, 0):
            return ((),)
        if h < 1 or p < 2 or f < 0:
            return ()
        out = []
        for k, a in enumerate(self.bank):
            H, P, F = a.sigma[:3]
            if h % H or p % P:
                continue
            pre_h, pre_p = h // H, p // P
            numerator = f - pre_p * F
            if numerator < 0 or numerator % H:
                continue
            out.extend(word + (k,) for word in
                       self.right_prefix(pre_h, pre_p, numerator // H))
        return tuple(out)

    def words_right(self, f):
        out = []
        for h, p in self.hp_values(f):
            out.extend(word for word in self.right_prefix(h, p, f)
                       if word and self.bank[word[-1]].ports)
        return tuple(sorted(out))

    def sigma(self, word):
        value = self.bank[word[0]].sigma
        for k in word[1:]:
            value = product(value, self.bank[k].sigma)
        return value

    def graph(self, word):
        rows, ports = [], 0
        for k in word:
            a = self.bank[k]
            n = len(rows)
            new = ((1 << len(a.rows)) - 1) << n
            rows = [row | new for row in rows] + [row << n for row in a.rows]
            ports |= a.ports << n
        return tuple(rows), ports


def strong(rows):
    full = (1 << len(rows)) - 1
    for v in range(len(rows)):
        seen = 1 << v
        old = -1
        while old != seen:
            old = seen
            for w, row in enumerate(rows):
                if seen >> w & 1:
                    seen |= row
        if seen != full:
            return False
    return True


def attach(rows, ports, clock):
    n = len(rows)
    root = clock + n
    suffix = ((1 << n) - 1) << clock
    out = [(((1 << clock) - 1) ^ ((1 << (v + 1)) - 1)) | suffix
           for v in range(clock)]
    out += [(row << clock) | ((1 << root) if ports >> v & 1 else 0)
            for v, row in enumerate(rows)]
    out += [((1 << clock) - 1) | (suffix ^ (ports << clock))]
    return tuple(out)


def targets(level):
    for digits in itertools.product((2, 4), repeat=level):
        intercept = -5
        for d in digits:
            intercept = 3 * intercept + d
        yield ''.join(map(str, digits)), intercept


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--levels', type=int, default=4)
    parser.add_argument('--max-atom-order', type=int, choices=(4, 5), default=4)
    parser.add_argument('--literal-limit', type=int, default=12)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('component_inverse.json'))
    args = parser.parse_args()
    require(args.levels >= 1, 'positive levels required')
    inv = Inverse(make_bank_five() if args.max_atom_order == 5 else make_bank())
    small = Inverse(make_bank(False))
    checks = {'inverse_equalities': 0, 'literal_signatures': 0, 'attachments': 0}
    # Inherited exact acceptance is a positive/hostile regression control.
    expected = {1: ['2', '4'], 2: ['22', '42'], 3: ['422']}
    for level, accepted in expected.items():
        words = small.words(7 * 3 ** level)
        require(words == small.words_right(7 * 3 ** level), 'old bank inverse disagreement')
        checks['inverse_equalities'] += 1
        intercepts = {small.sigma(w)[0] - small.sigma(w)[4] for w in words}
        require([d for d, a in targets(level) if a in intercepts] == accepted, 'old bank regression')
    # Complete comparison on a dense small universe, not just target f.
    for f in range(1, 101):
        require(inv.words(f) == inv.words_right(f), f'inverse disagreement at f={f}')
        checks['inverse_equalities'] += 1
    levels = []
    for level in range(args.levels + 1):
        f = 7 * 3 ** level
        words = inv.words(f)
        require(words == inv.words_right(f), f'target inverse disagreement at {f}')
        checks['inverse_equalities'] += 1
        records = []
        for word in words:
            sigma = inv.sigma(word)
            rows, ports = inv.graph(word)
            require(sigma[2] == f and sigma[0] <= f, 'wrong marginal or h bound')
            require(sigma[1] <= inv.C * f, 'wrong p bound')
            literal = len(rows) <= args.literal_limit
            if literal:
                require(direct(rows, ports) == sigma, 'literal signature disagreement')
                checks['literal_signatures'] += 1
            records.append({'word': [inv.bank[k].name for k in word],
                            'rows': rows, 'ports': ports, 'sigma': sigma,
                            'intercept': sigma[0] - sigma[4], 'literal_checked': literal})
        found = {}
        for digits, a in targets(level):
            witnesses = [i for i, rec in enumerate(records) if rec['intercept'] == a]
            found[digits] = {'intercept': a, 'witness_indices': witnesses}
            if witnesses:
                # Recount the smallest representative at two independent clocks.
                best = min(witnesses, key=lambda i: len(records[i]['rows']))
                rec = records[best]
                if len(rec['rows']) <= args.literal_limit:
                    for clock in (1, 2):
                        full = attach(rec['rows'], rec['ports'], clock)
                        require(strong(full), 'accepted attachment is not strong')
                        observed = direct(full, 1 << (len(full) - 1))
                        require(observed[0] == a + (1 << clock) * f, 'attachment count disagreement')
                        require(observed[2] - observed[0] == rec['sigma'][4], 'root defect disagreement')
                        checks['attachments'] += 1
        level_record = {'level': level, 'f': f, 'word_count': len(words),
                        'intercepts': sorted({r['intercept'] for r in records}),
                        'targets': found, 'candidates': records}
        levels.append(level_record)
        accepted = [d for d, v in found.items() if v['witness_indices']]
        print(f'level={level} f={f} candidates={len(words)} accepted={accepted}', flush=True)
    receipt = {'scope': 'FINITE-EXACT complete fixed-f acceptance in the declared atom bank; no native return claimed',
               'max_atom_order': args.max_atom_order,
               'atoms': [{'name': a.name, 'rows': a.rows, 'ports': a.ports, 'sigma': a.sigma} for a in inv.bank],
               'C': str(inv.C), 'levels': levels, 'checks': checks}
    content = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    args.output.write_text(content)
    print(json.dumps({'checks': checks, 'sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
