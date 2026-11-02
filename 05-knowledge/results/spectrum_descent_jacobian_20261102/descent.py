"""Exact minimum-cost inverse in the inherited singleton/D language.

The state and predecessor formulas are inherited from packet 20261031b.
This consumer minimizes actual added vertices instead of accepting the first
inverse branch. It never replaces source membership by lattice feasibility.
"""
from functools import lru_cache
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
PRIOR = ROOT/'05-knowledge/results/square_sum_spectrum_overnight_20261031b'
sys.path.insert(0, str(PRIOR))
from diamond_profiles import profiles, apply, single

BLOCKS = profiles()
SEED = (1, 0, 0, 1)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def g_values(h):
    return tuple(i*(h//5) for i in range(6)) if h % 5 == 0 else (0, h)


def admissible(state):
    F, q, g, h = state
    if min(F, h) < 1 or q < 0 or F % 2 != 1 or F < h:
        return False
    r = h
    while r % 5 == 0:
        r //= 5
    return (r == 1 and g in g_values(h) and q >= F-h+g
            and (q-F+h-g) % 2 == 0)


def f(q):
    total = 1
    while q:
        parent = q//2
        if q % 2 == 0:
            total += parent+parent % 2
        q = parent
    return total


def step(state, letter):
    kind, cut = letter
    return single(state, cut) if kind == 'S' else apply(state, BLOCKS[cut])


def rank(state):
    F, q, g, h = state
    return F+q+h


def predecessors(state):
    if not admissible(state):
        return ()
    F, q, g, h = state
    candidates = []
    if g == 0 and q % 2 == 0:
        for oldg in g_values(h):
            candidates.append(((F-q//2-oldg, q//2, oldg, h), ('S', 0)))
    if g == h and q >= h and (q-h) % 2 == 0:
        for oldg in g_values(h):
            candidates.append(((F, (q-h)//2, oldg, h), ('S', 1)))
    if h > 1:
        oldh = h//5
        digit = (q*pow(oldh, -1, 32)) % 32
        for b in BLOCKS:
            if b['q'] != digit or g != (5-b['t'])*oldh:
                continue
            numerator = q-digit*oldh
            if numerator < 0 or numerator % 32:
                continue
            oldq = numerator//32
            for oldg in g_values(oldh):
                num = F-b['p']*oldq-(5-b['s'])*oldg-(b['F']-5)*oldh
                if num % 5 == 0:
                    candidates.append(((num//5, oldq, oldg, oldh), ('D', b['mask'])))
    out = []
    for previous, letter in candidates:
        if not admissible(previous) or previous == state:
            continue
        need(rank(previous) < rank(state), 'strict inherited inverse rank')
        need(step(previous, letter) == state, 'literal predecessor replay')
        out.append((previous, letter))
    return tuple(sorted(set(out)))


def height_exponent(h):
    r = 0
    while h > 1:
        need(h % 5 == 0, 'height is a power of five')
        h //= 5
        r += 1
    return r


def clock_lower_bound(state):
    _, q, _, h = state
    r = height_exponent(h)
    ceil_log = (q+h-1).bit_length()
    return 4*r+max(0, ceil_log-5*r)


@lru_cache(None)
def optimal(state):
    if not admissible(state):
        return None
    F, q, g, h = state
    if h == 1:
        if g == q % 2 and F == f(q):
            word = tuple(('S', int(b)) for b in bin(q)[2:]) if q else ()
            return len(word), word
        return None
    best = None
    for previous, letter in predecessors(state):
        cost = 1 if letter[0] == 'S' else 4
        if best is not None and clock_lower_bound(previous)+cost > best[0]:
            continue
        answer = optimal(previous)
        if answer is None:
            continue
        candidate = answer[0]+cost, answer[1]+(letter,)
        if best is None or candidate < best:
            best = candidate
    if best is not None:
        need(best[0] >= clock_lower_bound(state), 'clock lower certificate')
    return best


@lru_cache(None)
def greedy(state):
    """Precisely specified hostile: choose the smallest predecessor rank."""
    if not admissible(state):
        return None
    if state[3] == 1:
        return optimal(state)
    options = predecessors(state)
    if not options:
        return None
    previous, letter = min(options, key=lambda item: (rank(item[0]), item))
    answer = greedy(previous)
    if answer is None:
        return None
    return answer[0]+(1 if letter[0] == 'S' else 4), answer[1]+(letter,)


@lru_cache(None)
def jump_optimal(state):
    """Strip a whole singleton suffix, then one D atom, in one inverse step."""
    if not admissible(state):
        return None
    F, q, g, h = state
    if h == 1:
        if g != q % 2 or F != f(q):
            return None
        word = tuple(('S', int(b)) for b in bin(q)[2:]) if q else ()
        return len(word), word
    r = height_exponent(h)
    if q == 0:
        need(F == h and g == 0, 'zero slope is the canonical ordinal product')
        return 4*r, (('D', 0),)*r
    best = None
    oldh = h//5
    for length in range(q.bit_length()+1):
        modulus = 1 << length
        address = (q*pow(h, -1, modulus)) % modulus if length else 0
        remainder = q-h*address
        if remainder < 0:
            continue
        qD = remainder//modulus
        # Leading zero singletons at slope zero are removable identities.
        if qD == 0 and length and address.bit_length() != length:
            continue
        digit = (qD*pow(oldh, -1, 32)) % 32
        matches = [b for b in BLOCKS if b['q'] == digit]
        need(len(matches) <= 1, 'inherited injective D cut digits')
        if not matches:
            continue
        b = matches[0]
        numerator = qD-digit*oldh
        if numerator < 0 or numerator % 32:
            continue
        oldq = numerator//32
        gD = (5-b['t'])*oldh
        suffix = tuple(('S', (address >> j) & 1) for j in reversed(range(length)))
        temporary = (0, qD, gD, h)
        for letter in suffix:
            temporary = single(temporary, letter[1])
        increase, finalq, finalg, _ = temporary
        if (finalq, finalg) != (q, g):
            continue
        FD = F-increase
        possible_g = (oldq % 2,) if oldh == 1 else g_values(oldh)
        for oldg in possible_g:
            num = FD-b['p']*oldq-(5-b['s'])*oldg-(b['F']-5)*oldh
            if num % 5:
                continue
            previous = (num//5, oldq, oldg, oldh)
            if not admissible(previous):
                continue
            need(rank(previous) < rank(state), 'jump inverse decreases rank')
            cost = length+4
            if best is not None and clock_lower_bound(previous)+cost > best[0]:
                continue
            answer = jump_optimal(previous)
            if answer is None:
                continue
            candidate = answer[0]+cost, answer[1]+(('D', b['mask']),)+suffix
            if best is None or candidate < best:
                best = candidate
    return best


def actual(word):
    """Direct graph constructor, starting with only the marked apex.

    The inherited constructor keeps one extra dominating padding vertex.
    Removing that vertex preserves all four responses and saves one vertex.
    """
    rows = [0]
    for kind, cut in word:
        block = (8, 1, 3, 6) if kind == 'D' else (0,)
        n, size = len(rows), len(block)
        rows.extend([0]*size)
        for v in range(1, n):
            rows[v] |= ((1 << size)-1) << n
        for i in range(size):
            rows[n+i] |= block[i] << n
            if cut >> i & 1:
                rows[0] |= 1 << (n+i)
            else:
                rows[n+i] |= 1
    return rows


def hamilton(rows):
    n = len(rows)
    table = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            ways = table[mask][v]
            if not ways:
                continue
            nxt = rows[v] & ~mask
            while nxt:
                bit = nxt & -nxt
                nxt -= bit
                table[mask | bit][bit.bit_length()-1] += ways
    return sum(table[-1])


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state', nargs=4, type=int, metavar=('F', 'q', 'g', 'h'))
    args = parser.parse_args()
    target = tuple(args.state)
    answer = jump_optimal(target)
    result = {'state': target, 'grammar_member': answer is not None,
              'scope': 'singleton/D grammar, not arbitrary-tournament minimum'}
    if answer is not None:
        result.update(minimum_vertices=answer[0]+1, word=answer[1],
                      rows=actual(answer[1]), clock_lower_bound=clock_lower_bound(target)+1)
    print(json.dumps(result, sort_keys=True))
