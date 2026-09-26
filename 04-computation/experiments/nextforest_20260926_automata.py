"""Exact residual-state, macro-carry, and zero-tail controls. No dependencies."""
from collections import deque
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def shortcut(n):
    return (3 * n + 1) // 2 if n & 1 else n // 2


def transition(state, bit):
    a, b = state
    v = b + 3 ** a * bit
    p = v & 1
    return (a + p, ((1 + 2 * p) * v + p) // 2), p


def direct_section(n, length):
    a = 0
    for _ in range(length):
        a += n & 1
        n = shortcut(n)
    return a, n


def section_scan(n, length):
    state = (0, 0)
    word = 0
    for j in range(length):
        state, bit = transition(state, (n >> j) & 1)
        word |= bit << j
    return state, word


def bits_output(m, c, u, width):
    value = 0
    for j in range(width):
        v = m * ((u >> j) & 1) + c
        value |= (v & 1) << j
        c = v // 2
    return value, c


def ceil_log2(n):
    return (n - 1).bit_length()


def main():
    print('nextforest_20260926_automata: exact controls')
    total = 0
    counts = []
    for length in range(13):
        states = set()
        words = set()
        for n in range(1 << length):
            state, word = section_scan(n, length)
            require(state == direct_section(n, length), ('section', n, length))
            a, b = state
            require(0 <= b < 3 ** a, ('bounds', state))
            require((a, b) == (0, 0) or (a >= 1 and b % 3 != 0), ('units', state))
            for q in (0, 1, 2, 7):
                aa, endpoint = direct_section(n + (1 << length) * q, length)
                require(aa == a and endpoint == 3 ** a * q + b, ('cylinder', n, length, q))
            states.add(state)
            words.add(word)
            total += 1
        require(len(words) == 1 << length, ('parity permutation', length))
        counts.append(len(states))
    print('source-prefix controls:', total, 'at every length 0..12; four tail values each')
    print('distinct residual states at exact depth 0..12:', counts)

    reached = 0
    for a in range(1, 8):
        m = 3 ** a
        state, _ = section_scan((1 << a) - 1, a)
        require(state == (a, m - 1), ('ones', a))
        seen = set()
        for _ in range(2 * 3 ** (a - 1)):
            require(state[0] == a and state not in seen, ('cycle repeat', a, state))
            seen.add(state)
            state, p = transition(state, state[1] & 1)
            require(p == 0, ('horizontal output', a))
        require(state == (a, m - 1), ('horizontal closure', a))
        require(seen == {(a, b) for b in range(m) if b % 3}, ('level saturation', a))
        for state in seen:
            nxt, p = transition(state, 1 - (state[1] & 1))
            require(p == 1 and nxt[0] == a + 1 and nxt[1] % 3 == 2, ('lift', state))
        reached += len(seen)
    require(reached + 1 == 3 ** 7, 'finite ternary bank count')
    print('all 2186 non-root states through a=7 reached; each horizontal cycle is complete')

    machine_cases = reset_cases = 0
    for a in range(7):
        m = 3 ** a
        h = ceil_log2(m)
        # Every possible initial carry, including carries not in the parity encoder.
        for b in range(m):
            states = {b}
            todo = deque([b])
            while todo:
                c = todo.popleft()
                for bit in (0, 1):
                    nxt = (m * bit + c) // 2
                    if nxt not in states:
                        states.add(nxt)
                        todo.append(nxt)
            require(states == set(range(m)), ('macro reachability', a, b))
            for u in (0, 1, (1 << (h + 2)) - 1, 27):
                width = max(h + 2, u.bit_length())
                out, c = bits_output(m, b, u, width)
                require(out + (1 << width) * c == m * u + b, ('macro value', a, b, u))
            machine_cases += 1
        # Zero continuation separates every state: its h-bit outputs are c itself.
        require({bits_output(m, c, 0, h)[0] for c in range(m)} == set(range(m)), ('minimality', a))
        for width in range(max(0, h - 1), h + 3):
            B = 1 << width
            actual = 0
            for u in range(B):
                low = (m * u) // B
                high = (m * u + m - 1) // B
                actual += low == high
            expected = max(0, B - m + 1)
            require(actual == expected, ('reset census', a, width, actual, expected))
            require(all(bits_output(m, c, 0, h)[1] == 0 for c in range(m)), ('zero flush', a))
            reset_cases += 1
    print('macro minimal-state controls:', machine_cases, 'initial machines for a=0..6')
    print('reset-word censuses:', reset_cases, '; formula max(0, 2^h-3^a+1) exact')

    contraction_steps = 0
    for n in range(1, 256):
        H = n.bit_length()
        state, _ = section_scan(n, H)
        b_actual = direct_section(n, H)[1]
        for _ in range(64):
            a, b = state
            beta = Fraction(b, 3 ** a)
            nxt, p = transition(state, 0)
            aa, bb = nxt
            beta_next = Fraction(bb, 3 ** aa)
            require(bb == shortcut(b_actual), ('fixed source', n, state))
            require(beta_next == beta / 2 + Fraction(p, 2 * 3 ** (a + p)), ('normalized', n))
            require(0 < beta_next <= 2 * beta / 3, ('contraction', n, state))
            require(Fraction(beta).denominator == 3 ** a, ('reduced denominator', n, state))
            b_actual = bb
            state = nxt
            contraction_steps += 1
    print('fixed zero-input tails:', contraction_steps, 'exact contractions; original n=1..255')

    # Same capacity explosion, opposite actual boundary behavior.
    for k in range(1, 17):
        root_state, _ = section_scan(1, 2 * k)
        require(root_state == (k, 1), ('terminal bank control', k))
        rise_state, _ = section_scan((1 << k) - 1, k)
        require(rise_state == (k, 3 ** k - 1), ('rising control', k))
    print('controls: n=1 at depth 2k has (a,b)=(k,1); n=2^k-1 at depth k has (k,3^k-1), k=1..16')
    print('PASS: all checks are explicit exceptions, active under python -O')


if __name__ == '__main__':
    main()
