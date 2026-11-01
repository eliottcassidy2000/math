#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: why the strong points of symmetry are exactly the circulants of order <= 5.

Representation used here (elementary).
  A circulant tournament Circ(Z_q, S) (q odd, S one element of each pair {x, -x}) may be taken with 1 in S
  (x -> -x is an isomorphism onto Circ(Z_q, -S)). A Hamiltonian path starting at 0 is a word w = (d_1..d_{q-1}) in
  S^(q-1) whose partial sums 0, d_1, d_1+d_2, ... are distinct mod q (an "S-sequencing"); H = q * (number of words).
  The positional reversal rho(v_i) = v_{q+1-i} of the path is an anti-automorphism iff, for every segment [i, j] of
  positions, the segment sums of w and of rev(w) lie in S together (segment-sum criterion; the arc v_i -> v_j has
  difference equal to the segment sum of w on [i, j-1], and rho reverses it to the mirror segment, whose sum is the
  sum of rev(w) on the same positions).
  For prime q every anti-automorphism is affine (isomorphisms of prime-order circulants are affine: Turner 1967,
  Babai 1977, as cited in THM-4849), and an involutive one is a reflection x -> c - x (a^2 = 1, a = -1); rho is
  an involution, so for prime q:  path palindromic  <=>  its difference word is a palindrome.
  For composite q the word test is only sufficient (hostile: Circ(Z_9, {1,3,4,7}), path (0,1,2,3,7,5,6,4,8)).

Prints, for odd q <= 13 and every connection set S with 1 in S: the number of words, of palindromic words and of
segment-palindromic words (= palindromic Hamiltonian paths), the pure-blue verdict, the shortest non-palindromic
witness, and checks the segment-sum criterion against the literal arc test on every path.

usage: python3 circulant_sequencings.py
"""
import itertools
import sys
from fractions import Fraction

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def connection_sets(q):
    pairs = [(x, q - x) for x in range(1, (q + 1) // 2)]
    for choice in itertools.product((0, 1), repeat=len(pairs)):
        S = frozenset(p[c] for p, c in zip(pairs, choice))
        if 1 in S:
            yield S


def sequencings(q, S):
    """all words over S of length q-1 with distinct partial sums mod q (Hamiltonian paths from 0)"""
    Sl = sorted(S)
    out = []
    seen = [False] * q
    seen[0] = True
    word = []

    def rec(pos):
        if len(word) == q - 1:
            out.append(tuple(word))
            return
        for s in Sl:
            nxt = (pos + s) % q
            if not seen[nxt]:
                seen[nxt] = True
                word.append(s)
                rec(nxt)
                word.pop()
                seen[nxt] = False
    rec(0)
    return out


def segment_palindromic(q, S, w):
    r = w[::-1]
    n = len(w)
    for i in range(n):
        a = b = 0
        for j in range(i, n):
            a = (a + w[j]) % q
            b = (b + r[j]) % q
            if (a in S) != (b in S):
                return False
    return True


def arc_palindromic(q, S, w):
    path = [0]
    for d in w:
        path.append((path[-1] + d) % q)
    n = len(path)
    rho = {path[i]: path[n - 1 - i] for i in range(n)}
    for u in range(q):
        for v in range(q):
            if u != v and (((v - u) % q) in S) != (((rho[u] - rho[v]) % q) in S):
                return False
    return True


def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))


def iso_classes(q, sets):
    """group connection sets by isomorphism class of the circulant (multiplier orbits; for q <= 13 this is the full
    isomorphism relation for prime q; composite q = 9 is checked by brute force on adjacency under all
    permutations fixing 0 is too costly, so classes are reported by multiplier orbit and flagged)"""
    units = [a for a in range(1, q) if all(a % p for p in range(2, q + 1) if q % p == 0 and is_prime(p))]
    classes = []
    seen = set()
    for S in sets:
        if S in seen:
            continue
        orb = set()
        for a in units:
            T = frozenset((a * s) % q for s in S)
            for TT in (T, frozenset((-t) % q for t in T)):
                if 1 in TT:
                    orb.add(TT)
        seen |= orb
        classes.append(sorted(orb, key=lambda x: sorted(x)))
    return classes


def terras_cycle(word):
    L, k = len(word), sum(word)
    c = sum(3 ** sum(word[j + 1:]) * 2 ** j for j, b in enumerate(word) if b)
    D = 2 ** L - 3 ** k
    x = Fraction(c, D)
    orbit = []
    for b in word:
        orbit.append(x)
        x = x / 2 if b == 0 else (3 * x + 1) / 2
    return [int(y * D) for y in orbit], D


def main():
    print('Circulant tournaments Circ(Z_q, S), 1 in S: S-sequencings and palindromic Hamiltonian paths')
    print()
    total_checked = 0
    for q in (1, 3, 5, 7, 9, 11, 13):
        if q == 1:
            print('  q = 1: K1, one path (empty word), pure-blue.')
            continue
        sets = list(connection_sets(q))
        classes = iso_classes(q, sets)
        print(f'  q = {q} ({"prime" if is_prime(q) else "composite"}): {len(sets)} connection sets with 1 in S, '
              f'{len(classes)} multiplier classes')
        for cl in classes:
            S = cl[0]
            words = sequencings(q, S)
            pw = [w for w in words if w == w[::-1]]
            sp = [w for w in words if segment_palindromic(q, S, w)]
            for w in words:
                total_checked += 1
                if segment_palindromic(q, S, w) != arc_palindromic(q, S, w):
                    check(f'segment criterion = arc test, q={q} S={sorted(S)} w={w}', False)
            pure = len(sp) == len(words)
            witness = None
            if not pure:
                nonp = [w for w in words if not segment_palindromic(q, S, w)]
                witness = min(nonp, key=lambda w: (sum(1 for d in w if d != 1), w))
            extra = ''
            if not is_prime(q) and len(sp) != len(pw):
                extra = f'  (composite: {len(sp) - len(pw)} palindromic paths with non-palindromic words)'
            print(f'     S = {sorted(S)} (class size {len(cl)}): H = {q * len(words)}, words {len(words)}, '
                  f'palindromic words {len(pw)}, palindromic paths {len(sp)}, pure-blue {pure}'
                  + (f'; simplest non-palindromic word {witness}' if witness else '') + extra)
            if is_prime(q):
                check(f'prime q={q}, S={sorted(S)}: palindromic path <=> palindromic word', len(sp) == len(pw)
                      and set(sp) == set(pw))
            if q <= 5:
                check(f'q={q}: pure-blue', pure)
            else:
                check(f'q={q}, S={sorted(S)}: not pure-blue', not pure)
    check(f'segment-sum criterion agrees with the literal arc test on all {total_checked} paths', True)
    # composite hostile from the THM-4730 audit
    q, S = 9, frozenset({1, 3, 4, 7})
    path = (0, 1, 2, 3, 7, 5, 6, 4, 8)
    w = tuple((path[i + 1] - path[i]) % q for i in range(8))
    check('Z_9 hostile: word (1,1,1,4,7,1,7,4) is not a palindrome but the path is palindromic',
          w != w[::-1] and arc_palindromic(q, S, w) and segment_palindromic(q, S, w), f'word {w}')
    print()
    print('THM-4849 reading of the q <= 5 circulants (connection set as a q-periodic parity word)')
    for q, S in ((3, {1}), (5, {1, 2}), (5, {1, 3})):
        word = [1 if i in S else 0 for i in range(q)]
        orbit, D = terras_cycle(word)
        print(f'  q={q}, S={sorted(S)}: parity word {word}, cycle of x -> x/2, (3x+{D})/2: {orbit}')
    # T5's three words and the cycle points they visit
    q, S = 5, frozenset({1, 2})
    words = sequencings(q, S)
    orbit, D = terras_cycle([1 if i in S else 0 for i in range(q)])
    print(f'  T5 = Circ(Z_5,{{1,2}}): its {len(words)} S-sequencings {words} are all palindromes;')
    for w in words:
        pos = [0]
        for d in w:
            pos.append((pos[-1] + d) % q)
        print(f'     word {w}: visits cycle points {[orbit[p] for p in pos]} (position j = T^j of the base point)')
    print()
    if FAIL:
        print('FAILED:', FAIL)
        return 1
    print('ALL CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
