#!/usr/bin/env python3
"""mac-mini-2026-11-01-fractranTM: explicit non-palindromic Hamiltonian paths for every circulant tournament of odd
order q >= 7 (connection set S with 1 in S), checked with the segment-sum criterion.

Case rule (k = (q-1)/2; chi = indicator of S on 1..q-1; a "rise" is u in [1, q-2] with u not in S, u+1 in S):
  A  if 2 in S:  pick s in S, s >= 3;          w_A(s) = 1^(s-2) . 2 . 1^(q-1-s) . s
  B  if 2 not in S and some rise u != k:       w_B(u) = 1^(u-1) . (q-u) . (u+1) . 1^(q-2-u)
  D  otherwise (then S = {1} u [k+1, q-2]):    w_D    = 1 . ((q-3) . 1)^(k-1) . (q-2)
Each word is an S-sequencing (partial sums 0 .. q-1 distinct mod q) and is not a palindrome (first letter differs
from last, or for B the block lengths differ because u != k). For prime q this already proves non-palindromicity
of the path (anti-automorphisms of prime-order circulants are affine; involutive ones are reflections). For
composite q the script checks the full segment-sum criterion, which is the literal arc test.

PRE-REGISTERED (before the first run of this file): for every odd q with 7 <= q <= 35 and every S with 1 in S, the
case rule produces a valid S-sequencing that fails the segment-sum criterion (so the path is not palindromic), and
for case A some choice of s works with the fixed segment "the suffix after the letter 2" whenever s - 1 is not in S.

usage: python3 witness_families.py [QMAX]
"""
import itertools
import sys

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def sets_with_one(q):
    pairs = [(x, q - x) for x in range(2, (q + 1) // 2)]
    for choice in itertools.product((0, 1), repeat=len(pairs)):
        yield frozenset([1] + [p[c] for p, c in zip(pairs, choice)])


def is_sequencing(q, S, w):
    if len(w) != q - 1 or any(d not in S for d in w):
        return False
    seen = {0}
    s = 0
    for d in w:
        s = (s + d) % q
        if s in seen:
            return False
        seen.add(s)
    return True


def first_failing_segment(q, S, w):
    r = w[::-1]
    n = len(w)
    for length in range(1, n + 1):
        for i in range(0, n - length + 1):
            a = sum(w[i:i + length]) % q
            b = sum(r[i:i + length]) % q
            if (a in S) != (b in S):
                return (i + 1, i + length)
    return None


def witness(q, S):
    k = (q - 1) // 2
    if 2 in S:
        cands = sorted(s for s in S if s >= 3)
        # prefer s with s-1 not in S (the suffix segment then fails by construction), else the smallest
        pref = [s for s in cands if (s - 1) not in S and s <= (q + 1) // 2]
        s = pref[0] if pref else cands[0]
        return 'A', (1,) * (s - 2) + (2,) + (1,) * (q - 1 - s) + (s,)
    rises = [u for u in range(1, q - 1) if u not in S and (u + 1) in S]
    other = [u for u in rises if u != k]
    if other:
        u = other[0]
        return 'B', (1,) * (u - 1) + (q - u, u + 1) + (1,) * (q - 2 - u)
    assert S == frozenset([1] + list(range(k + 1, q - 1))), (q, sorted(S))
    return 'D', (1,) + (q - 3, 1) * (k - 1) + (q - 2,)


def main():
    qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 35
    totals = {}
    for q in range(7, qmax + 1, 2):
        cnt = {'A': 0, 'B': 0, 'D': 0}
        bad = 0
        for S in sets_with_one(q):
            case, w = witness(q, S)
            cnt[case] += 1
            ok_seq = is_sequencing(q, S, w)
            ok_np = w != w[::-1]
            seg = first_failing_segment(q, S, w)
            if not (ok_seq and ok_np and seg is not None):
                bad += 1
                if bad <= 3:
                    print(f'    counterexample q={q} S={sorted(S)} case {case} w={w} seq={ok_seq} '
                          f'nonpal-word={ok_np} failing segment={seg}')
        totals[q] = cnt
        prime = all(q % p for p in range(2, int(q ** 0.5) + 1))
        check(f'q={q} ({"prime" if prime else "composite"}): all {sum(cnt.values())} connection sets get a valid '
              f'non-palindromic Hamiltonian path', bad == 0, f'cases {cnt}')
    print()
    print('case D occurs exactly once per q (S = {1} u [k+1, q-2}):', all(c['D'] == 1 for c in totals.values()))
    if FAIL:
        print('FAILED:', FAIL)
        return 1
    print('ALL CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
