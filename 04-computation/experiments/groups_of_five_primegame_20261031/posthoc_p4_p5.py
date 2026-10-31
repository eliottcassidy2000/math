#!/usr/bin/env python3
"""mac-mini-2026-10-31-fivegroups: POST-HOC (not pre-registered) anatomy of the two failed predictions P4 and P5 of
correspondence_tests.py. Nothing here changes those verdicts; it explains them.

P4: the event "some fraction is (odd-point product of a T5 cycle)/23" holds for exactly 864 / 17280 = 1/20 of
    relabellings. Only the three fractions (register * state)/state, namely F = 77/29, G = 95/23, H = 77/19, can
    become (5*19)/23 or (7*11)/23; each does so with probability 2 * (1/6)(1/5)(1/4) = 1/60, and the three events
    are disjoint (each needs a different state to map to 23).  The full observed coincidence (the fraction over 23
    carries one odd-point product and another numerator carries the other) is computed here too.
P5: lists the matched calibration triples and the chance that at least one of the four targets tried is matched.
"""
import itertools
from fractions import Fraction

PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
REG, STA = PR[:4], PR[4:]
HEAD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11)]
OWNER = HEAD + [(15, 2), (1, 7), (55, 1)]


def relabel(prog, mr, ms):
    m = dict(zip(REG, mr))
    m.update(dict(zip(STA, ms)))

    def img(x):
        y = 1
        for p in PR:
            while x % p == 0:
                x //= p
                y *= m[p]
        return y
    return [(img(a), img(b)) for a, b in prog]


def main():
    print('P4 anatomy (post-hoc)')
    tot = 0
    by_src = {5: 0, 6: 0, 7: 0}  # indices of F, G, H in the list
    joint = 0
    for mr in itertools.permutations(REG):
        for ms in itertools.permutations(STA):
            pr = relabel(OWNER, mr, ms)
            tot += 1
            hit = [j for j, (a, b) in enumerate(pr) if b == 23 and a in (95, 77)]
            for j in hit:
                by_src[j] = by_src.get(j, 0) + 1
            if hit:
                a23 = pr[hit[0]][0]
                other = 77 if a23 == 95 else 95
                if any(a == other for k, (a, b) in enumerate(pr) if k != hit[0]):
                    joint += 1
    print(f'  relabellings {tot}; source fraction of the hit: '
          f'{ {OWNER[j]: c for j, c in by_src.items()} }; total {sum(by_src.values())} = {Fraction(sum(by_src.values()), tot)}')
    print(f'  full coincidence (fraction over 23 carries one odd-point product, another numerator the other): {joint} '
          f'= {Fraction(joint, tot)} of relabellings')
    print('  reading: 1/60 per (register*state)/state fraction, three disjoint chances, so exactly 1/20; the full '
          'coincidence needs two more label matches.')
    print()
    print('P5 anatomy (post-hoc)')
    rate = Fraction(13, 136)
    p_any = 1 - (1 - rate) ** 4
    print(f'  calibration rate per target 13/136 = {float(rate):.4f}; chance that at least one of the four targets '
          f'(F, S, U, tau) is matched under independence: {float(p_any):.3f}')
    print('  the additive split F = S + U has no counterpart: S = (0,1,7) and U = (0,1,3) are matched by no invariant')


if __name__ == '__main__':
    main()
