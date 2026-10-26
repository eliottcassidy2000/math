#!/usr/bin/env python3
"""A4 = 1RB0LF_1LC1RA_0RE0RD_---1LE_1LF1RC_1LC1LA: exact reduction to a one-variable map with an Antihydra core (S56).

Rules (Dos Santos et al., A.4), C(a,b) = 0^inf E> 10^a 1^b 0^inf:
  C(a,b) -> C(2a+2, b-2a-2) if b > 2a+2;  C(a,1) -> C(2,2a+2);  C(a,2) -> Halt;
  C(a,2k) -> C((6a+2k+2)//4, (6a+2k+2)%4+1) (a>=1, 2<=k<=a+1);  C(a,2k+1) -> C((6a+2k+6)//4, (6a+2k+6)%4+1) (a>=1, 1<=k<=a).
Claims checked here.
 (1) Every mod rule yields b' in {1,3} (6a+2k+2 and 6a+2k+6 are even).
 (2) From C(a,3) (a>=1): C(a,3) -> C(f(a), 1 if a even else 3) with f(n) = floor(3n/2) + 2, Antihydra's map.
 (3) Visit map on C(a,1): B = 2a+2; cascade a_i = 2^(i+2)-2 subtracting s_i = 2^(i+3)-2 while the remainder exceeds
     2a_i+2; with final (m, r): r = 1 -> next visit a_m; r = 2 -> HALT; otherwise x = 6a_m + r + (2 if r even else 5),
     a' = x//4, and b' = 1 (visit a') or b' = 3 (Hydra run: repeat old=a; a=f(a) until old was even; then visit a).
 (4) Halting criterion: A4 halts iff some visit has a = 2^(m+2) - m - 4 (m >= 0).
Statistics: the visit sequence, growth per visit, the fraction of visits entering a Hydra run, run-length law, the
distance of each visit to the halting set, and the heuristic residual halting probability.
Usage: python a4_hydra_reduction_20261025.py [n_visits]
"""
import sys, json, os, math
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, HALT


def f(n):
    return 3 * n // 2 + 2


def visit_map(a):
    """Return ('HALT', None) or (next_visit_a, info)."""
    B = 2 * a + 2
    # the rule system applies the cascade from C(2, B)
    aa, b = 2, B
    m = 0
    while b > 2 * aa + 2:
        b -= 2 * aa + 2
        aa = 2 * aa + 2
        m += 1
    r = b
    if r == 1:
        return aa, {'m': m, 'r': r, 'run': 0}
    if r == 2:
        return HALT, {'m': m, 'r': r}
    x = 6 * aa + r + (2 if r % 2 == 0 else 5)
    a2, bp = x // 4, x % 4 + 1
    assert bp in (1, 3)
    if bp == 1:
        return a2, {'m': m, 'r': r, 'run': 0}
    run = 0
    while True:
        old = a2
        a2 = f(a2)
        run += 1
        if old % 2 == 0:
            break
    return a2, {'m': m, 'r': r, 'run': run}


def halting_set_member(a):
    # a = 2^(m+2) - m - 4 ?
    m = 0
    while (1 << (m + 2)) - m - 4 <= a:
        if (1 << (m + 2)) - m - 4 == a:
            return True
        m += 1
    return False


if __name__ == '__main__':
    nv = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    # (A) cross-check against the rule system: list of C(a,1) visits from the start C(0,3)
    p = (0, 3)
    visits_rule = []
    for _ in range(200000):
        if p == HALT:
            break
        if p[1] == 1:
            visits_rule.append(p[0])
            if len(visits_rule) >= 60:
                break
        p = step('A4', p, jump=True)
    a = visits_rule[0]
    visits_map = [a]
    for _ in range(len(visits_rule) - 1):
        a, info = visit_map(a)
        visits_map.append(a)
    ok = visits_map == visits_rule
    print(f'visit map reproduces the rule system on the first {len(visits_rule)} C(a,1) visits: {ok}')
    # also check claim (4) on small a by brute force through the rule system
    bad = []
    for a0 in range(0, 3000):
        q = (a0, 1)
        q = step('A4', q, jump=True)  # -> C(2, 2a0+2)
        # follow until next visit or halt
        for _ in range(100000):
            if q == HALT or q[1] == 1:
                break
            q = step('A4', q, jump=True)
        halts = (q == HALT)
        if halts != halting_set_member(a0):
            bad.append(a0)
    print(f'halting criterion "next step from C(a,1) halts iff a = 2^(m+2)-m-4" exact for a < 3000: {not bad} {bad[:5]}')
    # (B) long statistics
    a = visits_rule[0]
    runs = Counter()
    bits = []
    dmin = []
    hydra_frac = 0
    halted = False
    for i in range(nv):
        nxt, info = visit_map(a)
        if nxt == HALT:
            halted = True
            break
        runs[info['run']] += 1
        hydra_frac += info['run'] > 0
        a = nxt
        bits.append(a.bit_length())
    L = len(bits)
    print(f'{L} visits, halted: {halted}; final a has {bits[-1]} bits; mean bits gained per visit '
          f'{(bits[-1] - bits[0]) / max(1, L - 1):.4f}')
    print(f'fraction of visits entering a Hydra run: {hydra_frac / L:.4f}; run-length histogram (first 10): '
          f'{sorted(runs.items())[:10]}')
    # heuristic residual halting probability: one target per dyadic scale; P(hit at a visit with value ~2^k) ~ c 2^-k
    resid = sum(2.0 ** (-b) for b in bits)
    print(f'heuristic residual halting probability after the first visit (sum 2^-bits over visits): {resid:.3e}')
    json.dump({'reproduces_rule_system': ok, 'criterion_exact_below_3000': not bad, 'visits': L, 'halted': halted,
               'final_bits': bits[-1], 'bits_per_visit': (bits[-1] - bits[0]) / max(1, L - 1),
               'hydra_fraction': hydra_frac / L, 'run_hist': sorted(runs.items())[:20], 'residual': resid},
              open(os.path.join(HERE, 'a4_hydra_reduction.json'), 'w'), indent=1)
