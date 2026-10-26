#!/usr/bin/env python3
"""Double parity protection in A2, A5, A10 (S56).

A2/A5: the fuel w keeps its parity under every rule except the end cases; HALT needs w = 0 at an a = 0, even-b
state, so it needs an even-w episode, which only the o_small end case (a = 0, b odd, w <= 1) creates.
A10: HALT needs c = 3 at b = 0; c is odd only in episodes created by b0_small (b = 0, c <= 2).
We count cascade ends (anchor steps), dangerous episodes, the scale at their start, the final halting margin of
each episode (A2/A5: w at the last a=0 even-b state; A10: c - 3 at the last b=0 state), and compare with the
heuristic rates  P(episode per cascade) ~ alpha/L,  P(halt | episode) ~ beta/L.
Usage: python parity_episodes_20261025.py NAME macro_steps
"""
import sys, os, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rs_accel_20261025 as R
from deepmind_cryptids_20261025 import step, START as DSTART


def scale(nm, p):
    if nm in ('A2', 'A5'):
        a, b, w = p
        return 2 * a + b + 2 * w
    a, b, c = p
    return a + 2 * b + c


def dangerous(nm, p):
    if nm in ('A2', 'A5'):
        return p[2] % 2 == 0
    return p[2] % 2 == 1


def check(nm, p):
    """Halting-check state and its margin (0 = halts)."""
    if nm in ('A2', 'A5'):
        a, b, w = p
        if a == 0 and b % 2 == 0:
            return w
        return None
    a, b, c = p
    if b == 0 and c >= 1 and c % 2 == 1:
        return c - 3
    return None


if __name__ == '__main__':
    nm, N = sys.argv[1], int(float(sys.argv[2]))
    ls = R.LSTEP[nm]
    p = DSTART[nm]
    t0 = time.time()
    n_anchor = 0
    episodes = []
    cur = None
    halted = False
    anchor = R.ANCHOR[nm]
    for i in range(N):
        q, lab = ls(p)
        if q == 'HALT':
            halted = True
            break
        if anchor(q):
            n_anchor += 1
        dq = dangerous(nm, q)
        if dq and cur is None:
            cur = {'start_step': i, 'L': scale(nm, q), 'label': str(lab), 'last_margin': None, 'min_margin': None,
                   'checks': 0}
        if cur is not None:
            m = check(nm, q)
            if m is not None:
                cur['checks'] += 1
                cur['last_margin'] = m
                if cur['min_margin'] is None or abs(m) < abs(cur['min_margin']):
                    cur['min_margin'] = m
            if not dq:
                cur['len'] = i - cur['start_step']
                episodes.append(cur)
                cur = None
        p = q
    L_end = scale(nm, p)
    # heuristic: episode creation rate per anchor step vs 1/L, and halting-margin smallness
    tot_inv = sum(1.0 / e['L'] for e in episodes)
    small = sorted((abs(e['min_margin']), e['L']) for e in episodes if e['min_margin'] is not None)[:10]
    byL = {}
    for e in episodes:
        k = int(math.log2(max(2, e['L'])))
        byL[k] = byL.get(k, 0) + 1
    print(f'{nm}: {N:.2e} macro steps ({time.time() - t0:.0f}s), halted {halted}, final scale L = {L_end} '
          f'(2^{math.log2(max(2, L_end)):.1f}), anchor steps {n_anchor}, dangerous episodes {len(episodes)}', flush=True)
    print('  episodes by log2(L) at start:', sorted(byL.items()))
    print('  smallest |halting margin| seen in episodes (margin, L):', small)
    print('  sum over episodes of 1/L (expected halts if P(halt|episode) = c/L, per unit c):', f'{tot_inv:.4e}')
    lab_counts = {}
    for e in episodes:
        lab_counts[e['label']] = lab_counts.get(e['label'], 0) + 1
    print('  episode-creating labels:', lab_counts)
    json.dump({'machine': nm, 'steps': N, 'halted': halted, 'final_scale': L_end, 'anchor_steps': n_anchor,
               'episodes': len(episodes), 'by_log2L': sorted(byL.items()), 'smallest_margins': small,
               'sum_inv_L': tot_inv, 'labels': lab_counts,
               'episode_sample': episodes[:200]},
              open(os.path.join(HERE, f'parity_episodes_{nm}.json'), 'w'), indent=1, default=str)
