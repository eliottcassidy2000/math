#!/usr/bin/env python3
"""Exact unconditioned affine pair-chain survival; meet-in-the-middle.

Python 3, stdlib. No convergence or asymptotic-rate assumption. State (k,e)
means v=3**k*u+e, after swapping to k>=0 (and e>=0 when k=0).
Absorbed (0,0) has TWO self edges; edge multiplicity is never discarded.
"""
import argparse
import json
from collections import Counter, defaultdict
from time import monotonic


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def children(s):
    k, e = s
    if e % 2 == 0:
        out = ((k, e // 2), (k, (3 * e + 1 - 3**k) // 2))
    elif k:
        out = ((k + 1, (3 * e + 1) // 2),
               (k - 1, (e - 3**(k - 1)) // 2))
    else:
        out = ((1, (3 * e + 1) // 2), (1, -(3 * e - 1) // 2))
    return tuple((K, abs(E) if K == 0 else E) for K, E in out)


def parents(s):
    k, e = s
    if k:
        out = [(k, 2 * e), (k + 1, 2 * e + 3**k)]
        if (2 * e - 1 + 3**k) % 3 == 0:
            out.append((k, (2 * e - 1 + 3**k) // 3))
        if k > 1:
            if (2 * e - 1) % 3 == 0:
                out.append((k - 1, (2 * e - 1) // 3))
        else:
            for E in ((2 * e - 1) // 3, (1 - 2 * e) // 3):
                if E > 0 and E % 2 and (
                    e == (3 * E + 1) // 2 or e == -(3 * E - 1) // 2
                ):
                    out.append((0, E))
    else:
        out = [(0, 2 * e), (1, 2 * e + 1)]
        if e % 3 == 0:
            out.append((0, 2 * e // 3))
        if e:
            out.append((1, -2 * e + 1))
    return out


def step(D, adjacency):
    out = defaultdict(int)
    for s, weight in D.items():
        for t in adjacency(s):
            out[t] += weight
    return dict(out)


def audit():
    checks = 0
    # Independent edge inversion, with multiplicities and sign normalisation.
    for k in range(8):
        for e in range(0 if k == 0 else -100, 101):
            s = (k, e)
            for t, mult in Counter(children(s)).items():
                require(Counter(parents(t))[s] == mult, "forward/back edge")
                checks += 1
            for t, mult in Counter(parents(s)).items():
                require(Counter(children(t))[s] == mult, "back/forward edge")
                checks += 1
    F = {(0, 1): 1}
    forwards = [F]
    for t in range(1, 17):
        F = step(F, children)
        require(sum(F.values()) == 2**t, "probability mass")
        forwards.append(F)
    B = {(0, 0): 1}
    for r in range(13):
        for t in range(17 - r):
            a = sum(w * B.get(s, 0) for s, w in forwards[t].items())
            require(a == forwards[t+r].get((0, 0), 0), "split independence")
            checks += 1
        B = step(B, parents)
    # Independent actual integer orbits for every residue at depth 16.
    # Equal odd counts exclude isolated coincidences that are not affine merges.
    K = 16
    direct = [0] * (K + 1)
    for r in range(2**K):
        u, v = 2**(K + 40) + r, 2**(K + 40) + r + 1
        au = av = 0
        for j in range(1, K + 1):
            bu, bv = u % 2, v % 2
            au += bu
            av += bv
            u = (3*u+1)//2 if bu else u//2
            v = (3*v+1)//2 if bv else v//2
            if u == v and au == av:
                direct[j] += 1
    for j in range(K + 1):
        require(direct[j] == forwards[j].get((0, 0), 0) * 2**(K-j),
                "independent integer cylinder control")
        checks += 1
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--back', type=int, default=17)
    ap.add_argument('--forward', type=int, default=31)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    checks = audit()
    start = monotonic()
    B = {(0, 0): 1}
    rows = {0: {'depth': 0, 'unmerged_numerator': 1, 'denominator': 1}}
    back_sizes = [1]
    for r in range(1, args.back + 1):
        B = step(B, parents)
        back_sizes.append(len(B))
        rows[r] = {'depth': r, 'unmerged_numerator': 2**r - B.get((0, 1), 0),
                   'denominator': 2**r}
        print('back', r, len(B), round(monotonic()-start, 2), flush=True)
    F = {(0, 1): 1}
    forward_sizes = [1]
    for t in range(args.forward + 1):
        if t:
            F = step(F, children)
            forward_sizes.append(len(F))
        require(sum(F.values()) == 2**t, 'forward mass')
        direct = 2**t - F.get((0, 0), 0)
        if t in rows:
            require(direct == rows[t]['unmerged_numerator'], 'direct/split agreement')
            checks += 1
        absorbed = sum(w * B.get(s, 0) for s, w in F.items())
        K = t + args.back
        rows[K] = {'depth': K, 'unmerged_numerator': 2**K - absorbed,
                   'denominator': 2**K}
        print('forward', t, len(F), 'depth', K, 'live', 2**K - absorbed,
              round(monotonic()-start, 2), flush=True)
    result = []
    for K in sorted(rows):
        row = rows[K]
        if K:
            row['first_merge_numerator_at_this_depth'] = (
                2*rows[K-1]['unmerged_numerator'] - row['unmerged_numerator'])
            require(row['first_merge_numerator_at_this_depth'] >= 0, 'monotone survival')
        result.append(row)
    require(rows[20]['unmerged_numerator'] == 626933, 'inherited exact control')
    data = {'scope': 'unconditioned n,n+1 affine-cylinder pair chain; no rate theorem',
            'back_depth': args.back, 'forward_depth': args.forward,
            'audit_checks': checks, 'independent_integer_residues': 2**16,
            'back_state_counts': back_sizes, 'forward_state_counts': forward_sizes,
            'rows': result}
    with open(args.output, 'w') as f:
        json.dump(data, f, indent=2)
        f.write('\n')


if __name__ == '__main__':
    main()
