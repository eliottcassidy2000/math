#!/usr/bin/env python3
"""Independent full-permutation witness audit for the order-nine depth hostile."""
import argparse
import itertools
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail', action='store_true')
    args = parser.parse_args()
    need(not args.fail, 'intentional failure control')
    n, code = 9, 26485
    choices = [0]*n
    c = code
    for v in range(n, 1, -1):
        choices[v-1], c = c % v, c // v
    need(c == 0, 'valid factorial code')
    a, depths = [0]*n, [0]*n
    for v in range(1, n):
        q = choices[v]
        depths[v] = 1+depths[q-1] if q else 0
        for u in range(v):
            if q == u+1:
                a[v] |= 1 << u
            else:
                a[u] |= 1 << v
    need(a == [502,508,424,465,484,328,132,288,64], 'literal graph')
    children = [0]*(n+1)
    valid_orderings, minimum_depth, attaining_order = 0, n, None
    for p in itertools.permutations(range(n)):
        # Removing a newborn from any child Hamilton path gives one old
        # permutation. At most one bad old adjacency can be repaired, by
        # inserting the newborn precisely there. No transfer law is imported.
        bad = []
        for i, (u, v) in enumerate(zip(p, p[1:])):
            if not (a[u] >> v & 1):
                bad.append(i)
                if len(bad) == 2:
                    break
        if not bad:
            children[0] += 1  # Sink appended to this valid old path.
            # For a newborn beating only r, insert immediately before r.
            # It may also be appended unless r was the old last vertex.
            for r in range(n):
                children[r+1] += 1 + (r != p[-1])
        elif len(bad) == 1:
            children[p[bad[0]+1]+1] += 1
        # Check all possible chronological labelings of this very graph.
        older_mask, d, current_depth = 0, [0]*n, 0
        for v in p:
            older_out = a[v]&older_mask
            if older_out.bit_count() > 1:
                break
            d[v] = 1+d[older_out.bit_length()-1] if older_out else 0
            current_depth = max(current_depth, d[v])
            older_mask |= 1 << v
        else:
            valid_orderings += 1
            if current_depth < minimum_depth:
                minimum_depth, attaining_order = current_depth, p
    need(children == [217,2041,2241,653,769,1197,539,437,501,373], 'all literal child counts')
    need(valid_orderings == 16 and minimum_depth == 3, 'all-labeling depth obstruction')
    print(json.dumps({'status': 'PASS', 'n': n, 'code': code, 'outgoing_masks': a,
        'choices_v1_through_v8': choices[1:], 'history_depth': max(depths),
        'H': children[0], 'all_one_birth_child_counts_by_choice': children,
        'valid_chronological_orderings': valid_orderings,
        'minimum_depth_over_all_vertex_orderings': minimum_depth,
        'attaining_order': attaining_order,
        'method': 'all 9! old permutations, exact newborn insertion, all graph vertex orderings; no producer imports',
        'scope': 'Literal graph and all children; does not by itself exclude different shallow graphs sharing H or the child signature.'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
