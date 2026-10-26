#!/usr/bin/env python3
"""Independent actual-edge audit of the graceful overwrite quotient.
No producer imports. Explicit exceptions remain active under optimization.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import factorial
import json

CHECKS = Counter()


def need(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise RuntimeError(label)


def forest(n, edges):
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a,b in edges:
        a,b = find(a),find(b)
        if a == b:
            return False
        parent[a] = b
    return True


def run_family(k, extra):
    m = k*k+2*k+extra
    leaves = [i*(k+1) for i in range(1,k+1)]
    special = {leaf:i for i,leaf in enumerate(leaves,1)}
    fixed = [(special[d],d+special[d]) if d in special else (0,d)
             for d in range(k+1,m+1)]
    need(forest(m+1,fixed), 'fixed actual forest')
    need(all(not any(b in e for e in fixed) for b in leaves), 'exact isolated leaf vertices')
    core = set(range(m+1))-set(leaves)
    need(len(fixed) == len(core)-1 and all(set(e) <= core for e in fixed),
         'remaining core is one tree')
    options = [[(a,a+d) for a in range(m+1-d)] for d in range(1,k+1)]
    def leaf_of(e):
        hits = [i for i,b in enumerate(leaves) if b in e]
        return hits[0] if len(hits) == 1 else None
    contacts = [[e for e in op if leaf_of(e) is not None] for op in options]
    need(all(Counter(leaf_of(e) for e in op) == Counter({i:2 for i in range(k)})
             for op in contacts), 'exactly two actual edges per color leaf')
    states = []
    for pi in permutations(range(k)):
        for signs in product((-1,1),repeat=k):
            state = tuple(tuple(sorted((leaves[i],leaves[i]+s*d)))
                          for d,(i,s) in enumerate(zip(pi,signs),1))
            need(forest(m+1,fixed+list(state)), 'formula state actually graceful')
            states.append((pi,state))
    need(len(set(state for pi,state in states)) == factorial(k)*2**k,
         'literal state count')
    if k <= 3:
        brute = {tuple(s) for s in product(*options) if forest(m+1,fixed+list(s))}
        need(brute == {state for pi,state in states}, 'exhaustive all actual choices match state formula')
    profiles = defaultdict(set)
    requests = 0
    for pi,state in states:
        profile = []
        for width in range(1,(2 if k <= 3 else 1)+1):
            for colors in combinations(range(k),width):
                for replacements in product(*(options[d] for d in colors)):
                    after = list(state)
                    predicted = list(pi)
                    for d,e in zip(colors,replacements):
                        after[d] = e
                        predicted[d] = leaf_of(e)
                    expected = None not in predicted and set(predicted) == set(range(k))
                    actual = forest(m+1,fixed+after)
                    need(actual == expected, 'actual overwrite legality equals assignment test')
                    if actual:
                        need(tuple(leaf_of(e) for e in after) == tuple(predicted),
                             'accepted successor assignment exact')
                    if width == 1:
                        profile.append(actual)
                    requests += 1
        profiles[tuple(profile)].add(pi)
    need(len(profiles) == factorial(k) and all(len(v) == 1 for v in profiles.values()),
         'one-color profiles reconstruct precisely pi')
    # Sign visibility changes if the observer sees whether a request changed an edge.
    pi,state = states[0]
    other = next(s for p,s in states if p == pi and s != state)
    d = next(i for i in range(k) if state[i] != other[i])
    after_other = list(other)
    after_other[d] = state[d]
    need(forest(m+1,fixed+list(state)) and forest(m+1,fixed+after_other) and
         other[d] != state[d] and tuple(after_other) != other,
         'no-op/change-indicator hostile distinguishes signs')
    return {'k':k,'m':m,'literal_states':len(states),'guard_classes':len(profiles),
            'sign_fiber':2**k,'actual_overwrite_requests':requests,
            'max_batch_width':2 if k <= 3 else 1,
            'loop_requests_included':True}


def graceful_nonplanar():
    left,right = (0,1,2),(3,6,9)
    differences = sorted(b-a for a in left for b in right)
    need(differences == list(range(1,10)), 'K33 explicit graceful labels')
    odd_counts = sorted({a*(5-a) for a in range(6)})
    need(5 not in odd_counts, 'K5 odd difference obstruction')
    normalized = []
    for middle in combinations(range(1,10),3):
        labels = (0,)+middle+(10,)
        ds = [b-a for a,b in combinations(labels,2)]
        need(sorted(ds) != list(range(1,11)), 'all normalized K5 labels rejected')
        normalized.append(labels)
    need(len(normalized) == 84, 'normalized K5 complete search cardinality')
    return {'K33_labels':[left,right],'K33_differences':differences,
            'possible_K5_odd_difference_counts':odd_counts,'normalized_K5_sets':len(normalized)}


def main():
    rows = [run_family(k,extra) for k,extra in [(1,0),(2,0),(3,0),(4,0),(1,2),(2,2),(3,2)]]
    result = {'status':'PASS','producer_imports':False,'graceful_overwrite':rows,
              'nonplanar_gracefulness':graceful_nonplanar()}
    result['checks'] = dict(sorted(CHECKS.items()))
    result['total_checks'] = sum(CHECKS.values())
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
