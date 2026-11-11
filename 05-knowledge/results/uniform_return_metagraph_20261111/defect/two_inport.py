"""Positive two-inport defect automaton; exact finite-language experiments.

Letter bits say backbone vertex -> a,b; a -> b, both -> root,
root -> every backbone vertex. No inherited counting implementation.
"""
from collections import deque
from itertools import product
from pathlib import Path
import argparse
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows_for(word):
    # names r=0,a=1,b=2, then increasing backbone
    n = len(word) + 3
    rows = [sum(1 << v for v in range(3, n)), 5, 1] + [0]*len(word)
    for j, letter in enumerate(word):
        v = j+3
        rows[v] |= sum(1 << w for w in range(v+1, n))
        for port in range(2):
            if letter >> port & 1:
                rows[v] |= 1 << (port+1)
            else:
                rows[port+1] |= 1 << v
    return tuple(rows)


def path_tables(rows):
    n = len(rows)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            avail = rows[v] & ~mask
            while avail:
                bit = avail & -avail
                avail -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return [1]+[sum(x) for x in dp[1:]], dp


def actual_defect(word):
    rows = rows_for(word)
    hp, ends = path_tables(rows)
    full = (1 << len(rows))-1
    f = sum(ends[m][0]*hp[full ^ m] for m in range(1, full+1) if m & 1)
    return f-hp[full], hp[full]


# Four disjoint skeletons for q-s + delta.
# q-s: paths starting at a or b, with the two ports in either order.
# delta: first path ends at a, second starts at b, or conversely.
# Each skeleton has two backbone runs. A run is (leftport,rightport),
# with -1 an absent boundary. Its state 0=empty,1=last right-good,2=bad.
SKELETONS = (
    (((0,1),(1,-1)), True),
    (((1,0),(0,-1)), False),
    (((-1,0),(1,-1)), True),
    (((-1,1),(0,-1)), True),
)
STATES = tuple(product(range(3), repeat=2))
STATE_ID = {state:i for i,state in enumerate(STATES)}


def compile_transitions():
    table = []
    accept = []
    for runs, empty_bridge_legal in SKELETONS:
        finals = []
        for st in STATES:
            good = 2 not in st and (st[0] != 0 or empty_bridge_legal)
            finals.append(good)
        accept.append(tuple(finals))
        letters = []
        for letter in range(4):
            origins = []
            for st in STATES:
                targets = []
                for j, (left, right) in enumerate(runs):
                    # An empty run may begin only if left -> backbone.
                    if st[j] == 0 and left >= 0 and letter >> left & 1:
                        continue
                    ns = list(st)
                    ns[j] = 1 if right < 0 or letter >> right & 1 else 2
                    targets.append(STATE_ID[tuple(ns)])
                origins.append(tuple(targets))
            letters.append(tuple(origins))
        table.append(tuple(letters))
    return tuple(table), tuple(accept)


TRANS, ACCEPT = compile_transitions()
START = tuple((1,)+(0,)*8 for _ in SKELETONS)


def advance(weights, letter, cap=None):
    result = []
    for sk, old in enumerate(weights):
        new = [0]*9
        for origin, count in enumerate(old):
            for target in TRANS[sk][letter][origin]:
                new[target] += count
        if cap is not None:
            new = [min(cap, x) for x in new]
        result.append(tuple(new))
    return tuple(result)


def value(weights):
    return sum(x for sk, weights_sk in enumerate(weights)
               for x, good in zip(weights_sk, ACCEPT[sk]) if good)


def automaton_value(word):
    weights = START
    for letter in word:
        weights = advance(weights, letter)
    return value(weights)


def controls(max_length):
    queries = 0
    by_length = []
    for n in range(max_length+1):
        histogram = {}
        for word in product(range(4), repeat=n):
            expected, _ = actual_defect(word)
            found = automaton_value(word)
            require(found == expected, f"defect mismatch {word}: {found}!={expected}")
            queries += 1
            if word and word[-1] and any(x & 1 for x in word):
                histogram[found] = histogram.get(found, 0)+1
        by_length.append({"backbone_length":n, "strong_defects":histogram})
    return {"literal_DP_controls":queries, "by_length":by_length}


def build_dfa(budget, limit=1000000):
    # Keep all reachable capped states; no unjustified monotone-cost pruning.
    # flags: any a-contact, last letter nonzero. Empty word never accepts.
    first = (START, False, False)
    states, index, edges = [first], {first:0}, []
    for current in states:
        weights, seen_a, _ = current
        row = []
        for letter in range(4):
            nxt = (advance(weights, letter, budget+1), seen_a or bool(letter&1), bool(letter))
            if nxt not in index:
                require(len(states)<limit, f"state limit {limit} reached, no complete DFA claimed")
                index[nxt] = len(states)
                states.append(nxt)
            row.append(index[nxt])
        edges.append(row)
    accepted = [i for i,(w,a,last) in enumerate(states) if a and last and value(w)==budget]
    reverse = [[] for _ in states]
    for i, row in enumerate(edges):
        for j in row:
            reverse[j].append(i)
    productive = set(accepted)
    todo = deque(accepted)
    while todo:
        for i in reverse[todo.popleft()]:
            if i not in productive:
                productive.add(i)
                todo.append(i)
    accepted_set = set(accepted)
    # A complete structural certificate of the accepted language: remove
    # self loops, topologically sort, and enumerate the remaining paths.
    # This deliberately fails if a more complicated productive cycle occurs.
    loops = [(i,letter) for i in sorted(productive)
             for letter,j in enumerate(edges[i]) if i==j]
    indegree = {i:0 for i in productive}
    for i in productive:
        for j in edges[i]:
            if j in productive and i!=j:
                indegree[j] += 1
    zero = deque(i for i,v in indegree.items() if not v)
    topo = []
    while zero:
        i = zero.popleft()
        topo.append(i)
        for j in edges[i]:
            if j in productive and i!=j:
                indegree[j] -= 1
                if not indegree[j]:
                    zero.append(j)
    regex = []
    if len(topo)==len(productive):
        def walk(i, parts):
            own_loops = [letter for node,letter in loops if node==i]
            require(len(own_loops)<=1,"multiple loop letters require richer regex")
            if own_loops:
                parts = parts+[f"{own_loops[0]}*"]
            if i in accepted_set:
                regex.append(parts)
            for letter,j in enumerate(edges[i]):
                if j in productive and i!=j:
                    walk(j,parts+[str(letter)])
        if 0 in productive:
            walk(0,[])
    # Exact class counts by backbone length, not a claimed minimal DFA.
    counts = []
    distribution = {0:1}
    for n in range(41):
        counts.append(sum(count for state,count in distribution.items() if state in accepted_set))
        newer = {}
        for i, count in distribution.items():
            for j in edges[i]:
                if j in productive:
                    newer[j] = newer.get(j, 0)+count
        distribution = newer
    return {"budget":budget,"reachable_states":len(states),"productive_states":len(productive),
            "accepting_states":len(accepted),"counts_by_backbone_length_0_through_40":counts,
            "productive_edges":[[i,letter,j] for i in sorted(productive)
                                for letter,j in enumerate(edges[i]) if j in productive],
            "productive_accepting":accepted,"initial":0,
            "self_loops":loops,"acyclic_after_removing_self_loops":len(topo)==len(productive),
            "exact_disjoint_path_regex":regex}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--controls', type=int, default=4)
    parser.add_argument('--budget', type=int, nargs='*', default=[2,3,4,5,6])
    parser.add_argument('--limit', type=int, default=1000000)
    args = parser.parse_args()
    report = {"status":"FINITE-EXACT controls and complete capped automata when construction terminates",
              "controls":controls(args.controls),"automata":[]}
    for budget in args.budget:
        answer = build_dfa(budget,args.limit)
        report['automata'].append(answer)
        print(json.dumps({k:v for k,v in answer.items() if not k.startswith('productive_')}),flush=True)
    Path(__file__).with_name('two_inport.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
