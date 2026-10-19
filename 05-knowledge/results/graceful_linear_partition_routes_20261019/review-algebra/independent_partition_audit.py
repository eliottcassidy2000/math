"""Independent finite audit: exhaustive support systems, routing, and actual-state BFS.

No producer imports. Finite controls only; universal proof is reviewed separately.
"""
from collections import deque
from itertools import combinations, permutations, product
import json

CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

def bound(n):
    return 0 if n == 0 else 1 if n <= 2 else 2*n-3

def assignments(classes):
    return [s for s in product(*(range(len(c)) for c in classes))
            if None not in [classes[i][s[i]] for i in range(len(classes))]
            and len({classes[i][s[i]] for i in range(len(classes))}) == len(classes)]

def independent_route(classes, source, target):
    n = len(classes)
    supports = [set(c)-{None} for c in classes]
    suffix = [set().union(*supports[i:]) for i in range(n)] + [set()]
    tight = [i for i in range(n+1) if len(suffix[i]) == n-i]
    require(tight[0] == 0 and tight[-1] == n, 'fully tight endpoints')
    current = list(source)
    path = [tuple(source)]
    total_budget = 0

    def emit(update):
        new = current[:]
        for i, part in update.items():
            new[i] = target[i] if classes[i][target[i]] == part else classes[i].index(part)
        if new != current:
            require(sum(x != y for x, y in zip(new, current)) <= 2, 'move width')
            occupied = [classes[i][new[i]] for i in range(n)]
            require(None not in occupied and len(set(occupied)) == n, 'actual basis')
            current[:] = new
            path.append(tuple(new))

    for left, right in zip(tight, tight[1:]):
        block = list(range(left, right))
        b = len(block)
        total_budget += bound(b)
        actual_part = [classes[i][target[i]] for i in block]
        part_label = {p:j for j,p in enumerate(actual_part)}
        require(set(actual_part) == suffix[left]-suffix[right], 'target block parts')
        avail = [{part_label[p] for p in supports[i] if p in part_label} for i in block]
        require(avail[0] == set(range(b)), 'full buffer')
        for size in range(1,b):
            for subset in combinations(range(1,b),size):
                require(len(set().union(*(avail[j] for j in subset))) > size, 'strict Hall')
        before = len(path)
        sigma = [part_label[classes[i][current[i]]] for i in block]
        if sigma == list(range(b)):
            wrong = [j for j in range(b) if current[left+j] != target[left+j]]
            for k in range(0,len(wrong),2):
                emit({left+j:actual_part[j] for j in wrong[k:k+2]})
        else:
            cycles = []
            unseen = set(range(b))
            while unseen:
                start = min(unseen)
                cyc = [start]
                unseen.remove(start)
                j = sigma[start]
                while j != start:
                    cyc.append(j)
                    unseen.remove(j)
                    j = sigma[j]
                cycles.append(cyc)
            which = {j:k for k,cyc in enumerate(cycles) for j in cyc}
            root = which[0]
            seen = {root}
            children = {j:[] for j in range(b)}
            # Independent closure/arborescence discovery, not producer BFS.
            while len(seen) < len(cycles):
                candidate = next(((a,j) for a in reversed(range(b)) if which[a] in seen
                                  for j in reversed(range(1,b))
                                  if which[j] not in seen and a in avail[j]), None)
                require(candidate is not None, 'cycle reachability')
                a,j = candidate
                seen.add(which[j])
                children[a].append(j)

            def swap(j):
                holder = classes[left][current[left]]
                displaced = classes[left+j][current[left+j]]
                require(part_label[holder] in avail[j], 'available temporary/target part')
                emit({left:displaced,left+j:holder})

            def visit(entry, is_root=False):
                incoming = classes[left][current[left]]
                if not is_root:
                    swap(entry)
                cursor = sigma[entry]
                while True:
                    require(classes[left][current[left]] == actual_part[cursor], 'cycle cursor')
                    for child in children[cursor]:
                        visit(child)
                        require(classes[left][current[left]] == actual_part[cursor], 'child exact return')
                    if cursor == entry:
                        if not is_root:
                            swap(entry)
                        break
                    swap(cursor)
                    cursor = sigma[cursor]
                if not is_root:
                    require(classes[left][current[left]] == incoming, 'incoming part preserved')
                    require(all(current[left+j] == target[left+j] for j in cycles[which[entry]]), 'cycle exact final elements')

            visit(0, True)
            require(len(path)-before == b+len(cycles)-2, 'exact cycle count')
        require(current[left:right] == list(target[left:right]), 'block exact target')
        require(len(path)-before <= bound(b), 'block budget')
    require(tuple(current) == tuple(target), 'global exact target')
    require(len(path)-1 <= total_budget, 'sum of block budgets')
    # Separate direct width-two move gives rank-two global exception.
    if n >= 3:
        require(total_budget <= bound(n), 'global linear budget')
    return path

def audit_bfs(states):
    adjacency = [[j for j,t in enumerate(states) if i != j and sum(a != b for a,b in zip(s,t)) <= 2]
                 for i,s in enumerate(states)]
    diameter = 0
    for start in range(len(states)):
        distances = {start:0}
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j in adjacency[i]:
                if j not in distances:
                    distances[j] = distances[i]+1
                    queue.append(j)
        require(len(distances) == len(states), 'independent state graph connected')
        diameter = max(diameter,max(distances.values(),default=0))
    return diameter

def main():
    rows = []
    route_count = 0
    moves = 0
    max_diameter = 0
    for n in range(1,5):
        options = [[frozenset(s) for size in range(n-i,n+1) for s in combinations(range(n),size)]
                   for i in range(n)]
        systems = 0
        endpoint_pairs = 0
        max_d = 0
        for supports in product(*options):
            classes = [tuple(sorted(s)) for s in supports]
            states = assignments(classes)
            require(bool(states), 'staircase existence')
            d = audit_bfs(states)
            require(d <= bound(n), 'BFS diameter theorem')
            max_d = max(max_d,d)
            for x in states:
                for y in states:
                    path = independent_route(classes,x,y)
                    route_count += 1
                    endpoint_pairs += 1
                    moves += len(path)-1
            systems += 1
        rows.append({'rank':n,'all_support_systems':systems,'all_endpoint_pairs':endpoint_pairs,'maximum_exact_diameter':max_d})
        max_diameter = max(max_diameter,max_d)
    realization_rows = []
    for n in range(1,4):
        options = [[frozenset(s) for size in range(n-i,n+1) for s in combinations(range(n),size)]
                   for i in range(n)]
        systems = endpoint_pairs = max_d = 0
        for supports in product(*options):
            # Every available part has TWO actual realizations of each color;
            # each color also has an unavailable actual loop.
            classes = [(None,)+tuple(p for p in sorted(s) for _ in range(2)) for s in supports]
            states = assignments(classes)
            d = audit_bfs(states)
            require(d <= bound(n), 'actual-element BFS diameter')
            max_d = max(max_d,d)
            for x in states:
                for y in states:
                    path = independent_route(classes,x,y)
                    route_count += 1
                    endpoint_pairs += 1
                    moves += len(path)-1
            systems += 1
        realization_rows.append({'rank':n,'all_support_systems':systems,'all_actual_endpoint_pairs':endpoint_pairs,'maximum_exact_diameter':max_d})
    print(json.dumps({'status':'PASS','checks':CHECKS,'support_universe':'All labelled support systems satisfying |A_i| >= r-i for zero-indexed i; rank at most four.','pure_support_controls':rows,'parallel_universe':'All such support systems of rank at most three, with two distinguishable elements per available color/part pair and one loop per color.','actual_realization_controls':realization_rows,'independently_constructed_routes':route_count,'validated_moves':moves,'scope':'Independent implementation and exact BFS; no producer imports. Finite controls supplement the manuscript proof.'},indent=2))

if __name__ == '__main__':
    main()
