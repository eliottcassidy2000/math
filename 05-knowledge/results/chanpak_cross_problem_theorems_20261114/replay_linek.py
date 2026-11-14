from itertools import combinations, permutations

def rooted_orbit(m):
    """Physical closure modulo odd m, using invertible rerooting."""
    points = [(1 % m, 3 % m)]
    seen = set(points)
    def add(x, y):
        pair = x % m, y % m
        if pair not in seen:
            seen.add(pair)
            points.append(pair)
    for a, b in points:
        add(b, b + 2*a)
        for c, d in points:
            add(a*c, b*d)
    return seen

for m in (3, 5, 7, 9, 11, 15, 21, 25):
    orbit = rooted_orbit(m)
    assert {b for a, b in orbit} == set(range(m))
assert rooted_orbit(3) == {(0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(2,2)}
for power in range(3, 13):
    m = 2**power
    assert {(pow(3,d,m)+2) % m for d in range(m//4)} == {
        j for j in range(m) if j%8 in (3,5)
    }

def strong(a):
    for u in range(len(a)):
        todo = [u]
        seen = {u}
        for v in todo:
            for w, edge in enumerate(a[v]):
                if edge and w not in seen:
                    todo.append(w)
                    seen.add(w)
        if len(seen) != len(a):
            return False
    return True

def odd_cycles(a):
    n = len(a)
    for j in range(3, n+1, 2):
        for first, *rest in combinations(range(n), j):
            for suffix in permutations(rest):
                cycle = (first,*suffix,first)
                if all(a[v][w] for v,w in zip(cycle,cycle[1:])):
                    yield frozenset(cycle)

checked = 0
for m in range(3,6):
    edges = list(combinations(range(m),2))
    for mask in range(1 << len(edges)):
        a = [[0]*m for _ in range(m)]
        for j,(v,w) in enumerate(edges):
            a[v][w] = mask>>j&1
            a[w][v] = 1-a[v][w]
        if not strong(a):
            continue
        checked += 1
        cycles = list(odd_cycles(a))
        if m <= 4:
            assert len(cycles) == m-2
        else:
            assert any(all(x&y for x,y in combinations(four,2))
                       for four in combinations(cycles,4))
for m in range(3,10):
    a = [[int(w==v+1 or w<=v-2) for w in range(m)] for v in range(m)]
    cycles = list(odd_cycles(a))
    expected = {frozenset(range(i,i+j)) for j in range(3,m+1,2)
                for i in range(m-j+1)}
    assert set(cycles) == expected and len(cycles)==len(expected)

print(f"weighted residue controls and {checked} small strong hosts passed")
