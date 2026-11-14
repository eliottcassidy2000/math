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

def F(a,b,mod):
    return b%mod,(b+2*a)%mod

def period(m):
    # Return an even period of the full F action modulo m.
    a,b,c,d = 1%m,0,0,1%m
    initial = a,b,c,d
    L = 0
    while True:
        a,b=F(a,b,m)
        c,d=F(c,d,m)
        L+=1
        if (a,b,c,d)==initial and L%2==0:
            return L

def test_gluing(m,e):
    mod=m*2**e
    L=period(m)
    def extension_power(a,b,k):
        for _ in range(k):
            a,b=F(a,b,mod)
        return a,b
    cols=[extension_power(1,0,L),extension_power(0,1,L)]
    reps={}
    points=[]
    def add(a,b,n):
        label=(a%m,b%m,n%2)
        if label not in reps:
            reps[label]=(a%mod,b%mod,n%2)
            points.append(reps[label])
    def product(a,b,c,d):
        return tuple((u*a*c+v*b*d)%mod for u,v in zip(*cols))
    add(1,3,1)
    for a,b,n in points:
        add(*F(a,b,mod),n+1)
        for c,d,s in points:
            add(*product(a,b,c,d),n+s)

    A,B,_=reps[1%m,1%m,1]
    a,b=extension_power(1,1,L*((e+2+L-1)//L))
    assert (a%m,b%m)==(1%m,1%m)
    assert (a%(2**e),b%(2**e))==(pow(3,-1,2**e),-pow(3,-1,2**e)%(2**e))
    got=set()
    for x,y,n in points:
        for X,Y in [(x,y),(B*x,B*x+A*(y-x))]:
            for _ in range(2**(e-3)):
                got.add(Y%mod)
                X,Y=b*b*X,b*b*X+a*a*(Y-X)
                X%=mod
                Y%=mod
    assert got=={r for r in range(mod) if r%8 in (3,5)}

for m in (1,3,5,9,15):
    for e in (3,4,5):
        test_gluing(m,e)
print("mixed-modulus physical tree compiler controls passed")
