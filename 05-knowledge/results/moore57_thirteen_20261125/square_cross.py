"""Square-sum pentagons, the Moore cut inventory, local moments and 13 phases."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import ceil, floor, isqrt


def dihedral(word):
    return [s[i:]+s[:i] for s in (word, word[::-1]) for i in range(len(s))]


def square_window(n):
    adjacency = {x: {y for y in range(1, n+1) if x != y and isqrt(x+y)**2 == x+y}
                 for x in range(1, n+1)}
    open_paths = []
    closed_paths = []
    pentagons = set()
    for x, neighbors in adjacency.items():
        for a in neighbors:
            for b in adjacency[a] - {x}:
                for y in adjacency[b] - {x, a}:
                    if x >= y:
                        continue
                    closers = (neighbors & adjacency[y]) - {a,b}
                    path = x, a, b, y
                    open_paths.append(path)
                    if closers:
                        closed_paths.append(path)
                    for z in closers:
                        pentagons.add(min(dihedral(path+(z,))))
    words = Counter(min(dihedral(tuple(v%2 for v in cycle))) for cycle in pentagons)
    assert words == {(0, 0, 1, 1, 1): len(pentagons)}
    gaps = {d: (sum(p[3]-p[0] == d for p in open_paths),
                sum(p[3]-p[0] == d for p in closed_paths)) for d in (7,13,21)}
    return len(open_paths), len(closed_paths), sorted(pentagons), gaps


def hamilton_crossings(n):
    rows = [sum(1 << j for j in range(n) if i!=j and isqrt(i+j+2)**2 == i+j+2)
            for i in range(n)]
    full = (1 << n)-1
    paths, closed, ports = Counter(), Counter(), Counter()

    def extend(start, end, used, crossings, at_thirteen):
        if used == full:
            if start < end:
                paths[crossings] += 1
                ports[at_thirteen] += 1
                if rows[end] & (1 << start):
                    closed[crossings+((end-start)%2)] += 1
            return
        choices = rows[end] & ~used
        while choices:
            b = choices & -choices
            choices -= b
            j = b.bit_length()-1
            crossing = (end-j)%2
            extend(start, j, used|b, crossings+crossing,
                   at_thirteen+crossing*int(end==12 or j==12))
    for i in range(n):
        extend(i,i,1 << i,0,0)
    assert all(count % n == 0 for count in closed.values())
    return dict(sorted(paths.items())), {k:v//n for k,v in closed.items()}, ports


COLORS = tuple(sorted({min(dihedral(word)) for word in product(range(2), repeat=5)}))
WORDS = tuple(product(range(2), repeat=4))


def path_count(word, inside, cross, panel_size):
    a,b,c,d = word
    degree = lambda x,y: inside if x == y else cross
    return panel_size * degree(a,b) * (degree(b,c) - (a==c)) * (degree(c,d) - (b==d))


def solve(matrix, vector):
    """Exact elimination, with all residual equations checked."""
    rows = [[Fraction(z) for z in row+[rhs]] for row,rhs in zip(matrix, vector)]
    width = len(matrix[0])
    for i in range(width):
        pivot = next(j for j in range(i,len(rows)) if rows[j][i])
        rows[i], rows[pivot] = rows[pivot], rows[i]
        scale = rows[i][i]
        rows[i] = [v/scale for v in rows[i]]
        for j,row in enumerate(rows):
            if i != j and row[i]:
                rows[j] = [x-row[i]*y for x,y in zip(row,rows[i])]
    assert all(not any(row) for row in rows[width:])
    answer = [row[-1] for row in rows[:width]]
    assert all(x.denominator == 1 for x in answer)
    return tuple(map(int,answer))


def whole_census(inside, cross, panel_size):
    entries = [[sum(w[:4] == subword for w in dihedral(color)) for color in COLORS]
               for subword in WORDS]
    rhs = [path_count(w, inside, cross, panel_size) for w in WORDS]
    entries.append([0]*7+[1])
    base = solve(entries, rhs+[0])
    one = solve(entries, rhs+[1])
    slope = tuple(y-x for x,y in zip(base,one))
    assert slope == (-1,5,-5,-5,5,5,-5,1)
    return base, slope


def fixed_vertex_census(internal_pentagons):
    """Oriented cycles starting at a fixed vertex of colour zero."""
    fullwords = tuple(product(range(2), repeat=4))
    rows, rhs = [], []
    for i in range(4):
        for observed in product(range(2), repeat=3):
            # After deleting slot i+1 from the circular word, start after its gap.
            retained = (0,)+observed
            path = retained[i+1:]+retained[:i+1]
            rows.append([int(w[:i]+w[i+1:] == observed) for w in fullwords])
            rhs.append(path_count(path,32,25,1))
    return solve(rows+[[1]+[0]*15],rhs+[2*internal_pentagons])


def determinant(h):
    return (h[0][0]*(h[1][1]*h[2][2]-h[1][2]*h[2][1])
            - h[0][1]*(h[1][0]*h[2][2]-h[1][2]*h[2][0])
            + h[0][2]*(h[1][0]*h[2][1]-h[1][1]*h[2][0]))


def local_bounds():
    def minors(t):
        moments = [1,0,32,0,32*63,2*t]
        s = [value - Fraction(32**i,1625) for i,value in enumerate(moments)]
        return [determinant([[c*s[i+j]+sign*s[i+j+1] for j in range(3)]
                             for i in range(3)]) for c,sign in ((7,-1),(8,1))]
    y0, y1 = minors(0), minors(1)
    maximum, minimum = [-b/(y-b) for y,b in zip(y1,y0)]
    assert (minimum,maximum) == (Fraction(94736,11), Fraction(70432,7))
    return ceil(minimum), floor(maximum)


def phase_orbits():
    tuples = [word for tail in product(range(1,13), repeat=4)
              for a in [(-sum(tail))%13] if a for word in [(a,)+tail]]
    count = len({min(dihedral(w)) for w in tuples})
    assert (len(tuples), count) == (19140,1980)
    return count


if __name__ == "__main__":
    for n, expected in [(25,(99,10,2)),(32,(188,20,4))]:
        total, closed, cycles, gaps = square_window(n)
        assert (total,closed,len(cycles)) == expected
        print("Square-sum", n, "pentagons", cycles, "gap (paths, closed)", gaps)
        path_hist, cycle_hist, port_hist = hamilton_crossings(n)
        print("Hamilton paths/cycles by crossings; use of 12--13", path_hist, cycle_hist, port_hist)
        assert (sum(path_hist.values()),cycle_hist) == ((10,{}) if n==25 else (392,{16:1}))
    for a,b,N in [(2,1,5),(2,5,25),(32,25,1625)]:
        base, slope = whole_census(a,b,N)
        print("Degree", a+b, "constant/slope",base,slope)
    assert base == (6094400,-5486000,25636000,21086000,-4836000,-9386000,24986000,0)
    v0,v1 = fixed_vertex_census(0),fixed_vertex_census(1)
    assert tuple(y-x for x,y in zip(v0,v1)) == tuple(2*(-1)**sum(w) for w in WORDS)
    print("Rooted constants",v0,"local bounds",local_bounds())
    print("Dihedral zero-sum nonzero 13-frequency orbits",phase_orbits())
