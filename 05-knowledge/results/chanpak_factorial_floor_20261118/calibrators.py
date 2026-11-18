"""# Strong monomial calibrators

The companion `calibrators.py` gives literal strong local pairs on orders
5, 7 and 9, differing in one arc, whose whole optional-slot polynomials
differ by exactly k q for k=2,4,6,8,10. Outgoing rows encode arcs by bits,
least significant bit indexed 0. All row pairs, direct slot counts and
strong connectivity are independently checked by the script (run from
this directory).

Why finitely many slot reads certify the entire identity: odd-cycle
packing makes the polynomial odd for these odd-order blocks. After
subtracting k q and dividing by q, it is a polynomial in q^2 of degree
at most floor(n/2). The script uses floor(n/2)+1 distinct positive reads.

Section 4 of the main proof therefore plugs any of the five signed moves
into an arbitrary source at a single replaced vertex. The k q difference
affects only the one-run coordinate, so the two physical supertournaments
differ in one arc and their Hamiltonian counts differ by k H(source).
All optional-slot counts change by k P_s(source) in parallel. These are
strong hosts for any nontrivial strong source. This is a multiplicative
scale selector, not safe addition of independently calibrated differences:
different calibrators need not have the same baseline or higher profile.
"""

from verify import h

PAIRS = [
    (2,(6,12,8,17,7),(2,12,9,17,7)),
    (4,(6,24,42,81,37,75,23),(14,24,42,80,37,75,23)),
    (6,(6,12,8,113,39,71,23),(2,12,9,113,39,71,23)),
    (8,(480,277,401,87,193,158,38,330,120),(480,277,409,83,193,158,38,330,120)),
    (10,(46,376,386,436,197,20,429,291,49),(46,376,386,436,197,276,429,291,17)),
]


def slot_values(a):
    n = len(a)
    size = 1 << n
    # Every subset path count, with the empty path counted once.
    hp = [1]
    for S in range(1, size):
        vs = [i for i in range(n) if S >> i & 1]
        hp.append(h([sum(1 << k for k,j in enumerate(vs) if a[i] >> j & 1)
                     for i in vs]))
    prev = [1]+[0]*(size-1)
    out = []
    for _ in range(n//2+1):
        post = [0]*size
        for S in range(size):
            T = S
            while True:
                post[S] += prev[S^T]*hp[T]
                if T == 0:
                    break
                T = (T-1)&S
        out.append(post[-1])
        prev = post
    return out


def valid(a):
    n = len(a)
    if any((row >> i & 1) or row >= 1 << n for i,row in enumerate(a)):
        return False
    if any((a[i] >> j & 1)+(a[j] >> i & 1) != 1 for i in range(n) for j in range(i)):
        return False
    for i in range(n):
        reach = 1 << i
        for _ in a:
            for j in range(n):
                if reach >> j & 1:
                    reach |= a[j]
        if reach != (1 << n)-1:
            return False
    return True


if __name__=='__main__':
    for k,a,b in PAIRS:
        x,y = slot_values(a),slot_values(b)
        if not (valid(a) and valid(b) and
                sum((i^j).bit_count() for i,j in zip(a,b)) == 2 and
                [j-i for i,j in zip(x,y)] == [k*s for s in range(1,len(x)+1)]):
            raise RuntimeError((k,x,y))
        print(k,len(a),x,y)
