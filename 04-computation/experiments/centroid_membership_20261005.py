"""Finite centroid addresses versus rational point cycles, with exact local11 algebra.

No numerical tolerance, root-convergence oracle, or import-time census.
Run normally or with -O; all checks survive optimization.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import gcd, lcm, comb

import small_port_refinement_20261004 as inherited


CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def rational_point(point):
    """Canonical primitive numerators and positive common denominator."""
    if not isinstance(point, (tuple, list)) or len(point) != 3:
        raise ValueError("three exact barycentric coordinates required")
    if any(type(x) not in (int, F) for x in point):
        raise ValueError("coordinates must be exact integers or Fractions")
    point = tuple(map(F, point))
    if any(x < 0 for x in point) or sum(point) != 1:
        raise ValueError("point must lie in the closed simplex")
    denominator = lcm(*(x.denominator for x in point))
    numerators = tuple(int(x*denominator) for x in point)
    need(gcd(*numerators) == 1, "minimal common denominator")
    return numerators, denominator


def point_from_numerators(numerators):
    return tuple(F(x, sum(numerators)) for x in numerators)


def strip(numerators):
    """Strip the unique-minimum chart, without normalizing its denominator."""
    minimum = min(numerators)
    if minimum <= 0 or numerators.count(minimum) != 1:
        raise ValueError("strict interior fan chart required")
    c = numerators.index(minimum)
    a, b = (c+1) % 3, (c+2) % 3
    return c, (numerators[a]-minimum, numerators[b]-minimum, 3*minimum)


def encode(word):
    if type(word) is not tuple or any(type(c) is not int or c not in (0, 1, 2) for c in word):
        raise ValueError("finite fan word must be a tuple of exact slots0,1,2")
    numerators, denominator = (1, 1, 1), 3
    # The outside chart is the first letter, so construct from the inside.
    for c in reversed(word):
        a, b = (c+1) % 3, (c+2) % 3
        out = [0, 0, 0]
        out[c] = numerators[2]
        out[a] = 3*numerators[0]+numerators[2]
        out[b] = 3*numerators[1]+numerators[2]
        numerators, denominator = tuple(out), 3*denominator
    need(sum(numerators) == denominator and gcd(*numerators) == 1,
         "encoded point remains reduced")
    need(all(x % 3 == 1 for x in numerators), "finite-centroid numerator residue")
    return tuple(F(x, denominator) for x in numerators)


def decode(point):
    """Return the finite centroid word, or None for a well-typed nonmember."""
    numerators, denominator = rational_point(point)
    if min(numerators) == 0:
        return None
    power = 0
    remaining = denominator
    while remaining % 3 == 0:
        remaining //= 3
        power += 1
    if remaining != 1 or power == 0:
        return None
    word = []
    for _ in range(power-1):
        if numerators.count(min(numerators)) != 1:
            return None
        c, out = strip(numerators)
        common = gcd(*out)
        if common != 3:
            return None
        numerators = tuple(x//3 for x in out)
        word.append(c)
    return tuple(word) if numerators == (1, 1, 1) else None


@dataclass(frozen=True)
class Orbit:
    status: str
    denominator: int
    states: tuple
    letters: tuple
    cycle_start: int | None = None

    @property
    def period(self):
        return 0 if self.cycle_start is None else len(self.states)-self.cycle_start


def orbit(point):
    """Total exact rational-point classifier: center, boundary, tie, or cycle.

    Keep the initial reduced denominator fixed, so the state space is finite.
    Center is checked before ties. A noncentral tie is a stopping boundary,
    never a fabricated orientation choice.
    """
    numerators, denominator = rational_point(point)
    if min(numerators) == 0:
        return Orbit("boundary", denominator, (numerators,), ())
    seen, states, letters = {}, [], []
    while numerators not in seen:
        if len(set(numerators)) == 1:
            return Orbit("center", denominator, tuple(states)+(numerators,), tuple(letters))
        if numerators.count(min(numerators)) != 1:
            return Orbit("tie", denominator, tuple(states)+(numerators,), tuple(letters))
        seen[numerators] = len(states)
        states.append(numerators)
        c, numerators = strip(numerators)
        letters.append(c)
        need(sum(numerators) == denominator and min(numerators) > 0,
             "fixed-denominator interior state")
        need(len(states) <= comb(denominator-1, 2), "finite state bound")
    return Orbit("cycle", denominator, tuple(states), tuple(letters), seen[numerators])


def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def ma(a, b):
    return tuple(tuple(x+y for x, y in zip(r, s)) for r, s in zip(a, b))


def scale(c, a):
    return tuple(tuple(c*x for x in row) for row in a)


def matrix_power(a, k):
    need(type(k) is int and k >= 0, "nonnegative exact matrix exponent")
    out = tuple(tuple(int(i == j) for j in range(len(a))) for i in range(len(a)))
    while k:
        if k % 2:
            out = mm(out, a)
        a, k = mm(a, a), k//2
    return out


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def determinant2(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def transpose(a):
    return tuple(zip(*a))


def word_controls():
    count = frame_controls = 0
    points = set()
    for depth in range(7):
        for word in product(range(3), repeat=depth):
            point = encode(word)
            need(point not in points, "distinct finite words give distinct centroid points")
            points.add(point)
            need(decode(point) == word, "point word decoder")
            num, den = rational_point(point)
            need(den == 3**(depth+1), "denominator determines exact finite depth")
            classified = orbit(point)
            need(classified.status == "center" and classified.letters == word,
                 "finite-state and denominator decoders agree")
            if depth <= 4:
                matrix = inherited.fan_word_matrix(word)
                need(inherited.apply(matrix, (F(1, 3),)*3) == point,
                     "independent inherited matrix construction")
                for p in inherited.PERMS:
                    framed = inherited.mul(matrix, inherited.permutation_matrix(p))
                    need(inherited.apply(framed, (F(1, 3),)*3) == point,
                         "centroid loses terminal port permutation")
                    need(inherited.decode_fan(framed) == (word, p),
                         "retained matrix recovers the missing frame")
                    frame_controls += 1
            count += 1
    for k in range(1, 5):
        D = 3**k
        primitive = accepted = 0
        for a in range(1, D-1):
            for b in range(1, D-a):
                nums = (a, b, D-a-b)
                if gcd(*nums) != 1:
                    continue
                primitive += 1
                accepted += decode(point_from_numerators(nums)) is not None
        need(primitive == 4*3**(2*k-2)-3**k,
             "primitive exact denominator census")
        need(accepted == D//3, "all finite centroids at this denominator")
        need(F(accepted, primitive) == F(3, 4*D-9), "exact language fraction")
    return count, frame_controls


def rational_controls():
    cycle1 = ((1, 2, 6), (1, 5, 3), (4, 2, 3))
    cycle2 = ((2, 1, 6), (5, 1, 3), (2, 4, 3))
    for states, letters in ((cycle1, (0, 0, 1)), (cycle2, (1, 1, 0))):
        o = orbit(point_from_numerators(states[0]))
        need(o.status == "cycle" and o.denominator == 9 and o.states == states and
             o.letters == letters and o.cycle_start == 0, "primitive denominator9 period3")
        need(decode(point_from_numerators(states[0])) is None, "power-of3 denominator is insufficient")
    point = point_from_numerators(cycle1[0])
    for j in range(5):
        matrix = inherited.fan_word_matrix((0, 0, 1)*j)
        need(inherited.apply(matrix, point) == point, "fixed noncentral marker loses every repeated block")
        need(abs(inherited.determinant(matrix)) == F(1, 3**(3*j)),
             "the corresponding marked matrices remain distinct")
    totals = Counter()
    smallest_period3 = None
    total = 0
    for D in range(3, 19):
        for a in range(1, D-1):
            for b in range(1, D-a):
                nums = (a, b, D-a-b)
                if gcd(*nums) != 1:
                    continue
                o = orbit(point_from_numerators(nums))
                need((decode(point_from_numerators(nums)) is not None) == (o.status == "center"),
                     "two independent membership decisions")
                # An integer branch determinant3 can lose no other denominator prime.
                if nums.count(min(nums)) == 1:
                    _, out = strip(nums)
                    need(gcd(*out) in (1, 3), "primitive denominator loses at most one factor3")
                if o.period == 3 and smallest_period3 is None:
                    smallest_period3 = D
                if D == 9:
                    totals[o.status] += 1
                total += 1
    need(smallest_period3 == 9, "smallest denominator of a strict period3 cycle in complete range")
    need(totals == {"center": 3, "tie": 6, "cycle": 18}, "complete primitive denominator9 graph")
    need(orbit((F(1, 6), F(2, 6), F(3, 6))).period == 1,
         "period3 is not the first nonterminating rational point")
    for point in ((1, 0, 0), (F(1, 2), F(1, 2), 0), (0, F(1, 3), F(2, 3))):
        need(orbit(point).status == "boundary" and decode(point) is None, "boundary rejected immediately")
    need(orbit((F(1, 3),)*3).status == "center" and decode((F(1, 3),)*3) == (),
         "center stopped before minimum tie")
    need(orbit((F(1, 9), F(1, 9), F(7, 9))).status == "tie", "noncentral tie rejected")
    return total, dict(totals), cycle1, cycle2


def local11_controls():
    T0 = ((-1, 1, 0), (-1, 0, 1), (3, 0, 0))
    T1 = ((0, -1, 1), (1, -1, 0), (0, 3, 0))
    M = mm(T1, mm(T0, T0))
    need(M == ((-7, 4, 0), (-4, 0, 1), (12, -3, 0)), "return001 chronological order")
    need(tuple(sum(M[i][j]*x for j, x in enumerate((1, 2, 6))) for i in range(3)) == (1, 2, 6),
         "return fixes the rational point numerator")
    U = ((1, 0), (0, 1), (-1, -1))
    L = ((-7, 4), (-5, -1))
    need(mm(M, U) == mm(U, L), "sum-zero tangent-plane restriction")
    need(trace(L) == -8 and determinant2(L) == 27, "transverse polynomial z2+8z+27")
    need(trace(L)**2-4*determinant2(L) == -44, "quadratic discriminant minus44")
    Qmetric = ((5, -3), (-3, 4))
    need(Qmetric[0][0] > 0 and determinant2(Qmetric) == 11 and
         mm(transpose(L), mm(Qmetric, L)) == scale(27, Qmetric),
         "positive quadratic metric and exact block contraction")
    inverseL = ((F(-1, 27), F(-4, 27)), (F(5, 27), F(-7, 27)))
    need(mm(L, inverseL) == ((1, 0), (0, 1)) and trace(inverseL) == F(-8, 27) and
         determinant2(inverseL) == F(1, 27), "transverse inverse factor")
    target = (F(1, 9), F(2, 9), F(6, 9))
    for j in range(7):
        point = encode((0, 0, 1)*j)
        _, denominator = rational_point(point)
        need(denominator == 3**(3*j+1), "finite periodic-prefix denominator")
        x, y = point[0]-target[0], point[1]-target[1]
        need(5*x*x-6*x*y+4*y*y == F(4, 27**(j+1)),
             "exact center-to-cycle contraction, not a Collatz rank")
        repeated2 = inherited.fan_word_matrix((2,)*j)
        for vertex in ((1, 0, 0), (0, 1, 0)):
            need(inherited.apply(repeated2, vertex) == vertex,
                 "arbitrary infinite fan words need not shrink their cells")
    I = ((1, 0), (0, 1))
    C = scale(F(-1, 2), ma(L, scale(3, I)))
    need(C == ((F(2), F(-2)), (F(5, 2), F(-1))), "rational cubic root")
    need(ma(ma(mm(C, C), scale(-1, C)), scale(3, I)) == ((0, 0), (0, 0)),
         "cube root has polynomial z2-z+3")
    need(matrix_power(C, 3) == L, "cubic return exactly recovered")
    S, Sinv = ((2, 0), (0, 1)), ((F(1, 2), 0), (0, 1))
    C2 = mm(Sinv, mm(C, S))
    L2 = mm(Sinv, mm(L, S))
    need(C2 == ((2, -1), (5, -1)) and L2 == ((-7, 2), (-10, -1)),
         "index2 integral lattice sidecar")
    D = scale(-1, C2)
    need(trace(D) == -1 and determinant2(D) == 3 and matrix_power(D, 3) == scale(-1, L2),
         "local degree3 Euler-factor candidate and signed cubic match")
    Scomp = ((0, 1), (1, 2))
    companion = ((-1, -3), (1, 0))
    need(determinant2(Scomp) == -1 and mm(D, Scomp) == mm(Scomp, companion),
         "integral Hecke-companion conjugacy after the index2 lattice change")
    # Formal product q*prod(1-q^n)^2*(1-q^(11n))^2 through degree3:
    # only (1-q)^2*(1-q^2)^2 contribute through degree2 after removing q.
    first = (1, -2, 1)
    second = (1, 0, -2)
    product2 = tuple(sum(first[j]*second[k-j] for j in range(k+1)) for k in range(3))
    need(product2 == (1, -2, -1), "formal level11 eta-product coefficients a1,a2,a3")
    need(trace(D) == product2[2], "exact a3 local polynomial match")
    truncated = [1]+[0]*26
    for exponent in tuple(range(1, 27))+tuple(range(11, 27, 11)):
        for _ in range(2):
            for degree in range(26, exponent-1, -1):
                truncated[degree] -= truncated[degree-exponent]
    need((truncated[1], truncated[2], truncated[26]) == (-2, -1, 5),
         "independent formal product a2,a3,a27")
    mod3 = tuple(tuple(int(x) % 3 for x in row) for row in D)
    need(mod3 == ((1, 1), (1, 1)), "local3 reduction is singular, unlike the golden unit clock")
    # Independent trace recurrence, with determinant retained.
    previous, current = 2, -1
    traces = [previous, current]
    for k in range(2, 17):
        previous, current = current, -current-3*previous
        need(trace(matrix_power(D, k)) == current, "local trace recurrence")
        traces.append(current)
    need(traces[:4] == [2, -1, -5, 8], "square and cube traces have different conventions")
    reciprocal = [1, -1]
    for k in range(2, 17):
        reciprocal.append(-reciprocal[-1]-3*reciprocal[-2])
        need(traces[k] == reciprocal[k]-3*reciprocal[k-2],
             "trace versus reciprocal-factor coefficient sidecar")
    need(reciprocal[:4] == [1, -1, -2, 5], "local prime-power coefficients are not traces")
    return M, L, C2, D, traces[:9]


def hostiles():
    for call in (lambda: encode((True,)), lambda: encode((3,)),
                 lambda: decode((0.2, 0.3, 0.5)),
                 lambda: decode((1, 1, 1)), lambda: decode((True, 0, 0))):
        try:
            call()
        except ValueError:
            continue
        raise ValueError("invalid exact input accepted")


def main():
    words, frames = word_controls()
    points, counts, cycle1, cycle2 = rational_controls()
    M, L, C, D, traces = local11_controls()
    hostiles()
    print("PROVED: finite fan centroid points recover the full word; exact denominator=3^(depth+1).")
    print("MISSING SIDECAR: terminal port permutation is lost; arbitrary rational points need not terminate.")
    print("FINITE WORDS:", words, "through depth6; independent terminal-frame controls:", frames)
    print("RATIONAL UNIVERSE:", points, "primitive interior triples with denominator3..18.")
    print("DENOMINATOR9: primitive counts", counts, "; strict cycles", cycle1, "and", cycle2)
    print("CYCLE WORDS:001 and110, primitive period3; denominator6 already has fixed nonterminating points.")
    print("CENTER MEMBERS AT denominator3^k:3^(k-1); primitive total=4*3^(2k-2)-3^k; fraction3/(4*3^k-9).")
    print("RETURN001:", M, "; transverse:", L, "; characteristic z^2+8z+27, discriminant-44.")
    print("INDEX2 CUBIC ROOT:", C, "; signed local operator:", D, "; det(I-tD)=1+t+3t^2.")
    print("COMPANION: D*S=S*[[-1,-3],[1,0]], S=[[0,1],[1,2]], detS=-1; formal eta a3=-1.")
    print("LOCAL TRACES:", traces, "; signed local cube is minus the transverse return.")
    print("CONTROLLED LIMIT: E_Q(A001^j*center-p)=4/27^(j+1), Q=[[5,-3],[-3,4]], j0..6 checked.")
    print("INFINITE-WORD HOSTILE: A2^j fixes both endpoints of the same full edge for everyj.")
    print("SCOPE: exact local polynomial match only, not a conductor derivation/global dynamical equivalence.")
    print("BOUNDARIES: center before ties; noncentral ties/zero coordinates reject; denominator3-power alone fails.")
    print("CHECKS:", CHECKS)


if __name__ == "__main__":
    main()
