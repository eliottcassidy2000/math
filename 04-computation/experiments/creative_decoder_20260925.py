"""Exact two-direction Collatz proof compiler; no universal termination claim."""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


def check(p, msg="check failed"):
    if not p:
        raise RuntimeError(msg)


def odd_step(n, sign=1):
    check(n > 0 and n % 2 == 1, "odd positive domain")
    z = 3*n + sign
    k = (z & -z).bit_length()-1
    return z >> k, k


def iterate(n, count, sign=1):
    for _ in range(count):
        n = odd_step(n, sign)[0]
    return n


def rho(n, sign=1):
    while n % 8 == (5 if sign == 1 else 3):
        n = (n-sign)//4
    return n


@dataclass(frozen=True)
class Word:
    ks: tuple
    K: int
    C: int

    def __post_init__(self):
        K = C = 0
        for k in self.ks:
            check(type(k) is int and k >= 1, "positive integer exponent required")
            C = 3*C+2**K
            K += k
        check(type(self.K) is int and type(self.C) is int and (self.K, self.C) == (K, C),
              "word metadata disagrees with exact exponents")

    @classmethod
    def make(cls, ks):
        ks = tuple(ks)
        K = C = 0
        for k in ks:
            check(k >= 1, "zero exponent is not a legal odd Collatz edge")
            C = 3*C + 2**K
            K += k
        return cls(tuple(ks), K, C)

    @property
    def L(self):
        return len(self.ks)

    def inverse(self, n, sign=1):
        check(type(n) is int and n > 0 and n % 2 == 1, "positive odd target required")
        check(sign in (-1, 1), "invalid sign")
        z = 2**self.K*n-sign*self.C
        den = 3**self.L
        if z % den:
            return None
        y = z//den
        if y <= 0 or y % 2 == 0:
            return None
        return y

    def forward(self, n, sign=1):
        check(type(n) is int and n > 0 and n % 2 == 1, "positive odd source required")
        check(sign in (-1, 1), "invalid sign")
        z = 3**self.L*n+sign*self.C
        den = 2**self.K
        if z % den:
            return None
        y = z//den
        if y <= 0 or y % 2 == 0:
            return None
        return y

    def replay(self, n, sign=1):
        for k in self.ks:
            n, actual = odd_step(n, sign)
            check(k == actual, "wrong exact exponent")
        return n


@dataclass(frozen=True)
class Join:
    x: int
    y: int
    a: int
    b: int
    label: str

    def verify(self, sign=1):
        check(all(type(v) is int for v in (self.x, self.y, self.a, self.b)), "integer join fields required")
        check(self.a >= 0 and self.b >= 0, "nonnegative clocks required")
        check(sign in (-1, 1), "invalid sign")
        check(0 < self.y < self.x and self.x % 2 == self.y % 2 == 1,
              "join must strictly decrease positive odd input")
        check(iterate(self.x, self.a, sign) == iterate(self.y, self.b, sign),
              "false common-future claim")


def sibling_chain(n, sign=1):
    yield n
    while n % 8 == (5 if sign == 1 else 3):
        n = (n-sign)//4
        yield n


def base_joins(n, sign=1):
    u = odd_step(n, sign)[0]
    candidates = [(u, 1, 0, "forward")]
    candidates += [(s, 1, 1, "sibling") for s in sibling_chain(n, sign) if s != n]
    candidates += [(s, 2, 1, "forward-sibling") for s in sibling_chain(u, sign)]
    if (2*n-sign) % 3 == 0:
        candidates.append(((2*n-sign)//3, 0, 1, "inverse-one"))
    for y, a, b, label in candidates:
        if 0 < y < n and y % 2:
            yield Join(n, y, a, b, label)


def compositions(total, length):
    for cuts in combinations(range(1, total), length-1):
        ends = (0,)+cuts+(total,)
        yield tuple(ends[i+1]-ends[i] for i in range(length))


def inverse_bank(depth):
    bank = []
    for L in range(1, depth+1):
        maxK = (3**L).bit_length()-1
        for K in range(L, maxK+1):
            for ks in compositions(K, L):
                w = Word.make(ks)
                check(2**w.K < 3**w.L)
                bank.append(w)
    return bank


def inverse_joins(n, bank, sign=1):
    for w in bank:
        y = w.inverse(n, sign)
        if y is not None and y < n:
            yield Join(n, y, 0, w.L, "inverse-word"), w


def root_closure(limit, bank, forward=True, quotient=True):
    known = {1}
    for n in range(3, limit+1, 2):
        moves = [j for j in base_joins(n)
                 if (forward or j.label != "forward")
                 and (quotient or j.label not in {"sibling", "forward-sibling"})]
        moves += [j for j, w in inverse_joins(n, bank)]
        if any(j.y in known for j in moves):
            known.add(n)
    return known


def compile_interval(limit, bank, complete_suffixes=False):
    heights = {1: 0}
    proof = {}
    learned = []
    uses = {}
    for n in range(3, limit+1, 2):
        moves = list(base_joins(n))
        moves += [j for j, w in inverse_joins(n, bank)]
        for w in learned:
            y = w.forward(n)
            if y is not None and y < n:
                moves.append(Join(n, y, w.L, 0, "learned-port"))
        if not moves:
            x, ks = n, []
            while x >= n:
                x, k = odd_step(x)
                ks.append(k)
                check(len(ks) <= 10000, "finite experiment cap exceeded")
            w = Word.make(ks)
            check(2**w.K > 3**w.L, "forward descent slope")
            check(w.forward(n) == x)
            additions = [w]
            if complete_suffixes:
                additions += [Word.make(w.ks[i:]) for i in range(1, w.L)]
            for port in additions:
                check(2**port.K > 3**port.L, "suffix slope")
                if port not in learned:
                    learned.append(port)
            moves.append(Join(n, x, w.L, 0, "new-port"))
        j = min(moves, key=lambda j: (j.a+max(heights[j.y]-j.b, 0), j.y, j.label))
        j.verify()
        heights[n] = j.a+max(heights[j.y]-j.b, 0)
        check(iterate(n, heights[n]) == 1, "compiled root height")
        proof[n] = j
        uses[j.label] = uses.get(j.label, 0)+1
    return proof, heights, learned, uses


def main():
    print("PROVED RULES + FINITE-EXACT REPLAY; global Collatz remains OPEN")
    bank = inverse_bank(8)
    print("Inverse contracting port bank: depth<=8, words=", len(bank))
    for depth in (1, 2, 4, 6, 8):
        sub = [w for w in bank if w.L <= depth]
        ternary = 3**depth
        covered = set()
        for w in sub:
            mod = 3**w.L
            r = (w.C*pow(2**w.K, -1, mod)) % mod
            covered.update(range(r, ternary, mod))
        combined = sum(1 for n in range(1, 16*ternary, 2)
                       if n % 4 == 1 or n % 16 == 3 or n % ternary in covered)
        print("Depth", depth, "inverse density", Fraction(len(covered), ternary),
              "combined local descent density", Fraction(combined, 8*ternary))
    w12 = Word.make((1, 2))
    check(w12.C == 5 and w12.K == 3)
    for n in range(1, 10000, 2):
        y = w12.inverse(n)
        check((y is not None) == (n % 18 == 13), "two-step ternary guard")
        if y is not None:
            check(y < n and w12.replay(y) == n)
    residual144 = [n for n in range(1, 144, 2)
                   if not (n % 4 == 1 or n % 16 == 3 or n % 6 == 5 or n % 18 == 13)]
    print("Residual modulo144 for base+two-step rule:", residual144)
    check(len(residual144) == 15)
    print("Witness31 <- 41 <- 27:", w12.replay(27), "; rank reduction31->27")
    limit = 10000
    baseline = root_closure(limit, [], quotient=False)
    dualbase = root_closure(limit, [])
    allbank = root_closure(limit, bank)
    print("Root closure on odds<=10000, forward+reverse1 / all base / base+inverse bank:",
          len(baseline), len(dualbase), len(allbank))
    check(baseline <= dualbase <= allbank)
    print("First ten new root-certified inputs from inverse bank:", sorted(allbank-dualbase)[:10])
    plain = compile_interval(limit, bank)
    proof, heights, learned, uses = compile_interval(limit, bank, complete_suffixes=True)
    no_inverse = compile_interval(limit, [], complete_suffixes=True)
    print("Learned seed obligations, raw / suffix-complete / suffix-complete without inverse bank:",
          plain[3].get("new-port", 0), uses.get("new-port", 0), no_inverse[3].get("new-port", 0))
    print("Adaptive compiler: certified odds=", len(heights), "learned infinite ports=", len(learned))
    print("Chosen joins:", dict(sorted(uses.items())))
    print("Maximum expanded odd steps:", max(heights.values()))
    maxL = max(w.L for w in learned)
    a = 2*((maxL+3)//2)
    hostile = 2**a-1
    check(a > maxL+1 and hostile % 3 == 0)
    check(not list(base_joins(hostile)))
    check(not list(inverse_joins(hostile, bank)))
    check(not any((y := w.forward(hostile)) is not None and y < hostile for w in learned))
    print("Finite-bank obstruction: max forward word length=", maxL,
          "; uncovered all-ones input2^a-1 at a=", a,
          "; all larger even exponents also uncovered")
    print("First learned exact-exponent words:", [w.ks for w in learned[:8]])
    # Independent direct ordinary-map root check, separate clock and implementation.
    peak = 0
    for n in range(1, limit+1):
        x, steps = n, 0
        while x != 4:
            x = 3*x+1 if x % 2 else x//2
            steps += 1
            check(steps < 10000, "ordinary root check cap")
        peak = max(peak, steps)
    print("Independent ordinary root4 check: all1..10000; maximumsteps=", peak)
    # Every learned port has an infinite positive source cylinder, and carries
    # are replayed at separated points without running the generating search.
    probes = 0
    for w in bank:
        mod = 3**w.L
        r = w.C*pow(2**w.K, -1, mod) % mod
        n0 = r if r % 2 else r+mod
        if n0 == 0:
            n0 = 2*mod
        for t in (0, 1, 7):
            n = n0+2*mod*t
            y = w.inverse(n)
            check(y is not None and 0 < y < n)
            check(w.replay(y) == n)
            probes += 1
    for w in learned:
        den = 2**(w.K+1)
        r = (2**w.K-w.C)*pow(3**w.L, -1, den) % den
        # Move beyond exact carry threshold; no floating-point arithmetic.
        threshold = w.C//(2**w.K-3**w.L)+1
        n0 = r+max(0, (threshold-r+den-1)//den)*den
        for t in (0, 1, 7):
            n = n0+den*t
            y = w.forward(n)
            check(y is not None and y < n and w.replay(n) == y)
            probes += 1
    print("Independent algebra/guard replay at port samples:", probes)
    # Opposite-sheet controls: local proof rules preserve distinct cycles.
    for seed, cycle in ((1, (1,)), (5, (5,7)), (17, (17,25,37,55,41,61,91))):
        check(iterate(seed, len(cycle), -1) == seed)
        moves = list(base_joins(seed, -1))+[j for j,w in inverse_joins(seed, bank, -1)]
        check(not moves, "cycle minimum acquired a false strict certificate")
    print("Minus-sheet cycle minima1,5,17: no strict local join; allcycles replayed")
    rejected = 0
    for j in (Join(5, 1, 0, 1, "wrong sign"), Join(7, 7, 0, 0, "circular"),
              Join(31, 26, 0, 2, "lost parity"), Join(5, 3, -1, 1, "negative clock")):
        try:
            j.verify()
        except RuntimeError:
            rejected += 1
    for thunk in (lambda: Word.make((1,)).inverse(2),
                  lambda: Word((1,), 2, 1)):
        try:
            thunk()
        except RuntimeError:
            rejected += 1
    check(rejected == 6)
    print("Hostile sign, circle, parity, clock, target-domain and metadata certificates rejected:", rejected)
    print("Residual obligation: every positive odd input needs a finite chain of verified joins to1.")


if __name__ == "__main__":
    main()
