#!/usr/bin/env python3
"""
polygonal_polyhedral_vandermonde_tail.py — klein-2026-07-15-S313

The two generalization axes of the triangular numbers and the Vandermonde tail.

  POLYHEDRAL axis (dimension):  S_k(j) = C(j+k-1, k)   (k-simplex numbers; Pascal columns)
  POLYGONAL axis (gonality):    P(s, j) = (s-2)*C(j,2) + j   ((s)-gonal numbers; s = k+1 for column k)

Master identity (Vandermonde):   S_k(j) = sum_{i>=0} C(k-1, i) * C(j, i+1)
Polygonal = modes i <= 1:        P(k+1, j) = C(j,1) + (k-1)*C(j,2)
Difference = VANDERMONDE TAIL:   Delta_k(j) = sum_{i>=2} C(k-1, i) * C(j, i+1)

Consequences verified here:
  1. Pascal triangle row sums = 2^n; shallow-diagonal sums = Fibonacci.
  2. Polygonal triangle row sums = 1 + C(n+1,2) + C(n+1,4) = sum_{j<=4} C(n,j) = A000127(n+1) (Moser).
  3. Polygonal shallow-diagonal sums D(d) = 1,1,2,3,5,8,13,21,33,51,76,111,157,218,...
     with closed form D(d) = (2d^4 - 12d^3 + 64d^2 + 6d + 159 + (-1)^d (33 - 6d)) / 192.
  4. Delta diagonals: j=3: C(k-1,2); j=4: 4C(k-1,2)+C(k-1,3) = (k-1)(k-2)(k+9)/6; general row-j
     coefficients on C(k-1,i) are the Pascal row tail C(j,3), C(j,4), ..., C(j,j).
  5. Delta table = Pascal kernel integrated 3x along j and 2x along k (bivariate GF
     x^3 y^3 / ((1-x)^3 (1-y)^2 (1-x-y))).
  6. Fibonacci-pole cancellation: (1-x)^3 (1-x^2)^2 - x^8 is divisible by 1-x-x^2
     (certificate: x^8 == 13 - 21x mod (1-x-x^2), Fibonacci 13, 21). Hence D is quartic.
  7. tail(d) - tail(d-1) - tail(d-2) >= 0 (GF x^8/((1-x^2)^2 (1-x)^3)), so
     D(d) <= D(d-1) + D(d-2) with equality iff d <= 7  =>  Brown criterion => D is COMPLETE.
  8. Mode/layer decompositions: 2^n = 1 + sum_{i>=0} C(n+1, 2i+2);
     F(d+1) = 1 + sum_{i>=0} L_i(d);  Moser = modes i<=1 (row);  D = modes i<=1 (diagonal).
  9. Locker parity law: tau(n) odd <=> n is a perfect square (divisor-pairing involution).
 10. No-holes completeness: the F3 exchange walk (F_k + F_{k+1} -> F_{k+2}; 2F_k -> F_{k+1} + F_{k-2})
     terminates (integer potential w(2)=8, w(c)=5c drops by >=1 per move) at the Zeckendorf
     normal form; existence + uniqueness verified; Brown transfer to D and Moser verified.
"""
from math import comb, isqrt

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

def S(k, j):            # k-simplex number, j-th (j>=1); k=0 -> all ones
    return comb(j + k - 1, k)

def P(s, j):            # s-gonal number, j-th; conventions: s=2 -> naturals, s=1 -> all ones
    if s == 1: return 1
    return (s - 2) * comb(j, 2) + j

def fib(n, memo={0: 0, 1: 1}):
    if n not in memo:
        memo[n] = fib(n - 1) + fib(n - 2)
    return memo[n]

N = 40

# ---------- 1+2: row sums ----------
pascal_rows = all(sum(comb(n, k) for k in range(n + 1)) == 2 ** n for n in range(N))
check("Pascal row sums = 2^n", pascal_rows)

def polygonal_entry(n, k):  # entry in column k of row n (0 <= k <= n): (k+1)-gonal, position j = n-k+1
    return P(k + 1, n - k + 1) if k >= 1 else 1

R = [sum(polygonal_entry(n, k) for k in range(n + 1)) for n in range(N)]
check("Polygonal row sums = 1 + C(n+1,2) + C(n+1,4)",
      all(R[n] == 1 + comb(n + 1, 2) + comb(n + 1, 4) for n in range(N)))
check("  ... = sum_{j<=4} C(n,j)  (Moser truncation of 2^n)",
      all(R[n] == sum(comb(n, j) for j in range(5)) for n in range(N)))
print("   polygonal row sums:", R[:10], "= A000127(n+1)")

# ---------- 1+3: shallow-diagonal (skip-row) sums ----------
check("Pascal diagonal sums = Fibonacci",
      all(sum(comb(d - k, k) for k in range(d // 2 + 1)) == fib(d + 1) for d in range(N)))

def D(d):
    return sum(polygonal_entry(d - k, k) for k in range(d // 2 + 1))

Dlist = [D(d) for d in range(N)]
user_claim = [1, 1, 2, 3, 5, 8, 13, 21, 33, 51, 76, 111, 157, 218]
check("Polygonal diagonal sums start 1,1,2,3,5,8,13,21,33,51,76,111,157,218", Dlist[:14] == user_claim)

closed = lambda d: (2*d**4 - 12*d**3 + 64*d**2 + 6*d + 159 + (-1)**d * (33 - 6*d)) // 192
check("D(d) = (2d^4-12d^3+64d^2+6d+159+(-1)^d(33-6d))/192 exactly",
      all(D(d) == closed(d) and (2*d**4 - 12*d**3 + 64*d**2 + 6*d + 159 + (-1)**d * (33 - 6*d)) % 192 == 0
          for d in range(N)))
check("D(2m) = (m^4-3m^3+8m^2+6)/6 and D(2m+1) = (m^4-m^3+5m^2+7m+6)/6",
      all(D(2*m) == (m**4 - 3*m**3 + 8*m**2 + 6) // 6 and
          D(2*m + 1) == (m**4 - m**3 + 5*m**2 + 7*m + 6) // 6 for m in range(N // 2 - 1)))
print("   D continued:", Dlist[14:24])

# ---------- master identity + tail ----------
JK = 26
master = all(S(k, j) == sum(comb(k - 1, i) * comb(j, i + 1) for i in range(k))
             for k in range(1, JK) for j in range(1, JK))
check("Vandermonde master: S_k(j) = sum_i C(k-1,i) C(j,i+1)", master)
check("Polygonal = modes i<=1: P(k+1,j) = C(j,1) + (k-1)C(j,2)",
      all(P(k + 1, j) == comb(j, 1) + (k - 1) * comb(j, 2) for k in range(1, JK) for j in range(1, JK)))

def Delta(k, j):
    return S(k, j) - P(k + 1, j)

check("Tail: Delta_k(j) = sum_{i>=2} C(k-1,i) C(j,i+1)",
      all(Delta(k, j) == sum(comb(k - 1, i) * comb(j, i + 1) for i in range(2, k))
          for k in range(1, JK) for j in range(1, JK)))
check("Delta = 0 exactly on columns k<=2 and rows j<=2 (agreement zone)",
      all((Delta(k, j) == 0) == (k <= 2 or j <= 2) for k in range(1, JK) for j in range(1, JK)))

# ---------- 4: diagonals of Delta ----------
d1 = [Delta(k, 3) for k in range(3, 10)]
check("1st diagonal (j=3): Delta = C(k-1,2) = 1,3,6,10,15,... (triangulars)",
      d1 == [comb(k - 1, 2) for k in range(3, 10)] and d1[:5] == [1, 3, 6, 10, 15])
d2 = [Delta(k, 4) for k in range(3, 12)]
check("2nd diagonal (j=4): Delta = 4C(k-1,2)+C(k-1,3) = (k-1)(k-2)(k+9)/6",
      d2 == [4 * comb(k - 1, 2) + comb(k - 1, 3) for k in range(3, 12)]
      and d2 == [(k - 1) * (k - 2) * (k + 9) // 6 for k in range(3, 12)]
      and d2[:5] == [4, 13, 28, 50, 80])
print("   2nd diagonal continued:", d2)
d3 = [Delta(k, 5) for k in range(3, 12)]
check("3rd diagonal (j=5): Delta = 10C(k-1,2)+5C(k-1,3)+C(k-1,4)",
      d3 == [10 * comb(k - 1, 2) + 5 * comb(k - 1, 3) + comb(k - 1, 4) for k in range(3, 12)])
print("   3rd diagonal:", d3)
check("General: j-th diagonal coefficients on C(k-1,i) are C(j,i+1), i=2..j-1 (Pascal row tails)",
      all(Delta(k, j) == sum(comb(j, i + 1) * comb(k - 1, i) for i in range(2, j))
          for j in range(3, JK) for k in range(3, JK)))
check("Fixed column k: Delta_k(j) = sum_i C(k-1,i) C(j,i+1); e.g. Delta_3 = C(j,3) (tet - square = tet)",
      all(Delta(3, j) == comb(j, 3) for j in range(1, JK)) and
      all(Delta(4, j) == 3 * comb(j, 3) + comb(j, 4) for j in range(1, JK)))

# ---------- 5: Delta = iterated partial sums of Pascal ----------
M = 22
pas = [[comb(a + b, a) for b in range(M)] for a in range(M)]  # pas[a][b] = C(a+b, a)
def cumsum_axis0(T):
    U = [row[:] for row in T]
    for a in range(1, len(U)):
        for b in range(len(U[0])):
            U[a][b] += U[a - 1][b]
    return U
def cumsum_axis1(T):
    U = [row[:] for row in T]
    for a in range(len(U)):
        for b in range(1, len(U[0])):
            U[a][b] += U[a][b - 1]
    return U
T = pas
for _ in range(3): T = cumsum_axis0(T)   # three integrations along j-direction
for _ in range(2): T = cumsum_axis1(T)   # two integrations along k-direction
check("Delta_k(j) = triple-j/double-k partial sums of Pascal, shifted to (3,3)",
      all(Delta(k, j) == T[j - 3][k - 3] for j in range(3, M) for k in range(3, M)))

# ---------- 6: Fibonacci-pole cancellation ----------
def polymul(A, B):
    C_ = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        for j2, b in enumerate(B):
            C_[i + j2] += a * b
    return C_
one_x = [1, -1]
num = polymul(polymul(one_x, polymul(one_x, one_x)), polymul([1, 0, -1], [1, 0, -1]))  # (1-x)^3 (1-x^2)^2
num = num + [0] * (9 - len(num))
num[8] -= 1                                                                            # minus x^8
# divide by 1 - x - x^2
q, r = [], num[:]
for i in range(len(r) - 2):
    q.append(r[i]); r[i + 1] += r[i]; r[i + 2] += r[i]; r[i] = 0
check("(1-x)^3(1-x^2)^2 - x^8 divisible by 1-x-x^2 (phi-pole cancels => D quartic)",
      all(c == 0 for c in r))
print("   quotient Q(x) coeffs (D = Q/((1-x)^5(1+x)^2)):", q)

# ---------- 7: tail super-Fibonacci => D sub-Fibonacci => Brown ----------
tail = [fib(d + 1) - Dlist[d] for d in range(N)]
# coefficients of x^8/((1-x^2)^2 (1-x)^3): convolve series
def series_inv_pow(period, power, length):   # coeffs of 1/(1-x^period)^power
    out = [0] * length
    for t in range(0, length, period):
        out[t] = comb(t // period + power - 1, power - 1)
    return out
L_ = N
s1 = series_inv_pow(2, 2, L_); s2 = series_inv_pow(1, 3, L_)
conv = [sum(s1[a] * s2[t - a] for a in range(t + 1)) for t in range(L_)]
rhs = [(conv[d - 8] if d >= 8 else 0) for d in range(N)]
check("tail(d) - tail(d-1) - tail(d-2) = [x^d] x^8/((1-x^2)^2(1-x)^3) >= 0",
      all(tail[d] - tail[d - 1] - tail[d - 2] == rhs[d] for d in range(2, N)))
check("D(d) <= D(d-1)+D(d-2), equality iff d <= 7",
      all((Dlist[d] <= Dlist[d - 1] + Dlist[d - 2]) and
          ((Dlist[d] == Dlist[d - 1] + Dlist[d - 2]) == (d <= 7)) for d in range(2, N)))

def brown_complete(a):    # Brown criterion: a sorted nondecreasing, a[0]=1, a[k+1] <= 1 + sum(a[:k+1])
    s = 0
    for i, t in enumerate(a):
        if i == 0 and t != 1: return False
        if t > 1 + s: return False
        s += t
    return True
check("Brown criterion holds for D  => D is a complete (no-holes) sequence", brown_complete(Dlist))
check("Brown criterion holds for Moser A000127 row sums", brown_complete(R))

def greedy_rep(Nval, terms):
    rem, used = Nval, []
    for t in reversed(terms):
        if t <= rem:
            rem -= t; used.append(t)
    return rem == 0
check("every N in 1..5000 is a sum of DISTINCT D-terms (greedy, holes: none)",
      all(greedy_rep(n, sorted(set(Dlist))) or
          # greedy with duplicates allowed once for the repeated leading 1,1
          greedy_rep(n - 1, sorted(set(Dlist) - {1})) for n in range(1, 5001)))

# ---------- 8: mode/layer decompositions ----------
check("2^n = 1 + sum_{i>=0} C(n+1, 2i+2)  (even-binomial mode decomposition)",
      all(2 ** n == 1 + sum(comb(n + 1, 2 * i + 2) for i in range(n + 1)) for n in range(N)))
def layer(i, d):
    return sum(comb(k - 1, i) * comb(d - 2 * k + 1, i + 1) for k in range(1, d // 2 + 1))
check("F(d+1) = 1 + sum_{i>=0} L_i(d)  (Fibonacci Vandermonde-layer decomposition)",
      all(fib(d + 1) == 1 + sum(layer(i, d) for i in range(d)) for d in range(N)))
check("D(d) = 1 + L_0(d) + L_1(d)  (modes i<=1) and L_0(d) = floor(d^2/4)",
      all(Dlist[d] == 1 + layer(0, d) + layer(1, d) and layer(0, d) == d * d // 4 for d in range(N)))
check("Moser row sum = modes i<=1: 1 + C(n+1,2) + C(n+1,4)",
      all(R[n] == 1 + comb(n + 1, 2) + comb(n + 1, 4) for n in range(N)))

# ---------- 9: locker parity law ----------
LOCK = 10000
state = [False] * (LOCK + 1)
for step in range(1, LOCK + 1):
    for loc in range(step, LOCK + 1, step):
        state[loc] = not state[loc]
check("locker n open <=> n is a perfect square (simulation to 10000)",
      all(state[n] == (isqrt(n) ** 2 == n) for n in range(1, LOCK + 1)))
def tau_parity_by_involution(n):
    fixed = 0; pairs = 0
    for dd in range(1, isqrt(n) + 1):
        if n % dd == 0:
            if dd * dd == n: fixed += 1
            else: pairs += 1
    return fixed % 2, 2 * pairs + fixed
check("tau(n) = 2*(pairs of involution d<->n/d) + [n square]; tau odd <=> square",
      all(tau_parity_by_involution(n)[1] == sum(1 for dd in range(1, n + 1) if n % dd == 0)
          and (tau_parity_by_involution(n)[1] % 2 == 1) == (isqrt(n) ** 2 == n)
          for n in range(1, 2000)))

# ---------- 10: F3 exchange walk => Zeckendorf (termination + confluence) ----------
def weight(c): return 8 if c == 2 else 5 * c

def exchange_walk(ms):
    """Multiset of Fibonacci indices >= 2 (F_2=1, F_3=2, ...) as a Counter. Applies F3
    exchanges until no duplicates and no adjacent indices remain. Asserts the integer
    potential Phi = sum w(c) (w(2)=8, w(c)=5c) strictly decreases at every step, and that
    total steps <= initial Phi (the termination certificate)."""
    from collections import Counter
    cnt = Counter(ms)
    phi = sum(weight(c) * m for c, m in cnt.items())
    phi0, steps = phi, 0
    def add(c, delta):
        nonlocal phi
        cnt[c] += delta
        phi += weight(c) * delta
        if cnt[c] == 0: del cnt[c]
    while True:
        before = phi
        moved = False
        for k in sorted(cnt):
            if cnt[k] >= 2:                       # duplicate rules
                add(k, -2)
                if k == 2:   add(3, 1)
                elif k == 3: add(4, 1); add(2, 1)
                else:        add(k + 1, 1); add(k - 2, 1)
                moved = True; break
            if cnt.get(k + 1, 0) >= 1:            # merge rule {k, k+1} -> {k+2}
                add(k, -1); add(k + 1, -1); add(k + 2, 1)
                moved = True; break
        if not moved:
            assert steps <= phi0, "termination bound violated"
            return sorted(cnt.elements()), steps
        assert phi < before, "potential failed to decrease"
        steps += 1

def zeckendorf_greedy(n):
    out = []
    k = 2
    while fib(k + 1) <= n: k += 1
    while n:
        while fib(k) > n: k -= 1
        out.append(k); n -= fib(k); k -= 2
    return sorted(out)

ok10 = True
for n in range(1, 601):
    terminal, _ = exchange_walk([2] * n)              # start from n ones
    vals = sum(fib(c) for c in terminal)
    zeck = zeckendorf_greedy(n)
    legal = all(terminal[i + 1] - terminal[i] >= 2 for i in range(len(terminal) - 1))
    if not (vals == n and legal and terminal == zeck):
        ok10 = False; break
check("F3 exchange walk from n*[F_2]: terminates, value-preserving, lands on the unique "
      "Zeckendorf form (n=1..600)", ok10)

# also from random non-canonical representations
import itertools
ok10b = True
for n in range(1, 300):
    # a deliberately messy representation: binary-of-n over Fibonacci with repeats
    rep = []
    rem = n
    while rem:
        k = 2
        while fib(k + 1) <= rem: k += 1
        take = min(2, rem // fib(k)) or 1
        rep += [k] * take; rem -= take * fib(k)
    terminal, _ = exchange_walk(rep)
    if terminal != zeckendorf_greedy(n): ok10b = False; break
check("F3 exchange walk is confluent: messy start reps normalize to the same Zeckendorf form", ok10b)

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed, {len(fails)} failed ===")
if fails:
    for f in fails: print("FAILED:", f)
print()
print("OEIS-ready sequences:")
print("  D(d), d>=0:        ", Dlist[:20])
print("  2nd Delta diagonal:", [Delta(k, 4) for k in range(3, 13)])
print("  3rd Delta diagonal:", [Delta(k, 5) for k in range(3, 13)])
