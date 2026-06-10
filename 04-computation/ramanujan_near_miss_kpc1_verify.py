# ramanujan_near_miss_kpc1_verify.py
# ADVERSARIAL VERIFIER for session kind-pasteur-2026-06-10-S1, thread D
# (claims D7, D8). FRESH code; method deliberately DIFFERENT from worker:
#   - coefficients computed by direct power-series division (convolution
#     against D(x)), NOT via a companion matrix;
#   - the all-n-in-Z proof is routed through FACTORING the characteristic
#     cubic as (t+1)(t^2-83t+1) and an exponential-polynomial/Vandermonde
#     argument with at most 14 distinct bases -- NOT the worker's 27x27
#     Kronecker-cube Cayley-Hamilton route.
# Pure python integers. No floats (isqrt only, exact).

from math import isqrt

fails = 0
def check(name, ok, detail=""):
    global fails
    tag = "PASS" if ok else "FAIL"
    if not ok:
        fails += 1
    print(f"[{tag}] {name}" + (f"  {detail}" if detail else ""))

# ---------------------------------------------------------------- PART 1
print("=" * 72)
print("PART 1: series coefficients by direct division against D(x)")
print("=" * 72)
# D(x) = 1 - 82x - 82x^2 + x^3, low-to-high coefficients:
Dc = [1, -82, -82, 1]
NUM = {
    "a": [1, 53, 9],
    "b": [2, -26, -12],
    "c": [2, 8, -10],
}
NFWD = 41   # n = 0..40

def series_div(P, N):
    """coefficients u_0..u_{N-1} of P(x)/D(x) via exact convolution."""
    u = []
    for n in range(N):
        s = P[n] if n < len(P) else 0
        for k in range(1, 4):
            if n - k >= 0:
                s -= Dc[k] * u[n - k]
        u.append(s)   # Dc[0] == 1
    return u

seq = {name: series_div(P, NFWD) for name, P in NUM.items()}

print(f"  a_0..a_3 = {seq['a'][:4]}")
print(f"  b_0..b_3 = {seq['b'][:4]}")
print(f"  c_0..c_3 = {seq['c'][:4]}")
check("first triples match claim (1,2,2), (135,138,172), (11161,11468,14258)",
      (seq['a'][0], seq['b'][0], seq['c'][0]) == (1, 2, 2) and
      (seq['a'][1], seq['b'][1], seq['c'][1]) == (135, 138, 172) and
      (seq['a'][2], seq['b'][2], seq['c'][2]) == (11161, 11468, 14258))

# sanity: multiply back  D(x) * U(x) and confirm it reproduces P(x) exactly
def mul_trunc(A, B, N):
    out = [0] * N
    for i, ai in enumerate(A):
        if i >= N:
            break
        for j, bj in enumerate(B):
            if i + j >= N:
                break
            out[i + j] += ai * bj
    return out

for name, P in NUM.items():
    back = mul_trunc(Dc, seq[name], NFWD)
    want = P + [0] * (NFWD - len(P))
    check(f"D(x)*{name}(x) == numerator (back-multiplication, n<41)", back == want)

# homogeneous recurrence holds for n>=3 (numerator degree 2 < 3) -- verify
for name in NUM:
    u = seq[name]
    ok = all(u[n] == 82 * u[n - 1] + 82 * u[n - 2] - u[n - 3]
             for n in range(3, NFWD))
    check(f"{name}: homogeneous recurrence u_n=82u_(n-1)+82u_(n-2)-u_(n-3) for 3<=n<=40", ok)

# ---------------------------------------------------------------- PART 2
print()
print("=" * 72)
print("PART 2: canonical backward extension (trailing coefficient of D is 1)")
print("=" * 72)
NBWD = 12   # n = -1..-12
back = {}
for name in NUM:
    u = dict(enumerate(seq[name]))
    # from u_n = 82u_(n-1)+82u_(n-2)-u_(n-3):  u_(n-3) = 82u_(n-1)+82u_(n-2)-u_n
    for m in range(-1, -NBWD - 1, -1):       # m = n-3
        n = m + 3
        u[m] = 82 * u[n - 1] + 82 * u[n - 2] - u[n]
    back[name] = u

t_m1 = (back['a'][-1], back['b'][-1], back['c'][-1])
t_m2 = (back['a'][-2], back['b'][-2], back['c'][-2])
t_m3 = (back['a'][-3], back['b'][-3], back['c'][-3])
print(f"  n=-1 triple: {t_m1}")
print(f"  n=-2 triple: {t_m2}")
print(f"  n=-3 triple: {t_m3}")
check("D8: n=-1 triple is (-9, 12, 10)", t_m1 == (-9, 12, 10))
check("D8: n=-2 triple is (-791, 1010, 812)", t_m2 == (-791, 1010, 812))
check("D8: n=-3 triple is (-65601, 83802, 67402)", t_m3 == (-65601, 83802, 67402))
check("D8: (-9)^3 + 12^3 == 10^3 + (-1)^(-1)",
      (-9) ** 3 + 12 ** 3 == 10 ** 3 - 1)
check("D8: 9^3 + 10^3 == 12^3 + 1 == 1729 (taxicab reading)",
      9 ** 3 + 10 ** 3 == 1729 and 12 ** 3 + 1 == 1729)

# ---------------------------------------------------------------- PART 3
print()
print("=" * 72)
print("PART 3: d_n = a_n^3 + b_n^3 - c_n^3 - (-1)^n on n = -12..40 (53 values)")
print("=" * 72)
all_zero = True
worst = None
for n in range(-NBWD, NFWD):
    an = back['a'][n]
    bn = back['b'][n]
    cn = back['c'][n]
    sgn = 1 if n % 2 == 0 else -1
    d = an ** 3 + bn ** 3 - cn ** 3 - sgn
    if d != 0:
        all_zero = False
        worst = (n, d)
    if -3 <= n <= 3:
        print(f"  n={n:3d}: (a,b,c)=({an},{bn},{cn})  d_n={d}")
check("d_n == 0 exactly for ALL n in [-12, 40] (53 consecutive integers)",
      all_zero, f"first failure {worst}" if worst else "")
print(f"  digit lengths at n=40: a has {len(str(back['a'][39+1]))} digits")

# ---------------------------------------------------------------- PART 4
print()
print("=" * 72)
print("PART 4: independent PROOF that d_n == 0 for all n in Z")
print("=" * 72)
# Characteristic polynomial of the recurrence: f(t) = t^3 - 82t^2 - 82t + 1.
# Step 1 (exact): f(t) = (t+1)(t^2 - 83t + 1).  Verify by integer poly mult.
lhs = [1, -82, -82, 1][::-1]          # high-to-low: t^3 -82t^2 -82t +1
prod = [0, 0, 0, 0]
# (t+1)*(t^2-83t+1): coefficients high-to-low
p1 = [1, 1]            # t + 1
p2 = [1, -83, 1]       # t^2 - 83t + 1
conv = [0] * 4
for i, ai in enumerate(p1):
    for j, bj in enumerate(p2):
        conv[i + j] += ai * bj
check("f(t) == (t+1)(t^2-83t+1) exactly", conv == [1, -82, -82, 1],
      f"product coeffs (high->low) = {conv}")
# Step 2: the quadratic factor has discriminant 83^2-4 = 6885, NOT a square,
# and 6885 > 0: two distinct REAL irrational roots mu, 1/mu (product = 1),
# both positive (sum 83 > 0, product 1 > 0), mu != 1 since f(1) = -162... check:
disc = 83 * 83 - 4
r = isqrt(disc)
check("disc(t^2-83t+1) = 6885 > 0 and not a perfect square",
      disc == 6885 and r * r != disc, f"isqrt(6885)={r}, {r}^2={r*r}")
val1 = 1 - 83 + 1
check("t=1 is not a root of the quadratic (so mu != 1, mu > 1 > 1/mu > 0)",
      val1 != 0, f"q(1) = {val1}")
val_m1 = 1 + 83 + 1
check("t=-1 is not a root of the quadratic (eigenvalue -1 is simple)",
      val_m1 != 0, f"q(-1) = {val_m1}")
print("""
PROOF (independent of the worker's 27x27 Cayley-Hamilton route):
  f has three DISTINCT nonzero roots: -1, mu, 1/mu  (mu = (83+9*sqrt(85))/2,
  real, mu > 1; note 6885 = 81*85).  Each extended sequence a,b,c is defined
  on ALL of Z by the homogeneous order-3 recurrence (leading AND trailing
  coefficients of f are 1, so it runs both directions over Z), hence each is
  an exponential polynomial  A1(-1)^n + A2 mu^n + A3 mu^(-n)  on all of Z
  (distinct roots => no polynomial-in-n terms).  Cubing, every term of
  a_n^3 + b_n^3 - c_n^3 has base a product of three roots, i.e. of the form
  (-1)^e * mu^j with j in {-3..3}; together with the lone (-1)^n term
  (base -1 = -mu^0, already in the list), d_n is an exponential polynomial
  with at most 14 DISTINCT NONZERO real bases  {+-mu^j : j = -3..3}
  (distinct because mu > 1 => the mu^j are 7 distinct positives, and a
  positive never equals a negative).  An exponential polynomial with N
  distinct nonzero bases vanishing at N CONSECUTIVE integers is identically
  zero: the N x N matrix (lambda_i^(n0+k)) = diag(lambda_i^n0) * Vandermonde
  is invertible.  PART 3 exhibits 53 >= 14 consecutive exact zeros
  (n = -12..40), so d_n = 0 for ALL n in Z.   QED (claims D7, D8).
""")
check("53 consecutive verified zeros >= 14 required by the Vandermonde bound",
      53 >= 14)

# consistency notes for the worker's route (not needed for this proof):
# det C = -(-1)^3*f(0)... for monic cubic t^3+c2t^2+c1t+c0, product of roots
# = -c0 = -1, so det C = -1 and (-1)^n = (product of roots)^n is a root power
# of the Kronecker cube; f(-1) = 0 also directly. Both consistent with above.
fm1 = (-1) ** 3 - 82 * 1 - 82 * (-1) + 1
check("f(-1) == 0 (so (-1)^n already satisfies the order-3 recurrence)",
      fm1 == 0, f"f(-1) = {fm1}")
f0 = 1
check("f(0) == 1 (trailing coefficient 1: backward extension is integral and "
      "canonical)", f0 == 1)

print()
print("=" * 72)
print(f"VERIFIER SUMMARY: {'ALL CHECKS PASS' if fails == 0 else str(fails) + ' FAILURES'}")
print("=" * 72)
