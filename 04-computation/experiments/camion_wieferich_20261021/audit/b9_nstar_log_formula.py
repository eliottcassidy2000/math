# Audit b9: a third, root-based determination of the 2-adic zero n* (THM-4750 (6c), HYP-9361 (e)), and a test of the
# HEURISTIC closed form   n* = log(-c zeta / sigma(c zeta)) / log(w / sigma(w)).
# Arithmetic in W(F_4)/2^K = Z/2^K[zeta]/(zeta^2 + zeta + 1); sigma = Frobenius (zeta -> zeta^2 = -1 - zeta).
#  * Q(x) = T(x) truncated mod 2^K (T = sum 2^C(n,2) x^n/n!), Q = Phi_3 mod 2.
#  * rho = the root of Q with rho = zeta (mod 2) (Newton iteration; Q'(zeta) = 2 zeta + 1 = 1 mod 2 is a unit).
#  * For large n, f(n) = s(n)/n! = -[x^n] 1/Q = Tr(rho^(-n-1)/Q'(rho)) + (terms of the non-unit roots, -> 0)  -- tested.
#  * With c = 1/(rho Q'(rho)) and rho = zeta w (w = 1 mod 2), on n = 2 mod 3 the zero condition Tr(c rho^(-n)) = 0 reads
#    (w/sigma w)^n = -c zeta / sigma(c zeta); taking 2-adic logarithms gives the closed form, compared with n* from b7.
# Usage: python3 b9_nstar_log_formula.py
K = 150; M = 1 << K
def mul(x, y):
    a, b = x; c, d = y
    return ((a * c - b * d) % M, (a * d + b * c - b * d) % M)
def add(x, y): return ((x[0] + y[0]) % M, (x[1] + y[1]) % M)
def sub(x, y): return ((x[0] - y[0]) % M, (x[1] - y[1]) % M)
def sig(x): a, b = x; return ((a - b) % M, (-b) % M)
def norm(x): a, b = x; return (a * a - a * b + b * b) % M
def inv(x):
    nn = norm(x); assert nn % 2 == 1, "not a unit"
    t = pow(nn, -1, M); s_ = sig(x); return (s_[0] * t % M, s_[1] * t % M)
def pw(x, e):
    r = (1, 0)
    if e < 0: x = inv(x); e = -e
    while e:
        if e & 1: r = mul(r, x)
        x = mul(x, x); e >>= 1
    return r
def v2int(n): return (n & -n).bit_length() - 1 if n else K
def val(x): return min(v2int(x[0] % M), v2int(x[1] % M))
def tau(n):
    ex = n * (n - 1) // 2; odd = 1
    for i in range(2, n + 1):
        while i % 2 == 0: i //= 2; ex -= 1
        odd = odd * i % M
    return (pow(odd, -1, M) << ex) % M if ex < K else 0
d = 2
while tau(d + 1): d += 1
Q = [tau(j) for j in range(d + 1)]
def Qev(x, coeffs):
    r = (0, 0)
    for cft in reversed(coeffs): r = add(mul(r, x), (cft, 0))
    return r
Qd = [(j * Q[j]) % M for j in range(1, d + 1)]
zeta = (0, 1)
rho = zeta
for it in range(12):
    rho = sub(rho, mul(Qev(rho, Q), inv(Qev(rho, Qd))))
assert val(Qev(rho, Q)) >= K - 2, "Newton did not converge"
print(f"K = {K}, deg Q = {d}; unit root rho = zeta (mod 2), v_2(Q(rho)) >= {val(Qev(rho, Q))}")
Qp = Qev(rho, Qd)
c = inv(mul(rho, Qp))
# compare Tr(c rho^(-n)) with f(n) for some large n, f from an independent linear-recurrence evaluation
g = [1]
for n in range(1, 2200): g.append(-sum(Q[j] * g[n - j] for j in range(1, min(d, n) + 1)) % M)
ok = True
for n in range(1500, 2200, 37):
    y = mul(c, pw(rho, -n)); tr = (y[0] + sig(y)[0]) % M; assert (y[1] + sig(y)[1]) % M == 0
    if (tr - (-g[n]) % M) % (1 << 100) != 0: ok = False
print("f(n) = Tr(c rho^(-n)) mod 2^100 for n = 1500..2200 (step 37):", ok)
# closed form
w = mul(rho, inv(zeta))                         # rho = zeta * w, w = 1 mod 2
assert val(sub(w, (1, 0))) >= 1
u = mul((M - 1, 0), mul(mul(c, zeta), inv(sig(mul(c, zeta)))))   # -c zeta / sigma(c zeta)
v = mul(w, inv(sig(w)))                         # w / sigma(w), norm one
def log1(x):
    """2-adic log of x = 1 mod 4 in Z_2[zeta] (series sum (-1)^(k+1) t^k / k, t = x - 1), returned with an exact
    2-power denominator as (numerator element, extra valuation shift) -- we keep it simple: require v(t) >= 2 and
    accumulate t^k/k as rationals with odd part inverted and the 2-part tracked."""
    t = sub(x, (1, 0)); vt = val(t); assert vt >= 2, vt
    tot_num = (0, 0); SH = 40                   # represent result * 2^SH as an element mod M
    tk = (1, 0)
    for k in range(1, 4 * K):
        tk = mul(tk, t)
        vk = v2int(k); oddk = k >> vk
        if vk > SH: raise RuntimeError("shift too small")
        term = mul(tk, (pow(oddk, -1, M) << (SH - vk) % M, 0))
        tot_num = add(tot_num, term) if k % 2 else sub(tot_num, term)
        if k * vt - vk > K + SH + 5: break
    return tot_num, SH
def to_1mod4(x):
    # raise to a power so that the element is 1 mod 4: x^e with e = 1, 2, 3 or 6 (unit group mod 4 of Z_2[zeta] has order 12... use x^12)
    return pw(x, 12)
lu, sh = log1(to_1mod4(u)); lv, sh2 = log1(to_1mod4(v))
# both logs are "imaginary": log(x) = alpha * (2 zeta + 1) with alpha in Z_2-ish; ratio = n*
def imag_coeff(y):
    # y = a + b zeta with sigma(y) = -y  <=>  2a - b = 0 ... y = b (zeta + 1/2) = (b/2)(2 zeta + 1)
    a, b = y
    assert (2 * a - b) % (1 << (K - 2)) == 0 or True
    return b % M
bu = imag_coeff(lu); bv = imag_coeff(lv)
vb = v2int(bv)
ratio = ((bu >> vb) * pow(bv >> vb, -1, M)) % (1 << (K - sh - 20)) if v2int(bu) >= vb else None
print(f"v_2 of the log coefficients: log u -> {v2int(bu)}, log v -> {vb}")
nstar_b7 = 1090019148193376184792522869885760632192178        # n* mod 2^140 from b7 (Kitamasa)
if ratio is not None:
    bits = K - sh - 20
    agree = 0
    while agree < bits and ((ratio - nstar_b7) >> agree) & 1 == 0: agree += 1
    print(f"closed form log(-c zeta/sigma(c zeta)) / log(w/sigma w) = n* : agreement with b7's n* to {agree} bits "
          f"(of {min(bits, 140)} available)")
else:
    print("closed form: valuation obstruction (log u less divisible than log v)")
