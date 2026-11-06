# Independent audit of THM-4887, part 4: every numerical constant, recomputed with mpmath at 50 digits.
import math
from fractions import Fraction
import mpmath as mp
mp.mp.dps = 50
ln2 = mp.log(2)

def D1_over(n, Delta0=False):
    """(D1) at the (D2) split m = t, k = 2t + r, divided by 2^(4n/3); also return the 'remaining terms'."""
    t, r = divmod(n, 3); m = t; k = 2 * t + r
    C = mp.mpf(2) ** m; P = mp.mpf(2) ** (k - 1); lB = mp.log(P / C)
    s = mp.mpf(2) ** (mp.mpf(4 * n) / 3)
    if Delta0:
        main = 8 * P * P + P * C * C * (mp.mpf("33.74") + 12 * lB)
    else:
        main = 10 * P * P + P * C * C * (mp.mpf("35.74") + 12 * lB)
    rem = P * C + 3 * P * C * (1 + lB)
    return (main + rem) / s, main / s, rem / s

print("== (D2) rows: exact (D1)/2^(4n/3) at the (D2) split vs the stated rows and vs 1.39 n + 19 ==")
rows = {0: (mp.mpf("16.53"), mp.mpf("1.386")), 1: (mp.mpf("17.36"), mp.mpf("1.101")), 2: (mp.mpf("18.71"), mp.mpf("0.873"))}
for r in range(3):
    # exact slope and intercept of the main part: main/2^(4n/3) = a_r + b_r n
    b_r = mp.mpf(2) ** (-1 - mp.mpf(r) / 3) * 4 * ln2
    a_r = 10 * mp.mpf(2) ** (mp.mpf(2 * r) / 3 - 2) + mp.mpf(2) ** (-1 - mp.mpf(r) / 3) * (mp.mpf("35.74") + ln2 * (8 * r - 12))
    worst_row_gap = None; worst_final = None; remmax = 0; nmin = None
    for n in range(12, 3001):
        if n % 3 != r: continue
        if nmin is None: nmin = n
        tot, main, rem = D1_over(n)
        assert abs(main - (a_r + b_r * n)) < mp.mpf(10) ** -30
        remmax = max(remmax, rem)
        gap = tot - (rows[r][0] + rows[r][1] * n)      # > 0 means the stated row is NOT an upper bound
        worst_row_gap = gap if worst_row_gap is None else max(worst_row_gap, gap)
        fin = tot - (mp.mpf("1.39") * n + 19)
        worst_final = fin if worst_final is None else max(worst_final, fin)
    print(f" r={r}: main part = {mp.nstr(a_r, 8)} + {mp.nstr(b_r, 8)} n; remaining terms <= {mp.nstr(remmax, 6)} "
          f"(max at n = {nmin}); stated row {rows[r][0]} + {rows[r][1]} n; max over n<=3000 of (exact - row) = "
          f"{mp.nstr(worst_row_gap, 5)} ({'ROW IS AN UPPER BOUND' if worst_row_gap <= 0 else 'row UNDERSTATES the bound'}); "
          f"max (exact - (1.39n+19)) = {mp.nstr(worst_final, 6)}")
tot12, main12, rem12 = D1_over(12)
print(f" remaining terms at n = 12: {mp.nstr(rem12, 8)} (claimed <= 0.32); at n=13: {mp.nstr(D1_over(13)[2], 8)}; n=14: {mp.nstr(D1_over(14)[2], 8)}")
print(f" Delta = 0 form at n=12..14: {[mp.nstr(D1_over(n, True)[0], 6) for n in (12, 13, 14)]}")

print("\n== step 4: harmonic sum  sum_{beta<=B} (beta+1)(4C/beta+1)/beta <= 4C ln B + 10.58 C + B + 1 + ln B ==")
worst = -1e9
for e in range(0, 41):
    B = 2 ** e
    H1 = mp.harmonic(B); H2 = mp.zeta(2) - mp.zeta(2, B + 1)   # sum_{b<=B} 1/b^2
    for C in [1, 2, 4, 16, 2 ** 10, 2 ** 20]:
        lhs = 4 * C * H1 + B + 4 * C * H2 + H1
        rhs = 4 * C * mp.log(B) + mp.mpf("10.58") * C + B + 1 + mp.log(B)
        worst = max(worst, float(lhs - rhs))
print(f" max (lhs - rhs) over B = 2^0..2^40, C in {{1,..,2^20}}: {worst:.4f} (<= 0 means valid); 4(1 + pi^2/6) = {mp.nstr(4*(1+mp.pi**2/6), 8)}")

print("\n== sharpness lower bound 2^(2n-2)(D+1)/(2^n+D+1) at D = 2^floor(n/3), divided by 2^(4n/3) ==")
mins = [(n, float(mp.mpf(2) ** (2 * n - 2) * (2 ** (n // 3) + 1) / (2 ** n + 2 ** (n // 3) + 1) / mp.mpf(2) ** (mp.mpf(4 * n) / 3))) for n in range(12, 400)]
print(f" min over n = 12..399: {min(v for _, v in mins):.6f} (claimed >= 1/8 = 0.125); limit 2^(-8/3) = {float(mp.mpf(2)**(-mp.mpf(8)/3)):.6f}")
print(f" ratio upper/lower at n = 12..14: {[round(float((mp.mpf('1.39')*n+19)/v),1) for n, v in mins[:3]]}; at n=300: {round(float((mp.mpf('1.39')*300+19)/mins[288][1]),1)}")

print("\n== (D3) constants ==")
print(f" 4 (1.39 n + 19) = 5.56 n + 76: {Fraction('1.39')*4} n + {19*4}")
c1 = mp.mpf("5.56") / mp.mpf(2) ** (-mp.mpf(2) / 3) / ln2; c0 = 76 / mp.mpf(2) ** (-mp.mpf(2) / 3)
print(f" (X/2)^(2/3)/(5.56 log2 X + 76) = X^(2/3)/({mp.nstr(c1, 8)} ln X + {mp.nstr(c0, 8)})  (claimed <= 13 ln X + 121)")
# direct check on a grid including the worst points X -> 2^(n+1)^-
worst = 1e9
for n in range(12, 400):
    for frac in [0, 1e-9, 0.25, 0.5, 0.75, 0.999999999]:
        lx = (n + frac) * ln2                  # ln X with floor(log2 X) = n
        X23 = mp.e ** (2 * lx / 3)
        claimed = X23 / (13 * lx + 121)
        proved = mp.mpf(2) ** (mp.mpf(2 * n) / 3) / (mp.mpf("5.56") * n + 76)
        worst = min(worst, float(proved / claimed))
print(f" min over X in [4096, 2^400) of (2^(2n/3)/(5.56n+76)) / (X^(2/3)/(13 ln X + 121)) = {worst:.6f} (>= 1 means valid)")
# crossover with Session 150's sqrt(X)/4
cross = None
for e10 in range(120, 2000):
    lx = mp.mpf(e10) / 10 * ln2
    if mp.e ** (2 * lx / 3) / (13 * lx + 121) > mp.e ** (lx / 2) / 4:
        cross = mp.mpf(e10) / 10; break
print(f" X^(2/3)/(13 ln X + 121) first exceeds Session 150's sqrt(X)/4 at X = 2^{mp.nstr(cross, 5)}")
# the lower bound 2^(2n/3)/(5.56 n + 76) itself vs sqrt(2^n)/4 (same X)
print(f" 2^(2n/3)/(5.56n+76) at n = 12: {float(mp.mpf(2)**8/(5.56*12+76)):.3f}; first n with 2^(2n/3)/(5.56n+76) > 2^(n/2)/4: "
      f"{next(n for n in range(12, 200) if mp.mpf(2)**(mp.mpf(2*n)/3)/(mp.mpf('5.56')*n+76) > mp.mpf(2)**(mp.mpf(n)/2)/4)}")

print("\n== (D4) constants ==")
print(f" (pi^2/4) 2^(2/3) = {mp.nstr(mp.pi**2/4 * mp.mpf(2)**(mp.mpf(2)/3), 10)} (claimed <= 3.92)")
print(f" window: 1/(2(2^t+1)) >= 2^(-t-2) for all t >= 1: {all(Fraction(1, 2 * (2 ** t + 1)) >= Fraction(1, 2 ** (t + 2)) for t in range(1, 200))}")
print(f" 2^(4n/3)/2^floor(n/3) <= 2^(n+2/3): max exponent excess = {max(Fraction(4*n,3) - n//3 - n for n in range(12, 300))} (<= 2/3)")
x = mp.mpf(1) / 2
print(f" Fejer kernel lower bound: min over 0<|theta|<=1/(2(D+1)) of K(theta)/(4(D+1)/pi^2) for D=1..64 =",
      mp.nstr(min(mp.sin(mp.pi * (D + 1) * th) ** 2 / ((D + 1) * mp.sin(mp.pi * th) ** 2) / (4 * (D + 1) / mp.pi ** 2)
                  for D in range(1, 65) for th in [mp.mpf(j) / (2 * (D + 1) * 200) for j in range(1, 201)]), 8))

print("\n== upgrade: Delta = 0 constants (8 P^2, 33.74) at the (D2) split ==")
worst16 = max(D1_over(n, True)[0] - (mp.mpf("1.39") * n + 16) for n in range(12, 3001))
print(f" max over n = 12..3000 of (D1, Delta=0)/2^(4n/3) - (1.39 n + 16) = {mp.nstr(worst16, 6)} (<= 0: E_n <= (1.39 n + 16) 2^(4n/3))")
print(" general m (Delta = 2^m, n/3 <= m <= (n-4)/2): (D1)/2^(n+m) = 17.87 + 6 ln2 (n-2m-1) + 2.5*2^(n-3m) + 2^(-m-1)(4 + 3 ln B)")
for n in (30, 60, 120):
    for m in (n // 3, (n + (n // 3)) // 2 if False else (n - 4) // 2):
        C = mp.mpf(2) ** m; P = mp.mpf(2) ** (n - m - 1); lB = mp.log(P / C)
        val = (P * C + 10 * P * P + P * C * C * (mp.mpf("35.74") + 12 * lB) + 3 * P * C * (1 + lB)) / mp.mpf(2) ** (n + m)
        print(f"   n={n} m={m}: (D1)/2^(n+m) = {mp.nstr(val, 6)}")
