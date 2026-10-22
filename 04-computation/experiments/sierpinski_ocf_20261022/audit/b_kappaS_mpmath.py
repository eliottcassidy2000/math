# Independent high-precision value of kappa_S = prod_{i>=0} ((3^(2i+1)+1)/(3^(2i+1)-1))^(3^i)
from mpmath import mp, mpf, log, exp, atanh, cosh, cos, e
mp.dps = 50
s = mpf(0)
for i in range(0, 60):
    q = mpf(3) ** (2 * i + 1)
    s += mpf(3) ** i * log((q + 1) / (q - 1))
print("kappa_S (log-sum, 60 levels)       =", exp(s))
s2 = sum(2 * mpf(3) ** i * atanh(mpf(3) ** (-(2 * i + 1))) for i in range(60))
print("kappa_S (2 sum 3^i artanh 3^-(2i+1)) =", exp(s2))
print("THM-4755 printed value                = 2.79150988819834634766")
print("rho(S_k), k=1..8:", [mp.nstr(exp(sum(mpf(3) ** i * log((mpf(3) ** (2*i+1) + 1) / (mpf(3) ** (2*i+1) - 1)) for i in range(k))), 15) for k in range(1, 9)])
print("cosh1/cos1 =", cosh(1) / cos(1), " e =", e, " ratio cosh1/(e cos1) =", cosh(1) / (e * cos(1)))
