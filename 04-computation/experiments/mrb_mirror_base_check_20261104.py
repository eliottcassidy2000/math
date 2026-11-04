"""Re-check of the numbers quoted in the owner's MRB PDFs (read as untrusted data; nothing executed from them).
r_n = exp(-W0(-ln n / n)) is the mirror of n under x^(1/x) (r_n^(1/r_n) = n^(1/n), r_n < e < n for n >= 3; r_2 = 4).
B = sum_{n>=1} (-1)^n (r_n - 1) with r_1 = 1, r_2 = 4. Also checks B = 3 - length(E), E = U_{k>=1} [r_(2k+2), r_(2k+1)],
the MRB constant C = sum (-1)^n (n^(1/n) - 1), and r_n = e^(T(ln n / n)) with the tree function T(z) = -W0(-z) and its
Cayley series sum_{m>=1} m^(m-1) z^m / m!.
kind-pasteur-2026-11-04-S10."""
import mpmath as mp
mp.mp.dps = 40

def r(n):
    if n == 1: return mp.mpf(1)
    if n == 2: return mp.mpf(4)
    return mp.exp(-mp.lambertw(-mp.log(n) / n, 0).real)

for n in (3, 5, 10):
    x = r(n)
    print(f"r_{n} = {mp.nstr(x, 12)}   check r^(1/r) - n^(1/n) = {mp.nstr(x ** (1 / x) - mp.mpf(n) ** (mp.mpf(1) / n), 3)}")
print("r_4 =", mp.nstr(r(4), 15), "(rational mirror; r_2 = 4 by definition)")
B = mp.nsum(lambda n: (-1) ** int(n) * (r(int(n)) - 1), [1, mp.inf], method="alternating")
C = mp.nsum(lambda n: (-1) ** int(n) * (mp.mpf(int(n)) ** (mp.mpf(1) / int(n)) - 1), [1, mp.inf], method="alternating")
print("B =", mp.nstr(B, 25), " (PDF: 2.096791124391526125725030675)")
print("C =", mp.nstr(C, 25), " (MRB constant: 0.1878596424620671202485179)")
# length(E) is a slowly convergent non-alternating sum: use Euler-Maclaurin on the real-argument mirror
# (the default nsum method gives a wrong 4th decimal here).
rr = lambda x: mp.exp(-mp.lambertw(-mp.log(x) / x, 0).real)
lenE = mp.nsum(lambda k: rr(2 * k + 1) - rr(2 * k + 2), [1, mp.inf], method="euler-maclaurin")
print("3 - length(E) =", mp.nstr(3 - lenE, 22), " (Euler-Maclaurin; equals B)")
z = mp.log(3) / 3
T = -mp.lambertw(-z, 0).real
cay = mp.nsum(lambda m: mp.mpf(int(m)) ** (int(m) - 1) * z ** int(m) / mp.factorial(int(m)), [1, mp.inf])
print("tree function at ln3/3: -W0 =", mp.nstr(T, 20), " Cayley series =", mp.nstr(cay, 20), " e^T =", mp.nstr(mp.exp(T), 15), "= r_3")
