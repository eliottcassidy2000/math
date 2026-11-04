"""aud_mrb.py - independent numerical audit of section 3 of the S76 note (the two MRB notes).

Checks, with mpmath (own code):
 1. C_MRB = sum_{n>=1} (-1)^n (n^(1/n) - 1) to 60 digits (alternating-series acceleration, two methods);
 2. the Abel-Plana integral C = int_0^oo Im((1+it)^(1/(1+it))) / sinh(pi t) dt;
 3. M_k = int_0^oo t^(2k+1)/sinh(pi t) dt = 2 (2k+1)! (1 - 2^(-2k-2)) zeta(2k+2) / pi^(2k+2): M_0 = 1/4, M_1 = 1/8;
 4. exact Taylor coefficients of h(t) = Im exp(Log(1+it)/(1+it)) (sympy, rational);
 5. R = C - 3/16 and the remainder integral int (h - t + t^3/2)/sinh(pi t);
 6. residues Res G = -a_k M_k / 4 at z = k + 3/2 (k = 0, 1) from the stated prefactor identity;
 7. the W_(-1) inversion N(a_n) = n and the W_0 mirrors r_n (r_n^(1/r_n) = n^(1/n), r_n < e);
 8. the Euler-transformed interval series at M = 50, 51, 100, 101, 150, 151, 200, 201 terms, against the note's table;
 9. mpmath nsum on the paired series with its default and with explicit methods (the S76 'cautionary echo');
10. 3/16 = 2 - 29/16 and the 3-cycle of x^2 - 29/16.
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 80

# 1. the constant
f = lambda n: mp.power(n, mp.mpf(1) / n) - 1
C_alt = mp.nsum(lambda n: (-1) ** int(n) * f(n), [1, mp.inf], method='alternating')
# Levin transform as an independent second method
C_lev = mp.nsum(lambda n: (-1) ** int(n) * f(n), [2, mp.inf], method='levin')   # n = 1 term is 0
print("1. C_MRB (alternating)    =", mp.nstr(C_alt, 60))
print("   C_MRB (levin)          =", mp.nstr(C_lev, 60))
C = C_alt
ref_note = mp.mpf('0.1878596424620671202485179340542732300559')
print("   attachment value 0.1878596424620671202485179340542732300559, |diff| =", mp.nstr(abs(C - ref_note), 5))

# 2. Abel-Plana integral
mp.mp.dps = 45


def h(t):
    w = 1 + 1j * t
    return mp.im(mp.exp(mp.log(w) / w))


I = mp.quad(lambda t: h(t) / mp.sinh(mp.pi * t), [0, 1, 4, 16, 64, mp.inf])
print("2. Abel-Plana integral    =", mp.nstr(I, 40), " |I - C| =", mp.nstr(abs(I - C), 5))

# 3. moments
for k in range(3):
    Mk_int = mp.quad(lambda t: t ** (2 * k + 1) / mp.sinh(mp.pi * t), [0, 1, 4, 16, mp.inf])
    Mk_formula = 2 * mp.factorial(2 * k + 1) * (1 - mp.power(2, -2 * k - 2)) * mp.zeta(2 * k + 2) / mp.pi ** (2 * k + 2)
    print("3. M_%d integral = %s  formula = %s  identify = %s" % (k, mp.nstr(Mk_int, 30), mp.nstr(Mk_formula, 30), mp.identify(Mk_formula)))

# 4. exact Taylor coefficients of h
t = sp.symbols('t', real=True)
w = 1 + sp.I * t
expr = sp.exp(sp.log(w) / w)
ser = sp.series(expr, t, 0, 12).removeO()
ser = sp.expand(ser)
coeffs = [sp.nsimplify(sp.im(sp.expand(ser).coeff(t, k))) for k in range(12)]
print("4. Im-part Taylor coefficients of h, t^0..t^11:", coeffs)

# 5. R
mp.mp.dps = 45
R = C - mp.mpf(3) / 16
Rint = mp.quad(lambda t: (h(t) - t + t ** 3 / 2) / mp.sinh(mp.pi * t), [0, 1, 4, 16, 64, mp.inf])
print("5. R = C - 3/16           =", mp.nstr(R, 40))
print("   remainder integral     =", mp.nstr(Rint, 40), " |diff| =", mp.nstr(abs(R - Rint), 5))
print("   S76 note value R = 0.000359642462067120248517934..., matches to", mp.nstr(abs(R - mp.mpf('0.000359642462067120248517934')), 5))

# 6. residues from Res G = -a_k M_k / 4
a = [1, sp.Rational(-1, 2)]
M = [sp.Rational(1, 4), sp.Rational(1, 8)]
print("6. Res G at z = 3/2, 5/2:", [-a[k] * M[k] / 4 for k in range(2)], " contributions -4*Res:", [4 * a[k] * M[k] / 4 for k in range(2)], " sum", sum(a[k] * M[k] for k in range(2)))

# 7. Lambert W branches
mp.mp.dps = 40
print("7. W_(-1) inversion N(a_n) and W_0 mirrors r_n:")
for n in [2, 3, 4, 5, 6, 8, 10, 20]:
    z = mp.log(n) / n
    Nn = -mp.lambertw(-z, -1) / z
    rn = -mp.lambertw(-z, 0) / z
    print("   n=%2d  N(a_n) = %s   r_n = %s   check r_n^(1/r_n) - n^(1/n) = %s" % (n, mp.nstr(mp.re(Nn), 20), mp.nstr(mp.re(rn), 12), mp.nstr(mp.re(rn) ** (1 / mp.re(rn)) - mp.power(n, mp.mpf(1) / n), 3)))

# 8. Euler transform of the interval areas, as in the attachment's code (wp = 110 guard digits)
for wp in (110, 200):
    mp.mp.dps = wp
    s2 = mp.sqrt(2) - 1
    print("8. Euler transform, working precision %d digits:" % wp)
    for Mterms in (50, 51, 100, 101, 150, 151, 200, 201):
        d = [mp.exp(mp.log(n) / n) - 1 for n in range(3, Mterms + 3)]
        s = mp.mpf(0)
        for j in range(Mterms):
            s += d[0] / mp.power(2, j + 1)
            if j < Mterms - 1:
                d = [d[i] - d[i + 1] for i in range(len(d) - 1)]
        val = s2 - s
        mp.mp.dps = 80
        err = abs(val - C)
        mp.mp.dps = wp
        print("   M=%3d  %s   |error| = %s" % (Mterms, mp.nstr(val, 58), mp.nstr(err, 3)))
table = {51: '0.187859642462067120238962868492', 101: '0.187859642462067120248517934054273229264',
         151: '0.187859642462067120248517934054273230055903094900138656', 201: '0.187859642462067120248517934054273230055903094900138786'}
mp.mp.dps = 80
print("   attachment table values minus C:", {k: mp.nstr(mp.mpf(v) - C, 3) for k, v in table.items()})

# 9. the paired series and mpmath's defaults
mp.mp.dps = 30
paired = lambda k: f(2 * k + 1) - f(2 * k + 2)
for meth in (None, 'r+s', 'e', 'a', 'l'):
    try:
        if meth is None:
            v = mp.nsum(paired, [1, mp.inf])
        else:
            v = mp.nsum(paired, [1, mp.inf], method=meth)
        print("9. paired series sum_(k>=1) (a_(2k+1) - a_(2k+2)), method %-5s -> sqrt2 - 1 - sum = %s" % (meth, mp.nstr(mp.sqrt(2) - 1 - v, 15)))
    except Exception as ex:
        print("9. method", meth, "failed:", ex)

# 10. numerology
c = sp.Rational(-29, 16)
x = sp.Rational(-1, 4)
orb = [x]
for _ in range(3):
    x = x ** 2 + c
    orb.append(x)
print("10. 2 - 29/16 =", 2 - sp.Rational(29, 16), "; orbit of -1/4 under x^2 - 29/16:", orb)
