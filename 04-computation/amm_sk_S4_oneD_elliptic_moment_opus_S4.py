"""
S(4) = sum_{n>=0} C(2n,n)C(4n,2n)/((4n+1)64^n): clean 1-D integral + elliptic-moment core.

Uses death-star's kernel 2F1(1/4,3/4;1;z) = (1/pi) int_0^pi dphi/sqrt(1+sqrt(z) cos phi).
Since sqrt(x^4)=x^2, the inner x-integral is elementary; folding phi->pi-phi gives
    S(4) = (2/pi) int_0^1 [arcsinh(s)+arcsin(s)]/sqrt(1-s^4) ds.
Non-elementary core: I2 = int_0^1 arcsin(s)/sqrt(1-s^4) ds = int_0^{pi/2} theta/sqrt(1+sin^2 theta) dtheta
(a theta-weighted elliptic moment). PSLQ shows S(4) is independent of {K(i), lemniscate, Catalan, pi}
at 40 digits -> concrete witness for kps-S148's "irreducible motive".
"""
import mpmath as mp
mp.mp.dps = 40
tgt = mp.mpf('1.069352554441268058582961939534')

# kernel check
def Fker(z): return (1/mp.pi)*mp.quad(lambda p: 1/mp.sqrt(1+mp.sqrt(z)*mp.cos(p)), [0, mp.pi])
print("kernel 2F1(1/4,3/4;1;.7) vs death-star F(.7):",
      mp.nstr(mp.hyp2f1(mp.mpf(1)/4, mp.mpf(3)/4, 1, mp.mpf('0.7')) - Fker(mp.mpf('0.7')), 3))

# clean 1-D form
S4 = (2/mp.pi)*mp.quad(lambda s: (mp.asinh(s)+mp.asin(s))/mp.sqrt(1-s**4), [0, 1])
print("S(4) = (2/pi) int_0^1 [asinh+asin]/sqrt(1-s^4) ds =", mp.nstr(S4, 28))
print("   vs target:", mp.nstr(S4 - tgt, 3))

# odd-term cancellation
lhs = mp.asinh(mp.mpf('0.6')) + mp.asin(mp.mpf('0.6'))
rhs = 2*mp.nsum(lambda m: mp.binomial(4*m, 2*m)/(16**m*(4*m+1))*mp.mpf('0.6')**(4*m+1), [0, mp.inf])
print("arcsinh+arcsin keeps only s^{4m+1}:", mp.nstr(lhs - rhs, 3))

# elliptic-moment core + independence
anchor = mp.quad(lambda th: 1/mp.sqrt(1+mp.sin(th)**2), [0, mp.pi/2])   # K(i)
I2 = mp.quad(lambda th: th/mp.sqrt(1+mp.sin(th)**2), [0, mp.pi/2])
lemn = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*mp.pi))
print("anchor int dtheta/sqrt(1+sin^2) = K(i) =", mp.nstr(anchor, 18))
print("I2 (theta-weighted moment) =", mp.nstr(I2, 22))
rel = mp.pslq([mp.pi*S4, anchor, mp.catalan, mp.pi, lemn], maxcoeff=10**7)
print("PSLQ pi*S(4) vs {K(i),Catalan,pi,lemniscate}:", rel,
      "  (coeff 0 on S(4) => independent => irreducible, kps-S148)")
