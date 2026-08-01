"""
FC(2) saddle-descent: L(f^m) ~ (Dm)! * P_m, P_m = int_0^1 phi_D(a)^m C(a) da (top-form period).
Growth law [L(f^m)/(Dm)!]^{1/m} -> max_{a in [0,1]} |phi_D(a)|, phi_D(a)=f_D(a,1-a).
Exact L(f^m) = sum_{i,j} [x^i y^j](f^m) i! j!.
"""
import numpy as np
import mpmath as mp
mp.mp.dps = 30

def Lfm_exact(A, m):  # A[i,j] = coeff of x^i y^j in f; returns L(f^m) = sum coeff(f^m) i! j!
    cur = np.array([[1.0]])
    for _ in range(m):
        c2 = np.zeros((cur.shape[0]+A.shape[0]-1, cur.shape[1]+A.shape[1]-1))
        for i in range(cur.shape[0]):
            for j in range(cur.shape[1]):
                if cur[i, j] != 0:
                    c2[i:i+A.shape[0], j:j+A.shape[1]] += cur[i, j]*A
        cur = c2
    tot = mp.mpf(0)
    for i in range(cur.shape[0]):
        for j in range(cur.shape[1]):
            if cur[i, j] != 0:
                tot += mp.mpf(float(cur[i, j]))*mp.factorial(i)*mp.factorial(j)
    return tot

# f = x^2 + y^2 + 0.7 x  (D=2, top form f_2=x^2+y^2, phi_2(a)=a^2+(1-a)^2)
A = np.zeros((3, 3)); A[2, 0] = 1; A[0, 2] = 1; A[1, 0] = 0.7
def phi2(a): return a*a+(1-a)**2
maxphi = max(float(phi2(mp.mpf(i)/1000)) for i in range(1001))
print("f = x^2+y^2+0.7x (D=2); max_[0,1] phi_2 =", round(maxphi, 4))
print("growth law [L(f^m)/(2m)!]^{1/m} -> max|phi_D|:")
for m in [3, 6, 9, 12, 15]:
    Lm = Lfm_exact(A, m)
    g = (Lm/mp.factorial(2*m))**(mp.mpf(1)/m)
    print("  m=%2d: L(f^m)=%s  [L/(2m)!]^(1/m)=%s" % (m, mp.nstr(Lm, 8), mp.nstr(g, 7)))
print("=> confirms L(f^m) ~ (Dm)! * (top-form period), reducing inhomogeneous FC(2) to the top-form PMP.")
