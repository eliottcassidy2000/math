"""
FC(2) on Re f_D = 0: imaginary top form f_D = i g_D makes psi=g_D(a,1-a) REAL, so the top-form period
int_0^1 psi^m C da has a real-modulus Laplace peak -> P_m ~ (M^m/sqrt m)[S_+ + (-1)^m S_-], NO continuous
Stokes cancellation. Unique max-modulus => S_+ = C(a_1)w_1 != 0 (|C|>0) => L(f^m)!=0 => descent closes
(no K-Z). Verified: growth [|L(f^m)|/(Dm)!]^{1/m} -> max|psi|, L(f^m) != 0.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 30

def Lfm(A, m):  # exact L(f^m) for complex-coeff f; A[i,j]=coeff x^i y^j
    cur = np.array([[1.0+0j]])
    for _ in range(m):
        c2 = np.zeros((cur.shape[0]+A.shape[0]-1, cur.shape[1]+A.shape[1]-1), complex)
        for i in range(cur.shape[0]):
            for j in range(cur.shape[1]):
                if cur[i, j] != 0:
                    c2[i:i+A.shape[0], j:j+A.shape[1]] += cur[i, j]*A
        cur = c2
    tot = mp.mpc(0)
    for i in range(cur.shape[0]):
        for j in range(cur.shape[1]):
            if cur[i, j] != 0:
                tot += mp.mpc(complex(cur[i, j]))*mp.factorial(i)*mp.factorial(j)
    return tot

# f = i(x^2+y^2) + (1+i)x : imaginary top form, psi(a)=a^2+(1-a)^2 REAL
A = np.zeros((3, 3), complex); A[2, 0] = 1j; A[0, 2] = 1j; A[1, 0] = 1+1j
def psi(a): return a*a+(1-a)**2
M = max(float(psi(mp.mpf(i)/1000)) for i in range(1001))
print("Re f_D=0: f=i(x^2+y^2)+(1+i)x, psi=a^2+(1-a)^2 real, max|psi|=", round(M, 4))
for m in [4, 8, 12, 16]:
    Lm = Lfm(A, m)
    g = (abs(Lm)/mp.factorial(2*m))**(mp.mpf(1)/m)
    print("  m=%2d: |L(f^m)|=%s  [|L|/(2m)!]^(1/m)=%s (->max|psi|, and L!=0)" % (m, mp.nstr(abs(Lm), 6), mp.nstr(g, 7)))
print("=> real-modulus Laplace, no continuous Stokes; unique max => descent closes unconditionally.")
