import mpmath as mp, numpy as np
mp.mp.dps = 25
b = mp.sqrt(3)/2
def Lchi(p): return mp.power(3, -p)*(mp.zeta(p, mp.mpf(1)/3) - mp.zeta(p, mp.mpf(2)/3))
# Riesz/Epstein: E_{t^-p}(A) = sum_{a != 0} |a|^{-2p} = sum_D r(D) (b/D)^p = 6 b^p zeta(p) L(p, chi_-3)
c = np.sqrt(2/np.sqrt(3)); w = np.exp(1j*np.pi/3)
K = 1500
j = np.arange(-K, K+1)
J, L = np.meshgrid(j, j)
Z = c*(J + L*w); R2 = np.abs(Z)**2; R2[K, K] = np.inf
for p in (2, 3):
    direct = np.sum(R2**(-p))
    # tail beyond the hexagon of radius ~ K*c*sqrt(3)/2 : approximate by integral of 2*pi*r*r^{-2p} from Rin
    Rin = K*c*np.sqrt(3)/2
    mask_out = R2 > Rin**2
    direct_in = np.sum(np.where(mask_out, 0.0, R2**(-p)))
    tail = 2*np.pi*Rin**(2-2*p)/(2*p-2)
    formula = 6*b**p*mp.zeta(p)*Lchi(p)
    print(f"p={p}: lattice sum (inner disk + tail integral) = {direct_in + tail:.12f};  6 b^p zeta(p) L(p,chi_-3) = {mp.nstr(formula, 14)}")
# Gaussian at alpha=1 (self-dual point): theta_A(1) - 1 = sum_D r(D) exp(-2 pi D / sqrt3)
def chi3(d): return 0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
def r(D): return 6*sum(chi3(d) for d in range(1, D+1) if D % d == 0)
th = mp.fsum(r(D)*mp.exp(-2*mp.pi*D/mp.sqrt(3)) for D in range(1, 200))
print("theta_A(1)-1 via r(D) =", mp.nstr(th, 15))
# BHS constant
Clog = 2*mp.log(2) + mp.log(mp.mpf(2)/3)/2 + 3*mp.log(mp.sqrt(mp.pi)/mp.gamma(mp.mpf(1)/3))
print("BHS C_log = 2log2 + (1/2)log(2/3) + 3log(sqrt(pi)/Gamma(1/3)) =", mp.nstr(Clog, 15))
# Chowla-Selberg check: |eta(rho)| with rho = e^{2 pi i/3}; known |eta(rho)|^... = 3^{1/8} Gamma(1/3)^{3/2} / (2 pi)
rho = mp.e**(2j*mp.pi/3)
q = mp.e**(2j*mp.pi*rho)
eta = mp.e**(2j*mp.pi*rho/24)*mp.nprod(lambda n: 1 - q**n, [1, mp.inf])
print("|eta(rho)| =", mp.nstr(abs(eta), 15), " vs 3^(1/8) Gamma(1/3)^(3/2)/(2 pi) =", mp.nstr(mp.power(3, mp.mpf(1)/8)*mp.gamma(mp.mpf(1)/3)**1.5/(2*mp.pi), 15))
