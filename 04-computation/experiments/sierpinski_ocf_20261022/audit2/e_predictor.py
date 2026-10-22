# Audit2 (D): independent implementation of the HYP-9365 two-loop predictor
#   y(T) = n (1 - kappa/rho),  kappa = H 2^(n-1)/n!,  rho = det(nI - S)/|det(inI - S)|,
#   y_pred = -n X_1 - (D_2 + D_1^2 - D_1)/2,  X_1 = F8 + Theta - d^T G d  (THM-4757 (4), dense matrices),
#   D_1 = sum_k 2x^2/(1-x^4),  D_2 = sum_k [2x^2/(1-x^2)^2 + 2x^2/(1+x^2)^2],  x_k = Im eig(S)/n.
# Dense version for any regular tournament; FFT version for circulants (cross-checked against the dense one).
import math
import numpy as np
from mpmath import mp, mpf, factorial as mfac
mp.dps = 50

def two_loop_dense(S):
    S = np.asarray(S, dtype=float); n = S.shape[0]
    Sh = S / n; I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2)
    H = S2 @ np.linalg.inv(I - S2)
    O = Sh @ np.linalg.inv(I - S2)
    g, h = np.diag(G), np.diag(H)
    F8 = float(np.sum(2 * g * h - 1.5 * g * g))
    TH = float(np.sum(G * (H * H - O * O - (2.0 / 3.0) * G * G)))
    d = h - g
    DB = float(-d @ G @ d)
    ev = np.linalg.eigvals(Sh)
    x = ev.imag
    D1 = float(np.sum(2 * x ** 2 / (1 - x ** 4)))
    D2 = float(np.sum(2 * x ** 2 / (1 - x ** 2) ** 2 + 2 * x ** 2 / (1 + x ** 2) ** 2))
    E0 = (D2 + D1 * D1 - D1) / 2
    # rho exactly as a determinant ratio (float, via log-det)
    s1, l1 = np.linalg.slogdet(n * I - S)
    l2 = np.log(np.abs(np.linalg.det(1j * n * I - S)))
    rho = float(np.exp(l1 - l2)) if s1 > 0 else float('nan')
    yp = -n * (F8 + TH + DB) - E0
    return dict(rho=rho, y_pred=yp, nF8=n * F8, nTH=n * TH, nDB=n * DB, E0=E0, D1=D1, D2=D2)

def two_loop_circ(n, conn):
    conn = sorted({c % n for c in conn})
    k = np.arange(n)
    if n <= 64:
        x = np.zeros(n)
        for s in conn: x += 2 * np.sin(2 * np.pi * k * s / n) / n
    else:                                   # x_k = (2/n) sum_s sin(2 pi k s/n) = -(2/n) Im fft(indicator)_k
        a = np.zeros(n); a[conn] = 1
        x = -2 * np.imag(np.fft.fft(a)) / n
    g = x ** 2 / (1 - x ** 2); h = -x ** 2 / (1 + x ** 2); o = 1j * x / (1 + x ** 2)
    # first rows: M_(0,w) = (1/n) sum_k m_k omega^(-k w)  = fft(m)/n
    G0 = np.fft.fft(g) / n; H0 = np.fft.fft(h) / n; O0 = np.fft.fft(o) / n
    assert np.max(np.abs(G0.imag)) < 1e-12 and np.max(np.abs(O0.imag)) < 1e-12
    G0, H0, O0 = G0.real, H0.real, O0.real
    F8 = n * (2 * G0[0] * H0[0] - 1.5 * G0[0] ** 2)
    TH = n * np.sum(G0 * (H0 * H0 - O0 * O0 - (2.0 / 3.0) * G0 * G0))
    D1 = np.sum(2 * x ** 2 / (1 - x ** 4))
    D2 = np.sum(2 * x ** 2 / (1 - x ** 2) ** 2 + 2 * x ** 2 / (1 + x ** 2) ** 2)
    E0 = (D2 + D1 * D1 - D1) / 2
    logrho = 0.5 * np.sum(np.log((1 + x ** 2) / (1 - x ** 2)))
    return dict(rho=float(np.exp(logrho)), y_pred=float(-n * (F8 + TH) - E0), nF8=float(n * F8), nTH=float(n * TH), nDB=0.0,
                E0=float(E0), D1=float(D1), D2=float(D2))

def rho_circ_exact(n, conn):
    """rho via mpmath product over the spectrum (high precision)"""
    conn = sorted({c % n for c in conn})
    r = mpf(1)
    for k in range(n):
        x = sum(2 * mp.sin(2 * mp.pi * k * s / n) for s in conn) / n
        r *= (1 + x * x) / (1 - x * x)
    return mp.sqrt(r)

def y_exact(H, n, rho):
    kap = mpf(H) * mpf(2) ** (n - 1) / mfac(n)
    return float(n * (1 - kap / mpf(rho))), kap

def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])
