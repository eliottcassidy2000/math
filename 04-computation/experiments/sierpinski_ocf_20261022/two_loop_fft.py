# FFT version of the circulant two-loop predictor (two_loop_predictor.predict) for large n: towers S_k to k = 9,
# R_n to n = 10^5, and the closed-form limits.  mac-mini-2026-10-22-sierp.
import numpy as np, math

def predict_fft(n, conn):
    a = np.zeros(n); a[list(conn)] = 1
    x = -2 * np.imag(np.fft.fft(a)) / n
    g = x ** 2 / (1 - x ** 2); h = -x ** 2 / (1 + x ** 2); o = 1j * x / (1 + x ** 2)
    gr = np.real(np.fft.fft(g)) / n; hr = np.real(np.fft.fft(h)) / n; orr = np.real(np.fft.fft(o)) / n
    F8 = n * (2 * gr[0] * hr[0] - 1.5 * gr[0] ** 2)
    TH = n * np.sum(gr * (hr ** 2 - orr ** 2 - (2 / 3) * gr ** 2))
    D1 = np.sum(2 * x ** 2 / (1 - x ** 4)); D2 = np.sum(4 * x ** 2 * (1 + x ** 4) / (1 - x ** 4) ** 2)
    EL = (D2 + D1 ** 2 - D1) / 2
    return -n * (F8 + TH) - EL, n * F8, n * TH, EL

def tower_conn(n):
    out = []
    for z in range(1, n):
        t = z
        while t % 3 == 0: t //= 3
        if t % 3 == 1: out.append(z)
    return out

if __name__ == '__main__':
    print('towers S_k: k, n, y_pred, n F8, n Theta, E0')
    ys = []
    for k in range(2, 10):
        n = 3 ** k; y, f8, th, el = predict_fft(n, tower_conn(n)); ys.append((n, y))
        print(f"  k={k} n={n}: y_pred={y:.7f}  nF8={f8:+.6f}  nTheta={th:+.6f}  E0={el:.6f}", flush=True)
    # Richardson in 1/n (geometric n = 3^k)
    (n1, y1), (n2, y2), (n3, y3) = ys[-3:]
    r1 = (3 * y3 - y2) / 2; r0 = (3 * y2 - y1) / 2
    print(f"  Richardson (1/n): {r0:.7f}, {r1:.7f};  second order: {(9 * r1 - r0) / 8:.7f}")
    print('rotational R_n: n, y_pred')
    for n in (1001, 10001, 100001):
        y, f8, th, el = predict_fft(n, range(1, (n - 1) // 2 + 1))
        print(f"  n={n}: y_pred={y:.7f}  nF8={f8:+.6f}  nTheta={th:+.7f}  E0={el:.6f}")
    t, tp = math.tan(1), math.tanh(1)
    print(f"closed form for R_n: tan1 tanh1 + tan^2(1)/2 - 1 = {t * tp + t * t / 2 - 1:.7f};  nF8 -> {-2 * t * tp - 1.5 * t * t:.7f}, E0 -> {(t + 1 / math.cos(1) ** 2 + tp + 1 / math.cosh(1) ** 2 + (t + tp) ** 2 - t - tp) / 2:.7f}")
