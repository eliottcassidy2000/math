"""Exact J-series of the resonant coefficient of the 3-adic Syracuse law (session collatz-necklace-20260929, part 2).

mu_hat_h(2^s) = sum_J c_J with (Lemma R'', exact)  c_J = mu_hat_J(1) * W_J,  W_J = sum_{kappa>=0} 2^-kappa Zhat_{J+1}(kappa),
Zhat_j(kappa) = sum over top words (a_j..a_h), s - T_j = kappa >= 0, of 2^-T_j prod_{i>=j} omega_i(s - T_i);
c_0 = sum_kappa Zhat_1(kappa); c_h = 2^-s mu_hat_h(1).  |W_J| <= mass_h(J) = 2^-s C(s, h-J).
So the bottom factor of every c_J is the frequency-1 Fourier coefficient mu_hat_J(1) = E e(Y_J/3^J) (hypothesis H1),
not the weighted negative-power norm of Lemma R'.  The window recursion f_n(k) = sum_a 2^-a omega_n(k-a) f_{n-1}(k-a)
(S22's closed family, truncation a <= A) gives mu_hat_n(1) = f_n(0) for every n <= h and f_h(s) as an independent check.
Usage: python3 jseries.py h [A] [delta]   (s = round(h log2 3) - delta, default delta = 6)."""
import sys, math, time
import numpy as np
LOG23 = math.log2(3.0); Qc = 1 - 1/LOG23; THETA = math.log(2*Qc); IRATE = THETA*LOG23 - math.log(Qc/(1-Qc))

def phases_level(n, kmin, kmax):
    """omega_n(k) = exp(2 pi i (2^k mod 3^n)/3^n), k in [kmin,kmax] (negative k = modular inverse), complex array."""
    mod = 3**n; L = 2*3**(n-1); K = kmax - kmin + 1
    out = np.empty(K, dtype=np.complex128)
    if L <= 2*K + 4096:                      # periodic table over one period
        tab = np.empty(L); r = 1
        for k in range(L):
            tab[k] = r/mod; r = (r << 1) % mod
        idx = np.arange(kmin, kmax+1) % L
        return np.exp(2j*np.pi*tab[idx])
    if kmax >= 0:
        k0 = max(kmin, 0); r = pow(2, k0, mod); m = kmax - k0 + 1
        vals = np.empty(m)
        for i in range(m):
            vals[i] = r/mod; r = (r << 1) % mod
        out[k0-kmin:] = np.exp(2j*np.pi*vals)
    if kmin < 0:
        M = -kmin                                # need D_n(Mp) for Mp = 1..M
        R = (-pow(3, -n, 1 << (M+1))) % (1 << (M+1))   # -3^-n mod 2^(M+1): 2-adic digits b_1..b_(M+1)
        bstr = bin(R)[2:].zfill(M+1)[::-1]        # bstr[i] = b_(i+1)
        bits = (np.frombuffer(bstr.encode(), dtype=np.uint8) - 48).astype(np.float64)[:M]
        kern = 2.0**(-np.arange(1, 54))            # D(Mp) ~ sum_{t=0}^{52} b_(Mp-t) 2^-(t+1)
        D = np.convolve(bits, kern)[:M]            # D[Mp-1] for Mp = 1..M
        Mp = np.arange(1, M+1, dtype=np.float64)
        corr = np.exp2(-Mp) * (3.0**(-n) if n < 640 else 0.0)   # the real part 2^-Mp 3^-n of (2^-Mp mod 3^n)/3^n
        ph = np.exp(2j*np.pi*(D + corr))          # index Mp-1 <-> k = -Mp
        out[:M] = ph[::-1]                        # k = kmin..-1  <-> Mp = M..1
    return out

def run(h, A, delta, verbose=True):
    s = int(round(h*LOG23)) - delta
    t0 = time.time()
    w = 2.0**(-np.arange(1, A+1))
    LO = -h*A
    # ---- window recursion on [lo_n, s], lo_n = -(h-n)A; f_0 = 1
    prev_lo = LO - A; prev = np.ones(s - prev_lo + 1, dtype=np.complex128)   # f_0 on [LO-A, s]
    mu1 = np.zeros(h+1, dtype=np.complex128); mu1[0] = 1.0
    PH = {}   # phases on [0, s] per level for the top DP
    for n in range(1, h+1):
        lo = -(h-n)*A
        ph = phases_level(n, lo - A, s)          # omega_n on [lo-A, s]
        PH[n] = ph[(0-(lo-A)):]                   # slice for k in [0, s]
        P = ph * prev[(lo - A) - prev_lo : (s - prev_lo) + 1]   # omega_n(k) f_{n-1}(k) on [lo-A, s]
        cur = np.zeros(s - lo + 1, dtype=np.complex128)
        for a in range(1, A+1):
            cur += w[a-1] * P[A - a : A - a + (s - lo + 1)]      # index k-lo  <-  (k-a) - (lo-A) = k - lo + A - a
        prev, prev_lo = cur, lo
        mu1[n] = cur[0 - lo]
    f_hs = prev[s - prev_lo]
    t1 = time.time()
    # ---- top DP on kappa in [0, s]: Zhat_h(kappa) = 2^-(s-kappa) omega_h(kappa), kappa <= s-1
    Z = np.zeros(s+1, dtype=np.complex128)
    kap = np.arange(0, s+1)
    Z[:s] = 2.0**(-(s - kap[:s])) * PH[h][:s]
    W = np.zeros(h+1, dtype=np.complex128)        # W_J for J = 0..h-1 ; W[h] handled separately
    tw = 2.0**(-kap)
    W[h-1] = np.sum(tw * Z)                        # J = h-1: Zhat_{h}
    for j in range(h-1, 0, -1):                    # compute Zhat_j from Zhat_{j+1}
        G = np.zeros(s+1, dtype=np.complex128)
        for a in range(1, A+1):
            if a > s: break
            G[:s+1-a] += w[a-1] * Z[a:]
        Z = PH[j] * G
        if j >= 2: W[j-1] = np.sum(tw * Z)         # W_{J} with J = j-1 uses Zhat_{J+1} = Zhat_j
        else: c0 = np.sum(Z)                       # Zhat_1: c_0
    c = np.zeros(h+1, dtype=np.complex128)
    c[0] = c0
    for J in range(1, h): c[J] = mu1[J] * W[J]
    c[h] = mu1[h] * 2.0**(-s)
    tot = c.sum()
    scale = h**(-1.5) * math.exp(-h*IRATE)
    # masses 2^-s C(s, h-J)
    def lC(nn, kk): return math.lgamma(nn+1) - math.lgamma(kk+1) - math.lgamma(nn-kk+1)
    mass = np.array([math.exp(-s*math.log(2) + lC(s, h-J)) if h-J <= s else 0.0 for J in range(h+1)])
    if verbose:
        print(f"h={h} A={A} s={s} (delta={h*LOG23-s:.2f})  recursion {t1-t0:.1f}s, DP {time.time()-t1:.1f}s")
        print(f"  f_h(s) direct = {f_hs:.6e}  |.|={abs(f_hs):.6e};  sum_J c_J = {tot:.6e}  |.|={abs(tot):.6e};  rel diff {abs(tot-f_hs)/abs(f_hs):.2e}")
        print(f"  |f_h(s)|/(h^-3/2 e^-hI) = {abs(f_hs)/scale:.4f}   arg = {np.angle(f_hs):+.3f}")
        print("  J: |c_J|, mass_h(J), |mu_hat_J(1)|, |W_J|/mass (top coherence), gamma_J = c_J/scale (mod, arg)")
        for J in list(range(0, 26)) + [30, 40, 50, 60, 80, 100, 150, 200]:
            if J > h: break
            coh = abs(W[J])/mass[J] if (J >= 1 and J < h and mass[J] > 0) else float('nan')
            print(f"   {J:3d}: {abs(c[J]):.4e} {mass[J]:.4e} {abs(mu1[J]):.4e} {coh:.4f}   {abs(c[J])/scale:.4f} {np.angle(c[J]):+.2f}")
        ps = np.cumsum(c)
        print("  partial sums |sum_{J'<=J} c_J'|/scale at J = 5,10,20,40,80,160,320:", " ".join(f"{abs(ps[J])/scale:.4f}" for J in (5,10,20,40,80,160,320) if J <= h))
        # H1 data
        v = np.abs(mu1[1:]) * 3.0**(np.arange(1, h+1)/2)
        print(f"  |mu_hat_n(1)| 3^(n/2): max {v.max():.3f} at n={v.argmax()+1}; n=1..12: " + " ".join(f"{x:.3f}" for x in v[:12]))
        for (a_, b_) in ((20, h), (max(20, h//2), h), (max(20, 3*h//4), h)):
            if b_ - a_ >= 10:
                ns = np.arange(a_, b_+1); ys = np.log(np.abs(mu1[a_:b_+1]))
                sl = np.polyfit(ns, ys, 1)[0]
                print(f"  least-squares rate of |mu_hat_n(1)| over n={a_}..{b_}: {math.exp(sl):.4f} per level (3^-1/2 = 0.5774, log2 3 - 1 = 0.5850)")
        print(f"  sup_n (log2 3 - 1)^-n |mu_hat_n(1)| over n<=h = {max(abs(mu1[n])/(LOG23-1)**n for n in range(1,h+1)):.3f}; over 20<=n<=h = {max(abs(mu1[n])/(LOG23-1)**n for n in range(20,h+1)):.3f}")
    return dict(h=h, s=s, f=f_hs, tot=tot, c=c, mu1=mu1, mass=mass, W=W, scale=scale)

if __name__ == "__main__":
    h = int(sys.argv[1]); A = int(sys.argv[2]) if len(sys.argv) > 2 else 60; delta = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    run(h, A, delta)
