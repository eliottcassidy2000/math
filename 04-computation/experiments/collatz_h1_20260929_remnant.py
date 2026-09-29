"""The growing remnant: the negative family v_n(m) = mu_hat_n(2^-m) for n <= NMAX, 1 <= m <= MMAX (window recursion, a <= A),
the ridge path (argmax of 3^(n/2)|v_n(m)| over m in [1, MMAX]), and along the path the exact one-step energy budget
  |v_n(m)|^2 = sum_a 4^-a |v_{n-1}(m+a)|^2  +  2 Re sum_{a<a'} 2^-a-a' omega_n(-m-a) conj(omega_n(-m-a')) v_{n-1}(m+a) conj(v_{n-1}(m+a'))
             = direct (incoherent, random-phase prediction) + cross (coherence).
Random phases give E|v_n|^2 = direct: energy x 1/3 per level (the Parseval scale) transported with the step law 3 4^-a
(mean 4/3 = the observed ridge slope -1.3, std 2/3).  Growth against the Parseval scale = positive cross terms."""
import sys, math, numpy as np
sys.path.insert(0, '.')
from jseries import phases_level
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 140
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 160
A = 60
# window at level n: [-(NMAX-n)A - MMAX, 0]
prev_lo = -NMAX*A - MMAX - A; prev = np.ones(-prev_lo + 1, dtype=np.complex128)  # f_0 on [prev_lo, 0]
V = {}   # V[n] = array of v_n(m), m = 0..MMAX  (index m)
OM = {}  # OM[n] = omega_n(-M) for M = 0..MMAX+A (index M)
w = 2.0**(-np.arange(1, A+1))
for n in range(1, NMAX+1):
    lo = -(NMAX-n)*A - MMAX
    ph = phases_level(n, lo - A, 0)                # omega_n(k), k in [lo-A, 0]
    P = ph * prev[(lo-A) - prev_lo : (0 - prev_lo) + 1]
    cur = np.zeros(-lo + 1, dtype=np.complex128)   # k in [lo, 0]
    for a in range(1, A+1):
        cur += w[a-1] * P[A-a : A-a + (-lo+1)]
    prev, prev_lo = cur, lo
    V[n] = cur[(0-lo) - MMAX : (0-lo) + 1][::-1].copy()     # v_n(m) for m = MMAX..0 reversed -> index m
    OM[n] = ph[(0-(lo-A)) - (MMAX+A) : (0-(lo-A)) + 1][::-1].copy()  # omega_n(-M), M = 0..MMAX+A
print(f"negative family to n={NMAX}, m<=MMAX={MMAX}, A={A}")
print("ridge path: n, argmax m in [1,MMAX], v = 3^(n/2)|v_n(m)|, v_2(n), and the one-step energy budget along the path:")
print("   n   m*    v     v2(n)   direct/prev  cross/direct   gain=3|v_n|^2/|v_{n-1}(m_{n-1})|^2  top-3 term phases (deg, relative)")
path = {}
for n in range(1, NMAX+1):
    vv = 3.0**(n/2) * np.abs(V[n][1:MMAX+1]); m = int(vv.argmax()) + 1; path[n] = (m, float(vv[m-1]))
prevpk = None
for n in range(60, NMAX+1):
    m, vn = path[n]
    terms = np.array([w[a-1] * OM[n][m+a] * V[n-1][m+a] for a in range(1, A+1) if m+a <= MMAX])
    direct = float(np.sum(np.abs(terms)**2)); tot = abs(np.sum(terms))**2
    cross = tot - direct
    v2 = 0; t = n
    while t % 2 == 0: t //= 2; v2 += 1
    # relative phases of the three largest terms
    idx = np.argsort(-np.abs(terms))[:3]
    ph = np.degrees(np.angle(terms[idx] * np.conj(np.sum(terms))))
    gain = 3*abs(V[n][m])**2 / (abs(V[n-1][path[n-1][0]])**2) if prevpk else float('nan')
    print(f"  {n:3d}  {m:3d}  {vn:6.2f}   {v2}     {direct/ (abs(V[n-1][path[n-1][0]])**2):.3f}        {cross/direct:+.3f}         {gain:6.3f}      a={idx+1} {np.round(ph,0)}")
    prevpk = True
# cumulative: product of per-level energy gains vs observed
print("\nsummary over the wave n=84..134: geometric-mean per-level factor of v:", (path[134][1]/path[84][1])**(1/50) if path[84][1] > 0 else None)
# the frequency-1 coefficient itself (m=0) around the arrival
print("m=0 (frequency 1), 3^(n/2)|mu_hat_n(1)| for n=120..140:", " ".join(f"{3.0**(n/2)*abs(V[n][0]):.2f}" for n in range(120, 141)))
np.save('remnant_V.npy', np.array([V[n] for n in range(1, NMAX+1)]))
