"""
mac-mini-2026-07-07-S41 (HYP-4817) -- PZ on the window-avoidance profile U_theta
+ the gap-moment ladder, on the k=13 density-floor leg.

FRAME (new pieces, building on the S130-S134 burst + monad/boxeph audits):
  For a 13-set E and x ~ U[0,1), let g_1..g_13 be the circular gaps of {frac(e x)}.
  U_theta(x) = sum_j (g_j - theta)_+   (opus-S131's E[U] object; window-avoidance integral:
               U_theta(x) = meas{s : arc (s, s+theta) contains no point}).
  IDENTITIES used:
    (1) E[maxgap] = int_0^1 mu_theta dtheta            (mean = integral of tails)
    (2) sum_j g_j^p = p(p-1) int_0^1 u^{p-2} U_u du    (all gap moments = transforms of U)
    (3) maxgap >= (sum_j g_j^p)^{1/(p-1)} pointwise    (the LADDER; p=2 is kps-S58's
        failed length-biased bound sum g^2; p -> infty recovers maxgap)
    (4) PALEY-ZYGMUND on U:  mu_theta = P(maxgap>theta) = P(U_theta>0) >= E[U]^2/E[U^2]
        -- both moments are BALANCED-lattice sums (s-averaging forces sum m_i = 0,
        killing all pair relations; 3-term APs of E lead). Sharper than P(U>0)>=E[U]
        (opus-S131) whenever E[U^2] < E[U], and continuous-anchor vs monad-S1's 14-window CE.
  TARGETS: m_P = 14249/252252 ~ 0.05649 (THM-530 floor, the k=13 hlarge bar);
           1/7 for the mean; T* = 1/7 + (6/7)m_P ~ 0.191275 (monad HYP-4787 bar).

  QUESTIONS:
   Q1: does PZ-on-U (E[U]^2/E[U^2] at theta=1/7) clear m_P uniformly on the adversarial
       bank AND under descent minimizing the PZ ratio itself?
   Q2: at which p does the ladder E[(sum g^p)^{1/(p-1)}] clear 1/7 (and T*) at the
       E[maxgap]-record families? (p=2 known to FAIL at 0.135.)
   Q3: how much of E[U_{1/7}] is captured by main term (1-u)^13 + 3-AP balanced
       relations only? (the provability question)
NUMERIC (grid); exactness deferred. Grid convergence checked at the end.
"""
import numpy as np
from math import gcd
from functools import reduce
import random as rnd
rnd.seed(41)

GRID = 200_000
xs = (np.arange(GRID) + 0.5) / GRID
THETA = 1.0 / 7.0
MP = 14249 / 252252
TSTAR = 1/7 + (6/7) * MP

def gaps(E):
    """gaps[x_index, j] for the 13 circular gaps of {frac(e x)}."""
    ph = np.mod(np.outer(xs, np.array(E, float)), 1.0)
    ph.sort(axis=1)
    return np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)

def stats(E, theta=THETA):
    g = gaps(E)
    mg = g.max(axis=1)
    U = np.clip(g - theta, 0, None).sum(axis=1)
    out = {
        'Emax': mg.mean(),
        'mu': (mg > theta).mean(),
        'EU': U.mean(),
        'EU2': (U * U).mean(),
    }
    out['PZ'] = out['EU'] ** 2 / out['EU2'] if out['EU2'] > 0 else 0.0
    # ladder p=2..6: E[(sum g^p)^{1/(p-1)}]
    for p in (2, 3, 4, 5, 6):
        Y = (g ** p).sum(axis=1)
        out[f'lad{p}'] = (Y ** (1.0 / (p - 1))).mean()
    return out

def norm(E):
    E = sorted(set(E)); E = [e - E[0] for e in E]
    g = reduce(gcd, E[1:]) if len(E) > 1 else 1
    return tuple(e // g for e in E) if g else tuple(E)

BANK = {
    'AP {1..13}':                 list(range(1, 14)),
    'GW {1..11,13,24}':           list(range(1, 12)) + [13, 24],
    'death-star 2*{1..12}u{13}':  [2*i for i in range(1, 13)] + [13],
    'monad record evens+{11,13}': [2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
    'opus stretched {0,2..12,17,28}': [0] + list(range(2, 13)) + [17, 28],
    'monad S1-min [6,8,..61]':    [6, 8, 10, 11, 12, 13, 14, 16, 18, 25, 36, 43, 61],
    'two-scale {1..6}u50*{1..7}': list(range(1, 7)) + [50*i for i in range(1, 8)],
    'random big':                 sorted(rnd.sample(range(1, 2000), 13)),
    'primes13':                   [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41],
}

print("=== S41 Q1+Q2: PZ-on-U and the gap-moment ladder at theta=1/7 ===")
print(f"bars: m_P={MP:.5f} (tail), 1/7={THETA:.5f} (mean), T*={TSTAR:.5f} (mean, monad)")
print(f"{'family':34s} {'E[mg]':>7s} {'mu':>6s} {'E[U]':>6s} {'PZ=EU^2/EU2':>11s} "
      f"{'lad2':>6s} {'lad3':>6s} {'lad4':>6s} {'lad5':>6s} {'lad6':>6s}")
for name, E in BANK.items():
    s = stats(E)
    flag = ' PZ<m_P!' if s['PZ'] < MP else ''
    print(f"{name:34s} {s['Emax']:7.4f} {s['mu']:6.3f} {s['EU']:6.3f} {s['PZ']:11.4f} "
          f"{s['lad2']:6.4f} {s['lad3']:6.4f} {s['lad4']:6.4f} {s['lad5']:6.4f} {s['lad6']:6.4f}{flag}")

print("\nladder verdicts at the E[maxgap]-record families: which p clears 1/7? T*?")
for name in ('monad record evens+{11,13}', 'death-star 2*{1..12}u{13}', 'opus stretched {0,2..12,17,28}'):
    s = stats(BANK[name])
    clears17 = [p for p in (2, 3, 4, 5, 6) if s[f'lad{p}'] > THETA]
    clearsT = [p for p in (2, 3, 4, 5, 6) if s[f'lad{p}'] > TSTAR]
    print(f"  {name:34s} clears 1/7 at p={clears17}; clears T* at p={clearsT}")

# ---- Q1 descent: minimize the PZ ratio itself (adversarial correctness discipline) ----
print("\n=== Q1 descent: adversarially MINIMIZE PZ = E[U]^2/E[U^2] (coarse grid, refine finalists) ===")
GRID_C = 30_000
xs_c = (np.arange(GRID_C) + 0.5) / GRID_C
def pz_coarse(E, theta=THETA):
    ph = np.mod(np.outer(xs_c, np.array(E, float)), 1.0)
    ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    U = np.clip(g - theta, 0, None).sum(axis=1)
    m2 = (U * U).mean()
    return (U.mean() ** 2 / m2) if m2 > 0 else 0.0

best_overall = []
for start_name in ('AP {1..13}', 'monad record evens+{11,13}', 'monad S1-min [6,8,..61]', 'random big'):
    E = list(BANK[start_name]); cur = pz_coarse(E)
    for it in range(400):
        i = rnd.randrange(13)
        cand = E.copy()
        cand[i] = max(0, cand[i] + rnd.choice([-3, -2, -1, 1, 2, 3, rnd.randint(-40, 40)]))
        if len(set(cand)) < 13:
            continue
        cv = pz_coarse(cand)
        if cv < cur:
            E, cur = cand, cv
    best_overall.append((cur, norm(E), start_name))
best_overall.sort()
print(f"{'start':30s} {'coarse min PZ':>13s}  shape")
for cur, E, sn in best_overall:
    print(f"{sn:30s} {cur:13.4f}  {E}")
worstPZ, worstE, _ = best_overall[0]
sw = stats(list(worstE))
print(f"\nfinalist re-eval at fine grid: PZ={sw['PZ']:.4f} mu={sw['mu']:.3f} E[U]={sw['EU']:.4f} "
      f"(m_P={MP:.4f}; PZ>m_P? {sw['PZ'] > MP})")

# ---- Q3: balanced-lattice structure of E[U_u]: main term + 3-AP relations ----
print("\n=== Q3: E[U_{1/7}] vs main term (1-u)^13 and the balanced 3-AP relations ===")
u = THETA
main = (1 - u) ** 13
def count_3aps(E):
    S = set(E); c = 0
    for a in E:
        for cc in E:
            if a < cc and (a + cc) % 2 == 0 and (a + cc) // 2 in S and (a + cc)//2 not in (a, cc):
                c += 1
    return c
print(f"main term (1-1/7)^13 = {main:.5f}")
print(f"{'family':34s} {'E[U]':>7s} {'E[U]-main':>9s} {'#3APs':>6s} {'deficit/3AP':>11s}")
for name, E in BANK.items():
    s = stats(E)
    n3 = count_3aps(E)
    d = s['EU'] - main
    per = d / n3 if n3 else float('nan')
    print(f"{name:34s} {s['EU']:7.4f} {d:9.4f} {n3:6d} {per:11.5f}")

# grid convergence sanity: AP E[maxgap] should be 93/440 = 0.211364
s_ap = stats(BANK['AP {1..13}'])
print(f"\ngrid sanity: E[maxgap](AP) = {s_ap['Emax']:.6f} vs 93/440 = {93/440:.6f} "
      f"(|err| = {abs(s_ap['Emax'] - 93/440):.2e})")
print(f"            mu_1/7(AP) = {s_ap['mu']:.6f} vs 477/1078 = {477/1078:.6f} "
      f"(|err| = {abs(s_ap['mu'] - 477/1078):.2e})")
