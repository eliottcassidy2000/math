"""THM-3012 addendum 5 (R5): the CERTIFIED bounded exclusion for the residual
constant Lambda = 2G + D of the exact reduction  S(4) = (2 sqrt2 K(1/r2)/pi^2) Lambda.

Instrument calibration lives in
  04-computation/sk_S4_lemniscatic_eisenstein_reduction_thm3012.py  (checks C6a-C6c5):
the SAME PSLQ scanner blind-recovers pi*S(3) (with its irrational sqrt3 multiplier),
the weight-2 Eisenstein Eichler value sum 1/(n(e^{2 pi n}-1)) at tau=i in terms of
Gamma(1/4), the Ramanujan coth sums, and the Gamma(1/4)^8 K-moments.  A null here is
therefore a bounded statement, not evidence of nonexistence.
"""
import mpmath as mp, itertools, time, sys
from fractions import Fraction
DPS = 210
mp.mp.dps = DPS
r = mp.sqrt; pi = mp.pi

def beta_series(dps):
    with mp.workdps(dps+30):
        b = Fraction(0); tot = mp.mpf(0); m = 0
        while True:
            m += 1
            b = (Fraction((2*m-1)**2)*b + 1)/Fraction(4*m*m)
            t = mp.mpf(b.numerator)/mp.mpf(b.denominator)/mp.mpf(2)**m
            tot += t
            if m > 60 and t < mp.mpf(10)**(-(dps+25)):
                break
        return +(tot/2)

K = mp.ellipk(mp.mpf(1)/2); G = mp.catalan; z3 = mp.zeta(3)
g14 = mp.gamma(mp.mpf(1)/4); g13 = mp.gamma(mp.mpf(1)/3)
g18 = mp.gamma(mp.mpf(1)/8); g38 = mp.gamma(mp.mpf(3)/8)
L2 = mp.log(2); Ls = mp.log(1+r(2))
D = pi*beta_series(DPS)/K
Lam = 2*G + D
U = 5*pi**2/12 - 2*Lam
S4 = 2*r(2)*K*Lam/pi**2

ATOM = {
 'one': mp.mpf(1), 'pi': pi, 'pi^2': pi**2, 'pi^3': pi**3, 'pi^4': pi**4,
 'G': G, 'G*pi': G*pi, 'G/pi': G/pi, 'G^2': G**2,
 'log2': L2, 'Ls': Ls, 'log2^2': L2**2, 'Ls^2': Ls**2, 'log2*Ls': L2*Ls,
 'pi*log2': pi*L2, 'pi*Ls': pi*Ls,
 'K': K, 'K^2': K**2, 'K^2/pi': K**2/pi, 'K^2/pi^2': K**2/pi**2,
 'K^2/pi^3': K**2/pi**3, 'K/pi': K/pi, 'K/pi^2': K/pi**2,
 'pi/K': pi/K, 'pi^2/K': pi**2/K, 'pi^2/K^2': pi**2/K**2,
 'pi^3/K^2': pi**3/K**2, 'pi^4/K^2': pi**4/K**2, 'pi^4/K^4': pi**4/K**4,
 'K^4/pi^2': K**4/pi**2, 'K^4/pi^3': K**4/pi**3, 'K^4/pi^4': K**4/pi**4,
 'K*G': K*G, 'G/K': G/K, 'K*G/pi': K*G/pi, 'K^2*G/pi^2': K**2*G/pi**2,
 'K*Ls': K*Ls, 'K*log2': K*L2, 'K^2*Ls/pi': K**2*Ls/pi,
 'zeta3': z3, 'zeta3/pi': z3/pi, 'zeta3/pi^2': z3/pi**2, 'zeta3/K': z3/K,
 'g13^3/pi': g13**3/pi, 'g13^6/pi^4': g13**6/pi**4,
 'g18g38': g18*g38, '(g18g38)^2/pi': (g18*g38)**2/pi,
 'log(g14)': mp.log(g14), 'logpi': mp.log(pi),
}
MULT = {'1': mp.mpf(1), 'r2': r(2), '1/r2': 1/r(2), 'r3': r(3), 'r5': r(5),
        '2^(1/4)': mp.mpf(2)**(mp.mpf(1)/4)}
PROD = [(f"{m}*{a}", MULT[m]*ATOM[a]) for m in MULT for a in ATOM]
NP = len(PROD)
print(f"products {NP}  pairs {NP*(NP-1)//2}  dps {DPS}")
sys.stdout.flush()
TOL = mp.mpf(10)**-135
H = 10**5

def pairscan(target, name):
    t0 = time.time(); hits = []
    for c in itertools.combinations(PROD, 2):
        rel = mp.pslq([target, c[0][1], c[1][1]], maxcoeff=H, maxsteps=30000, tol=TOL)
        if rel and rel[0] != 0:
            hits.append((rel, c[0][0], c[1][0]))
    print(f"{name}: PAIR hits {len(hits)}  [{time.time()-t0:.0f}s]")
    for h in hits[:20]:
        print("    ", h)
    sys.stdout.flush()

# ---- tier 1: very large coefficient bound on the SIMPLEST shapes --------------
SIMPLE = {k: ATOM[k] for k in
          ['one', 'pi', 'pi^2', 'pi^3', 'G', 'G*pi', 'K', 'K^2/pi', 'K^2/pi^2',
           'pi^2/K^2', 'pi^3/K^2', 'K*G/pi', 'zeta3', 'log2', 'Ls', 'K^4/pi^3']}
SIMPLE['r2'] = r(2); SIMPLE['r2*pi^2'] = r(2)*pi**2
print("\nTIER 1  size<=2 over %d simple atoms, 190-digit tol" % len(SIMPLE))
for tname, tv in [('Lambda', Lam), ('U', U), ('S(4)', S4)]:
    for HH in (10**8, 10**12, 10**20):
        found = []
        for size in (1, 2):
            for c in itertools.combinations(list(SIMPLE.items()), size):
                rel = mp.pslq([tv]+[x[1] for x in c], maxcoeff=HH, maxsteps=200000,
                              tol=mp.mpf(10)**-190)
                if rel and rel[0] != 0:
                    found.append((rel, [x[0] for x in c]))
        print(f"  {tname:7s} H={HH:<22d} hits {len(found)} {found[:2]}")
        sys.stdout.flush()

# ---- tier 2: all size-3 subsets of the independent atom pool ------------------
print("\nTIER 2  all size-3 subsets of the %d-atom pool, H=1e5, 135-digit tol" % len(ATOM))
for tname, tv in [('Lambda', Lam), ('S(4)', S4)]:
    t0 = time.time(); hits = []
    for c in itertools.combinations(list(ATOM.items()), 3):
        rel = mp.pslq([tv]+[x[1] for x in c], maxcoeff=10**5, maxsteps=60000, tol=TOL)
        if rel and rel[0] != 0:
            hits.append((rel, [x[0] for x in c]))
    print(f"  {tname}: hits {len(hits)}  [{time.time()-t0:.0f}s]  {hits[:3]}")
    sys.stdout.flush()

# ---- tier 3: all pairs of products alpha*atom ---------------------------------
print("\nTIER 3  all %d pairs of %d products alpha*atom, H=1e5, 135-digit tol"
      % (NP*(NP-1)//2, NP))
# tier 3 is run on Lambda only: S(4) = (2 sqrt2 K/pi^2) Lambda and
# U = 5 pi^2/12 - 2 Lambda are exact rescalings, already covered at tiers 1-2.
for tgt, nm in [(Lam, 'Lambda')]:
    pairscan(tgt, nm)
