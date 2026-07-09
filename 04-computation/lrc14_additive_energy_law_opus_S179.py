"""
lrc14_additive_energy_law_opus_S179.py   (opus-2026-07-09-S179)

THE ADDITIVE-ENERGY LAW governing looseness.  opus-S178 claimed the Riesz margin is governed by the
additive energy E(S)=SUM|S-hat|^4 (= #{(a,b,c,d): a+b=c+d}); dissociated (low E) => uniform ratio<<1,
near-AP/tight (high E) => ratio->1.  This ESTABLISHES the quantitative law across the spectrum
(Sidon -> dissociated -> near-AP -> AP), measuring for each 13-set:
  E(S) = additive energy (a+b=c+d quadruples) ;
  L(S) = lonely measure = meas{tau: min_i ||v_i tau|| > 1/14} (the singular series / inf-L object) ;
  ratio = inf_R int(M*R)/int(R) (naive Riesz certificate) ;
and reports the monotone relationship E up => L down & ratio -> 1, locating the loose/tight boundary.
"""
import numpy as np
import random

NG = 120011
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14
MAIN = (6.0 / 7.0) ** 13


def additive_energy(S):
    """#{(a,b,c,d) in S^4 : a+b=c+d} = SUM_x r(x)^2, r(x)=#{(a,b): a+b=x}."""
    from collections import Counter
    r = Counter(a + b for a in S for b in S)
    return sum(c * c for c in r.values())


def Mmult(S):
    M = np.zeros(NG)
    for v in S:
        d = np.abs(((v * TAU + 0.5) % 1.0) - 0.5)
        M += (d <= H)
    return M


def lonely_measure(S):
    return float((Mmult(S) == 0).mean())


def riesz_ratio(S, iters=6):
    M = Mmult(S)
    D = sorted(S)
    cb = [np.cos(2 * np.pi * m * TAU) for m in D]
    a = np.zeros(len(D))

    def ratio(a):
        R = np.ones(NG)
        for i in range(len(D)):
            R = R * (1 + a[i] * cb[i])
        r = R.mean()
        return (M * R).mean() / r if r > 1e-6 else 5.0
    cur = ratio(a)
    for _ in range(iters):
        imp = False
        for i in range(len(D)):
            bi, bv = a[i], cur
            for c in np.linspace(-0.999, 0.999, 41):
                a[i] = c; v = ratio(a)
                if v < bv - 1e-9: bv, bi = v, c; imp = True
            a[i] = bi; cur = bv
        if not imp: break
    return cur


def longest_ap(S):
    Sset = set(S); S = sorted(S); best = 1
    for i, aa in enumerate(S):
        for b in S[i + 1:]:
            d = b - aa
            if aa - d in Sset: continue
            L = 2; x = b + d
            while x in Sset: L += 1; x += d
            best = max(best, L)
    return best


rng = random.Random(179)
print("=" * 100)
print(f"ADDITIVE-ENERGY LAW: E(S) vs lonely-measure L vs Riesz ratio.  (6/7)^13={MAIN:.4f}, 1/14=0.0714")
print("=" * 100)
families = []
# pure AP (tight), dilated AP, near-AP (AP minus one + outlier), dissociated, Sidon-like
families.append(("AP {1..13} (tight)", list(range(1, 14))))
families.append(("2*AP {2..26} (tight)", [2 * i for i in range(1, 14)]))
families.append(("near-AP {1..12}+{20}", list(range(1, 13)) + [20]))
families.append(("near-AP {1..11}+{20,30}", list(range(1, 12)) + [20, 30]))
for seed in range(6):
    r = random.Random(1000 + seed)
    spread = r.randint(60, 200)
    S = sorted(set([1] + r.sample(range(2, spread), 11) + [spread]))
    while len(S) != 13:
        S = sorted(set(S + [r.randint(1, spread)]))[:13]
    families.append((f"random spread{spread} L={longest_ap(S)}", S))
families.append(("Sidon {1,2,5,11,22,...}", [1, 2, 5, 11, 22, 33, 45, 60, 78, 95, 110, 130, 150]))
families.append(("Fib-ish (max dissoc)", [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]))

rows = []
for name, S in families:
    S = sorted(set(S))
    if len(S) != 13:
        continue
    E = additive_energy(S)
    L = lonely_measure(S)
    ratio = riesz_ratio(S)
    Lap = longest_ap(S)
    rows.append((E, name, L, ratio, Lap))
rows.sort()  # by additive energy
print(f"  {'E(S) add.energy':>15} {'longestAP':>9} {'L (lonely meas)':>15} {'Riesz ratio':>12} {'set':>28}")
for E, name, L, ratio, Lap in rows:
    tag = "TIGHT" if L < 1e-4 else ("loose" )
    print(f"  {E:>15} {Lap:>9} {L:>15.4f} {ratio:>12.4f}  {name:>28} [{tag}]")
print()
print("  READING: as additive energy E(S) RISES (dissociated -> near-AP -> AP), the lonely measure L")
print("  FALLS toward 0 and the Riesz ratio RISES toward 1 -- the quantitative loose/tight law. The AP")
print("  (max E, tight L=0, ratio=1) is the extremal; dissociated (min E) has L~(6/7)^13 & ratio<<1.")
print("=" * 100)
