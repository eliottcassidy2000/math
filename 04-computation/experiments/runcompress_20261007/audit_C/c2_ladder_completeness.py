#!/usr/bin/env python3
"""audit_C item 2: THM-4603 (2) ladder completeness, tested on merges from MANY states (independent code).

For a pair-chain state (k0, e0) (u = 3^k0 v + e0), take random v, run the two Terras orbits in lockstep with the debt
(k += par(u) - par(v)) until the first coincidence u_s = v_s (an absorption iff k_s = 0) or a horizon.
At an absorption, check every assertion of THM-4603 (2):
  (a) z_u, z_v (last odd values strictly before s) exist, z_u != z_v, and U(z_u) = U(z_v) = Z (first common odd value);
  (b) ladder: z_v = 4^i z_u + (4^i-1)/3 (child) or the mirror (source), i >= 1, letters c and c+2i;
  (c) heads: |h_v| = |h_u| + k0, and the equal-time condition t_u = t_v + 2i (child) / t_v = t_u + 2i (source), where
      t_u, t_v are the Terras times of z_u, z_v (= g + Sum h, g = initial halvings);
  (d) merge time s = max(t_u, t_v) + 1;
  (e) the head identity is an identity of affine maps of v (slopes and constants), hence holds at the anchor
      c* = e0/(1 - 3^k0) when k0 != 0 (anchored form F_hv(c*) = 4^i F_hu(c*) + (4^i-1)/3, generalised to even starts).
States: universal (3,-26); re-anchored (-5, 3^-5 - 1) [negative debt]; (4,10); (6,-728); (1,2) [D=1 chain]; (2,8);
(-2, 1/9 - 1); (5,1) [the (J,1) collapse state]; (1,1) [reset]; (0,1) [n+1 vs n]; (-3, random) and random (k, e).
v is a random integer of 60-400 bits; u = A(v) may be a rational with 3-power denominator when k0 < 0 (2-adic chain).
Edge cases recorded: an orbit with no odd value before the merge; inadmissible starts (THM-4581 3'), where the letter
difference b can be odd (half-ladders) - searched separately.
"""
import random
from fractions import Fraction as Fr
from collections import Counter

rnd = random.Random(4603)


def par(q):
    return q.numerator & 1


def T(q):
    return (3 * q + 1) / 2 if par(q) else q / 2


def oddpart(q):
    while par(q) == 0:
        q = q / 2
    return q


def U(q):
    return oddpart(3 * q + 1)


def v2q(q):
    n = q.numerator
    return (n & -n).bit_length() - 1


def affine_of_prefix(x0, s):
    """the Terras orbit of x0 for s steps as an affine map x -> (alpha x + beta) valid on x0's parity class: returns
    list of (alpha, beta) for times 0..s"""
    maps = [(Fr(1), Fr(0))]
    x = x0
    a, b = Fr(1), Fr(0)
    for t in range(s):
        if par(x):
            a, b = 3 * a / 2, (3 * b + 1) / 2
        else:
            a, b = a / 2, b / 2
        x = T(x)
        maps.append((a, b))
    return maps


def run(k0, e0, v0, H):
    """lockstep with integer numerators over a common odd denominator d (the 3x+d form of T on (1/d)Z):
    x = a/d, T: a -> (3a + d)/2 (a odd), a/2 (a even). Returns (s, k_s, hist_u, hist_v) as Fractions at the first
    coincidence, or None (horizon, or both orbits on the trivial cycle {1,2})."""
    u0 = Fr(3) ** k0 * v0 + e0
    d = u0.denominator
    au, av, k = u0.numerator, v0 * d, k0
    hu, hv = [au], [av]
    one, two = d, 2 * d
    for s in range(1, H + 1):
        pu, pv = au & 1, av & 1
        k += pu - pv
        au = (3 * au + d) >> 1 if pu else au >> 1
        av = (3 * av + d) >> 1 if pv else av >> 1
        hu.append(au)
        hv.append(av)
        if au == av:
            return s, k, [Fr(a, d) for a in hu], [Fr(a, d) for a in hv]
        if (au == one or au == two) and (av == one or av == two):
            return None
    return None


stats = Counter()
ladder_kinds = Counter()


def analyse(k0, e0, s, hu, hv):
    U0, V0 = hu[0], hv[0]
    oddu = [t for t in range(s) if par(hu[t])]
    oddv = [t for t in range(s) if par(hv[t])]
    if not oddu or not oddv:
        stats['edge: an orbit has no odd value before the merge'] += 1
        return
    tu, tv = oddu[-1], oddv[-1]
    zu, zv = hu[tu], hv[tv]
    Z = oddpart(hu[s])
    assert zu != zv and U(zu) == Z and U(zv) == Z
    cu, cv = v2q(3 * zu + 1), v2q(3 * zv + 1)
    b = cv - cu
    if b % 2:
        stats['odd letter difference (half-ladder)'] += 1
        return
    i = abs(b) // 2
    assert i >= 1
    if b > 0:
        assert zv == 4 ** i * zu + Fr(4 ** i - 1, 3)
        kind = 'child'
        assert tu == tv + 2 * i
    else:
        assert zu == 4 ** i * zv + Fr(4 ** i - 1, 3)
        kind = 'source'
        assert tv == tu + 2 * i
    # heads: odd values before z
    assert len(oddv) - 1 == len(oddu) - 1 + k0
    # merge time
    assert s == max(tu, tv) + 1
    # affine identity of the head maps as maps of v: z_v(v) = 4^i z_u(A(v)) + (4^i-1)/3 (child) or mirror
    mu = affine_of_prefix(U0, tu)[tu]
    mv = affine_of_prefix(V0, tv)[tv]
    # z_u as a map of v: mu o A, A(v) = 3^k0 v + e0
    au, bu = mu[0] * Fr(3) ** k0, mu[0] * e0 + mu[1]
    av, bv = mv
    if kind == 'child':
        assert av == 4 ** i * au and bv == 4 ** i * bu + Fr(4 ** i - 1, 3)
    else:
        assert au == 4 ** i * av and bu == 4 ** i * bv + Fr(4 ** i - 1, 3)
    if k0 != 0:
        cs = e0 / (1 - Fr(3) ** k0)
        Fu, Fv = mu[0] * cs + mu[1], av * cs + bv      # head maps evaluated at the anchor
        if kind == 'child':
            assert Fv == 4 ** i * Fu + Fr(4 ** i - 1, 3)
        else:
            assert Fu == 4 ** i * Fv + Fr(4 ** i - 1, 3)
        stats['anchored identity at c* checked'] += 1
    ladder_kinds[(kind, i)] += 1
    stats['merges fully verified'] += 1
    mv = hu[s]
    if mv.denominator == 1 and mv.numerator in (1, 2):
        stats['merge on the trivial cycle (value 1 or 2)'] += 1
    elif par(mv):
        stats['merge at an odd value (final letter c = 1)'] += 1
    else:
        stats['merge at an even value'] += 1
    if hu[s].denominator != 1:
        stats['merge at a non-integer (2-adic) value'] += 1


states = {
    'universal (3,-26)': (3, Fr(-26)),
    're-anchored (-5, 3^-5 - 1)': (-5, Fr(1, 243) - 1),
    '(4,10)': (4, Fr(10)),
    '(6,-728)': (6, Fr(-728)),
    'D=1 chain (1,2)': (1, Fr(2)),
    'D=2 chain (2,8)': (2, Fr(8)),
    '(-2, 1/9 - 1)': (-2, Fr(1, 9) - 1),
    '(J,1) J=5': (5, Fr(1)),
    'reset (1,1)': (1, Fr(1)),
    'translation (0,1)': (0, Fr(1)),
}
per_state = {}
for name, (k0, e0) in states.items():
    before = stats['merges fully verified']
    before_edge = stats['edge: an orbit has no odd value before the merge']
    tried = 0
    for trial in range(500):
        bits = rnd.choice((60, 120, 250, 400))
        v0 = rnd.getrandbits(bits) | (1 << (bits - 1))
        if rnd.random() < 0.8:
            v0 |= 1
        r = run(k0, e0, v0, 3000)
        tried += 1
        if r is None:
            continue
        s, k, hu, hv = r
        if k != 0:
            stats['coincidence with nonzero debt ' + ('at a value <= 2 (trivial cycle)' if (hu[s].denominator == 1 and hu[s].numerator in (1, 2)) else 'at a value > 2 (sporadic)')] += 1
            continue
        analyse(k0, e0, s, hu, hv)
    per_state[name] = (stats['merges fully verified'] - before,
                       stats['edge: an orbit has no odd value before the merge'] - before_edge, tried)

# random admissible states, including negative debt
for trial in range(1500):
    k0 = rnd.randint(-6, 6)
    if k0 >= 0:
        e0 = Fr(rnd.randint(-10 ** 6, 10 ** 6))
    else:
        e0 = Fr(rnd.randint(-10 ** 6, 10 ** 6), 3 ** (-k0))
    if k0 == 0 and e0 == 0:
        continue
    bits = rnd.choice((60, 150, 300))
    v0 = rnd.getrandbits(bits) | (1 << (bits - 1))
    r = run(k0, e0, v0, 3000)
    if r is None:
        continue
    s, k, hu, hv = r
    if k != 0:
        stats['coincidence with nonzero debt ' + ('at a value <= 2 (trivial cycle)' if (hu[s].denominator == 1 and hu[s].numerator in (1, 2)) else 'at a value > 2 (sporadic)')] += 1
        continue
    analyse(k0, e0, s, hu, hv)
per_state['random admissible (k in [-6,6])'] = None

print("per state: (fully verified ladder merges, edge merges with an orbit lacking an odd value before the merge, trials)")
for name, val in per_state.items():
    print(f"  {name:32s} {val}")
print("stats:", dict(stats))
print("ladder kinds (kind, i):", dict(sorted(ladder_kinds.items())))

# ---- explicit edge case: an orbit without an odd value before the merge, admissible state ----
# state (-1, -1/3): u = v/3 - 1/3. v = 16 gives u = 5: T(5) = 8 = T(16), debt -1 + 1 - 0 = 0.
k0, e0, v0 = -1, Fr(-1, 3), 16
r = run(k0, e0, v0, 10)
print("edge example (k0,e0) = (-1,-1/3), v = 16, u = 5:", "merge at s =", r[0], "debt", r[1], "value", r[2][r[0]],
      "; v has no odd value before the merge (z_v undefined)")

# ---- inadmissible starts (3-adic excess, THM-4581 3'): search for odd letter differences ----
found = []
for trial in range(4000):
    k0 = rnd.randint(-3, 3)
    m = max(0, -k0) + rnd.randint(1, 2)          # excess denominator
    e0 = Fr(rnd.randint(-3000, 3000), 3 ** m)
    if e0.denominator <= 3 ** max(0, -k0):
        continue
    v0 = rnd.getrandbits(40) | 1
    r = run(k0, e0, v0, 400)
    if r is None or r[1] != 0:
        continue
    s, k, hu, hv = r
    oddu = [t for t in range(s) if par(hu[t])]
    oddv = [t for t in range(s) if par(hv[t])]
    if not oddu or not oddv:
        continue
    zu, zv = hu[oddu[-1]], hv[oddv[-1]]
    b = v2q(3 * zv + 1) - v2q(3 * zu + 1)
    if b % 2:
        found.append((k0, e0, v0, s, b, zu, zv))
        if len(found) >= 3:
            break
print("inadmissible starts with an ODD letter difference at the merge (not a 4z+1 ladder):", len(found))
for f in found:
    k0, e0, v0, s, b, zu, zv = f
    print(f"   state ({k0}, {e0}), v = {v0}: merge at s = {s}, letters differ by b = {b}, z_u = {zu}, z_v = {zv}, "
          f"z_v = 2^b z_u + (2^b-1)/3: {zv == 2**b * zu + Fr(2**b - 1, 3) or zu == 2**(-b) * zv + Fr(2**(-b) - 1, 3)}")
print("ALL ASSERTIONS PASSED")
