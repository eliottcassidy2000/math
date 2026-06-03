#!/usr/bin/env python3
"""
claudebox-2026-06-03-S585 : Lonely Runner through two lenses — abstract functions & rep theory.

Two-lens exploration of the Lemma A (randomness / circuit-free) vs Lemma B (structure / 3-term
fold) dichotomy.

REP-THEORY DICTIONARY (the lens we mostly use):
  speed set S ⊂ ℤ_{>0}          ↔  finite subset of the character lattice ℤ̂ of ℝ/ℤ
  runner v at time t            ↔  character χ_v(t) = e^{2πi v t}
  3-term relation v_c=v_a+v_b   ↔  FUSION  χ_a·χ_b = χ_c
  danger B_i = {t:‖v_i t‖<δ}    ↔  arc indicator; Fourier coeffs b̂(m) = sin(2π m δ)/(π m)
  sieve overlap ∫Π_{i∈T}1_{B_i} ↔  Σ_{m: Σ m_i v_i = 0} Π_i b̂(m_i)   (sum over RELATIONS)
  circuit-free                  ↔  no nontrivial small relations ⇒ only m=0 ⇒ INDEPENDENCE
  independence baseline         ↔  good-measure ≈ Π(1-2δ) = (1-2δ)^k → e^{-2} > 0   (Lemma A)

ABSTRACT-FUNCTION (Haskell) DICTIONARY (the inspiration lens):
  ev : ℤ → (ℝ/ℤ → ℝ/ℤ),  ev v t = {v t}        -- curry on v ⇒ the character χ_v
  curry ev at witness t* ⇒ position : Fin k → ℝ/ℤ
  tournament arc a→b := sign( Im χ_{v_a - v_b}(t*) )   -- arcs indexed by the DIFFERENCE set S−S

This file tests the rep-theory proof skeleton for Lemma A and probes where it needs work.
"""

import itertools, math, random
from math import gcd

# --------------------------------------------------------------------------- #
def frac_dist(x):
    """‖x‖ = distance to nearest integer."""
    y = x - math.floor(x)
    return min(y, 1 - y)

def gap(S, N=20000):
    """G(S) = max_t min_{v∈S} ‖v t‖, fine grid (t = a/N)."""
    best = 0.0
    bestt = 0.0
    for a in range(1, N):
        t = a / N
        m = min(frac_dist(v * t) for v in S)
        if m > best:
            best, bestt = m, t
    return best, bestt

def good_measure(S, delta, N=20000):
    """measure{ t∈[0,1) : min_v ‖v t‖ ≥ δ } — the all-lonely set, by grid."""
    cnt = 0
    for a in range(N):
        t = a / N
        if all(frac_dist(v * t) >= delta for v in S):
            cnt += 1
    return cnt / N

# --- relation structure (the rep-theory content) -------------------------- #
def three_term_count(S):
    """#{(a,b,c)∈S³ : a+b=c}  — Schur/fusion triples inside S (a,b unordered)."""
    Sset = set(S)
    c = 0
    for a, b in itertools.combinations(S, 2):
        if a + b in Sset:
            c += 1
        # also a+a=c doublings
    for a in S:
        if 2 * a in Sset:
            c += 1
    return c

def additive_energy(S):
    """E(S) = #{(a,b,c,d)∈S⁴ : a+b=c+d}  (the |p_S|⁴ / 4-term richness)."""
    from collections import Counter
    sums = Counter()
    for a in S:
        for b in S:
            sums[a + b] += 1
    return sum(v * v for v in sums.values())

def small_relations(S, mrange=1):
    """nontrivial integer relations Σ m_i v_i = 0 with coeffs in [-mrange,mrange], not all 0.
    mrange=1 catches 3-term (+1,+1,-1) and 4-term (+1,+1,-1,-1) circuits."""
    rels = []
    k = len(S)
    for m in itertools.product(range(-mrange, mrange + 1), repeat=k):
        if any(m) and sum(mi * vi for mi, vi in zip(m, S)) == 0:
            rels.append(m)
    return rels

# --- the Fourier / character identity for overlaps ------------------------ #
def bhat(m, delta):
    """Fourier coefficient of the arc indicator 1_{‖·‖<δ}: b̂(0)=2δ, b̂(m)=sin(2π m δ)/(π m)."""
    if m == 0:
        return 2 * delta
    return math.sin(2 * math.pi * m * delta) / (math.pi * m)

def overlap_via_relations(T, delta, mrange=6):
    """∫ Π_{i∈T} 1_{B_i} = Σ_{m:Σ m_i v_i=0} Π b̂(m_i), truncated to |m_i|≤mrange.
    The rep-theory identity: overlaps are sums over relations among the speeds in T."""
    total = 0.0
    for m in itertools.product(range(-mrange, mrange + 1), repeat=len(T)):
        if sum(mi * vi for mi, vi in zip(m, T)) == 0:
            prod = 1.0
            for mi in m:
                prod *= bhat(mi, delta)
            total += prod
    return total

def overlap_direct(T, delta, N=40000):
    """direct grid measure of ∩_{i∈T} B_i, to validate the relation-sum identity."""
    cnt = 0
    for a in range(N):
        t = a / N
        if all(frac_dist(v * t) < delta for v in T):
            cnt += 1
    return cnt / N

# --------------------------------------------------------------------------- #
def main():
    print("=" * 76)
    print("S585  Lonely Runner — character/sieve exploration (two lenses)")
    print("=" * 76)

    # --- [1] validate the rep-theory identity: overlaps = sums over relations -
    print("\n[1] REP-THEORY IDENTITY  ∫Π 1_{B_i} = Σ_{Σ m_i v_i=0} Π b̂(m_i)")
    print("    (validates that sieve overlaps are character/relation sums)")
    for T in [(3,), (3, 5), (2, 3, 5), (1, 2, 3)]:
        d = 1 / (len(T) + 1)  # use a δ; for a pair just pick δ=0.15
        d = 0.15
        rel = overlap_via_relations(T, d, mrange=8)
        dirv = overlap_direct(T, d)
        print(f"    T={str(T):14s}  relation-sum={rel:+.4f}  direct={dirv:.4f}  "
              f"indep=Π(2δ)={(2*d)**len(T):.4f}")

    # --- [2] independence baseline vs actual good-measure, circuit-free cases -
    print("\n[2] LEMMA A baseline:  circuit-free ⇒ good-measure ≈ (1-2δ)^k → e^{-2}")
    print("    k | 3term | #rel(|m|≤1) | G-δ (margin) | good-meas | (1-2δ)^k | e^-2")
    random.seed(585)
    rows = []
    for k in range(3, 9):
        delta = 1 / (k + 1)
        # find a circuit-free set (no 3-term, no |m|≤1 relations) of size k
        for _try in range(4000):
            S = sorted(random.sample(range(1, 6 * k + 4), k))
            if three_term_count(S) == 0 and len(small_relations(S, 1)) == 0:
                break
        tt = three_term_count(S)
        nrel = len(small_relations(S, 1))
        G, _ = gap(S)
        gm = good_measure(S, delta)
        base = (1 - 2 * delta) ** k
        rows.append((k, S, tt, nrel, G - delta, gm, base))
        print(f"    {k} |   {tt}   |     {nrel}      |   {G-delta:+.4f}   |  {gm:.4f}  "
              f"|  {base:.4f} | {math.exp(-2):.4f}")
    print(f"    e^-2 = {math.exp(-2):.4f}; circuit-free good-measure should hug the baseline, "
          f"margin G-δ > 0.")

    # --- [3] 3-term-rich vs circuit-free: where the margin collapses ----------
    print("\n[3] LEMMA B regime:  3-term relations shrink the margin toward 0 (the hard configs)")
    print("    config             | 3term | energy | G-δ     | good-meas")
    k = 6
    delta = 1 / (k + 1)
    tests = {
        "circuit-free":      None,   # filled below
        "AP 1..k (max 3term)": list(range(1, k + 1)),
        "geometric 1,2,4,..": [2 ** i for i in range(k)],
        "Sidon-ish":         [1, 2, 5, 11, 22, 33],
    }
    # circuit-free sample
    for _try in range(4000):
        S = sorted(random.sample(range(1, 40), k))
        if three_term_count(S) == 0 and len(small_relations(S, 1)) == 0:
            tests["circuit-free"] = S
            break
    for name, S in tests.items():
        if S is None:
            continue
        tt = three_term_count(S)
        en = additive_energy(S)
        G, _ = gap(S)
        gm = good_measure(S, delta)
        print(f"    {name:18s} |   {tt:2d}  |  {en:4d}  | {G-delta:+.4f} | {gm:.4f}")

    # --- [4] the 4-term-rich-but-3-term-free safety claim ---------------------
    print("\n[4] USER CLAIM: 4-term-rich (high energy) BUT 3-term-free ⇒ still safe (margin>0)")
    print("    searching for high-energy, 3-term-free configs and checking the margin...")
    k = 6
    delta = 1 / (k + 1)
    best = None
    for _try in range(60000):
        S = sorted(random.sample(range(1, 30), k))
        if three_term_count(S) != 0:
            continue
        en = additive_energy(S)
        if best is None or en > best[1]:
            G, _ = gap(S)
            best = (S, en, G - delta)
    S, en, margin = best
    print(f"    highest-energy 3-term-free found: S={S}  energy={en}  "
          f"3term={three_term_count(S)}  margin G-δ={margin:+.4f}")
    print(f"    => {'SAFE (margin>0), supports the claim' if margin > 0 else 'HARD — counterexample?!'}")

    # --- [5] ABSTRACT-FUNCTION lens: tournament arcs from difference characters
    print("\n[5] HASKELL LENS: curry ev at witness t* ⇒ tournament arcs = sign Im χ_{v_a-v_b}(t*)")
    S = [1, 2, 3, 5, 8]
    delta = 1 / (len(S) + 1)
    G, tstar = gap(S)
    print(f"    S={S}, witness t*={tstar:.5f}, G={G:.4f}, δ={delta:.4f}")
    print(f"    arcs (a→b if Im e^{{2πi(v_a−v_b)t*}}>0), indexed by difference set S−S:")
    arcs = []
    for a, b in itertools.permutations(S, 2):
        val = math.sin(2 * math.pi * (a - b) * tstar)
        if val > 0:
            arcs.append((a, b))
    # score sequence of the tournament
    score = {v: 0 for v in S}
    for (a, b) in arcs:
        score[a] += 1
    print(f"    out-degrees (tournament score seq): {sorted(score.values())}")
    print(f"    #arcs={len(arcs)} (of {len(S)*(len(S)-1)//2} pairs); "
          f"difference set size |S−S|={len(set(a-b for a in S for b in S))}")

    # --- [6] THETA-FUNCTION form: good-measure = Σ_{m∈Λ=ker(v)} Π ĝ(m_i) -------
    # ĝ(0)=1-2δ, ĝ(m)=-b̂(m). Constant term = (1-2δ)^k = independence; corrections
    # need a nonzero relation vector. Circuit-free ⇒ Λ has no short vectors ⇒ const dominates.
    def ghat(m, delta):
        return (1 - 2 * delta) if m == 0 else -bhat(m, delta)
    def theta_good(S, delta, M=3):
        """Σ_{m: Σ m_i v_i = 0, |m_i|≤M} Π ĝ(m_i)  — truncated theta sum over the relation lattice."""
        const = 0.0; corr = 0.0
        for m in itertools.product(range(-M, M + 1), repeat=len(S)):
            if sum(mi * vi for mi, vi in zip(m, S)) != 0:
                continue
            term = 1.0
            for mi in m:
                term *= ghat(mi, delta)
            if any(m):
                corr += term
            else:
                const += term
        return const, corr
    print("\n[6] THETA-FUNCTION form  good-measure = Σ_{m∈Λ=ker(v)} Π ĝ(m_i)")
    print("    Lemma A = 'the relation-lattice theta sum is dominated by its constant term'")
    print("    config        | (1-2δ)^k const | corrections | const+corr | actual good | λ1(Λ)=min relation")
    k = 6; delta = 1 / (k + 1)
    probes = {
        "circuit-free": tests["circuit-free"],
        "{1..6} (AP)":  list(range(1, 7)),
        "{7..12} (AP↑)":list(range(7, 13)),
        "geometric":    [2 ** i for i in range(6)],
    }
    for name, S in probes.items():
        const, corr = theta_good(S, delta, M=3)
        gm = good_measure(S, delta)
        rels = small_relations(S, 2)
        lam1 = min((sum(abs(x) for x in r) for r in rels), default=99)  # ℓ1 length of shortest relation
        print(f"    {name:14s}| {const:.4f}         | {corr:+.4f}     | {const+corr:.4f}     | "
              f"{gm:.4f}      | {lam1}")
    print("    => circuit-free: corrections ≈ 0, theta ≈ const > 0 (Lemma A); AP {1..6}: large negative")
    print("       corrections cancel the const → good=0 (tight). λ1(Λ) small = short 3-term relation.")

    # --- [7] THE AP-TRANSLATION FLIP: energy-invariant, 3-term-destroying -----
    print("\n[7] AP-TRANSLATION FLIP: {m+1..m+k} keeps additive energy but kills 3-term relations")
    print("    m | set          | energy | 3term | G-δ      | regime")
    k = 6
    for m in [0, 2, 4, 6, 8]:
        S = list(range(m + 1, m + k + 1))
        en = additive_energy(S)
        tt = three_term_count(S)
        delta = 1 / (k + 1)
        G, _ = gap(S)
        reg = "TIGHT (Lemma B)" if abs(G - delta) < 1e-3 else "safe (Lemma A)"
        print(f"    {m} | {str(S):13s}| {en:4d}   |   {tt}   | {G-delta:+.4f} | {reg}")
    print("    => energy constant (translation-invariant); 3-term count drops to 0 once m≥k;")
    print("       the gap flips TIGHT→safe exactly when the 3-term relations vanish.")

if __name__ == "__main__":
    main()
