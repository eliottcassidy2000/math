#!/usr/bin/env python3
"""HYP-7960 referee — the measure-dominance functor Phi (mac-mini-2026-07-19-S123).

Phi(W) = directed graph on the speeds of W with

    i -> j   iff   P(w_i, w_j) := meas{ t in [0,1) : ||w_i t|| > ||w_j t|| } > 1/2,
    i ~ j    iff   P(w_i, w_j) = 1/2  (exact tie).

All arithmetic exact (Fraction).  Facts this referee checks / measures:

  (T1) dilation invariance: P(ga, gb) = P(a, b).
  (T2) the PARITY-TIE LEMMA (proved by hand, verified here):
       for a, b both odd, the involution t -> 1/2 - t reverses dominance
       pointwise (||a(1/2-t)|| = 1/2 - ||at|| for odd a), hence P(a,b) = 1/2.
       Equivalently: after gcd reduction, ties happen EXACTLY at
       nu_2(a) = nu_2(b) ... this last "exactly" is EMPIRICAL (mixed pairs
       could in principle tie accidentally); the referee reports any
       mixed-parity exact tie as a counterexample to the sharp form.
  (T3) the crossing lattice: sign changes of ||at|| - ||bt|| occur only at
       multiples of 1/(a+b) and 1/(a-b) — the (D,s) pair-sum moduli.
  (S)  the upset census: primitive mixed pairs a < b with P(a,b) > 1/2
       ("slower dominates faster").  Null hypothesis: none (faster always
       weakly dominates).  Any upset is arithmetic signal.
  (F)  family analysis: for each named family, the exact matrix, tie
       structure, upset arcs, cyclic triples among decided triples,
       acyclicity of the decided part, and the recursive dyadic word.

Question under test (boxeph-S110 follow-up (b), opus-S399 lead §7.3): does
Phi carry the metagraph's transitive-pole isolation to the M-floor — i.e. do
the tight families sit in a distinguished (transitive-like) fiber of Phi?
"""

from fractions import Fraction as F
from math import gcd, floor
from itertools import combinations

# ---------------------------------------------------------------- exact core

def dist(x: F) -> F:
    """||x|| = distance to nearest integer, exact."""
    fr = x - floor(x)
    return min(fr, 1 - fr)


def P(a: int, b: int) -> F:
    """meas{t in [0,1): ||at|| > ||bt||}, exact rational."""
    assert a != b and a > 0 and b > 0
    g = gcd(a, b)
    a, b = a // g, b // g          # (T1) dilation invariance by substitution
    pts = set()
    for m in (2 * a, 2 * b, a + b, abs(a - b)):
        if m:
            for k in range(m + 1):
                pts.add(F(k, m))
    pts = sorted(pts)
    tot = F(0)
    for u, v in zip(pts, pts[1:]):
        mid = (u + v) / 2
        if dist(a * mid) > dist(b * mid):
            tot += v - u
    return tot


def nu2(x: int) -> int:
    n = 0
    while x % 2 == 0:
        x //= 2
        n += 1
    return n


# ------------------------------------------------------- referee: T1, T2, S

def referee_pair_theory(limit: int = 40):
    print("== (T1) dilation invariance spot checks ==")
    for (a, b, g) in [(1, 2, 3), (2, 3, 5), (3, 4, 7), (1, 3, 2)]:
        lhs, rhs = P(g * a, g * b), P(a, b)
        assert lhs == rhs, (a, b, g, lhs, rhs)
        print(f"  P({g*a},{g*b}) = P({a},{b}) = {lhs}  OK")

    print(f"\n== (T2)+(S) primitive pair sweep a<b<= {limit} ==")
    odd_ties = mixed_ties = upsets = decided = 0
    upset_list = []
    for b in range(2, limit + 1):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            p = P(a, b)
            both_odd = (a % 2 == 1 and b % 2 == 1)
            if p == F(1, 2):
                if both_odd:
                    odd_ties += 1
                else:
                    mixed_ties += 1
                    print(f"  !! MIXED-PARITY EXACT TIE: ({a},{b})")
            else:
                assert not both_odd, f"parity-tie lemma violated at ({a},{b}): P={p}"
                decided += 1
                if p > F(1, 2):
                    upsets += 1
                    upset_list.append((a, b, p))
    print(f"  both-odd ties: {odd_ties} (ALL both-odd pairs tie: lemma verified)")
    print(f"  mixed-parity exact ties: {mixed_ties}")
    print(f"  decided mixed pairs: {decided}, upsets (a<b, a wins): {upsets}")
    if upset_list:
        print("  upset pairs (slower dominates faster):")
        for a, b, p in upset_list:
            print(f"    ({a:>3},{b:>3})  P={p}  nu2=({nu2(a)},{nu2(b)})")
    return upset_list


# --------------------------------------------------------- family machinery

def phi_matrix(W):
    n = len(W)
    M = {}
    for i in range(n):
        for j in range(i + 1, n):
            M[(i, j)] = P(W[i], W[j])
    return M


def analyze_family(name, W, m_note=""):
    W = sorted(W)
    n = len(W)
    M = phi_matrix(W)

    def rel(i, j):
        """+1 if w_i dominates w_j, -1 if dominated, 0 tie (i<j given)."""
        p = M[(i, j)]
        if p > F(1, 2):
            return 1
        if p < F(1, 2):
            return -1
        return 0

    ties = [(i, j) for (i, j) in M if rel(i, j) == 0]
    # upset arc: slower speed dominates faster speed (i<j sorted => w_i<w_j)
    upset_arcs = [(W[i], W[j], M[(i, j)]) for (i, j) in M if rel(i, j) == 1]

    # decided 3-cycles among fully decided triples
    dec_c3 = 0
    dec_triples = 0
    for i, j, k in combinations(range(n), 3):
        r_ij, r_ik, r_jk = rel(i, j), rel(i, k), rel(j, k)
        if 0 in (r_ij, r_ik, r_jk):
            continue
        dec_triples += 1
        # wins for each in the triple; cyclic iff each has exactly 1 win
        wij = 1 if r_ij == 1 else 0
        wik = 1 if r_ik == 1 else 0
        wjk = 1 if r_jk == 1 else 0
        wins_i = wij + wik
        wins_j = (1 - wij) + wjk
        wins_k = (1 - wik) + (1 - wjk)
        if wins_i == wins_j == wins_k == 1:
            dec_c3 += 1

    # acyclicity of the decided digraph (arcs: dominator -> dominated)
    adj = [[False] * n for _ in range(n)]
    for (i, j) in M:
        r = rel(i, j)
        if r == 1:
            adj[i][j] = True
        elif r == -1:
            adj[j][i] = True
    # Kahn on decided arcs
    indeg = [0] * n
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                indeg[j] += 1
    order, stack = [], [i for i in range(n) if indeg[i] == 0]
    while stack:
        v = stack.pop()
        order.append(v)
        for j in range(n):
            if adj[v][j]:
                indeg[j] -= 1
                if indeg[j] == 0:
                    stack.append(j)
    acyclic = len(order) == n

    # recursive dyadic word: parity pattern in increasing order, recursing
    # on the even sublayer (divided by 2) — the tie-clique skeleton of Phi
    def word(ws):
        if not ws:
            return ""
        s = "".join("O" if w % 2 else "E" for w in sorted(ws))
        ev = [w // 2 for w in sorted(ws) if w % 2 == 0]
        return s + ("|" + word(ev) if ev else "")

    print(f"\n-- {name}  W={W}")
    if m_note:
        print(f"   M (cited): {m_note}")
    print(f"   ties: {len(ties)}  (both-odd predicted: "
          f"{sum(1 for (i,j) in M if W[i]%2==1 and W[j]%2==1 and (W[i]|W[j])>=1)} check below)")
    both_odd_pairs = sum(
        1 for (i, j) in M
        if (W[i] // gcd(W[i], W[j])) % 2 == 1 and (W[j] // gcd(W[i], W[j])) % 2 == 1
    )
    print(f"   equal-nu2-after-reduction pairs: {both_odd_pairs} "
          f"(= tie count? {'YES' if both_odd_pairs == len(ties) else 'NO'})")
    print(f"   upset arcs (slower dominates faster): "
          f"{[(a, b, str(p)) for a, b, p in upset_arcs] if upset_arcs else 'NONE'}")
    print(f"   decided triples: {dec_triples}, cyclic among them: {dec_c3}")
    print(f"   decided digraph acyclic: {acyclic}")
    print(f"   dyadic word: {word(W)}")
    return {
        "ties": len(ties), "upsets": len(upset_arcs), "c3": dec_c3,
        "acyclic": acyclic,
    }


def referee_deficit_law(limit: int = 60):
    """THE DOMINANCE DEFICIT LAW (found + proved this session):

        for coprime a < b:   P(a,b) = 1/2 - delta(a,b),
        delta = 1/(2*(b^2-a^2))  if a+b odd (mixed parity),
        delta = 0                if a,b both odd.

    Proof (sketch, Fourier): P(a,b) = meas{gamma(t) in D} where gamma = (at,bt)
    and D = {||x||+||y|| < 1/2} is the L1-diamond (via the half-shift identity
    sigma: t -> t+1/2, valid mixed-parity; the torus is the 2-diamond
    checkerboard).  sign(||x||-||y||) has coefficients -4/(pi^2 (m^2-n^2)) on
    m+n odd, 0 else; restricting to the geodesic line (m,n) = k(b,-a) leaves
    k(b-a) odd, i.e. k odd AND d odd; sum over odd k of 1/k^2 = pi^2/8 gives
    integral = -[d odd]/(ds).  Fejer + dominated convergence make it airtight.
    This referee checks the law EXACTLY on every primitive pair to `limit`.
    """
    print(f"\n== (T4) THE DOMINANCE DEFICIT LAW, exact check a<b<= {limit} ==")
    checked = 0
    for b in range(2, limit + 1):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            p = P(a, b)
            d, s = b - a, a + b
            if (a % 2 == 1) and (b % 2 == 1):
                want = F(1, 2)
            else:
                want = F(1, 2) - F(1, 2 * d * s)
            assert p == want, (a, b, p, want)
            checked += 1
    print(f"  {checked} primitive pairs: P(a,b) == 1/2 - [mixed]/(2(b^2-a^2))  ALL EXACT")


def active_pair_parity_census():
    """Rung-lattice coordinate: slack-1 rungs s=14D-1 are odd => any active
    (straddling) pair u+v=s is MIXED parity => strict dominance deficit
    1/(2*d*s); floor rungs s=14D are even => both-odd balanced (tie, delta=0)
    active pairs are admissible.  Census over known named families."""
    print("\n== (T5) active-pair parity census on named (D,s) data ==")
    named = [
        ("AP {1..13}",            1, 14, (1, 13)),
        ("GW {1..11,13,24}",      1, 14, (1, 13)),
        ("deep well {1..12,182}", 14, 183, (1, 182)),
        ("{1..11,13,36}",         3, 41, (5, 36)),
        ("{1..12,26}",            2, 27, (1, 26)),
        ("{1..29,31,120} (N=32)", 4, 127, (7, 120)),
    ]
    for name, D, s, (u, v) in named:
        assert u + v == s, name
        mixed = (u + v) % 2 == 1
        d = v - u
        delta = F(1, 2 * d * s) if mixed else F(0)
        print(f"  {name:<24} M=D/s={D}/{s:<4} s parity: {'ODD  -> mixed pair, delta=1/'+str(2*d*s) if mixed else 'EVEN -> both-odd tie possible, delta=' + ('0 (pair ('+str(u)+','+str(v)+') both odd)' if u%2==1 and v%2==1 else str(delta))}")
    print("  slack-1 rungs in the gap (D=4..7): s = 55, 69, 83, 97 — ALL ODD")
    print("  => every slack-1 realizer's active pair is mixed-parity (nonzero deficit);")
    print("  the tight floor s=14D is EVEN (both-odd balanced straddle admissible).")


def main():
    upsets = referee_pair_theory(40)
    referee_deficit_law(60)
    active_pair_parity_census()

    print("\n== (F) family analysis ==")
    fams = [
        ("AP (tight #1)",            list(range(1, 14)),        "1/14 (tight)"),
        ("GW doubling (tight #2)",   list(range(1, 12)) + [13, 24], "1/14 (tight)"),
        ("deep well (covering min)", list(range(1, 13)) + [182], "14/183"),
        ("slack-1 D=3",              list(range(1, 12)) + [13, 36], "3/41"),
        ("slack-1 D=2",              list(range(1, 13)) + [26],  "2/27"),
        ("near-AP {1..12,14}",       list(range(1, 13)) + [14],  "1/13 (contains {1..12} core)"),
        ("translate {2..14}",        list(range(2, 15)),         "(runs up, THM-1225)"),
        ("translate {3..15}",        list(range(3, 16)),         "(runs up, THM-1225)"),
        ("k=10 multi-killer extr.",  list(range(1, 11)) + [13, 22, 84], "(covering, L=2227/105105)"),
        ("BF38 band family",         [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35], "1/8 (killed by spread)"),
        ("powers of 2",              [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096], "~0.33 (loose)"),
        ("odds only (non-covering)", [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25], "1/2 at t=1/2"),
    ]
    rows = []
    for name, W, note in fams:
        rows.append((name, analyze_family(name, W, note)))

    print("\n== summary table ==")
    print(f"{'family':<28} {'ties':>5} {'upsets':>7} {'dec-c3':>7} {'acyclic':>8}")
    for name, r in rows:
        print(f"{name:<28} {r['ties']:>5} {r['upsets']:>7} {r['c3']:>7} {str(r['acyclic']):>8}")

    print("\nDone. Interpretation belongs to HYP-7960 in the hypotheses INDEX.")


if __name__ == "__main__":
    main()
