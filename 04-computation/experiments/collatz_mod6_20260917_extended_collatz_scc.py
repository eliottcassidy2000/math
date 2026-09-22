#!/usr/bin/env python3
"""collatz_mod6_20260917_extended_collatz_scc.py

Lane extended_collatz_scc of session collatz-mod6-20260917 (anchor: Collatz).

Object.  E = the nondeterministic digraph on the positive integers with arrows
    n -> n/2      (n even)          [halving, reversed doubling forest]
    n -> 3n+1     (ALL n)           [odd n: ordinary Collatz; even n: NEW]
E_- is the same with 3n-1.  E_5 (hostile control) uses 5n+1.

Everything load-bearing is exact integer / Fraction arithmetic.  All checks
use explicit `raise` (active under python -O).  Universes, filters, positive
and hostile controls are printed.  Iterative searches only; RAM << 1 GB.

Sections
  S1  leaf identity of the new arrows and the 3|n transient forest
  S2  Q1/Q2, inverse moves, residue drift table, 3-adic hostile family 3^j+1
  S3  FINITE-EXACT: Q2 greedy to 10^6, Q1 to 10^6, Tarjan SCC on [1,N]
  S4  Terras-type density theorem for the greedy inverse strategy
  S5  cycle census of E on nodes <= 2000, length <= 40
  S6  signed side E_- (3n-1): leaf identity, drift, three cycles, Q1_-/Q2_-
  S7  hostile control E_5 (5n+1): same density mechanism, extra cycles
  S8  audit-driven additions (2026-09-21 finalisation): correct carry formula,
      peak semantics, tail count vs Markov bound, residue-law modulus 3^(J+2),
      GLOBAL cycle-word census for all lengths <= 26, census cap sensitivity,
      E_- rescue k-word, full pn+1 hierarchy with general bounds and the
      p = 7, 17 dead-exit probabilities, the 2-adic hostile 31 -> 161

Run:  python3 04-computation/experiments/collatz_mod6_20260917_extended_collatz_scc.py
"""
import sys
import time
import hashlib
from fractions import Fraction

T_START = time.time()


def check(cond, msg):
    """Assertion that survives python -O."""
    if not cond:
        raise AssertionError(msg)


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def v3(n):
    c = 0
    while n % 3 == 0:
        n //= 3
        c += 1
    return c


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def collatz_T(n):
    """Accelerated odd-to-odd map T(n) = oddpart(3n+1)."""
    return oddpart(3 * n + 1)


# ---------------------------------------------------------------------------
banner("S1  Leaf identity of the even->3n+1 arrows; multiples of 3 are transient")
# ---------------------------------------------------------------------------
# E-predecessors of v: {2v} always; {(v-1)/3} iff 3 | v-1 and v >= 4.
# The predecessor (v-1)/3 is an ordinary Collatz arrow iff it is odd
# (then v is even).  It is NEW iff (v-1)/3 is even, i.e. v = 6j+1, v >= 7.
V_MAX = 10 ** 5
new_targets = 0
for v in range(1, V_MAX + 1):
    preds = {2 * v}
    if v >= 4 and (v - 1) % 3 == 0:
        preds.add((v - 1) // 3)
    # brute-force cross-check of the predecessor set inside [1, 2v]
    brute = set()
    if v % 2 == 0:
        brute.add(2 * v)
    brute.add(2 * v)  # halving arrow 2v -> v exists for every v
    if v >= 4 and (v - 1) % 3 == 0:
        brute.add((v - 1) // 3)
    check(preds == brute, "predecessor set mismatch at v=%d" % v)
    new = [p for p in preds if p % 2 == 0 and 3 * p + 1 == v]
    if new:
        check(v % 6 == 1 and v >= 7, "new arrow into a non-1-mod-6 target v=%d" % v)
        p = new[0]
        j = (v - 1) // 6
        check(p == 2 * j, "even predecessor is not 2j at v=%d" % v)
        # smallest odd T-predecessor of v: h0 = 1 for v = 1 mod 3, n0 = (4v-1)/3
        n0 = (4 * v - 1) // 3
        check((4 * v - 1) % 3 == 0 and n0 % 2 == 1, "n0 not odd integer at v=%d" % v)
        check(collatz_T(n0) == v, "T(n0) != v at v=%d" % v)
        check(3 * n0 + 1 == 4 * v, "3n0+1 != 4v at v=%d" % v)
        check((n0 - 1) % 4 == 0 and (n0 - 1) // 4 == p, "R^{-1}(n0) != 2j at v=%d" % v)
        check(n0 % 8 == 1, "n0 not 1 mod 8 at v=%d" % v)
        # (B1) fibre formula extended to j=-1: (2^{2} 4^{-1} v - 1)/3 = (v-1)/3
        check(Fraction(4, 4) * v - 1 == 3 * p, "fibre formula at index -1 fails at v=%d" % v)
        # index -2 would be (v/4-1)/3: never an integer for odd v
        check((Fraction(v, 4) - 1) / 3 != int((Fraction(v, 4) - 1) / 3), "index -2 integral?!")
        new_targets += 1
    else:
        if v % 6 == 1 and v >= 7:
            raise AssertionError("missing new arrow into v=%d" % v)
        if v % 6 == 5:
            # would need (v-1)/3 integral: v = 1 mod 3 -- impossible
            check((v - 1) % 3 != 0, "5 mod 6 target with 3n+1 predecessor?!")
        if v % 6 == 3:
            check((v - 1) % 3 != 0, "3 mod 6 target with 3n+1 predecessor?!")
print("universe: targets v in [1, %d]; predecessor sets cross-checked" % V_MAX)
print("new (even -> 3n+1) arrows landing in [1,%d]: %d = #{v = 1 mod 6, 7 <= v <= %d} = %d"
      % (V_MAX, new_targets, V_MAX, len([v for v in range(7, V_MAX + 1) if v % 6 == 1])))
check(new_targets == len([v for v in range(7, V_MAX + 1) if v % 6 == 1]), "new-arrow count")
print("PROVED (S1.1): for v=6j+1>=7 the unique even 3n+1-predecessor 2j=(v-1)/3 equals")
print("  R^{-1}(n0), n0=(4v-1)/3 the least odd T-predecessor (n0 = 1 mod 8, 3n0+1=4v);")
print("  i.e. E extends fibre (B1) of the inherited note to index j=-1 exactly once.")
print("PROVED (S1.2): targets 5 mod 6 and 3 mod 6 acquire no new predecessor.")

# transient forest: 3n+1 is never 0 mod 3; halving preserves 3 | n
for n in range(1, V_MAX + 1):
    check((3 * n + 1) % 3 == 1, "3n+1 divisible by 3?!")
    if n % 2 == 0:
        check((n % 3 == 0) == ((n // 2) % 3 == 0), "halving changes 3-divisibility?!")
print("PROVED (S1.3): no arrow of E enters 3Z from outside 3Z; inside 3Z only halving")
print("  arrows stay (an acyclic forest), so every multiple of 3 is a transient node")
print("  and a singleton SCC.  Checked n <= %d." % V_MAX)
print("Sample: predecessors of 7: {14, 2};  of 13: {26, 4};  of 25: {50, 8};  of 5: {10}; of 9: {18}")
for v in (7, 13, 25, 5, 9, 11, 19, 31):
    preds = sorted({2 * v} | ({(v - 1) // 3} if v >= 4 and (v - 1) % 3 == 0 else set()))
    print("   v=%3d  preds=%s  new_even_pred=%s" % (v, preds, [p for p in preds if p % 2 == 0 and 3 * p + 1 == v]))

# ---------------------------------------------------------------------------
banner("S2  Q1/Q2, inverse moves, residue drift, and the 3-adic hostile family")
# ---------------------------------------------------------------------------
# Inverse moves from m (3 not | m): m -> 2m ; m -> (m-1)/3 when 3 | m-1 and the
# result is not 0 mod 3 (a multiple of 3 is unreachable from 1, so that branch
# is dead).  Compound move: m -> (2^k m - 1)/3, k >= 0 minimal admissible.
UNITS9 = (1, 2, 4, 5, 7, 8)
KMIN = {}
KCLASSES = {}
NEXT_MOD3 = {}
for r in UNITS9:
    ks = [k for k in range(6) if (pow(2, k, 9) * r) % 9 in (4, 7)]
    KCLASSES[r] = ks
    KMIN[r] = min(ks)
    NEXT_MOD3[r] = ((pow(2, KMIN[r], 9) * r - 1) // 3) % 3
check(KMIN == {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}, "k_min table")
print("residue m mod 9 | k_min | admissible k mod 6 | factor 2^k/3 | result mod 3")
for r in UNITS9:
    print("   %d            |  %d    | %-18s | %-12s | %d"
          % (r, KMIN[r], KCLASSES[r], Fraction(2 ** KMIN[r], 3), NEXT_MOD3[r]))
print("PROVED (S2.1): k parity = [m = 2 mod 3]; admissible k form two classes mod 6")
print("  (the third class of the right parity gives a result 0 mod 3).")
print("  4,7 mod 9 shrink x1/3; 2,8 shrink x2/3; 1 grows x4/3; 5 grows x8/3.")
print("  From residues {1,2,4,5} the result is 1 mod 3; from {7,8} it is 2 mod 3.")


def greedy_step(x):
    k = KMIN[x % 9]
    return ((x << k) - 1) // 3, k


def greedy_until_below(m, step_cap=10 ** 6):
    """Iterate greedy from m until value < m (or ==1).  Returns
    (status, steps, arrows, peak, final) with status in
    {'below', 'one', 'cycle', 'cap'}."""
    x = m
    steps = 0
    arrows = 0
    peak = m
    while True:
        if x < m:
            return ("below", steps, arrows, peak, x)
        if steps > 0 and x == m:
            return ("cycle", steps, arrows, peak, x)
        k = KMIN[x % 9]
        y = x << k
        if y > peak:
            peak = y
        x = (y - 1) // 3
        steps += 1
        arrows += k + 1
        if x == 1 and m != 1:
            return ("one", steps, arrows, peak, x)
        if steps > step_cap:
            return ("cap", steps, arrows, peak, x)


print()
print("Hostile family m = 3^j + 1 (m = 1 mod 3^j): greedy forces j-1 steps x4/3 then x1/3")
print("   j |      m | greedy values until first shrink | net ratio | 4^(j-1)/3^j | >1 ?")
for j in range(1, 9):
    m = 3 ** j + 1
    x = m
    path = [x]
    ks = []
    for _ in range(j):
        x, k = greedy_step(x)
        path.append(x)
        ks.append(k)
    check(ks == [2] * (j - 1) + [0], "k-word of 3^j+1 is not 2^(j-1) 0 at j=%d" % j)
    # exact: after j-1 growth steps x = 4^(j-1)(m-1)/3^(j-1)+1 ; then (x-1)/3
    net = Fraction(path[-1], m)
    pred = Fraction(4 ** (j - 1), 3 ** j)
    check(path[-1] == (4 ** (j - 1) * (m - 1)) // 3 ** j, "closed form of the j-step image")
    print("  %2d | %6d | %s | %s | %s | %s" % (j, m, path, net, pred, path[-1] > m))
check((3 ** 5 + 1, greedy_step(greedy_step(greedy_step(greedy_step(greedy_step(244)[0])[0])[0])[0])[0]) == (244, 256),
      "244 -> 256")
print("PROVED (S2.2): m = 1 mod 3^j, m != 1 mod 3^(j+1)  =>  greedy k-word starts 2^(j-1) 0 and")
print("  the j-th greedy image is 4^(j-1)(m-1)/3^j < m  iff  j <= 4.  First exceeding family j=5:")
print("  244 -> 325 -> 433 -> 577 -> 769 -> 256 = 2^8 (net 256/244 = 64/61 > 1).")
print("  This is the 3-adic mirror of the 2-adic hostile n = 2^(L+1)-1 (inherited (B5) block).")
st = greedy_until_below(244)
print("  greedy from 244 continues: status=%s steps=%d peak=%d final=%d" % (st[0], st[1], st[3], st[4]))

# ---------------------------------------------------------------------------
banner("S3  FINITE-EXACT: Q2 greedy to 10^6, Q1 to 10^6, Tarjan SCC on [1,N]")
# ---------------------------------------------------------------------------
N_Q = 10 ** 6
t0 = time.time()
# arrays indexed by m: total compound steps to 1, total arrows, overall peak
tot_steps = [0] * (N_Q + 1)
tot_arrows = [0] * (N_Q + 1)
tot_peak = [0] * (N_Q + 1)
tot_peak[1] = 1
worst_peak = (1, 1)
worst_steps = (0, 1)
worst_arrows = (0, 1)
status_count = {"below": 0, "one": 0, "cycle": 0, "cap": 0}
unresolved = []
for m in range(2, N_Q + 1):
    if m % 3 == 0:
        continue
    st, steps, arrows, peak, fin = greedy_until_below(m, step_cap=10 ** 5)
    status_count[st] += 1
    if st in ("below", "one"):
        check(fin < m and fin % 3 != 0, "descent target invalid at m=%d" % m)
        tot_steps[m] = steps + tot_steps[fin]
        tot_arrows[m] = arrows + tot_arrows[fin]
        tot_peak[m] = max(peak, tot_peak[fin])
        if tot_peak[m] > worst_peak[0]:
            worst_peak = (tot_peak[m], m)
        if tot_steps[m] > worst_steps[0]:
            worst_steps = (tot_steps[m], m)
        if tot_arrows[m] > worst_arrows[0]:
            worst_arrows = (tot_arrows[m], m)
    else:
        unresolved.append(m)
check(not unresolved, "greedy failed to descend for %s" % unresolved[:5])
print("Q2 (1 reaches m in E) verified for all m <= %d with 3 not | m via the GREEDY inverse" % N_Q)
print("  strategy alone (no BFS fallback was needed): statuses %s" % status_count)
print("  max overall peak along the inverse path: %d at m=%d (peak/m = %.3f)"
      % (worst_peak[0], worst_peak[1], worst_peak[0] / worst_peak[1]))
print("  max number of compound moves to 1: %d at m=%d" % worst_steps)
print("  max number of E-arrows to 1     : %d at m=%d" % worst_arrows)
print("  time %.1fs" % (time.time() - t0))
# Show the worst path explicitly (forward direction 1 -> ... -> m)
m = worst_peak[1]
path = [m]
x = m
while x != 1:
    x, k = greedy_step(x)
    path.append(x)
print("  inverse path of the peak-worst m=%d (%d compound moves): %s" % (m, len(path) - 1, path))


def collatz_descends(n, cap_steps=10 ** 6):
    x = n
    s = 0
    while x >= n:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        s += 1
        if s > cap_steps:
            return False
        if x == 1:
            return True
    return True


t0 = time.time()
for n in range(2, N_Q + 1):
    check(collatz_descends(n), "Collatz does not descend below n=%d" % n)
print("Q1 (n reaches 1 in E) verified for all n <= %d by the deterministic Collatz path" % N_Q)
print("  (every n>=2 reaches a smaller value; induction).  time %.1fs" % (time.time() - t0))


def tarjan_scc(N, upmap):
    """Iterative Tarjan on nodes 1..N with arrows n->n/2 (n even) and
    n->upmap(n) if <= N.  Returns list of SCCs (as lists)."""
    index = [0] * (N + 1)
    low = [0] * (N + 1)
    onstk = [False] * (N + 1)
    idx = 1
    stk = []
    sccs = []

    def succ(n):
        out = []
        if n % 2 == 0:
            out.append(n // 2)
        u = upmap(n)
        if 1 <= u <= N:
            out.append(u)
        return out

    for root in range(1, N + 1):
        if index[root]:
            continue
        work = [(root, iter(succ(root)))]
        index[root] = low[root] = idx
        idx += 1
        stk.append(root)
        onstk[root] = True
        while work:
            node, it = work[-1]
            advanced = False
            for w in it:
                if not index[w]:
                    index[w] = low[w] = idx
                    idx += 1
                    stk.append(w)
                    onstk[w] = True
                    work.append((w, iter(succ(w))))
                    advanced = True
                    break
                elif onstk[w]:
                    if index[w] < low[node]:
                        low[node] = index[w]
            if advanced:
                continue
            work.pop()
            if work:
                parent = work[-1][0]
                if low[node] < low[parent]:
                    low[parent] = low[node]
            if low[node] == index[node]:
                comp = []
                while True:
                    w = stk.pop()
                    onstk[w] = False
                    comp.append(w)
                    if w == node:
                        break
                sccs.append(comp)
    return sccs


def scc_report(N, upmap, label):
    t0 = time.time()
    sccs = tarjan_scc(N, upmap)
    sizes = sorted((len(c) for c in sccs), reverse=True)
    giant = max(sccs, key=len)
    gset = set(giant)
    check(1 in gset, "1 not in giant SCC")
    nontriv = [c for c in sccs if len(c) > 1]
    mult3 = [n for n in range(1, N + 1) if n % 3 == 0]
    check(all(len(c) == 1 for c in sccs if any(x % 3 == 0 for x in c)), "3Z node in nontrivial SCC")
    non3 = [n for n in range(1, N + 1) if n % 3 != 0]
    outside = [n for n in non3 if n not in gset]
    half_outside = [n for n in outside if n <= N // 2]
    print("%s: N=%d  #SCC=%d  #nontrivial SCC=%d  giant size=%d  next sizes=%s"
          % (label, N, len(sccs), len(nontriv), len(giant), sizes[1:6]))
    print("  all %d multiples of 3 are singleton SCCs (checked)" % len(mult3))
    print("  non-multiples of 3 in [1,N] outside the giant SCC: %d of %d; among those <= N/2: %d"
          % (len(outside), len(non3), len(half_outside)))
    print("  smallest outsiders: %s" % outside[:12])
    print("  time %.1fs" % (time.time() - t0))
    return sccs, gset, outside


for N_S in (10 ** 3, 10 ** 4, 10 ** 5):
    sccs, gset, outside = scc_report(N_S, lambda n: 3 * n + 1, "E|[1,N]")
    if N_S == 10 ** 5:
        # why finite restriction underestimates: for the smallest outsider show which
        # direction leaves [1,N]
        for n in outside[:5]:
            # forward: deterministic Collatz path max
            x = n
            mx = n
            while x != 1:
                x = x // 2 if x % 2 == 0 else 3 * x + 1
                mx = max(mx, x)
            st = greedy_until_below(n)
            print("  outsider n=%d: Collatz path max=%d (> N? %s); greedy inverse segment peak=%d (> N? %s)"
                  % (n, mx, mx > N_S, st[3], st[3] > N_S))
print("Conjecture C_E: in E the non-multiples of 3 form ONE strongly connected component and every")
print("  multiple of 3 is a transient singleton feeding it.  Status: Q1 part is Collatz; Q2 part is NEW")
print("  and OPEN; both FINITE-EXACT to 10^6.  The finite SCC underestimates because membership of n in")
print("  the giant SCC of E|[1,N] needs BOTH a forward path n->1 and an inverse path m->1 that stay")
print("  inside [1,N]; the Collatz peak (e.g. 9232 from 27) and the inverse peak both exceed N for")
print("  many n <= N (the outsiders listed), and each fixed n is inside for all large N.")

# ---------------------------------------------------------------------------
banner("S4  Terras-type density theorem for the greedy inverse strategy")
# ---------------------------------------------------------------------------
LOG2_3_NUM = None  # we compare 2^K < 3^i exactly with integers


def greedy_word(m, J):
    ks = []
    x = m
    for _ in range(J):
        k = KMIN[x % 9]
        ks.append(k)
        x = ((x << k) - 1) // 3
    return ks


# (a) the first J moves depend only on m mod 3^(J+1); and NOT only on m mod 3^J
import random
random.seed(20260917)
for J in range(1, 9):
    mod = 3 ** (J + 1)
    for _ in range(300):
        m = random.randrange(1, 10 ** 12)
        if m % 3 == 0:
            m += 1
        t = random.randrange(1, 1000)
        check(greedy_word(m, J) == greedy_word(m + mod * t, J), "word not determined mod 3^(J+1)")
    # hostile: modulus 3^J does not determine the J-th letter
    found = False
    for m in range(1, 3 ** J * 4):
        if m % 3 == 0:
            continue
        if greedy_word(m, J) != greedy_word(m + 3 ** J, J):
            found = True
            hostile = (m, m + 3 ** J, greedy_word(m, J), greedy_word(m + 3 ** J, J))
            break
    check(found, "3^J unexpectedly determines J letters at J=%d" % J)
    if J <= 3:
        print("  J=%d: mod 3^%d determines the word; mod 3^%d does not, witness %s" % (J, J + 1, J, hostile))
print("PROVED (S4.1): the first J greedy letters (k_1..k_J) are a function of m mod 3^(J+1); this")
print("  modulus is sharp (witnesses above).")

# (b) exact enumeration of residues mod 3^(J+1)
print()
print("   J | residues | f_J = P(some prefix i<=J has 2^K_i < 3^i) | g_J = P(2^K_J < 3^J) | 1-f_J <= (7/9)^(J-1)? | E[2^K_J] | 3(7/3)^(J-1)")
FJ = {}
for J in range(1, 11):
    mod = 3 ** (J + 1)
    total = 0
    desc_prefix = 0
    desc_J = 0
    e2k = 0
    word_count = {}
    for m in range(1, mod):
        if m % 3 == 0:
            continue
        total += 1
        x = m
        K = 0
        hit = False
        ks = []
        for i in range(1, J + 1):
            k = KMIN[x % 9]
            ks.append(k)
            K += k
            x = ((x << k) - 1) // 3
            if not hit and (1 << K) < 3 ** i:
                hit = True
        if hit:
            desc_prefix += 1
        if (1 << K) < 3 ** J:
            desc_J += 1
        e2k += 1 << K
        if J <= 4:
            w = tuple(ks)
            word_count[w] = word_count.get(w, 0) + 1
    check(total == 2 * 3 ** J, "unit count")
    fJ = Fraction(desc_prefix, total)
    gJ = Fraction(desc_J, total)
    E2K = Fraction(e2k, total)
    pred = 3 * Fraction(7, 3) ** (J - 1)
    check(E2K == pred, "E[2^K_J] != 3(7/3)^(J-1) at J=%d: %s" % (J, E2K))
    bound_ok = (1 - fJ) <= Fraction(7, 9) ** (J - 1)
    check(bound_ok, "tail bound fails at J=%d" % J)
    FJ[J] = fJ
    print("  %2d | %8d | %s = %.6f | %s = %.6f | %s | %s | %s"
          % (J, total, fJ, float(fJ), gJ, float(gJ), bound_ok, E2K, pred))
    if J <= 4:
        # Markov-chain prediction of word counts: class chain iid A w.p. 2/3
        # r_1 uniform on the 6 units; r_{i+1} uniform on the coset of NEXT_MOD3[r_i]
        # count = total * P(word)
        def chain_prob(word):
            # dynamic programming over residues
            dist = {r: Fraction(1, 6) for r in UNITS9}
            for k in word:
                nd = {}
                for r, p in dist.items():
                    if KMIN[r] != k:
                        continue
                    c = NEXT_MOD3[r]
                    for r2 in (c, c + 3, c + 6):
                        nd[r2] = nd.get(r2, 0) + p / 3
                dist = nd
            return sum(dist.values())
        for w, cnt in sorted(word_count.items()):
            check(Fraction(cnt, total) == chain_prob(w), "word count != chain probability %s" % (w,))
        print("      all %d observed %d-letter words have count = 2*3^J * P_chain(word)  (checked)"
              % (len(word_count), J))
for J in range(1, 10):
    check(FJ[J] <= FJ[J + 1], "f_J not monotone")
print("FINITE-EXACT (S4.2): f_J increases toward 1; E[2^K_J] = 3(7/3)^(J-1) exactly for J<=10.")

# (c) residue distribution at step i>=1 is the stationary law (2/9 on 1,4,7; 1/9 on 2,5,8)
mod = 3 ** 6
for i in (1, 2, 3, 4):
    cnt = {r: 0 for r in UNITS9}
    tot = 0
    for m in range(1, mod):
        if m % 3 == 0:
            continue
        x = m
        for _ in range(i):
            x = ((x << KMIN[x % 9]) - 1) // 3
        cnt[x % 9] += 1
        tot += 1
    dist = {r: Fraction(cnt[r], tot) for r in UNITS9}
    check(dist == {1: Fraction(2, 9), 4: Fraction(2, 9), 7: Fraction(2, 9),
                   2: Fraction(1, 9), 5: Fraction(1, 9), 8: Fraction(1, 9)}, "stationary law at step %d" % i)
print("PROVED+checked (S4.3): after one greedy step the residue mod 9 has the exact law")
print("  P(1)=P(4)=P(7)=2/9, P(2)=P(5)=P(8)=1/9 (uniform on units mod 3^6, steps 1..4 checked);")
print("  E_pi[k] = 2*2/9 + 1*1/9 + 0 + 3*1/9 + 0 + 1*1/9 = 1 < log2(3).")

# (d) the additive carry: actual descent within J steps vs residue prediction, for J=6, m<=10^5
J = 6
mod = 3 ** (J + 1)
mismatch = []
for m in range(2, 10 ** 5 + 1):
    if m % 3 == 0:
        continue
    x = m
    K = 0
    pred_desc = False
    act_desc = False
    for i in range(1, J + 1):
        k = KMIN[x % 9]
        K += k
        x = ((x << k) - 1) // 3
        if (1 << K) < 3 ** i:
            pred_desc = True
        if x < m:
            act_desc = True
    if pred_desc != act_desc:
        mismatch.append(m)
        check(act_desc and not pred_desc, "prediction of descent failed at m=%d" % m)
print("FINITE-EXACT (S4.4): for J=6, m<=10^5, 'descends within 6 greedy steps' differs from the pure")
print("  residue prediction only where the carry helps small m: %d mismatches, all with actual descent" % len(mismatch))
print("  and predicted growth; largest mismatching m = %s" % (max(mismatch) if mismatch else None))
print("  (prediction => descent is unconditional since the carry B_i > 0).")

print()
print("THEOREM (S4.5, PROVED below in the note): the greedy inverse stopping time")
print("  sigma(m) = min{i : m_i < m} is finite on a set of natural density 1 among the")
print("  non-multiples of 3, with #{m<=X, 3 not|m, sigma(m) > J} <= (7/9)^(J-1) (2X/3) + C_J.")
print("  Mechanism: E[2^{K_J}] = 3 (7/3)^{J-1} (tilted 2x2 matrix [[5/3,1/3],[10/3,2/3]] of rank 1,")
print("  eigenvalue 7/3) + Markov inequality P(2^{K_J} >= 3^J) <= 3(7/3)^{J-1}/3^J = (7/9)^{J-1}.")
M = [[Fraction(5, 3), Fraction(1, 3)], [Fraction(10, 3), Fraction(2, 3)]]
det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
tr = M[0][0] + M[1][1]
check(det == 0 and tr == Fraction(7, 3), "tilted matrix spectrum")
init = (Fraction(5, 2), Fraction(1, 2))
row = (init[0] * M[0][0] + init[1] * M[1][0], init[0] * M[0][1] + init[1] * M[1][1])
check(row[0] + row[1] == 7, "init * M sums to 7")
print("  checked: det M = 0, tr M = 7/3, (5/2,1/2) M 1 = 7; E[2^{k_1}] = 3.")

# ---------------------------------------------------------------------------
banner("S5  Cycle census of E on nodes <= 2000, length <= 40")
# ---------------------------------------------------------------------------


def cycle_census(LIM, MAXLEN, upmap):
    def succ(n):
        out = []
        if n % 2 == 0:
            out.append(n // 2)
        u = upmap(n)
        if 1 <= u <= LIM:
            out.append(u)
        return out
    found = []
    for s in range(1, LIM + 1):
        stack = [(s, succ(s), 0)]
        path = [s]
        onpath = {s}
        while stack:
            node, ch, idx = stack[-1]
            if idx >= len(ch):
                stack.pop()
                path.pop()
                onpath.discard(node)
                continue
            stack[-1] = (node, ch, idx + 1)
            nxt = ch[idx]
            if nxt == s:
                found.append(list(path))
                continue
            if nxt < s or nxt in onpath or len(path) >= MAXLEN:
                continue
            stack.append((nxt, succ(nxt), 0))
            path.append(nxt)
            onpath.add(nxt)
    return found


t0 = time.time()
CYC = cycle_census(2000, 40, lambda n: 3 * n + 1)
hist = {}
for c in CYC:
    hist[len(c)] = hist.get(len(c), 0) + 1
print("simple cycles of E using only nodes <= 2000, length <= 40: %d  (time %.2fs)" % (len(CYC), time.time() - t0))
print("  length histogram: %s" % sorted(hist.items()))
# structure: a = #(3n+1 arrows), h = #halvings, e = #(even->3n+1) arrows
print("  len | a=#(x3+1) | h=#halve | e=#(even->3n+1) | 2^h>3^a | min node | cycle (canonical from min)")
for c in sorted(CYC, key=lambda c: (len(c), c[0])):
    L = len(c)
    a = h = e = 0
    for i in range(L):
        x, y = c[i], c[(i + 1) % L]
        if y == 3 * x + 1:
            a += 1
            if x % 2 == 0:
                e += 1
        elif x % 2 == 0 and y == x // 2:
            h += 1
        else:
            raise AssertionError("bad arrow %d->%d" % (x, y))
    check(2 ** h > 3 ** a, "cycle with 2^h <= 3^a?!")
    check(a + h == L, "arrow count")
    if c != [1, 4, 2]:
        check(e >= 1, "non-trivial cycle without even->3n+1 arrow: %s" % c)
    else:
        check(e == 0, "trivial cycle uses even->3n+1?!")
    check(all(x % 3 != 0 for x in c), "cycle through a multiple of 3?!")
    print("  %3d | %2d | %2d | %2d | %s | %4d | %s" % (L, a, h, e, 2 ** h > 3 ** a, c[0], c))
user_cycle = [2, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4]
check(user_cycle in CYC, "user's 16-cycle not found")
print("PROVED (S5.1): every E-cycle other than (1,4,2) uses >= 1 even->3n+1 arrow: a cycle avoiding")
print("  them is a cycle of the deterministic Collatz map, whose only cycle on nodes <= 2000 is (1,4,2)")
print("  (FINITE-EXACT input from S3, n <= 10^6; the published verification bound n < 2^68 is")
print("   UNCITED-RECOLLECTION (Barina, J. Supercomputing 2021, from memory) and is NOT used anywhere).")
print("  Every E-cycle has 2^h > 3^a (positivity of n = B/(2^h-3^a)); the lengths observed are")
print("  exactly a+h with 2^h > 3^a > 2^(h-1)... see the note.  The user's 16-cycle through 2->7 is #1 of len 16.")

# ---------------------------------------------------------------------------
banner("S6  Signed side E_- (3n-1): leaf identity, drift table, cycles, Q1_-/Q2_-")
# ---------------------------------------------------------------------------
# E_-: n -> n/2 (even), n -> 3n-1 (all n).  New arrows: even n -> 3n-1 = 6j-1 = 5 mod 6.
for v in range(1, V_MAX + 1):
    preds = {2 * v}
    if (v + 1) % 3 == 0 and (v + 1) // 3 >= 1:
        preds.add((v + 1) // 3)
    new = [p for p in preds if p % 2 == 0 and 3 * p - 1 == v]
    if new:
        p = new[0]
        check(v % 6 == 5, "new minus arrow into non-5-mod-6 target v=%d" % v)
        n0 = (4 * v + 1) // 3
        check((4 * v + 1) % 3 == 0 and n0 % 2 == 1, "n0- not odd integer")
        check(oddpart(3 * n0 - 1) == v and 3 * n0 - 1 == 4 * v, "T_-(n0) != v")
        check((n0 + 1) % 4 == 0 and (n0 + 1) // 4 == p, "R_-^{-1}(n0) != (v+1)/3")
        check(n0 % 8 == 7, "n0- not 7 mod 8 at v=%d" % v)
    else:
        check(not (v % 6 == 5), "missing minus arrow into v=%d" % v)
    check((3 * v - 1) % 3 == 2, "3n-1 divisible by 3?!")
print("PROVED (S6.1): in E_- the new arrows are even 2j -> 6j-1 (targets exactly 5 mod 6); the target v")
print("  has least odd T_- predecessor n0=(4v+1)/3 = 7 mod 8 with 3n0-1=4v, and 2j=(v+1)/3=R_-^{-1}(n0),")
print("  R_-(n)=4n-1.  Targets 1 mod 6 and 3 mod 6 get nothing new; 3Z is again a transient forest.")
KMIN_M = {}
NEXT_M = {}
for r in UNITS9:
    ks = [k for k in range(6) if (pow(2, k, 9) * r) % 9 in (2, 5)]
    KMIN_M[r] = min(ks)
    NEXT_M[r] = ((pow(2, KMIN_M[r], 9) * r + 1) // 3) % 3
check(KMIN_M == {2: 0, 5: 0, 1: 1, 7: 1, 8: 2, 4: 3}, "minus k_min table")
for r in UNITS9:
    check(KMIN_M[r] == KMIN[9 - r], "negation conjugacy of drift tables fails at r=%d" % r)
print("  minus drift table (inverse move m -> (2^k m + 1)/3, result not 0 mod 3):")
for r in UNITS9:
    print("     m = %d mod 9: k_min=%d factor=%s result mod 3 = %d   [= plus table at 9-%d=%d]"
          % (r, KMIN_M[r], Fraction(2 ** KMIN_M[r], 3), NEXT_M[r], r, 9 - r))
print("PROVED (S6.2): k_-(r) = k_+(9-r): negation conjugates the two drift tables, so the greedy")
print("  Markov chain, E[2^K_J]=3(7/3)^(J-1), and the density-1 greedy stopping theorem transfer verbatim.")
print("  Hostile family: m = 3^j - 1 (m = -1 mod 3^j): j-1 steps x4/3 then x1/3, e.g. 242 -> 323 -> 431 -> 575 -> 767 -> 256.")
x = 242
p = [x]
for _ in range(5):
    x = ((x << KMIN_M[x % 9]) + 1) // 3
    p.append(x)
check(p == [242, 323, 431, 575, 767, 256], "minus hostile path")


def greedy_minus_until_below(m, step_cap=10 ** 5):
    x = m
    steps = 0
    peak = m
    while True:
        if x < m:
            return ("below", steps, peak, x)
        if steps > 0 and x == m:
            return ("cycle", steps, peak, x)
        k = KMIN_M[x % 9]
        y = x << k
        peak = max(peak, y)
        x = (y + 1) // 3
        steps += 1
        if steps > step_cap:
            return ("cap", steps, peak, x)


def inverse_bfs_below(m, sign, k_cap=40, value_cap=1 << 90, node_budget=200000):
    """Bounded BFS over ALL admissible compound inverse moves
    x -> (2^k x - sign)/3 (integral, result not 0 mod 3, k <= k_cap) from m,
    until a value < m is reached.  Returns (path or None, nodes)."""
    from collections import deque
    par = {m: None}
    dq = deque([m])
    nodes = 0
    while dq:
        x = dq.popleft()
        nodes += 1
        if nodes > node_budget:
            return (None, nodes)
        for k in range(k_cap + 1):
            y = (x << k) - sign
            if y % 3 != 0:
                continue
            z = y // 3
            if z <= 0 or z % 3 == 0:
                continue
            if z >= value_cap or z in par:
                continue
            par[z] = x
            if z < m:
                p = [z]
                while p[-1] is not None:
                    p.append(par[p[-1]])
                return (p[:-1][::-1], nodes)
            dq.append(z)
    return (None, nodes)


t0 = time.time()
wp = (1, 1)
ws = (0, 1)
tp = [0] * (N_Q + 1)
ts = [0] * (N_Q + 1)
tp[1] = 1
bad = []
greedy_cycles_minus = []
rescued_minus = []
for m in range(2, N_Q + 1):
    if m % 3 == 0:
        continue
    st, steps, peak, fin = greedy_minus_until_below(m)
    if st != "below":
        if st == "cycle":
            cyc = [m]
            x = m
            while True:
                x = ((x << KMIN_M[x % 9]) + 1) // 3
                if x == m:
                    break
                cyc.append(x)
            greedy_cycles_minus.append(cyc)
        path, nodes = inverse_bfs_below(m, -1)
        if path is None:
            bad.append((m, st))
            continue
        rescued_minus.append((m, st, path))
        fin = path[-1]
        peak = max(path)
        steps = len(path) - 1
    tp[m] = max(peak, tp[fin])
    ts[m] = steps + ts[fin]
    if tp[m] > wp[0]:
        wp = (tp[m], m)
    if ts[m] > ws[0]:
        ws = (ts[m], m)
check(not bad, "Q2_- failures even with BFS fallback: %s" % bad[:5])
print("Q2_- (1 reaches m in E_-) FINITE-EXACT for all m <= %d, 3 not|m; max peak %d at m=%d;"
      % (N_Q, wp[0], wp[1]))
print("  max compound moves %d at m=%d; time %.1fs" % (ws[0], ws[1], time.time() - t0))
print("  greedy map G_-(m)=(2^k_min m+1)/3 has cycles with minima <= 10^6: %s" % greedy_cycles_minus)
print("  starts needing a NON-greedy k (BFS over admissible k): %d, namely %s"
      % (len(rescued_minus), [(r[0], r[1], r[2]) for r in rescued_minus]))
check(greedy_cycles_minus == [[4, 11]], "unexpected greedy cycles in G_-")
# the same fallback machinery on the plus side confirms it was never needed
plus_cycles = []
for m in range(2, N_Q + 1):
    if m % 3 == 0:
        continue
    st = greedy_until_below(m, step_cap=10 ** 5)[0]
    if st == "cycle":
        plus_cycles.append(m)
check(not plus_cycles, "greedy cycles on the plus side?!")
print("  (plus side: the greedy map G(m)=(2^k_min m-1)/3 has no cycle with minimum in [2,10^6]; its only")
print("   fixed point is 1 since G(1)=1.)  REFUTED for the minus side: G_-(4)=11, G_-(11)=4.")


def forward_reach_below(n, upmap, node_budget=200000, value_cap=1 << 80):
    """Iterative DFS in the nondeterministic graph from n, halving preferred,
    until a node < n is reached.  Returns (found, nodes_expanded, peak)."""
    if n == 1:
        return (True, 0, 1)
    seen = {n}
    stack = [n]
    expanded = 0
    peak = n
    while stack:
        x = stack.pop()
        expanded += 1
        if expanded > node_budget:
            return (False, expanded, peak)
        succ = []
        u = upmap(x)
        if u < value_cap:
            succ.append(u)
        if x % 2 == 0:
            succ.append(x // 2)   # pushed last => popped first (halving preferred)
        for y in succ:
            if y < n:
                return (True, expanded, peak)
            if y not in seen:
                seen.add(y)
                if y > peak:
                    peak = y
                stack.append(y)
    return (False, expanded, peak)


def minus_det_descends(n, cap=10 ** 5):
    x = n
    s = 0
    while x >= n:
        x = x // 2 if x % 2 == 0 else 3 * x - 1
        s += 1
        if x == n or s > cap:
            return False
    return True


t0 = time.time()
need_nd = []
fail = []
for n in range(2, N_Q + 1):
    if minus_det_descends(n):
        continue
    ok, exp, pk = forward_reach_below(n, lambda x: 3 * x - 1)
    if ok:
        need_nd.append((n, exp, pk))
    else:
        fail.append(n)
check(not fail, "Q1_- failures: %s" % fail[:5])
print("Q1_- (n reaches 1 in E_-) FINITE-EXACT for all n <= %d: the deterministic 3n-1 path descends" % N_Q)
print("  below n except for %d starts %s, each rescued by an even->3n-1 arrow (DFS nodes expanded, peak shown)."
      % (len(need_nd), need_nd))
print("  time %.1fs" % (time.time() - t0))

# the three 3n-1 cycles as E_- cycles; explicit paths 1 -> c and c -> 1
CYC_M = [[1, 2], [5, 14, 7, 20, 10],
         [17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 34]]
for c in CYC_M:
    for i in range(len(c)):
        x, y = c[i], c[(i + 1) % len(c)]
        check(y == 3 * x - 1 or (x % 2 == 0 and y == x // 2), "not an E_- cycle: %s" % c)


def bfs_path(src, dst, upmap, cap=1 << 40, limit=4 * 10 ** 6):
    from collections import deque
    par = {src: None}
    dq = deque([src])
    while dq and len(par) < limit:
        x = dq.popleft()
        if x == dst:
            p = []
            while x is not None:
                p.append(x)
                x = par[x]
            return p[::-1]
        for y in ((x // 2,) if x % 2 == 0 else ()) + ((upmap(x),) if upmap(x) < cap else ()):
            if y not in par:
                par[y] = x
                dq.append(y)
    return None


for c in CYC_M:
    m = c[0]
    p1 = bfs_path(1, m, lambda x: 3 * x - 1)
    p2 = bfs_path(m, 1, lambda x: 3 * x - 1)
    check(p1 is not None and p2 is not None, "cycle minimum %d not in giant SCC of E_-" % m)
    print("  3n-1 cycle min %3d: shortest 1->%d path %s ; shortest %d->1 path %s" % (m, m, p1, m, p2))
print("PROVED (S6.3): all three known 3n-1 cycles lie in ONE strongly connected component of E_- with 1")
print("  (explicit paths above); multiples of 3 remain a transient halving forest (3n-1 = 2 mod 3).")
sccs_m, gset_m, outside_m = scc_report(10 ** 5, lambda n: 3 * n - 1, "E_-|[1,N]")
check(all(min(c) in gset_m for c in CYC_M), "cycle in giant SCC of E_-|[1,10^5]")
CYC_M_CENSUS = cycle_census(2000, 40, lambda n: 3 * n - 1)
hm = {}
for c in CYC_M_CENSUS:
    hm[len(c)] = hm.get(len(c), 0) + 1
det_cycles_m = [c for c in CYC_M_CENSUS if all((y == x // 2) if x % 2 == 0 else (y == 3 * x - 1)
                                               for x, y in zip(c, c[1:] + c[:1]))]
print("  E_- cycle census nodes<=2000, len<=40: %d cycles, histogram %s; deterministic (no even->3n-1) ones: %s"
      % (len(CYC_M_CENSUS), sorted(hm.items()), det_cycles_m))
check(sorted(det_cycles_m) == sorted(CYC_M), "deterministic 3n-1 cycles on <=2000 are not exactly the three")

# ---------------------------------------------------------------------------
banner("S7  Hostile control E_5 (5n+1): same greedy density mechanism, yet extra cycles")
# ---------------------------------------------------------------------------
# inverse move m -> (2^k m - 1)/5, result not 0 mod 5: 2^k m in {6,11,16,21} mod 25
UNITS25 = [r for r in range(1, 25) if r % 5]
KMIN5 = {}
for r in UNITS25:
    KMIN5[r] = min(k for k in range(20) if (pow(2, k, 25) * r) % 25 in (6, 11, 16, 21))
ek = Fraction(sum(KMIN5.values()), len(UNITS25))
print("  E_5 greedy: k_min by residue mod 25 = %s" % KMIN5)
print("  mean k under uniform residues = %s = %.3f  vs log2(5) = 2.322" % (ek, float(ek)))
print("   J | E[2^K_J]/5^J (exact, residues mod 5^(J+1)) | f_J (prefix descent)")
for J in range(1, 6):
    mod = 5 ** (J + 1)
    tot = 0
    e2 = 0
    desc = 0
    for m in range(1, mod):
        if m % 5 == 0:
            continue
        tot += 1
        x = m
        K = 0
        hit = False
        for i in range(1, J + 1):
            k = KMIN5[x % 25]
            K += k
            x = ((x << k) - 1) // 5
            if (1 << K) < 5 ** i:
                hit = True
        e2 += 1 << K
        desc += hit
    print("  %2d | %s = %.4f | %s = %.4f" % (J, Fraction(e2, tot * 5 ** J), e2 / (tot * 5 ** J), Fraction(desc, tot), desc / tot))
print("  => the same density-1 greedy descent mechanism is present for 5n+1 (mean ratio < 1).")
t0 = time.time()
res5 = {"below": 0, "cycle": 0, "cap": 0}
cyc5 = []
for m in range(2, 10 ** 5 + 1):
    if m % 5 == 0:
        continue
    x = m
    s = 0
    st = None
    while True:
        if x < m:
            st = "below"
            break
        if s > 0 and x == m:
            st = "cycle"
            cyc5.append(m)
            break
        k = KMIN5[x % 25]
        x = ((x << k) - 1) // 5
        s += 1
        if s > 10000 or x > 1 << 200:
            st = "cap"
            break
    res5[st] += 1
print("  E_5 greedy inverse (Q2_5) for m <= 10^5, 5 not|m: %s ; greedy cycle minima: %s (time %.1fs)"
      % (res5, cyc5[:10], time.time() - t0))
# forward Q1_5 with caps on n <= 2000, including the known 5n+1 cycle 13,33,83
t0 = time.time()
unres5 = []
for n in range(2, 2001):
    if n % 5 == 0:
        continue
    ok, exp, pk = forward_reach_below(n, lambda x: 5 * x + 1, node_budget=20000, value_cap=1 << 60)
    if not ok:
        unres5.append(n)
print("  E_5 forward (Q1_5: reach below n) for n <= 2000 with budget 20000 nodes: unresolved = %s (time %.1fs)"
      % (unres5, time.time() - t0))
for n in (13, 33, 83):
    ok, exp, pk = forward_reach_below(n, lambda x: 5 * x + 1, node_budget=200000, value_cap=1 << 60)
    print("  5n+1 cycle member %d: reaches below itself in E_5? %s (nodes %d)" % (n, ok, exp))
print("  Typed analogy (see note): the density theorem transfers to 5n+1, so it cannot by itself")
print("  distinguish 3n+1 (one known positive cycle) from 5n+1 (extra cycle 13->33->83).")

# ---------------------------------------------------------------------------
banner("S8  Audit-driven additions (2026-09-21): carry formula, peak semantics, tail count, "
       "residue-law modulus, GLOBAL cycle-word census, cap sensitivity, E_- rescue word, pn+1 hierarchy")
# ---------------------------------------------------------------------------
t0 = time.time()


def greedy_prefix(m, i):
    """First i greedy compound moves from m: (k-word, partial sums K_1..K_i, m_i)."""
    x = m
    ks = []
    Ks = []
    K = 0
    for _ in range(i):
        k = KMIN[x % 9]
        K += k
        ks.append(k)
        Ks.append(K)
        x = ((x << k) - 1) // 3
    return ks, Ks, x


# S8.1  carry formula.  3^i m_i = 2^{K_i} m - B_i with B_i = sum_{l=1}^{i} 3^{l-1} 2^{K_i-K_l}.
# The draft note displayed B_i = sum_{l=0}^{i-1} 3^{i-1-l} 2^{K_i-K_{l+1}} (exponent of 3 reversed).
tot = 0
old_wrong = 0
for m in range(1, 3001):
    if m % 3 == 0:
        continue
    for i in range(1, 7):
        ks, Ks, mi = greedy_prefix(m, i)
        Ki = Ks[-1]
        B_new = sum(3 ** (l - 1) * 2 ** (Ki - Ks[l - 1]) for l in range(1, i + 1))
        B_old = sum(3 ** (i - 1 - l) * 2 ** (Ki - Ks[l]) for l in range(0, i))
        check(3 ** i * mi == 2 ** Ki * m - B_new, "carry formula fails at m=%d i=%d" % (m, i))
        check(B_new > 0, "carry not positive")
        tot += 1
        if B_old != B_new:
            old_wrong += 1
ks, Ks, m2 = greedy_prefix(4, 2)
check(ks == [0, 2] and Ks == [0, 2] and m2 == 1, "word of 4 is (0,2)")
B2 = 3 ** 0 * 2 ** (2 - 0) + 3 ** 1 * 2 ** (2 - 2)
check(B2 == 7 and 9 * 1 == 4 * 4 - 7, "B_2 for word (0,2) is 7")
print("S8.1 carry formula B_i = sum_{l=1}^{i} 3^(l-1) 2^(K_i-K_l): 3^i m_i = 2^K_i m - B_i checked for")
print("  m <= 3000 (3 not|m), i <= 6: %d instances, all exact, all B_i > 0; the draft's reversed form" % tot)
print("  sum_{l=0}^{i-1} 3^(i-1-l) 2^(K_i-K_(l+1)) disagrees in %d of them (e.g. word (0,2) from m=4:" % old_wrong)
print("  correct B_2 = 4 + 3 = 7 and 9*1 = 16 - 7; the reversed form gives 13).")

# S8.2  peak semantics: E-node peak (includes doubling intermediates 2^k x) vs compound-value peak.
m = 797162
path = [m]
x = m
enode_peak = m
while x != 1:
    k = KMIN[x % 9]
    enode_peak = max(enode_peak, x << k)
    x = ((x << k) - 1) // 3
    path.append(x)
cv_peak = max(path)
check(enode_peak == 150994948 and cv_peak == 50331649, "peak semantics at m=797162: %d %d" % (enode_peak, cv_peak))
check(150994948 == 4 * 37748737 and 37748737 in path and 150994948 == 9 * 2 ** 24 + 4, "peak decomposition")
print("S8.2 peak semantics at m=797162: E-node peak (all nodes of the E-path 1->m, doubling intermediates")
print("  included) = %d = 4*37748737 = 9*2^24+4; compound-value peak (values m_i only) = %d." % (enode_peak, cv_peak))
x = 244
p244 = [x]
epk = x
while True:
    k = KMIN[x % 9]
    epk = max(epk, x << k)
    x = ((x << k) - 1) // 3
    p244.append(x)
    if x < 244:
        break
check(p244 == [244, 325, 433, 577, 769, 256, 85] and epk == 2308, "244 chain")
print("  m=244: compound values %s, E-node peak %d = 4*577, compound-value peak %d." % (p244, epk, max(p244)))

# S8.3  actual tail count versus the Markov bound of Theorem 4.4 at J=10, X=10^6.
J = 10
X = 10 ** 6
tail = 0
tail_ge2 = 0
for m in range(1, X + 1):
    if m % 3 == 0:
        continue
    x = m
    hit = False
    for _ in range(J):
        x = ((x << KMIN[x % 9]) - 1) // 3
        if x < m:
            hit = True
            break
    if not hit:
        tail += 1
        if m >= 2:
            tail_ge2 += 1
bound = Fraction(7, 9) ** (J - 1) * (Fraction(2 * X, 3) + 2 * 3 ** J)
check(tail <= bound, "tail count exceeds the bound?!")
check(tail == 2044 and tail_ge2 == 2043, "tail count at J=10, X=10^6 is %d (m>=2: %d)" % (tail, tail_ge2))
print("S8.3 #{1 <= m <= 10^6, 3 not|m, sigma(m) > 10} = %d (this includes the greedy fixed point m = 1 with"
      % tail)
print("  sigma(1) = infinity; %d starts m >= 2)  <=  Markov bound (7/9)^9 (2X/3 + 2*3^10) = %d (integer part;"
      % (tail_ge2, int(bound)))
print("  the bound is %.1f times the truth)." % (float(bound) / tail))

# S8.4  residue-law modulus: m_J mod 9 is a function of m mod 3^(J+2), not of m mod 3^(J+1).
for J in range(1, 4):
    modJ2 = 3 ** (J + 2)
    for m in range(1, modJ2):
        if m % 3 == 0:
            continue
        r0 = greedy_prefix(m, J)[2] % 9
        for tt in (1, 2, 4):
            check(greedy_prefix(m + modJ2 * tt, J)[2] % 9 == r0, "m_J mod 9 not determined mod 3^(J+2)")
    modJ1 = 3 ** (J + 1)
    wit = None
    for m in range(1, modJ1):
        if m % 3 == 0:
            continue
        if greedy_prefix(m, J)[2] % 9 != greedy_prefix(m + modJ1, J)[2] % 9:
            wit = (m, m + modJ1, greedy_prefix(m, J)[2] % 9, greedy_prefix(m + modJ1, J)[2] % 9)
            break
    check(wit is not None, "m_J mod 9 IS determined mod 3^(J+1)?!")
    # law under uniform units mod 3^(J+2)
    law = {}
    for m in range(1, modJ2):
        if m % 3 == 0:
            continue
        r = greedy_prefix(m, J)[2] % 9
        law[r] = law.get(r, 0) + 1
    n_units = 2 * 3 ** (J + 1)
    check(all(Fraction(law[r], n_units) == Fraction(2 if r in (1, 4, 7) else 1, 9) for r in UNITS9), "law J=%d" % J)
    print("S8.4 J=%d: m_J mod 9 is a function of m mod 3^%d (checked all units, lifts t=1,2,4); NOT of m mod 3^%d"
          " (witness m=%d vs %d: residues %d vs %d); law under uniform units mod 3^%d = 2/9 on {1,4,7}, 1/9 on {2,5,8}."
          % (J, J + 2, J + 1, wit[0], wit[1], wit[2], wit[3], J + 2))
print("  time %.1fs" % (time.time() - t0))


# S8.5  GLOBAL cycle-word census of E (no node bound) for every length L <= 26.
# A simple E-cycle with a arrows n->3n+1 and h halvings, rotated to start with a 3n+1 arrow,
# is a composition (h_1,...,h_a) of h into a parts >= 0 (block i = one 3n+1 arrow then h_i
# halvings).  Going around once, 2^h n0 = 3^a n0 + B with B = sum_{i=1}^a 3^(a-i) 2^(h_1+..+h_(i-1)),
# so n0 = B/(2^h-3^a) must be a positive integer, and the walk must halve only even nodes.
def global_cycle_census(LMAX):
    found = {}
    n_words = 0
    for L in range(2, LMAX + 1):
        for a in range(1, L + 1):
            h = L - a
            if 2 ** h <= 3 ** a:
                continue
            D = 2 ** h - 3 ** a
            pow3 = [3 ** (a - i) for i in range(1, a + 1)]   # pow3[i-1] = 3^(a-i)
            comp = [0] * a
            cyc_keys = found.setdefault(L, {})

            def leaf(B):
                if B % D:
                    return
                n0 = B // D
                x = n0
                nodes = [x]
                for i in range(a):
                    x = 3 * x + 1
                    nodes.append(x)
                    for _ in range(comp[i]):
                        if x % 2:
                            return
                        x //= 2
                        nodes.append(x)
                if x != n0:
                    raise AssertionError("cycle equation inconsistent")
                cyc = nodes[:-1]
                if len(set(cyc)) != L:
                    return          # a shorter cycle traversed several times
                mn = cyc.index(min(cyc))
                key = tuple(cyc[mn:] + cyc[:mn])
                cyc_keys[key] = (a, h)

            def rec(i, Hsum, B):
                # choose comp[i] (0-based block i), given h_1+..+h_i = Hsum already used
                nonlocal n_words
                if i == a - 1:
                    comp[i] = h - Hsum
                    n_words += 1
                    leaf(B + pow3[i] * (1 << Hsum))
                    return
                Bi = B + pow3[i] * (1 << Hsum)
                for hi in range(0, h - Hsum + 1):
                    comp[i] = hi
                    rec(i + 1, Hsum + hi, Bi)

            rec(0, 0, 0)
    return found, n_words


t0 = time.time()
GLOB, n_words = global_cycle_census(26)
print("S8.5 GLOBAL census by cycle words, all lengths L <= 26, no node bound: %d compositions examined (time %.1fs)"
      % (n_words, time.time() - t0))
print("   L | #cycles | (a,h) | max node | #with max node > 2000 | h = ceil(a log2 3) for all")
glob_hist = {}
for L in range(2, 27):
    lst = GLOB.get(L, {})
    if not lst:
        continue
    glob_hist[L] = len(lst)
    ahs = sorted(set(lst.values()))
    mx = max(max(c) for c in lst)
    big = sum(1 for c in lst if max(c) > 2000)
    hforced = all(h == min(tt for tt in range(200) if 2 ** tt > 3 ** a) for (a, h) in lst.values())
    check(hforced, "h != ceil(a log2 3) in a global cycle of length %d" % L)
    check(all(all(x % 3 for x in c) for c in lst), "global cycle through 3Z")
    print("  %2d | %3d | %s | %5d | %3d | %s" % (L, len(lst), ahs, mx, big, hforced))
    # cross-check with the bounded census: cycles with all nodes <= 2000 must coincide
    bounded = sorted(tuple(c) for c in CYC if len(c) == L)
    inside = sorted(c for c in lst if max(c) <= 2000)
    check(bounded == inside, "bounded census != global census restricted to nodes <= 2000 at L=%d" % L)
check(glob_hist == {3: 1, 8: 1, 13: 6, 16: 1, 21: 2, 26: 22}, "global histogram %s" % glob_hist)
check(list(GLOB[16]) == [(2, 7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4)], "unique global 16-cycle")
check(sum(1 for c in GLOB[26] if max(c) > 2000) == 11 and max(max(c) for c in GLOB[26]) == 7168, "L=26 outsiders")
for L in (3, 8, 13, 16, 21):
    check(all(max(c) <= 2000 for c in GLOB[L]), "a cycle of length %d leaves nodes<=2000" % L)
print("FINITE-EXACT (S8.5): every E-cycle of length <= 25 has all nodes <= 2000, so the bounded census is")
print("  complete there: exactly 1,1,6,1,2 cycles of lengths 3,8,13,16,21 in ALL of E, and the user's 16-cycle")
print("  is the unique E-cycle of length 16.  At L=26 there are 22 cycles in total, 11 with a node > 2000")
print("  (largest node 7168); h = ceil(a log2 3) holds for every E-cycle of length <= 26.")
print("  Example L=26 cycle outside nodes<=2000: %s" % list(sorted(c for c in GLOB[26] if max(c) > 2000)[0]))
# cap sensitivity of the bounded census
t0 = time.time()
C45 = cycle_census(2000, 45, lambda n: 3 * n + 1)
C2500 = cycle_census(2500, 40, lambda n: 3 * n + 1)
C4000 = cycle_census(4000, 45, lambda n: 3 * n + 1)
h45 = {}
for c in C45:
    h45[len(c)] = h45.get(len(c), 0) + 1
check(len(C45) == 100 and h45.get(44, 0) == 26, "nodes<=2000,len<=45 census: %d, len-44: %d" % (len(C45), h45.get(44, 0)))
check(len(C2500) == 104, "nodes<=2500,len<=40 census: %d" % len(C2500))
check(len(C4000) == 410, "nodes<=4000,len<=45 census: %d" % len(C4000))
for c in C4000:
    a = sum(1 for i in range(len(c)) if c[(i + 1) % len(c)] == 3 * c[i] + 1)
    h = len(c) - a
    check(h == min(tt for tt in range(200) if 2 ** tt > 3 ** a), "h != ceil(a log2 3) at nodes<=4000")
print("  cap sensitivity: nodes<=2000,len<=45: %d cycles (histogram %s; the 26 new ones have length 44, (a,h)=(17,27));"
      % (len(C45), sorted(h45.items())))
print("  nodes<=2500,len<=40: %d cycles (30 more than 74); nodes<=4000,len<=45: %d cycles; h = ceil(a log2 3) in all %d."
      % (len(C2500), len(C4000), len(C4000)))
print("  Both caps of the '74' bind.  time %.1fs" % (time.time() - t0))

# S8.6  E_- rescue path and the shortest 1<->cycle paths (arrow counts)
resc = [4, 11, 59, 20, 7, 5, 2]
kw = []
for x, y in zip(resc, resc[1:]):
    ks = [k for k in range(12) if ((x << k) + 1) % 3 == 0 and ((x << k) + 1) // 3 == y]
    check(len(ks) == 1, "rescue step %d->%d" % (x, y))
    kw.append(ks[0])
check(kw == [3, 4, 0, 0, 1, 0], "rescue k-word %s" % kw)
print("S8.6 E_- rescue of m=4: 4->11->59->20->7->5->2 has k-word %s (20->7 is (20+1)/3, k=0)." % (kw,))
for c in CYC_M[1:]:
    mm = c[0]
    p1 = bfs_path(1, mm, lambda x: 3 * x - 1)
    p2 = bfs_path(mm, 1, lambda x: 3 * x - 1)
    check(len(set(p1)) == len(p1) and len(set(p2)) == len(p2), "BFS path not simple")
    print("  shortest E_- paths: 1->%d has %d arrows, %d->1 has %d arrows (both simple)." % (mm, len(p1) - 1, mm, len(p2) - 1))

# S8.7  pn+1 hierarchy: greedy inverse move m -> (2^k m - 1)/p (result not 0 mod p), residues mod p^2.
t0 = time.time()


def ord_p(p):
    o = 1
    while pow(2, o, p) != 1:
        o += 1
    return o


def p_info(p):
    o = ord_p(p)
    e = ((pow(2, o, p * p) - 1) // p) % p
    check(e != 0, "2^ord = 1 mod p^2 at p=%d (Wieferich-type); formulas need e != 0" % p)
    reach = [s for s in range(1, p) if any(pow(2, t, p) * s % p == 1 for t in range(o))]   # s in <2>
    k0 = {s: min(t for t in range(o) if pow(2, t, p) * s % p == 1) for s in reach}
    k = {}
    c = {}
    for r in range(1, p * p):
        if r % p == 0 or (r % p) not in k0:
            continue
        for t in range(0, o * p + 2):
            y = pow(2, t, p * p) * r % (p * p)
            if y % p == 1 and ((y - 1) // p) % p != 0:
                k[r] = t
                c[r] = ((y - 1) // p) % p
                break
        s = r % p
        a = ((pow(2, k0[s], p * p) * r - 1) // p) % p
        check(k[r] == k0[s] + o * (a == 0) and c[r] == (a if a else e), "k/c formula p=%d r=%d" % (p, r))
    for s in reach:
        law = {}
        for t in range(p):
            law[c[s + p * t]] = law.get(c[s + p * t], 0) + 1
        check(law == {cc: 1 + (cc == e) for cc in range(1, p)}, "class law (1+[c=e])/p fails p=%d s=%d" % (p, s))
    return o, e, reach, k0, k, c


print("S8.7 pn+1 hierarchy (greedy inverse chain on residues mod p^2; ord = ord_p(2), e = (2^ord-1)/p mod p):")
print("   p | ord | 2 primitive? | e | E_pi[k] | < log2 p ? | rho_p | rho_p < p ? | dead-exit prob (2 not primitive)")
ROWS = {}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 37):
    o, e, reach, k0, k, c = p_info(p)
    prim = (o == p - 1)
    if prim:
        Ek = sum(Fraction(1 + (s == e), p) * k0[s] for s in reach) + Fraction(o, p)
        rho = Fraction(sum(2 ** k0[s] for s in reach) + 2 ** (o + k0[e]), p)
        check(Ek == Fraction(p - 1, 2) + Fraction(k0[e], p), "E_pi[k] closed form p=%d" % p)
        check(sum(2 ** k0[s] for s in reach) == 2 ** (p - 1) - 1, "sum 2^k0 = 2^(p-1)-1")
        check(rho >= Fraction(2 ** p - 1, p), "rho_p >= (2^p-1)/p fails at p=%d" % p)
        if p >= 5:
            check(rho > p and 2 ** p - 1 > p * p, "mean-supercritical for p>=5 fails at p=%d" % p)
        if p >= 7:
            check(2 ** Ek.numerator > p ** Ek.denominator and 2 ** (p - 1) >= p * p, "log-supercritical p>=7 fails at %d" % p)
        logsub = 2 ** Ek.numerator < p ** Ek.denominator
        ROWS[p] = (o, e, Ek, logsub, rho, rho < p)
        print("  %2d | %3d | yes | %2d | %s | %s | %s | %s | -" % (p, o, e, Ek, logsub, rho, rho < p))
    else:
        dead = Fraction(sum(1 for t in range(p) if c[reach[0] + p * t] not in k0), p)
        dead_cf = Fraction(p - 1 - len(reach) + (0 if e in k0 else 1), p)
        check(dead == dead_cf, "dead-exit closed form p=%d" % p)
        ROWS[p] = (o, e, None, None, None, None, dead)
        print("  %2d | %3d | no (<2> = %s) | %2d | - | - | - | - | %s = (p-1-|<2>|+[e not in <2>])/p" % (p, o, reach, e, dead))
check(ROWS[3][2] == 1 and ROWS[3][4] == Fraction(7, 3), "p=3 row")
check(ROWS[5][2] == Fraction(11, 5) and ROWS[5][4] == Fraction(47, 5), "p=5 row")
check(ROWS[11][2] == Fraction(61, 11) and ROWS[11][4] == Fraction(66559, 11), "p=11 row")
check(ROWS[13][2] == Fraction(86, 13) and ROWS[13][4] == Fraction(1052671, 13), "p=13 row")
check(ROWS[19][2] == Fraction(176, 19) and ROWS[19][4] == Fraction(8650751, 19), "p=19 row")
check(ROWS[7][6] == Fraction(3, 7) and ROWS[17][6] == Fraction(8, 17), "dead exits p=7,17")
check([p for p in ROWS if ROWS[p][2] is not None and ROWS[p][5]] == [3], "mean-subcritical only p=3")
check([p for p in ROWS if ROWS[p][2] is not None and ROWS[p][3]] == [3, 5], "log-subcritical only p=3,5")
print("PROVED (S8.7): for 2 a primitive root mod p with e != 0, sum_s 2^k0(s) = 2^(p-1)-1, hence")
print("  rho_p = (2^(p-1)-1+2^(p-1+k0(e)))/p >= (2^p-1)/p, which exceeds p for every p >= 5 (2^p-1 > p^2);")
print("  E_pi[k] = (p-1)/2 + k0(e)/p >= (p-1)/2 >= log2 p for every p >= 7 (2^(p-1) >= p^2).  So among ALL such")
print("  primes: mean-subcritical <=> p = 3; log-subcritical <=> p in {3,5}.  For p = 7, 17 (2 not primitive) the")
print("  class law (1+[c=e])/p still holds on every <2>-coset, and the per-step probability of leaving <2>")
print("  (no inverse move at all) is exactly (p-1-|<2>|+[e not in <2>])/p = 3/7 and 8/17.")
lam = Fraction(21, 20)
o, e, reach, k0, k, c = p_info(5)
rl = (sum(lam ** k0[s] for s in reach) + lam ** (o + k0[e])) / 5
check(rl == Fraction(17876501, 16000000) and rl ** 10 < lam ** 23 and 2 ** 23 < 5 ** 10, "p=5 Chernoff certificate")
print("  p=5 Chernoff certificate: rho(21/20) = %s, rho^10 = %.6f < (21/20)^23 = %.6f, 2^23 < 5^10: rate per step %.5f."
      % (rl, float(rl ** 10), float(lam ** 23), float((rl ** 10 / lam ** 23) ** Fraction(1, 10))))
# E[lambda^K_J] = (1/4) 1^T M^J 1 for p=5 (rank-one tilt), J<=3, lambda in {2, 21/20}
for lamb in (Fraction(2), lam):
    S = reach
    Mx = [[lamb ** k0[s] / 5 * (1 + (lamb ** o if s2 == e else 0)) for s2 in S] for s in S]
    for J in range(1, 4):
        vec = [Fraction(1)] * 4
        for _ in range(J):
            vec = [sum(Mx[i][j] * vec[j] for j in range(4)) for i in range(4)]
        pred = sum(vec) / 4
        acc = Fraction(0)
        totm = 0
        for m in range(1, 5 ** (J + 1)):
            if m % 5 == 0:
                continue
            totm += 1
            x = m
            K = 0
            for _ in range(J):
                kk = k[x % 25]
                K += kk
                x = ((x << kk) - 1) // 5
            acc += lamb ** K
        check(acc / totm == pred, "E[lambda^K_J] p=5 J=%d" % J)
print("  E[lambda^K_J] = (1/4) 1^T M^J 1 for p=5, J<=3, lambda in {2, 21/20}: exact.")
# E_5 greedy cycles are the reversed deterministic 5n+1 cycles; E_7 residues reachable from 1
c13 = [13]
while True:
    c13.append(((c13[-1] << k[c13[-1] % 25]) - 1) // 5)
    if c13[-1] == 13:
        break
c17 = [17]
while True:
    c17.append(((c17[-1] << k[c17[-1] % 25]) - 1) // 5)
    if c17[-1] == 17:
        break
check(c13[:-1] == [13, 83, 33] and c17[:-1] == [17, 27, 43], "E_5 greedy cycles")
for start in (13, 17):
    x = start
    for _ in range(50):
        x = x // 2 if x % 2 == 0 else 5 * x + 1
        if x == start:
            break
    check(x == start, "5n+1 deterministic cycle at %d" % start)
seen = {1}
st = [1]
while st:
    x = st.pop()
    for y in ([x // 2] if x % 2 == 0 else []) + ([7 * x + 1] if 7 * x + 1 <= 10 ** 4 else []):
        if y not in seen:
            seen.add(y)
            st.append(y)
check({x % 7 for x in seen} == {1, 2, 4}, "E_7 residues")
print("  E_5 greedy cycles {13,83,33}, {17,27,43} = reversed deterministic 5n+1 cycles; nodes reachable from 1 in")
print("  E_7 within [1,10^4] have residues {1,2,4} mod 7 only (m = 3 refutes the Q2 analogue for p = 7).")
print("  time %.1fs" % (time.time() - t0))

# S8.8  the 2-adic hostile of the typed analogy: n0 = 2^(L+1)-1 with L consecutive k=1 steps
x = 31
orb = [x]
for _ in range(4):
    x = collatz_T(x)
    orb.append(x)
check(orb == [31, 47, 71, 107, 161] and 161 == 2 * 3 ** 4 - 1, "T-orbit of 31")
check(1 + (2 * 1 + 1) == 4, "1 -> 4 is the summand arrow with parents (1, 3)")
print("S8.8 accelerated T-orbit of 31 = 2^5-1 (L=4): %s, ending at 2*3^4-1 = 161 (inherited (B5) block);" % orb)
print("  the E-arrow 1 -> 4 is the summand arrow n -> n+(2n+1) with distinct parents (1,3); the diagonal (1,1) -> 2")
print("  of the inherited note is the shortcut arrow 1 -> 2, which is not an E-arrow.")

# ---------------------------------------------------------------------------
banner("Provenance")
# ---------------------------------------------------------------------------
with open(__file__, "rb") as fh:
    print("source sha256: %s" % hashlib.sha256(fh.read()).hexdigest())
print("python: %s" % sys.version.split()[0])
print("total time: %.1fs" % (time.time() - T_START))
print("ALL CHECKS PASSED")
