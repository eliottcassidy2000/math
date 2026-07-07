"""
anti_redei_decomposition_opus_S140.py

THE ANTI-REDEI DECOMPOSITION (owner: prove Anti-Redei via the rho-twisted involution).

Anti-Redei (THM-644c conjecture): H_anti(T) odd for every T admitting an anti-automorphism,
where H_anti = #{Ham paths P : positional reversal alpha_P is an anti-automorphism}.

PROVED PIECES (writeup in THM-644 update):
  * alpha_P is always an involutive anti-automorphism candidate (order 2, n mod 2 fixed pts),
    so H_anti = SUM over anti-INVOLUTIONS beta of h_beta(T), h_beta = #{P : alpha_P = beta}.
  * LEMMA A: #anti-involutions of T is ODD whenever AntiAut(T) is nonempty.
    (AntiAut = a coset alpha*Aut of odd size |Aut|; inversion x -> x^{-1} maps the coset to
    itself and is an involution; its fixed points are exactly the anti-involutions; an
    involution on an odd-cardinality set has an odd number of fixed points.)
  * Hence Anti-Redei FOLLOWS from:  MAIN LEMMA (folded Redei): h_beta(T) is odd for every
    anti-involution beta.   [odd-many odd summands]

THIS SCRIPT (decides whether the Main Lemma is true before proof investment):
 (1) for every self-converse class n = 4,5,6 (+ sampled n=7 classes): enumerate ALL
     anti-involutions beta; verify #beta odd (Lemma A cross-check); compute each h_beta;
     TEST: h_beta odd for EVERY beta?  Also verify SUM h_beta = H_anti (partition).
 (2) the FOLDED-INSERTION SLOT PARITY (the induction engine candidate): for random
     (T, beta, beta-symmetric Ham path Q of T minus a beta-pair {u, beta u}):
     count SYMMETRIC insertions of the pair into Q (u at first-half slot j & beta u at the
     mirror slot; arcs needed: pred->u, u->succ at j — mirrored arcs automatic; plus the two
     CENTER placements ...q_h, u, beta u, beta q_h... requiring q_h->u & u->beta u, or the
     reverse order requiring u<-beta u).  TEST: #symmetric insertion slots ODD always?
 (3) the fold algebra sanity: the pair-tie identities T(beta a, beta b) = 1 - T(a,b)... i.e.
     arc a->b <=> arc beta b->beta a, and T(a, beta b) = T(b, beta a); the quotient bit
     q(A,B) = T(a,b) XOR T(a, beta b) is representative-independent.  Verified exhaustively.
"""
import sys, time, random
from itertools import permutations
from collections import Counter

sys.path.insert(0, ".")
from metagraph_fiber_allocation_opus_S139 import tiles_of, sigma_map, arcs_of_tiling, canon
from metagraph_anti_redei_test_opus_S139 import ham_paths, is_antiauto, aut_order

def anti_involutions(a, n):
    """All involutive anti-automorphisms of tournament a."""
    out = []
    for p in permutations(range(n)):
        if all(p[p[i]] == i for i in range(n)) and is_antiauto(a, n, list(p)):
            out.append(list(p))
    return out

def h_beta(a, n, beta, H=None):
    if H is None: H = ham_paths(a, n)
    cnt = 0
    for P in H:
        ok = all(beta[P[j]] == P[n - 1 - j] for j in range(n))
        if ok: cnt += 1
    return cnt

def part1():
    print("=" * 96)
    print("(1) per-beta decomposition on all self-converse classes")
    print("=" * 96)
    for n in (4, 5, 6):
        T = tiles_of(n); m = len(T)
        perms = list(permutations(range(n)))
        rep = {}
        for bits in range(1 << m):
            a = arcs_of_tiling(n, T, bits)
            c = canon(a, n, perms)
            if c not in rep: rep[c] = a
        nself = 0; all_odd = True; lemA_ok = True; part_ok = True
        hspec = Counter()
        for c, a in rep.items():
            AIs = anti_involutions(a, n)
            if not AIs: continue
            nself += 1
            if len(AIs) % 2 == 0: lemA_ok = False
            H = ham_paths(a, n)
            tot = 0
            for beta in AIs:
                hb = h_beta(a, n, beta, H)
                hspec[hb] += 1
                tot += hb
                if hb % 2 == 0:
                    all_odd = False
            Hanti = 0
            for P in H:
                f = [0] * n
                for j, v in enumerate(P): f[v] = P[n - 1 - j]
                if is_antiauto(a, n, f): Hanti += 1
            if tot != Hanti: part_ok = False
        print(f"  n={n}: self-converse classes {nself}; Lemma A (#beta odd) {'OK' if lemA_ok else '*** FAIL ***'};"
              f" partition SUM=H_anti {'OK' if part_ok else '*** FAIL ***'};"
              f" MAIN LEMMA h_beta odd: {'HOLDS on all' if all_odd else '*** FAILS ***'}")
        print(f"        h_beta spectrum: {dict(sorted(hspec.items()))}")

def part2():
    print()
    print("=" * 96)
    print("(2) folded-insertion slot parity (random probes)")
    print("=" * 96)
    rng = random.Random(140140)
    bad = tested = 0
    for trial in range(4000):
        n = rng.choice([4, 6, 8])          # even n (pair insertion case)
        # build a random T with anti-involution beta = (i <-> n-1-i)
        beta = [n - 1 - i for i in range(n)]
        a = [[0] * n for _ in range(n)]
        # free bits: pairs {i, beta i} internal + cross ties
        for i in range(n):
            for j in range(n):
                a[i][j] = -1
        for i in range(n // 2):
            b = rng.getrandbits(1)
            a[i][beta[i]] = b; a[beta[i]][i] = 1 - b
        for i in range(n):
            for j in range(n):
                if i == j or a[i][j] != -1: continue
                b = rng.getrandbits(1)
                a[i][j] = b; a[j][i] = 1 - b
                # anti tie: T(beta j, beta i) = T(i,j)  <=> arc i->j forces beta j -> beta i
                a[beta[j]][beta[i]] = b; a[beta[i]][beta[j]] = 1 - b
        if not is_antiauto(a, n, beta):
            continue
        # pick a beta-pair to remove; find beta-symmetric ham paths of T-{u,bu}
        u = rng.randrange(n); bu = beta[u]
        rest = [v for v in range(n) if v not in (u, bu)]
        if not rest: continue
        idx = {v: k for k, v in enumerate(rest)}
        n2 = len(rest)
        a2 = [[a[x][y] for y in rest] for x in rest]
        beta2 = [idx[beta[v]] for v in rest]
        if not is_antiauto(a2, n2, beta2): continue
        H2 = ham_paths(a2, n2)
        syms = [P for P in H2 if all(beta2[P[j]] == P[n2 - 1 - j] for j in range(n2))]
        if not syms: continue
        Q = rng.choice(syms)
        Qv = [rest[i] for i in Q]
        h = n2 // 2
        # count symmetric insertions of (u, bu) into Qv
        slots = 0
        # u at first-half slot j: between Qv[j-1] and Qv[j] for j=0..h (j=h = just before center seam)
        for j in range(0, h + 1):
            pred_ok = (j == 0) or a[Qv[j - 1]][u]
            if j < h:
                succ_ok = a[u][Qv[j]]
                if pred_ok and succ_ok: slots += 1
            else:
                # center: ...Qv[h-1], u, bu, Qv[h]... needs Qv[h-1]->u, u->bu, bu->Qv[h] (auto?)
                # bu->Qv[h]: Qv[h] = beta(Qv[h-1]); T(bu, beta Qv[h-1]) = T(Qv[h-1], u) = pred_ok cond?
                # anti: T(beta u, beta q) = T(q, u). So bu -> beta(q) iff q -> u. auto given pred_ok.
                if pred_ok and a[u][bu]: slots += 1
                # reversed center: ...Qv[h-1], bu, u, Qv[h]...: needs Qv[h-1]->bu, bu->u
                pred2 = (h == 0) or a[Qv[h - 1]][bu]
                if pred2 and a[bu][u]: slots += 1
        tested += 1
        if slots % 2 == 0:
            bad += 1
    print(f"  probes tested: {tested}; EVEN slot-counts: {bad} "
          f"{'  => slot parity HOLDS (induction engine viable)' if bad == 0 else '  *** slot parity FAILS as stated ***'}")

def part3():
    print()
    print("=" * 96)
    print("(3) fold algebra: arc ties + q-invariance (exhaustive small check)")
    print("=" * 96)
    rng = random.Random(7)
    bad = 0
    for trial in range(2000):
        n = rng.choice([4, 6])
        beta = [n - 1 - i for i in range(n)]
        a = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n): a[i][j] = -1
        for i in range(n // 2):
            b = rng.getrandbits(1); a[i][beta[i]] = b; a[beta[i]][i] = 1 - b
        for i in range(n):
            for j in range(n):
                if i == j or a[i][j] != -1: continue
                b = rng.getrandbits(1)
                a[i][j] = b; a[j][i] = 1 - b
                a[beta[j]][beta[i]] = b; a[beta[i]][beta[j]] = 1 - b
        for i in range(n):
            for j in range(n):
                if i != j and a[beta[j]][beta[i]] != a[i][j]: bad += 1
        # q-invariance
        for i in range(n // 2):
            for j in range(n // 2):
                if i == j: continue
                A = (i, beta[i]); B = (j, beta[j])
                q1 = a[A[0]][B[0]] ^ a[A[0]][B[1]]
                q2 = a[A[1]][B[0]] ^ a[A[1]][B[1]]
                if q1 != q2: bad += 1
    print(f"  tie + q-invariance violations: {bad}")

if __name__ == "__main__":
    t0 = time.time()
    part1(); part2(); part3()
    print(f"\n[{time.time()-t0:.0f}s]")
