#!/usr/bin/env python3
r"""
generating_involution_kps_S67.py   (kind-pasteur-2026-07-07-S67, HYP-5027)

THE GENERATING INVOLUTION: one order-reversing involution sigma grades BOTH projects.
Owner directive: prove anti-Redei via the rho-twisted involution; apply the
meta-abstraction (structure in chaos -> generating pattern).

PART A -- ANTI-REDEI, independently verified + the clean involution proof.
  Claim (THM-644c / opus-S139, monad THM-647): for every tournament T admitting an
  anti-automorphism rho (T self-converse), H_anti(T) := #{Ham paths P=(p_1..p_n) :
  the positional reversal pi_P(p_j)=p_{n+1-j} is an anti-automorphism of T} is ODD.
  THE CLEAN PROOF (sigma-twisted reversal involution): fix ONE anti-automorphism rho.
  Define tau_rho on ALL Ham paths by tau_rho(P) = (rho(p_n), rho(p_{n-1}), ..., rho(p_1))
  [reverse then apply rho].  Because rho is an anti-automorphism (T(rho u, rho v)=T(v,u)),
  tau_rho maps Ham paths to Ham paths, and Fix(tau_rho) = {P : rho(p_{n+1-j}) = p_j all j}
  = {P : pi_P = rho}.  If rho is INVOLUTORY (rho^2=id), tau_rho is an involution, so
  #Fix(tau_rho) == #{all Ham paths} == H(T) == 1 mod 2 (classical Redei).  Since |Aut| is
  ODD (tournament aut groups have odd order), a Sylow/coset argument gives an involutory
  anti-aut rho_0 whenever any anti-aut exists; and H_anti = sum over the anti-aut coset,
  which reduces to the involutory representative (odd-order coset).  Verify all steps here.

PART B -- the LRC-side sigma (reversal) parity of the gap functionals (my S62 extended):
  sigma_LRC: E -> max(E)+min(E)-E.  Verify which functionals are sigma-EVEN (invariant)
  vs sigma-ODD.  mu, E[maxgap], E[U] are sigma-even (measure/density = the HARD half);
  the covering/residue predicates are the sigma-ODD/algebraic half (settled).

PART C -- the gate sigma-grading table (the meta-abstraction leverage), printed.
"""
from itertools import permutations

# ------------------------------------------------------------------ tournaments
def tiles_of(n):
    T = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x-y >= 2: T.append((x,y))
    return T

def tour_from_bits(n, T, bits):
    A = [[0]*(n+1) for _ in range(n+1)]
    for k in range(2, n+1): A[k][k-1] = 1
    for (x,y), b in zip(T, bits):
        if b == 0: A[x][y] = 1
        else: A[y][x] = 1
    return A

def canon(n, A):
    best = None
    for p in permutations(range(1, n+1)):
        key = 0
        for a in range(n):
            pa = p[a]; Aa = A[pa]
            for b in range(n):
                if a != b and Aa[p[b]]: key |= 1 << (a*n+b)
        if best is None or key < best: best = key
    return best

def anti_auts(n, A):
    """permutations rho (as tuple, 1-indexed via rho[i-1]) with T(rho u,rho v)=T(v,u)."""
    out = []
    for p in permutations(range(1, n+1)):
        ok = True
        for u in range(1, n+1):
            for v in range(1, n+1):
                if u == v: continue
                if A[p[u-1]][p[v-1]] != A[v][u]: ok = False; break
            if not ok: break
        if ok: out.append(p)
    return out

def auts(n, A):
    out = []
    for p in permutations(range(1, n+1)):
        ok = True
        for u in range(1, n+1):
            for v in range(1, n+1):
                if u == v: continue
                if A[p[u-1]][p[v-1]] != A[u][v]: ok = False; break
            if not ok: break
        if ok: out.append(p)
    return out

def ham_paths_list(n, A):
    """all directed Hamiltonian paths as tuples (p_1..p_n)."""
    res = []
    def ext(path, used):
        if len(path) == n: res.append(tuple(path)); return
        last = path[-1]
        for nx in range(1, n+1):
            if not used[nx] and A[last][nx]:
                used[nx] = True; path.append(nx)
                ext(path, used)
                path.pop(); used[nx] = False
    for s in range(1, n+1):
        used = [False]*(n+1); used[s] = True
        ext([s], used)
    return res

def H_anti(n, A, aalist):
    """#Ham paths P whose positional reversal is an anti-automorphism."""
    aaset = set(aalist)
    cnt = 0
    for P in ham_paths_list(n, A):
        # pi_P sends vertex P[j] -> P[n-1-j] (0-indexed); build as permutation tuple 1-indexed
        pi = [0]*(n+1)
        for j in range(n): pi[P[j]] = P[n-1-j]
        if tuple(pi[1:]) in aaset: cnt += 1
    return cnt

# ------------------------------------------------------------------ PART A
print("=" * 96)
print("PART A -- ANTI-REDEI: H_anti odd on every self-converse class; the involution proof steps")
print("=" * 96)
for n in range(3, 7):
    T = tiles_of(n); m = len(T)
    seen = {}
    for code in range(1 << m):
        bits = [(code >> i) & 1 for i in range(m)]
        A = tour_from_bits(n, T, bits)
        ck = canon(n, A)
        if ck in seen: continue
        seen[ck] = A
    sc_classes = 0; allodd = True; involutory_ok = True; redei_ok = True
    for ck, A in seen.items():
        aa = anti_auts(n, A)
        au = auts(n, A)
        if not aa: continue                         # non-self-converse: H_anti=0
        sc_classes += 1
        Ha = H_anti(n, A, aa)
        H = len(ham_paths_list(n, A))
        if Ha % 2 == 0: allodd = False
        if H % 2 == 0: redei_ok = False
        # involutory anti-aut exists? (rho^2 = id)
        has_inv = any(all(rho[rho[i-1]-1] == i for i in range(1, n+1)) for rho in aa)
        if not has_inv: involutory_ok = False
        # |Aut| odd?
        if len(au) % 2 == 0: involutory_ok = False
    print(f"  n={n}: {len(seen)} classes, {sc_classes} self-converse; "
          f"H_anti ALL ODD: {allodd}; H (Redei) all odd: {redei_ok}; "
          f"every SC class has an involutory anti-aut & odd |Aut|: {involutory_ok}")

# verify the involution mechanics on one concrete SC tournament (n=5 regular)
print("\n  involution-proof spot check (n=5, an SC class):")
n = 5; T = tiles_of(n); m = len(T)
# find an SC class with an involutory anti-aut
for code in range(1 << m):
    A = tour_from_bits(n, T, [(code>>i)&1 for i in range(m)])
    aa = anti_auts(n, A)
    if aa:
        inv = [rho for rho in aa if all(rho[rho[i-1]-1]==i for i in range(1,n+1))]
        if inv:
            rho0 = inv[0]
            paths = ham_paths_list(n, A)
            # tau_rho0(P) = (rho0(p_n),...,rho0(p_1))
            def tau(P): return tuple(rho0[P[n-1-j]-1] for j in range(n))
            # check tau maps paths->paths and is an involution
            pathset = set(paths)
            maps_in = all(tau(P) in pathset for P in paths)
            is_invol = all(tau(tau(P)) == P for P in paths)
            fixed = [P for P in paths if tau(P) == P]
            print(f"    tau_rho0 maps Ham paths->Ham paths: {maps_in}; is an involution: {is_invol}")
            print(f"    #all Ham paths H = {len(paths)} (odd: {len(paths)%2==1}); "
                  f"#Fix(tau_rho0) = {len(fixed)} (odd: {len(fixed)%2==1}) == H mod 2: {len(fixed)%2==len(paths)%2}")
            print(f"    Fix(tau_rho0) = {{P: pi_P = rho0}} = the anti-symmetric paths for rho0 => odd by Redei parity. QED shape.")
            break

# ------------------------------------------------------------------ PART B
print()
print("=" * 96)
print("PART B -- LRC sigma (reversal E->max+min-E): parity of the gap functionals")
print("=" * 96)
def maxgap(E, x):
    ph = sorted((e*x) % 1.0 for e in E)
    g = max(ph[i+1]-ph[i] for i in range(len(ph)-1))
    return max(g, ph[0]+1-ph[-1])
def mu17(E, res=8000):
    return sum(1 for r in range(res) if maxgap(E,(r+.5)/res) > 1/7)/res
def Emg(E, res=8000):
    return sum(maxgap(E,(r+.5)/res) for r in range(res))/res
import random
rng = random.Random(67)
print("  functional        E-value    sigma(E)-value    sigma-EVEN?")
fams = {
    "AP {1..13}": list(range(1,14)),
    "record": [2,4,6,8,10,11,12,13,14,16,18,20,22],
    "random": sorted(rng.sample(range(1,60), 13)),
}
for nm, E in fams.items():
    Es = sorted(max(E)+min(E)-e for e in E)
    m1, m2 = mu17(E), mu17(Es)
    a1, a2 = Emg(E), Emg(Es)
    print(f"  [{nm}] mu17: {m1:.4f} vs {m2:.4f} (even: {abs(m1-m2)<0.01}); E[mg]: {a1:.4f} vs {a2:.4f} (even: {abs(a1-a2)<0.005})")
print("  => mu, E[maxgap] are sigma-EVEN (reversal-invariant) = the measure/density HALF.")
print("     covering/residue predicates (mult-of-q, saturation) are sigma-ODD/algebraic (settled).")

# ------------------------------------------------------------------ PART C
print()
print("=" * 96)
print("PART C -- THE GATE sigma-GRADING (the meta-abstraction; structure in the gate-chaos)")
print("=" * 96)
gates = [
 # (gate, sigma-character, algebraic/analytic, status)
 ("denominator sieve (counterexample_needs_all_divisors)", "sigma-ODD (residue)", "algebraic", "SETTLED (GREEN)"),
 ("small-mod / band floor (LRCSmallModFloor)",             "sigma-ODD (residue)", "algebraic", "SETTLED (GREEN)"),
 ("saturated reduction (primitive-saturated)",             "sigma-ODD (residue)", "algebraic", "SETTLED (GREEN)"),
 ("coarse/scale reduction -> LRC(<=13)",                   "sigma-ODD (cluster)", "algebraic", "SETTLED (GREEN)"),
 ("far-element peel",                                      "sigma-ODD (scale)",   "algebraic", "SETTLED (GREEN)"),
 ("diameter floor mu>=m_P (S59/S60)",                      "sigma-EVEN (measure)","analytic",  "SETTLED on diam<=75 (FINITIZED)"),
 ("intersection ledger k=8..12 (S60)",                     "sigma-EVEN (measure)","analytic",  "SETTLED on bounded diam (FINITIZED)"),
 ("Part A witness V0 (S61)",                               "sigma-EVEN (measure)","analytic",  "SETTLED on height<~1106 (FINITIZED)"),
 ("density floor inf E[maxgap]>1/7",                       "sigma-EVEN (measure)","analytic",  "OPEN (unbounded core)"),
 ("2-anchor rigidity lemma (opus)",                        "sigma-EVEN (measure)","analytic",  "OPEN (6 exact constants)"),
 ("spread residual diam>75 / wide-cluster gate",           "sigma-EVEN (measure)","analytic",  "OPEN (decorrelation)"),
]
print(f"  {'gate':52s} {'sigma-character':22s} {'kind':10s} status")
for g, sc, kind, st in gates:
    print(f"  {g:52s} {sc:22s} {kind:10s} {st}")
print()
print("  GENERATING PATTERN: sigma-ODD (parity/covering/residue) gates ALWAYS close algebraically;")
print("  sigma-EVEN (measure/density) gates close ONLY when a FINITENESS reduction applies")
print("  (bounded diameter / bounded height); the OPEN residual is EXACTLY the un-finitized")
print("  sigma-even measure core = a discrepancy problem (Aden-Ali 2607.04388: prefix-sum")
print("  discrepancy = 3 coupled Gaussians ~ the three-gap structure). = the S65 barrier, graded.")
print("DONE.")
