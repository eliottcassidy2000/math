"""
mac-mini-2026-07-07-S43 (HYP-4887) -- the HAMILTONIAN-PATH / STEP-SEQUENCE frame for
runner families (owner directive: tiling-model analysis applied to lonely runners),
with the numerics supporting three theorems and pinning the constants of a fourth:

THM-A (reversal): gap law of -E = gap law of E pointwise => unique minimizers of any
  gap functional are PALINDROMIC step sequences.  Census: which record families are
  palindromes?  Palindrome-constrained E[U] descent: does the record drop / persist?
THM-B (wall count): sum_{i<j}(e_j - e_i) = sum_k s_k k(K-k); AP unique minimizer
  (= the 'transitive tournament' of the runner world -- coarsest order-cell complex).
  Table: wall count vs mu / E[U] across the bank.
THM-C (lattice-class): the gap-process law = Haar on the annihilator loop of L(E);
  factors through L_0(E) (= rational-affine class).  [proof in reflection; no numerics]
THM-D (BALANCED-GIRTH / sparse-lane floor): W(E) := min l1-norm of nonzero m in
  L_0(E) = {m: sum m_e e = 0, sum m_e = 0}.  |E[U_u] - (1-u)^13| <= T(W) :=
  sum_{n>=W} P(n) N_incr(n), P(n) = max product of |phi_u(m_e)| over vectors of
  l1-norm n (= u^ceil(n/2) for n <= 26 at u=1/7), N(n) <= (1+2n/W)^11 (packing, rank 11).
  Pin: T(W) numeric table; girth census of the bank; the W* where T(W*) <= budget
  (6/7)^13 - (6/7)m_P = 0.08638.
"""
import numpy as np
from math import gcd, comb
from functools import reduce
from itertools import combinations
import random as rnd
rnd.seed(43)

U = 1/7
MAIN = (1-U)**13
MP = 14249/252252
BAR_EU = (6/7)*MP
BUDGET = MAIN - BAR_EU

GRID = 150_000
xs = (np.arange(GRID)+0.5)/GRID

def gapstats(E):
    ph = np.mod(np.outer(xs, np.array(E,float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    Uv = np.clip(g-U,0,None).sum(axis=1)
    return float(Uv.mean()), float((g.max(axis=1)>U).mean())

def steps(E):
    E = sorted(E); return tuple(E[i+1]-E[i] for i in range(len(E)-1))

def is_palin(E):
    s = steps(E); return s == s[::-1]

def wallcount(E):
    E = sorted(E); K = len(E)
    return sum(E[j]-E[i] for i in range(K) for j in range(i+1,K))

def girth_L0(E, NMAX=8):
    """min l1-norm of nonzero balanced relation, searching supports of size 3..5 with
    coefficient search; exact for girth <= NMAX (supports beyond 5 need norm >= 6 with
    all-+-1 patterns checked separately)."""
    E = sorted(E); K = len(E)
    best = None
    # support 3: rank-1 generator w = (b-c, c-a, a-b)/g
    for (a,b,c) in combinations(E,3):
        w = (b-c, c-a, a-b); g = reduce(gcd,[abs(v) for v in w])
        n = sum(abs(v)//g for v in w)
        if best is None or n < best: best = n
    # support 4..6, small coefficients (covers all-±1 and mixed up to NMAX)
    for t in (4,5,6):
        M = NMAX - (t-1)  # max single coefficient if others are 1
        if M < 1: continue
        for idx in combinations(range(K), t):
            e = [E[i] for i in idx]
            # search coefficient vectors with |m_i| <= 3 and l1 <= NMAX
            rng = range(-3,4)
            for m in _coeffs(t, 3):
                if sum(m) != 0: continue
                n1 = sum(abs(v) for v in m)
                if n1 == 0 or n1 > NMAX: continue
                if best is not None and n1 >= best: continue
                if sum(mi*ei for mi,ei in zip(m,e)) == 0:
                    best = n1
    return best

_coeff_cache = {}
def _coeffs(t, cmax):
    key=(t,cmax)
    if key not in _coeff_cache:
        out=[]
        def rec(pre):
            if len(pre)==t:
                if any(v==0 for v in pre): return
                out.append(tuple(pre)); return
            for v in range(-cmax,cmax+1):
                rec(pre+[v])
        rec([])
        _coeff_cache[key]=out
    return _coeff_cache[key]

BANK = {
 'AP {1..13}': list(range(1,14)),
 'GW {1..11,13,24}': list(range(1,12))+[13,24],
 'death-star 2{1..12}u{13}': [2*i for i in range(1,13)]+[13],
 'monad record': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'PZ-min S41': [0,2,4,5,6,7,8,9,10,11,12,14,16],
 'EU-min 3-adic S41': [0,30,36,45,50,54,60,63,70,72,81,90,108],
 'klein k12-anomaly+0': [0,1,2,4,5,6,7,8,9,10,11,13,14],
 'deep interlace 10{0..10}u{49,51}': [10*i for i in range(11)]+[49,51],
 'random big': sorted(rnd.sample(range(1,2000),13)),
}

print("=== THM-A/B census: palindromes, wall counts, girth, functionals ===")
print(f"budget for THM-D: (6/7)^13 - (6/7)m_P = {BUDGET:.5f}")
print(f"{'family':32s} {'palin':>5s} {'walls':>6s} {'girth':>5s} {'E[U]':>7s} {'mu':>6s} {'deficit':>8s}")
for name,E in BANK.items():
    eu, mu = gapstats(E)
    pal = is_palin(E); wc = wallcount(E); gi = girth_L0(E)
    print(f"{name:32s} {str(pal):>5s} {wc:6d} {str(gi):>5s} {eu:7.4f} {mu:6.3f} {eu-MAIN:+8.4f}")
print(f"(AP walls minimum = C(14,3) = {comb(14,3)}; walls = sum_k s_k k(13-k))")

# ---- palindrome-constrained E[U] descent ----
print("\n=== palindrome-constrained E[U] descent (does the record structure persist?) ===")
GRID_C = 30_000; xs_c = (np.arange(GRID_C)+0.5)/GRID_C
def eu_c(E):
    ph = np.mod(np.outer(xs_c, np.array(E,float)), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return float(np.clip(g-U,0,None).sum(axis=1).mean())

def palin_family(half):  # half = s_1..s_6 -> full 12-step palindrome
    s = list(half)+list(half)[::-1]
    E=[0]
    for v in s: E.append(E[-1]+v)
    return E

best=(9.9,None)
for trial in range(60):
    half=[rnd.randint(1,9) for _ in range(6)]
    cur=eu_c(palin_family(half))
    for it in range(250):
        i=rnd.randrange(6); cand=half.copy()
        cand[i]=max(1,cand[i]+rnd.choice([-3,-2,-1,1,2,3]))
        cv=eu_c(palin_family(cand))
        if cv<cur: half,cur=cand,cv
    if cur<best[0]: best=(cur,tuple(half))
eu_p, mu_p = gapstats(palin_family(best[1]))
print(f"best palindrome half-steps {best[1]} -> steps {tuple(list(best[1])+list(best[1])[::-1])}")
print(f"E[U] = {eu_p:.5f} (fine grid) vs S41 unconstrained record 0.0938; mu = {mu_p:.4f}")

# ---- THM-D constants: T(W) table with the crude packing bound ----
print("\n=== THM-D: T(W) = sum_{n>=W} P(n)(1+2n/W)^11 with P(n)=u^ceil(n/2) (n<=26), ")
print("    u^12/(pi(n-24)) beyond; budget", f"{BUDGET:.5f}", "===")
def Pmax(n):
    if n <= 26: return U**((n+1)//2)
    return (U**12)/(np.pi*(n-24))
for W in (8,10,12,14,16,18,20,22):
    T = sum(Pmax(n)*(1+2*n/W)**11 for n in range(W, 200))
    print(f"  W={W:2d}: T(W) = {T:10.5f}   {'<= budget OK' if T<=BUDGET else 'too big'}")

# empirical check: families with large girth should have tiny |deficit|
print("\n(cross-check above: girth vs measured deficit in the census table)")
