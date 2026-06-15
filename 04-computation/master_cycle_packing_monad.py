#!/usr/bin/env python3
"""
master_cycle_packing_monad.py  (monad-explorer-2026-06-15-S5)   [pure-python, no numpy]

THE MASTER CYCLE-PACKING POLYNOMIAL Phi unifying the spectrum (char poly) and the OCF (H).

For a tournament T on n vertices with adjacency A, a *linear subdigraph* L is a set of
vertex-disjoint directed cycles (the empty set allowed). Define the master polynomial in
one variable y_k per cycle length k:

    Phi(T; {y_k}) = sum_{L linear subdigraph} prod_{C in L} y_{|C|}.

Two specializations:
  (1) SPECTRAL / Sachs Coefficient Theorem:  y_k = -x^{-k}  =>  x^n * Phi = det(x I - A).
      Equivalently  [coeff of x^{n-m} in det(xI-A)] = e_m^signed
        := sum_{L : |V(L)| = m} (-1)^{#cycles(L)}    (signed ALL-length packing count).
      These are (+-) the elementary symmetric functions of the eigenvalues => SPECTRAL.
  (2) OCF / Redei H:  y_k = 2*[k odd]  =>  Phi = sum_{odd packings} 2^{#cycles} = I(Omega,2) = H.

So H and the WHOLE characteristic polynomial are two evaluations of ONE packing polynomial.
The "non-spectral defect" (THM-505) is the gap between the all-length signed CURVE
(y_k=-x^{-k}) and the odd-only unsigned POINT (y_k=2[k odd]).

Char poly computed exactly via Faddeev-LeVerrier (integer arithmetic).
"""
import sys, itertools, random
from fractions import Fraction

# ---------------------------------------------------------------------------
def random_tournament(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.randint(0, 1):
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def matmul(A, B, n):
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for k in range(n):
            a = Ai[k]
            if a:
                Bk = B[k]
                for j in range(n):
                    Ci[j] += a*Bk[j]
    return C

def charpoly_int(A, n):
    """Faddeev-LeVerrier: coefficients of det(xI - A) = x^n + c1 x^{n-1} + ... + cn.
    Returns [1, c1, ..., cn] (length n+1); coeff of x^{n-m} is out[m]. Exact integers."""
    I = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    M = [row[:] for row in I]
    coeffs = [1]
    AM = None
    for k in range(1, n+1):
        AM = matmul(A, M, n)
        tr = sum(AM[i][i] for i in range(n))
        ck = Fraction(-tr, k)
        assert ck.denominator == 1, f"non-integer charpoly coeff {ck}"
        ck = ck.numerator
        coeffs.append(ck)
        if k < n:
            # M_{k+1} = A M_k + c_k I
            M = [[AM[i][j] + (ck if i==j else 0) for j in range(n)] for i in range(n)]
    return coeffs

# ---------------------------------------------------------------------------
def all_cycles(A, n):
    """All directed cycles (rotation classes): list of (length, frozenset(vertices)).
    Counted once: min vertex = start; all other vertices larger."""
    cycles = []
    for start in range(n):
        path = [start]; visited = {start}
        def dfs(u):
            for w in range(n):
                if w == start:
                    if len(path) >= 3 and A[u][start]:
                        cycles.append((len(path), frozenset(path)))
                elif w > start and w not in visited and A[u][w]:
                    visited.add(w); path.append(w)
                    dfs(w)
                    path.pop(); visited.discard(w)
        dfs(start)
    return cycles

def count_ham_paths(A, n):
    cnt = 0
    for perm in itertools.permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            cnt += 1
    return cnt

def enumerate_packings(cycles):
    """All linear subdigraphs as (length_multiset_tuple, num_cycles)."""
    n_c = len(cycles)
    vsets = [vs for (_, vs) in cycles]
    lens = [L for (L, _) in cycles]
    out = []
    def rec(start, used, chosen):
        out.append((tuple(sorted(chosen)), len(chosen)))
        for i in range(start, n_c):
            if not (vsets[i] & used):
                rec(i+1, used | vsets[i], chosen + [lens[i]])
    rec(0, frozenset(), [])
    return out

def packing_invariants(A, n):
    cycles = all_cycles(A, n)
    packings = enumerate_packings(cycles)
    e_signed = {}
    H_odd = H_even = H_all = 0
    sgn_odd = sgn_even = 0
    Nlam = {}
    for (lam, nc) in packings:
        m = sum(lam)
        e_signed[m] = e_signed.get(m, 0) + (-1)**nc
        Nlam[lam] = Nlam.get(lam, 0) + 1
        all_odd  = all(L % 2 == 1 for L in lam)
        all_even = all(L % 2 == 0 for L in lam)
        H_all += 2**nc
        if all_odd:
            H_odd += 2**nc; sgn_odd += (-1)**nc
        if all_even and lam != ():     # exclude vacuous empty from H_even (no even cycle)
            H_even += 2**nc; sgn_even += (-1)**nc
    return dict(e_signed=e_signed, H_odd=H_odd, H_even=H_even, H_all=H_all,
                sgn_odd=sgn_odd, sgn_even=sgn_even, Nlam=Nlam)

# ---------------------------------------------------------------------------
def verify_coefficient_theorem(n, n_samples, seed=1):
    rng = random.Random(seed); fails = checked = 0
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        inv = packing_invariants(A, n)
        cp = charpoly_int(A, n)  # coeff of x^{n-m} is cp[m]
        ok = True
        for m in range(0, n+1):
            sachs = 1 if m == 0 else inv['e_signed'].get(m, 0)
            if cp[m] != sachs:
                ok = False
                if fails < 3:
                    print(f"   [n={n}] COEFF-THM FAIL m={m}: charpoly={cp[m]} sachs={sachs}")
                break
        checked += 1; fails += (0 if ok else 1)
    print(f" Coefficient theorem (Sachs) n={n}: {checked-fails}/{checked} match")
    return fails == 0

def verify_OCF(n, n_samples, seed=2):
    rng = random.Random(seed); fails = checked = 0
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        inv = packing_invariants(A, n)
        Hp = count_ham_paths(A, n)
        if Hp != inv['H_odd']:
            fails += 1
            if fails <= 3:
                print(f"   [n={n}] OCF FAIL: hampaths={Hp} I(Om,2)={inv['H_odd']}")
        checked += 1
    print(f" OCF  H == I(Omega_odd, 2)  n={n}: {checked-fails}/{checked} match")
    return fails == 0

def verify_skeleton_sachs(n, n_samples, seed=3):
    rng = random.Random(seed); fails = checked = 0
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        inv = packing_invariants(A, n); Nl = inv['Nlam']
        c6 = Nl.get((6,), 0); D33 = Nl.get((3,3), 0)
        e6 = inv['e_signed'].get(6, 0); pred = -c6 + D33
        if e6 != pred:
            fails += 1
            if fails <= 3:
                print(f"   [n={n}] SKEL FAIL m=6: e6={e6} -c6+D33={pred} (c6={c6},D33={D33})")
        checked += 1
    print(f" Sachs skeleton  e6^signed == D33 - c6  n={n}: {checked-fails}/{checked} match")
    return fails == 0

def cospectral_face_probe(n, n_samples, seed=7):
    rng = random.Random(seed)
    buckets = {}
    invnames = ['H_odd', 'H_even', 'H_all', 'sgn_odd', 'sgn_even']
    for _ in range(n_samples):
        A = random_tournament(n, rng)
        cp = tuple(charpoly_int(A, n))
        inv = packing_invariants(A, n)
        b = buckets.setdefault(cp, {k: set() for k in invnames})
        for k in invnames:
            b[k].add(inv[k])
    any_split = sum(1 for cp,d in buckets.items() if max(len(d[k]) for k in invnames) > 1)
    print(f"\n cospectral probe n={n}: {len(buckets)} cospectral classes from {n_samples} samples "
          f"({any_split} split by some face)")
    for k in invnames:
        split = sum(1 for cp, d in buckets.items() if len(d[k]) > 1)
        verdict = "NON-SPECTRAL" if split > 0 else "spectral (constant on all classes)"
        print(f"   {k:9s}: split in {split:4d} / {len(buckets)} classes  -> {verdict}")

# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("="*78)
    print(" MASTER CYCLE-PACKING POLYNOMIAL  Phi  -- spectrum & OCF as two specializations")
    print("="*78)
    print("\n[1] Coefficient theorem  det(xI-A) coeff  ==  e_m^signed (Sachs all-length packing)")
    for n in (3,4,5,6,7):
        verify_coefficient_theorem(n, 200 if n < 7 else 120)

    print("\n[2] OCF   H(Hamiltonian paths)  ==  I(Omega_odd, 2)  (unsigned odd-only, fugacity 2)")
    for n in (3,4,5,6,7):
        verify_OCF(n, 150 if n < 7 else 80)

    print("\n[3] Sachs reading of the THM-505 defect lock  (e6^signed = D33 - c6)")
    for n in (6,7,8):
        verify_skeleton_sachs(n, 150 if n < 8 else 80)

    print("\n[4] NON-SPECTRAL BOUNDARY of each face of Phi (cospectral-class probe)")
    print("    H_odd = the OCF; H_even/H_all = even/all-cycle packing analogues; sgn_* = signed")
    for n in (4,5,6,7):
        cospectral_face_probe(n, 6000 if n <= 6 else 9000)
    print("\nDONE.")


# ---------------------------------------------------------------------------
# [5] The spectral SKELETON of H as a Z-linear combination of char-poly coeffs e_m
#     (cleaner than THM-505's power-sum W_k form; e_m ARE packing counts -> natural basis)
def verify_skeleton_in_e_basis():
    """n=7:  H = (1 - 2e3 - 2e5 + 4e6) + 4c6 + 2c7.
       n=8:  H = (1 - 2e3 - 2e5 + 4e6 + 4e8) + 4c6 + 2c7 + 4c8 - 4 D44.
       where e_m are char-poly coeffs (coeff of x^{n-m}), c_k = #k-cycles,
       D44 = # vertex-disjoint pairs of 4-cycles."""
    print("\n[5] Spectral skeleton as a Z-combination of char-poly coefficients e_m")
    # n=7
    rng = random.Random(11); f7 = c7tot = 0
    for _ in range(2000):
        A = random_tournament(7, rng); inv = packing_invariants(A, 7); Nl = inv['Nlam']
        cp = charpoly_int(A, 7)            # cp[m] = e_m^signed
        e3,e5,e6 = cp[3],cp[5],cp[6]
        c6 = Nl.get((6,),0); c7 = Nl.get((7,),0)
        H = inv['H_odd']
        skel = 1 - 2*e3 - 2*e5 + 4*e6
        if H != skel + 4*c6 + 2*c7: f7 += 1
        c7tot += 1
    print(f"   n=7  H == (1-2e3-2e5+4e6) + 4c6 + 2c7 :  {c7tot-f7}/{c7tot} match")
    # n=8
    rng = random.Random(12); f8 = c8tot = 0
    for _ in range(1500):
        A = random_tournament(8, rng); inv = packing_invariants(A, 8); Nl = inv['Nlam']
        cp = charpoly_int(A, 8)
        e3,e5,e6,e8 = cp[3],cp[5],cp[6],cp[8]
        c6 = Nl.get((6,),0); c7 = Nl.get((7,),0); c8 = Nl.get((8,),0)
        D44 = Nl.get((4,4),0)
        H = inv['H_odd']
        skel = 1 - 2*e3 - 2*e5 + 4*e6 + 4*e8
        if H != skel + 4*c6 + 2*c7 + 4*c8 - 4*D44: f8 += 1
        c8tot += 1
    print(f"   n=8  H == (1-2e3-2e5+4e6+4e8) + 4c6 + 2c7 + 4c8 - 4D44 :  {c8tot-f8}/{c8tot} match")
    # Confirm the e-basis skeleton is constant on cospectral classes (it must be: e_m spectral)
    rng = random.Random(13); buckets = {}
    for _ in range(8000):
        A = random_tournament(7, rng); cp = charpoly_int(A,7)
        skel = 1 - 2*cp[3] - 2*cp[5] + 4*cp[6]
        buckets.setdefault(tuple(cp), set()).add(skel)
    bad = sum(1 for s in buckets.values() if len(s) > 1)
    print(f"   n=7  e-basis skeleton constant on all {len(buckets)} cospectral classes: "
          f"{'YES' if bad==0 else f'NO ({bad} split)'}")

if __name__ == '__main__':
    verify_skeleton_in_e_basis()


def verify_skeleton_n9():
    """n=9 canonical (simple-cycle-carrier) Sachs-basis form:
       H = (1 - 2e3 - 2e5 + 4e6 + 4e8) + 4c6 + 2c7 + 4c8 + 2c9 - 4*D44 + 8*T333.
       Pattern conjecture:  skeleton = 1 - 2e3 - 2e5 + 4*sum_{even m>=6} e_m ;
       carriers = simple cycles (even weight 4, odd weight 2) + overlap configs."""
    import random
    rng = random.Random(21); f = tot = 0
    for _ in range(800):
        A = random_tournament(9, rng); inv = packing_invariants(A, 9); Nl = inv['Nlam']
        cp = charpoly_int(A, 9)
        e3,e5,e6,e8 = cp[3],cp[5],cp[6],cp[8]
        c6=Nl.get((6,),0); c7=Nl.get((7,),0); c8=Nl.get((8,),0); c9=Nl.get((9,),0)
        D44=Nl.get((4,4),0); T333=Nl.get((3,3,3),0)
        H = inv['H_odd']
        skel = 1 - 2*e3 - 2*e5 + 4*e6 + 4*e8
        pred = skel + 4*c6 + 2*c7 + 4*c8 + 2*c9 - 4*D44 + 8*T333
        if H != pred: f += 1
        tot += 1
    print(f"   n=9  H == (1-2e3-2e5+4e6+4e8) + 4c6+2c7+4c8+2c9 - 4D44 + 8T333 :  {tot-f}/{tot} match")


def verify_eulerchar_projection():
    """n=7: sgn_odd = I(Omega,-1) = -reduced-Euler-char of odd-cycle packing complex.
       Claim: sgn_odd = (1 + e3 + e5 + e6) + (c6 - c7)  -> non-spectral content is the
       1-D combination c6-c7, a projection of H's 2-D content (c6,c7).  And sgn_odd splits
       a cospectral class IFF (c6-c7) varies, while H splits IFF (2c6+c7) varies."""
    import random
    rng = random.Random(31); f = tot = 0
    for _ in range(3000):
        A = random_tournament(7, rng); inv = packing_invariants(A, 7); Nl = inv['Nlam']
        cp = charpoly_int(A, 7); e3,e5,e6 = cp[3],cp[5],cp[6]
        c6 = Nl.get((6,),0); c7 = Nl.get((7,),0)
        pred = (1 + e3 + e5 + e6) + (c6 - c7)
        if inv['sgn_odd'] != pred: f += 1
        tot += 1
    print(f"   n=7  sgn_odd == (1+e3+e5+e6) + (c6-c7) :  {tot-f}/{tot} match")
    # within cospectral classes: does sgn_odd split <=> Delta(c6-c7)!=0, H split <=> Delta(2c6+c7)!=0 ?
    rng = random.Random(32); buckets = {}
    for _ in range(12000):
        A = random_tournament(7, rng); inv = packing_invariants(A,7); Nl = inv['Nlam']
        cp = tuple(charpoly_int(A,7))
        c6 = Nl.get((6,),0); c7 = Nl.get((7,),0)
        b = buckets.setdefault(cp, {'sgn':set(),'H':set(),'diff':set(),'comb':set(),'pair':set()})
        b['sgn'].add(inv['sgn_odd']); b['H'].add(inv['H_odd'])
        b['diff'].add(c6-c7); b['comb'].add(2*c6+c7); b['pair'].add((c6,c7))
    sgn_split = sum(1 for d in buckets.values() if len(d['sgn'])>1)
    H_split   = sum(1 for d in buckets.values() if len(d['H'])>1)
    # consistency: sgn splits <=> diff varies ; H splits <=> comb varies
    bad_sgn = sum(1 for d in buckets.values() if (len(d['sgn'])>1) != (len(d['diff'])>1))
    bad_H   = sum(1 for d in buckets.values() if (len(d['H'])>1)   != (len(d['comb'])>1))
    # classes where the (c6,c7) pair varies but sgn does NOT split (covariation cancelled by x=-1)
    cov_cancel = sum(1 for d in buckets.values() if len(d['pair'])>1 and len(d['sgn'])==1)
    print(f"   n=7  cospectral classes: {len(buckets)};  sgn_odd splits {sgn_split}, H splits {H_split}")
    print(f"        sgn split <=> Delta(c6-c7)!=0 : {'CONSISTENT' if bad_sgn==0 else f'{bad_sgn} mismatch'}")
    print(f"        H   split <=> Delta(2c6+c7)!=0: {'CONSISTENT' if bad_H==0 else f'{bad_H} mismatch'}")
    print(f"        classes where (c6,c7) varies but sgn_odd does NOT (x=-1 cancels covariation): {cov_cancel}")
