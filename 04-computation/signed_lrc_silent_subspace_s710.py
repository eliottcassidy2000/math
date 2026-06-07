"""Signed-LRC: is each silent move's SILENT SET an affine F_2-subspace?  (monad-explorer-S710)

Goal: turn deficiency(C) into a NON-BRUTE count so we can reach C=63,75,81,99,105
(n=32,38,41,50,53) where 2^{n-2} brute force is infeasible (S708/S708b open frontier).

Setup (THM-417 / HYP-2273 convention, matches signed_lrc_subgroup_lattice_s708.py):
  C = 2n-1, runners (magnitudes) 1..n-1, gauge-fix sign of magnitude 1 to +1,
  free signs eps on magnitudes 2..n-1 (so n-2 free bits).
  S_eps = { eps_i * i  mod C }.  Two cuts COLLIDE (homometric) iff equal difference
  multiset over Z/C.  deficiency = 2^{n-2} - #classes.

Order-block lattice V = span_{F2}{ H_d : d|C, 1<d<C },  H_d = positive half of order-d
subgroup = { x in 1..(C-1)/2 : ord(x)=... } actually H_d := multiples of C/d in the half
system.  dim V = tau(C)-2.

A move D (subset of magnitudes 2..n-1, given as an F2 vector over the free coords) is
SILENT at eps iff eps^D ~ eps (homometric).  By the A.B lemma (S708b) <=> A_t B_t = 0 all t,
B_t = sum_{i in D} eps_i sin(2pi t i /C), A_t = Phi_t - B_t.

THIS SCRIPT (brute, small C): for each generator/element D of V, compute silent(D) exactly
and TEST whether it is an affine F_2 subspace of the free sign space; record its dimension.
If affine, deficiency becomes a linear-algebra computation reachable at any C.
"""
import sys
from itertools import product
from collections import defaultdict, Counter


def factor_str(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f)) if f else str(C)


def proper_divisors(C):
    return [d for d in range(2, C) if C % d == 0]


def order_d_halfsystem(C, d):
    """H_d = the order-d subgroup's positive half = { (C//d)*j mod C : j } cap [1,(C-1)/2]."""
    K = set(((C // d) * j) % C for j in range(d))
    return frozenset(x for x in K if 1 <= x <= (C - 1) // 2)


def diff_multiset(S, C):
    cnt = [0] * C
    for a in S:
        for b in S:
            if a != b:
                cnt[(a - b) % C] += 1
    return tuple(cnt)


def all_cuts(n):
    """Yield (free_bits_tuple, eps_array). magnitude 1 fixed +1; free signs on 2..n-1."""
    C = 2 * n - 1
    nfree = n - 2
    for bits in range(1 << nfree):
        eps = [1] * (n - 1)        # index i -> magnitude i+1
        for b in range(nfree):
            if (bits >> b) & 1:
                eps[b + 1] = -1
        yield bits, eps


def is_affine_subspace(points):
    """points: set of ints (bitmask over free coords). Affine subspace over F2 iff
    p ^ q ^ r in set for all p,q,r in set (equivalently a coset of a linear subspace).
    Return (True, dim) or (False, None)."""
    pts = list(points)
    if not pts:
        return False, None
    base = pts[0]
    # candidate linear subspace = { p ^ base }
    lin = set(p ^ base for p in pts)
    # must be closed under XOR and contain 0
    if 0 not in lin:
        return False, None
    lin_list = list(lin)
    basis = []
    for v in lin_list:
        cur = v
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    # size check
    if (1 << len(basis)) != len(lin):
        return False, None
    # closure check (cheap): every XOR of two elements in lin
    for a in lin_list[:64]:
        for b in lin_list[:64]:
            if (a ^ b) not in lin:
                return False, None
    return True, len(basis)


def analyze(n, max_brute=1 << 22, verbose=True):
    C = 2 * n - 1
    nfree = n - 2
    if (1 << nfree) > max_brute:
        if verbose:
            print(f"C={C} (n={n}): 2^{nfree} too big for brute, skipping brute analysis")
        return None
    runners = list(range(1, n))      # magnitudes; runners[k] = k+1
    mag_to_free = {runners[k]: k - 1 for k in range(1, n - 1)}  # magnitude -> free-bit index (mag>=2)

    # --- order-block lattice V ---
    divs = proper_divisors(C)
    halves = {d: order_d_halfsystem(C, d) for d in divs}
    # represent each H_d as an F2 vector over free coords (drop magnitude 1 if present-it never is)
    def half_to_mask(Hd):
        m = 0
        for x in Hd:
            if x == 1:
                # magnitude 1 is gauge-fixed; H_d should never contain it for proper d
                return None
            m |= 1 << mag_to_free[x]
        return m
    gen_masks = {}
    for d in divs:
        msk = half_to_mask(halves[d])
        if msk:
            gen_masks[d] = msk
    # F2 span of generator masks = V
    basisV = []
    for msk in gen_masks.values():
        cur = msk
        for b in basisV:
            cur = min(cur, cur ^ b)
        if cur:
            basisV.append(cur); basisV.sort(reverse=True)
    dimV = len(basisV)
    Velems = [0]
    for b in basisV:
        Velems = Velems + [x ^ b for x in Velems]
    Velems = sorted(set(Velems))

    # --- brute homometry classes ---
    sig2bits = defaultdict(list)
    bits2sig = {}
    for bits, eps in all_cuts(n):
        S = sorted((eps[i] * runners[i]) % C for i in range(n - 1))
        sig = diff_multiset(S, C)
        sig2bits[sig].append(bits)
        bits2sig[bits] = sig
    classes = list(sig2bits.values())
    sizes = Counter(len(c) for c in classes)
    defic = (1 << nfree) - len(classes)

    # --- for each move D in V, compute silent(D) and test affine ---
    silent_dim = {}
    silent_affine = {}
    for D in Velems:
        if D == 0:
            continue
        silentset = set()
        for bits in range(1 << nfree):
            if bits2sig[bits] == bits2sig[bits ^ D]:
                silentset.add(bits)
        aff, dim = is_affine_subspace(silentset)
        silent_dim[D] = (len(silentset), dim if aff else None)
        silent_affine[D] = aff

    if verbose:
        print(f"\n=== C={C}={factor_str(C)} (n={n}, free={nfree}, 2^free={1<<nfree}) ===")
        print(f"  deficiency={defic}  classes={len(classes)}  size-hist={dict(sorted(sizes.items()))}")
        print(f"  dim V = {dimV}  (tau-2 = {len(divs)})  |V|={len(Velems)}")
        all_aff = all(silent_affine.values())
        print(f"  ALL silent(D) affine F2-subspaces? {all_aff}")
        for D in sorted([d for d in Velems if d], key=lambda x: (bin(x).count('1'), x)):
            cnt, dim = silent_dim[D]
            mags = [k + 2 for k in range(nfree) if (D >> k) & 1]
            tag = ""
            for d, msk in gen_masks.items():
                if msk == D:
                    tag = f" =H_{d}"
            print(f"    D mags={mags}{tag}: |silent|={cnt} affine={silent_affine[D]} dim={dim}")
    return dict(C=C, n=n, defic=defic, sizes=dict(sizes), dimV=dimV,
                silent_dim={D: silent_dim[D] for D in silent_dim}, all_affine=all(silent_affine.values()))


if __name__ == "__main__":
    results = []
    # composite C, brute-forceable
    for n in [5, 8, 11, 13, 14, 17, 18, 20, 23]:   # C=9,15,21,25,27,33,35,39,45
        r = analyze(n)
        if r:
            results.append(r)
    print("\n\n=== SUMMARY (all silent sets affine?) ===")
    for r in results:
        print(f"  C={r['C']:3d}={factor_str(r['C']):8s} defic={r['defic']:6d} "
              f"sizes={r['sizes']} dimV={r['dimV']} all_affine={r['all_affine']}")
