#!/usr/bin/env python3
"""
legendre333_mim.py — search for a Legendre pair of length 333 where the two
sequences are invariant under DIFFERENT multiplier subgroups H_A, H_B of
(Z/333)^*, with small common intersection.

Background: a Legendre pair (A,B) of odd length L is a pair of {+1,-1}
sequences with PAF_A(s) + PAF_B(s) = -2 for all s != 0 (cyclic autocorrelation).
It yields a Hadamard matrix of order 2L+2; L=333 gives 668, the smallest open
Hadamard order.  arXiv:2607.20765 (July 2026) excludes pairs fixed by a COMMON
subgroup H with |H| >= 7, leaving |H| <= 6 and *asymmetric* structure open.

This program enumerates, for each pair of subgroups (H_A, H_B) of (Z/333)^*
with |H_A|,|H_B| >= LO and |H_A ∩ H_B| <= 6, all H_A-invariant candidate A and
H_B-invariant candidate B, computes exact integer PAFs (via orbit-collapsed
correlation), and hash-joins on the vector PAF_A restricted to full resolution:
match iff PAF_A(s) = -2 - PAF_B(s) for all s in 1..332.

Feasibility count: number of orbits of H on Z_333 is small for |H| large
(e.g. |H|=36 -> ~13 orbits -> 2^12 sign patterns per side).

Exact arithmetic throughout (numpy int64); any hit is re-verified from scratch.
"""
import numpy as np
from itertools import product
import sys, time

L = int(sys.argv[1]) if len(sys.argv) > 1 else 333
ALLOW_COMMON = "--allow-common" in sys.argv   # for positive controls (e.g. L=31)
LO = 7

def units(L):
    return [u for u in range(1, L) if np.gcd(u, L) == 1]

US = units(L)

def subgroup_generated(gens):
    S = {1}
    frontier = [1]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = (x * g) % L
            if y not in S:
                S.add(y); frontier.append(y)
    return frozenset(S)

def all_subgroups():
    """All subgroups of (Z/333)^* (order 216 = Z_6 x Z_36 as abstract group).
    Enumerate by closing over all subsets of a generating set is too big;
    instead: collect all cyclic subgroups, then close under joins until stable."""
    subs = set()
    for u in US:
        subs.add(subgroup_generated([u]))
    changed = True
    while changed:
        changed = False
        cur = list(subs)
        for i in range(len(cur)):
            for j in range(i + 1, len(cur)):
                J = subgroup_generated(list(cur[i]) + list(cur[j]))
                if J not in subs:
                    subs.add(J); changed = True
    return sorted(subs, key=len)

def orbits_of(H):
    """Orbits of H acting by multiplication on Z_333."""
    seen = np.zeros(L, dtype=bool)
    orbs = []
    for x in range(L):
        if not seen[x]:
            o = sorted({(h * x) % L for h in H})
            for y in o:
                seen[y] = True
            orbs.append(o)
    return orbs

def paf_all(seq):
    """Exact integer PAF vector for s=0..L-1 via FFT-free O(L^2/…) -> use numpy roll trick."""
    a = np.asarray(seq, dtype=np.int64)
    # correlation via FFT would be float; L=333 small: do exact with stride tricks
    out = np.empty(L, dtype=np.int64)
    for s in range(L):
        out[s] = int(np.dot(a, np.roll(a, -s)))
    return out

def paf_all_fast(seq):
    """Exact PAF via complex FFT with rounding check (L=333, entries ±1: safe)."""
    a = np.asarray(seq, dtype=np.float64)
    fa = np.fft.rfft(a, n=L)
    c = np.fft.irfft(fa * np.conj(fa), n=L)
    r = np.rint(c).astype(np.int64)
    return r

def enumerate_side(H, tag, cap_bits=24):
    """All H-invariant ±1 sequences on Z_333 with row sum = ±1.
    Returns list of (paf_vector_bytes, assignment) — but we only keep the dict
    from paf-key -> one representative assignment.  Row sum must be ±1 because
    sum condition: a^2 + b^2 = 2 with a=sum(A), b=sum(B)."""
    orbs = orbits_of(H)
    k = len(orbs)
    if k - 1 > cap_bits:
        return None, orbs
    sizes = np.array([len(o) for o in orbs], dtype=np.int64)
    # incidence matrix for building sequences quickly
    M = np.zeros((k, L), dtype=np.int64)
    for i, o in enumerate(orbs):
        M[i, o] = 1
    reps = {}
    t0 = time.time()
    total = 1 << k
    for bits in range(total):
        signs = np.array([1 if (bits >> i) & 1 else -1 for i in range(k)], dtype=np.int64)
        rs = int(np.dot(signs, sizes))
        if rs != 1 and rs != -1:
            continue
        seq = np.dot(signs, M)  # ±1 sequence
        paf = paf_all_fast(seq)
        # spot-verify one exact entry to guard FFT
        s_probe = 7
        exact = int(np.dot(seq, np.roll(seq, -s_probe)))
        assert exact == paf[s_probe], "FFT PAF mismatch"
        key = paf[1:].tobytes()  # s=1..332
        reps.setdefault(key, bits)
    return reps, orbs

def bits_to_seq(bits, orbs):
    seq = np.zeros(L, dtype=np.int64)
    for i, o in enumerate(orbs):
        v = 1 if (bits >> i) & 1 else -1
        for x in o:
            seq[x] = v
    return seq

def main():
    print(f"enumerating subgroups of (Z/{L})^* ...")
    subs = all_subgroups()
    print(f"  found {len(subs)} subgroups; orders: {sorted(set(len(s) for s in subs))}")
    big = [H for H in subs if len(H) >= LO]
    # precompute orbit counts
    info = []
    for H in big:
        orbs = orbits_of(H)
        info.append((H, len(orbs)))
    info.sort(key=lambda t: t[1])
    print(f"  {len(big)} subgroups with |H|>={LO}:")
    for H, k in info:
        print(f"    |H|={len(H):3d} orbits={k:3d} enumerable={'YES' if k<=25 else 'no'}")
    feas = [(H, k) for (H, k) in info if k <= 25]
    print(f"  {len(feas)} enumerable sides")

    # enumerate each side once, cache
    sides = {}
    for H, k in feas:
        t0 = time.time()
        reps, orbs = enumerate_side(H, tag=str(len(H)))
        sides[H] = (reps, orbs)
        print(f"  side |H|={len(H)} orbits={k}: {len(reps)} distinct PAF classes "
              f"({time.time()-t0:.1f}s)")
        sys.stdout.flush()

    # now try all pairs with small intersection
    target_hits = []
    hs = list(sides.keys())
    for i in range(len(hs)):
        for j in range(len(hs)):
            HA, HB = hs[i], hs[j]
            inter = len(HA & HB)
            if inter > 6 and not ALLOW_COMMON:
                continue  # excluded by common-multiplier obstruction... but note:
                          # obstruction is about the pair being fixed by the COMMON
                          # subgroup; HA∩HB fixes both, so |HA∩HB|<=6 required.
            repsA, orbsA = sides[HA]
            repsB, orbsB = sides[HB]
            # join: need PAF_A(s) = -2 - PAF_B(s) for all s>=1
            cnt = 0
            for keyB, bitsB in repsB.items():
                pafB = np.frombuffer(keyB, dtype=np.int64)
                want = (-2 - pafB).tobytes()
                if want in repsA:
                    bitsA = repsA[want]
                    A = bits_to_seq(bitsA, orbsA)
                    B = bits_to_seq(bitsB, orbsB)
                    # exact re-verification
                    pa = paf_all(A); pb = paf_all(B)
                    ok = all(pa[s] + pb[s] == -2 for s in range(1, L))
                    print(f"HIT?! |HA|={len(HA)} |HB|={len(HB)} inter={inter} verified={ok}")
                    if ok:
                        target_hits.append((A.tolist(), B.tolist()))
                        np.save("legendre333_A.npy", A)
                        np.save("legendre333_B.npy", B)
                cnt += 1
            print(f"pair |HA|={len(HA)},|HB|={len(HB)},∩={inter}: joined {cnt} classes, "
                  f"hits so far {len(target_hits)}")
            sys.stdout.flush()
    print(f"DONE. total verified Legendre pairs found: {len(target_hits)}")
    if not target_hits:
        print("NEGATIVE: no Legendre pair of length 333 exists with A,B invariant "
              "under (possibly different) subgroups of (Z/333)^* of order >= 7 "
              "with |H_A ∩ H_B| <= 6, among enumerable orbit counts <= 25.")

if __name__ == "__main__":
    main()
