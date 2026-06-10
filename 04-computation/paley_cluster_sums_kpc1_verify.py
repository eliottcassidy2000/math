#!/usr/bin/env python3
"""
paley_cluster_sums_kpc1_verify.py
ADVERSARIAL VERIFIER, kind-pasteur-2026-06-10-S1, thread E-R31-prediction.

Independent re-computation of the exact cluster/joint character integrals
claimed in E1 (and the zeros claimed in E1/E2), by a DIFFERENT method than the
worker's Moebius-over-set-partitions + einsum engine:

  METHOD: direct DFS enumeration over injective tuples (C helper), with
  translation reduction (fix x_0 = 0, factor p), scaling reduction for
  even-edge structures (fix x_1 = 1, factor p-1; the affine group acts freely
  on injective tuples and preserves the weight when #edges is even), and a
  complementary-sum trick at the last level (sum over unused c of chi(c-prev)
  = - sum over used points, since the full sum vanishes).

  VALIDATION of the C engine: pure-Python brute force with NO reduction at
  p=7 (all injective tuples) and translation-ONLY reduction at p=11 -- so the
  scaling reduction and the last-level trick are themselves checked.

Definitions (matching the worker's stated conventions exactly):
  A_L   : single run, L+1 distinct points, L edges
  J22   : runs (2,2)   6 distinct points
  J23   : runs (2,3)   7 distinct points  (claimed 0: odd run)
  J24   : runs (2,4)   8 distinct points
  J33   : runs (3,3)   8 distinct points  (claimed 0: odd runs)
  J222  : runs (2,2,2) 9 distinct points
  d22 = J22 - A2^2 ; d24 = J24 - A2*A4 ; d33 = J33 ; d3 = J222 - A2^3 - 3*d22*A2

Claimed (E1): A6(31)=3281040, A8(19)=10243584, A8(23)=34327040,
  A8(31)=202487040, d22(31)=-292020, d24(31)=-25935840, d3(31)=+241614000.
Claimed zeros: A5 = J23 = 0 at all p (negation), J33 = 0 at all p (E2 lemma).
Monad-stored anchors: A4{7:336,11:1760,19:10944,23:20240,31:52080},
  A6{7:1008,11:22880,19:361152,23:870320}, A8{7:0,11:154880}, A2=p(p-1).
"""
import json, os, subprocess, sys, time
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
T0 = time.time()

STRUCTS = {
    # name: prev array (-1 = run start), edges count implied
    'A2':   [-1, 0, 1],
    'A4':   [-1, 0, 1, 2, 3],
    'A5':   [-1, 0, 1, 2, 3, 4],
    'A6':   [-1, 0, 1, 2, 3, 4, 5],
    'A8':   [-1, 0, 1, 2, 3, 4, 5, 6, 7],
    'J22':  [-1, 0, 1, -1, 3, 4],
    'J23':  [-1, 0, 1, -1, 3, 4, 5],
    'J24':  [-1, 0, 1, -1, 3, 4, 5, 6],
    'J33':  [-1, 0, 1, 2, -1, 4, 5, 6],
    'J222': [-1, 0, 1, -1, 3, 4, -1, 6, 7],
}
PRIMES = [3, 7, 11, 19, 23, 31]


def legendre(p):
    qr = set((x * x) % p for x in range(1, p))
    return [0] + [1 if d in qr else -1 for d in range(1, p)]


def brute_full(p, prev):
    """Pure-Python: NO symmetry reduction at all. Only for p=7."""
    chi = legendre(p)
    n = len(prev)
    if n > p:
        return 0
    tot = 0
    for tup in permutations(range(p), n):
        w = 1
        for i in range(1, n):
            if prev[i] >= 0:
                w *= chi[(tup[i] - tup[prev[i]]) % p]
        tot += w
    return tot


def brute_trans(p, prev):
    """Pure-Python: translation reduction ONLY (fix x_0=0, multiply by p)."""
    chi = legendre(p)
    n = len(prev)
    if n > p:
        return 0
    tot = 0
    for rest in permutations(range(1, p), n - 1):
        tup = (0,) + rest
        w = 1
        for i in range(1, n):
            if prev[i] >= 0:
                w *= chi[(tup[i] - tup[prev[i]]) % p]
        tot += w
    return p * tot


def main():
    print("=" * 78)
    print("STEP 1: build + run the C enumerator (independent method)")
    print("=" * 78)
    csrc = os.path.join(HERE, "paley_cluster_sums_kpc1_verify_helper.c")
    cexe = os.path.join(HERE, "paley_cluster_sums_kpc1_verify_helper.exe")
    subprocess.run(["gcc", "-O3", "-o", cexe, csrc], check=True)
    out = subprocess.run([cexe], capture_output=True, text=True, check=True).stdout
    VAL = {p: {} for p in PRIMES}
    for line in out.strip().splitlines():
        name, ps, vs = line.split()
        VAL[int(ps)][name] = int(vs)
    print(f"  C run done ({time.time()-T0:.1f}s)")

    print()
    print("=" * 78)
    print("STEP 2: validate C engine vs pure-Python brute force")
    print("=" * 78)
    ok = True
    for name, prev in STRUCTS.items():
        bf = brute_full(7, prev)
        match = (bf == VAL[7][name])
        ok &= match
        print(f"  p=7  {name:>5}: brute(full)={bf:>10d}  C={VAL[7][name]:>10d}  "
              f"{'OK' if match else 'MISMATCH'}")
    for name in ['A2', 'A4', 'A5', 'A6', 'J22', 'J23', 'J24', 'J33', 'A8', 'J222']:
        bf = brute_trans(11, STRUCTS[name])
        match = (bf == VAL[11][name])
        ok &= match
        print(f"  p=11 {name:>5}: brute(trans)={bf:>10d}  C={VAL[11][name]:>10d}  "
              f"{'OK' if match else 'MISMATCH'}  ({time.time()-T0:.0f}s)")
    assert ok, "C ENGINE FAILED VALIDATION"
    print("  C engine validated (incl. scaling reduction + last-level trick).")

    print()
    print("=" * 78)
    print("STEP 3: full table of exact integrals (independent values)")
    print("=" * 78)
    names = ['A2', 'A4', 'A5', 'A6', 'A8', 'J22', 'J23', 'J24', 'J33', 'J222']
    print(f"{'p':>3} | " + " | ".join(f"{n:>13}" for n in names))
    for p in PRIMES:
        print(f"{p:>3} | " + " | ".join(f"{VAL[p][n]:>13d}" for n in names))

    print()
    print("deltas: d22=J22-A2^2, d24=J24-A2*A4, d33=J33, d3=J222-A2^3-3*d22*A2")
    DELTA = {}
    print(f"{'p':>3} | {'d22':>12} | {'d24':>15} | {'d33':>8} | {'d3':>15}")
    for p in PRIMES:
        V = VAL[p]
        d22 = V['J22'] - V['A2'] ** 2
        d24 = V['J24'] - V['A2'] * V['A4']
        d33 = V['J33']
        d3 = V['J222'] - V['A2'] ** 3 - 3 * d22 * V['A2']
        DELTA[p] = dict(d22=d22, d24=d24, d33=d33, d3=d3)
        print(f"{p:>3} | {d22:>12d} | {d24:>15d} | {d33:>8d} | {d3:>15d}")

    print()
    print("scaled coefficients (E4 claim: d22/p^3 -> ~-10, d24/p^4 -> ~-28, "
          "d3/p^4 -> ~+262):")
    for p in PRIMES[1:]:
        D = DELTA[p]
        print(f"  p={p:>2}: d22/p^3={D['d22']/p**3:>9.4f}  d24/p^4={D['d24']/p**4:>9.4f}"
              f"  d3/p^4={D['d3']/p**4:>9.4f}")

    print()
    print("=" * 78)
    print("STEP 4: check WORKER CLAIMS (computed independently above)")
    print("=" * 78)
    claims = [
        ("A2 = p(p-1) at all p", all(VAL[p]['A2'] == p * (p - 1) for p in PRIMES)),
        ("A6(31) = 3281040",  VAL[31]['A6'] == 3281040),
        ("A8(19) = 10243584", VAL[19]['A8'] == 10243584),
        ("A8(23) = 34327040", VAL[23]['A8'] == 34327040),
        ("A8(31) = 202487040", VAL[31]['A8'] == 202487040),
        ("d22(31) = -292020", DELTA[31]['d22'] == -292020),
        ("d24(31) = -25935840", DELTA[31]['d24'] == -25935840),
        ("d3(31) = +241614000", DELTA[31]['d3'] == 241614000),
        ("A5 = 0 at all p",   all(VAL[p]['A5'] == 0 for p in PRIMES)),
        ("J23 = 0 at all p",  all(VAL[p]['J23'] == 0 for p in PRIMES)),
        ("J33 = 0 at all p (E2 lemma)", all(VAL[p]['J33'] == 0 for p in PRIMES)),
        ("p=7 anchors A2=42,A4=336,J22=-336",
         VAL[7]['A2'] == 42 and VAL[7]['A4'] == 336 and VAL[7]['J22'] == -336),
        ("monad A4 {7:336,11:1760,19:10944,23:20240,31:52080}",
         all(VAL[p]['A4'] == v for p, v in
             {7: 336, 11: 1760, 19: 10944, 23: 20240, 31: 52080}.items())),
        ("monad A6 {7:1008,11:22880,19:361152,23:870320}",
         all(VAL[p]['A6'] == v for p, v in
             {7: 1008, 11: 22880, 19: 361152, 23: 870320}.items())),
        ("monad A8 {7:0,11:154880}",
         VAL[7]['A8'] == 0 and VAL[11]['A8'] == 154880),
    ]
    allok = True
    for txt, res in claims:
        allok &= res
        print(f"  [{'CONFIRMED' if res else 'REFUTED  '}] {txt}")

    dump = {str(p): {**VAL[p], **DELTA[p]} for p in PRIMES}
    jpath = os.path.join(ROOT, "05-knowledge", "results",
                         "paley_cluster_sums_kpc1_verify.json")
    with open(jpath, "w") as f:
        json.dump(dump, f, indent=1)
    print(f"\n  exact integers dumped to {jpath}")
    print(f"  VERDICT STEP 4: {'ALL CONFIRMED' if allok else 'SOME REFUTED'}")
    print(f"  [elapsed {time.time()-T0:.1f}s]")


if __name__ == "__main__":
    main()
