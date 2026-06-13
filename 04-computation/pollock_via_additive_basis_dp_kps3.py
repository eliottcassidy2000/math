#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pollock's conjectures via the repo's additive-basis minimal-summands DP.
kind-pasteur-2026-06-13-S3.

The repo's additive-basis machinery (HYP-1953/1962/1963, additive_basis_normal_forms_s494)
computes, for an atom family B, the EXACT minimal number of summands g_B(n) to write n
as a sum of atoms of B.  Pollock's conjectures are precisely the bounded-arity statement
"max_n g_B(n) = k" for B = tetrahedral / octahedral / icosahedral / dodecahedral / cubes.

We swap the atom generator to the 3D figurate (polytope) numbers and run the standard
coin-problem DP up to a bound, reporting:
  * the running max arity (the "Pollock bound" as seen up to N),
  * the EXCEPTIONAL SET (n requiring the full bound) for each family,
and cross-check against the literature (effective asymptotic + finite check; the
Basak-Dong-Saettone-Zaharescu IMRN 2025 refinements, Brady 2016, Watson 1952):
  - cubes (Waring g(3)=9): only 23, 239 need 9.
  - tetrahedral Te_k=C(k+2,3): 241 exceptions need a 5th; largest 343867 (OEIS A000797). OPEN.
  - octahedral O_k=k(2k^2+1)/3: Pollock bound 7 (Brady 2016, large N).
  - icosahedral I_k=k(5k^2-5k+2)/2: Pollock's 13 is WRONG; true bound 15; the 13-counterexamples
    are {47,83,94,95,119} (95 needs 15). Basak et al. 2025.
  - dodecahedral D_k=k(3k-1)(3k-2)/2 ... (use the standard formula); Pollock's 21 WRONG; true 22;
    sole counterexample 79. Basak et al. 2025.
Nothing hardcoded as "expected" (MISTAKE-062): we COMPUTE the exceptional sets and compare.
"""

import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def atoms_upto(gen, N):
    out = []
    k = 1
    while True:
        v = gen(k)
        if v > N:
            break
        if v >= 1:
            out.append(v)
        k += 1
    return sorted(set(out))


def min_summands_dp(atoms, N):
    """g[n] = min # atoms summing to n (n>=1); g[0]=0. coin-problem DP (positive atoms)."""
    INF = float('inf')
    g = [0] + [INF] * N
    for a in atoms:
        for n in range(a, N + 1):
            if g[n - a] + 1 < g[n]:
                g[n] = g[n - a] + 1
    return g


def report(name, gen, N, lit_bound=None, lit_note=""):
    atoms = atoms_upto(gen, N)
    g = min_summands_dp(atoms, N)
    # ignore n=0; require all n in [1,N] representable (figurate families incl. 1 do)
    reps = [g[n] for n in range(1, N + 1)]
    if any(r == float('inf') for r in reps):
        unr = [n for n in range(1, N + 1) if g[n] == float('inf')]
        print(f"   {name}: NOT all representable up to {N} (missing {len(unr)}, first {unr[:5]}) — atom set issue", flush=True)
        return
    mx = max(reps)
    # the exceptional set at each arity level
    from collections import Counter
    cnt = Counter(reps)
    need_max = sorted(n for n in range(1, N + 1) if g[n] == mx)
    # also the set needing >= (mx) and the threshold beyond which g<=mx-1
    last_at_max = need_max[-1] if need_max else None
    print(f"   {name}: max arity up to {N} = {mx}; "
          f"count needing {mx} = {cnt[mx]} (largest {last_at_max}); "
          f"arity histogram {dict(sorted(cnt.items()))}", flush=True)
    if lit_bound is not None:
        verdict = "MATCHES" if mx == lit_bound else f"DIFFERS (lit={lit_bound})"
        print(f"      vs literature bound {lit_bound}: {verdict}.  {lit_note}", flush=True)
    return mx, need_max, g


def main():
    print("=== Pollock's conjectures via the repo additive-basis min-summands DP ===", flush=True)
    print("   (atom generator swapped to 3D figurate numbers; nothing hardcoded)\n", flush=True)

    N = 400000
    print(f"--- cubes  k^3  (Waring g(3)=9; only 23,239 need 9) ---", flush=True)
    mx, nm, g = report("cubes", lambda k: k**3, N, 9, "expect exactly {23,239} at arity 9")
    print(f"      numbers needing 9 cubes: {nm}", flush=True)

    print(f"\n--- tetrahedral  C(k+2,3)  (Pollock 5; OPEN; 241 exceptions, largest 343867) ---", flush=True)
    mx, nm, g = report("tetrahedral", lambda k: k*(k+1)*(k+2)//6, N, 5,
                       "expect 241 numbers at arity 5, largest 343867 (A000797)")
    print(f"      #needing 5 tetrahedral = {len(nm)}; largest = {nm[-1]}; first few = {nm[:8]}", flush=True)

    print(f"\n--- octahedral  k(2k^2+1)/3  (Pollock 7; Brady 2016 large N) ---", flush=True)
    report("octahedral", lambda k: k*(2*k*k+1)//3, N, 7, "Brady: 7 suffices for large N")

    print(f"\n--- icosahedral  k(5k^2-5k+2)/2  (Pollock 13 WRONG; true 15; cex {47,83,94,95,119}) ---", flush=True)
    mx, nm, g = report("icosahedral", lambda k: k*(5*k*k-5*k+2)//2, 200000, 15,
                       "Basak et al 2025: true bound 15; 13 fails at 47,83,94,119; 95 needs 15")
    # report the numbers needing >13 (the counterexamples to Pollock's 13)
    over13 = sorted(n for n in range(1, 200001) if g[n] > 13)
    print(f"      numbers needing >13 icosahedral (counterexamples to Pollock's 13): {over13}", flush=True)

    print(f"\n--- dodecahedral  k(3k-1)^2... use k(3k^2-3k+1)? standard: (3k^3-3k^2+k)? ---", flush=True)
    # dodecahedral number D_k = k(3k-1)(3k-2)/2  (OEIS A006566)
    mx, nm, g = report("dodecahedral", lambda k: k*(3*k-1)*(3*k-2)//2, 200000, 22,
                       "Basak et al 2025: true bound 22; 21 fails only at 79")
    over21 = sorted(n for n in range(1, 200001) if g[n] > 21)
    print(f"      numbers needing >21 dodecahedral (counterexamples to Pollock's 21): {over21[:20]}", flush=True)


if __name__ == "__main__":
    main()
