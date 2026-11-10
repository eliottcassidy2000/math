# Session239: faster sequence functions that lead to actual proofs

**COMPLETE, 10 November 2026. Independently audited. Chan–Pak remains OPEN.**

The session produces six exact sequence functions, a measurable A000568
implementation improvement, uniform arithmetic refinements, and a sharp
all-order theorem excluding an infinite strong symmetry sector.

## Results at a glance

1. **Faster A000568 prefix.** The factorial-scaled future-response recurrence
   computes all terms through200 in7.48s, versus27.14s for the inherited
   compressed Fraction recurrence:3.63x on this machine. Its exact16-bit
   residue prefix takes1.12s. At80 the median speedup is4.13x. This is an
   implementation improvement to an inherited quotient, not a new Burnside
   identity or a fastest-known claim. Historical compiled GMP wins at some
   smaller orders. [Proof/code](burnside/PROOF.md), [benchmarks](burnside/scale.json),
   [THM-5097](../../../01-canon/theorems/THM-5097-factorial-gauge-integral-cycle-response-recurrence.md).
2. **Related sequences share the engine.** The same integral states produce
   graphs A000088, digraphs A000273, and oriented graphs A001174 via a reversal
   palette. Strong unlabeled A051337 and strong labeled A054946 use the
   inherited SCC inverse transforms. [Palette proof](burnside/PALETTE_PROOF.md)
   and [six-function CLI](generate.py). The speed ratio above is for A000568
   only; no timing advantage is asserted for the other sequences.
3. **Choose precision to obtain a proof.** For positive odd n put
   U_n=A000568(n)/2^((n-1)/2). Computing b bits requires only cycle counts r
   with (r-1)^2/2-v2(r!)<b. Classifying the surviving strata proves
   U_n==-n (mod4) for every odd n>=3, and a mod8 refinement involving coprime
   triples. [Proof](precision/PROOF.md),
   [THM-5098](../../../01-canon/theorems/THM-5098-precision-cycle-cutoff-and-normalized-tournament-congruences.md).
4. **Retain path values to reach Chan–Pak.** A value-refined fixed-permutation
   SCC recurrence works in the multiplicative algebra [a][b]=[ab]. For a
   requested target h it admits exact divisor-only pruning. Scalar A000568
   alone forgets these values, and the identity-profile input remains costly.
   [Proof](proof_bridges/PROOF.md).
5. **Sharp sevenfold symmetry gap.** A strong tournament with an automorphism
   of order7 either has exactly7 vertices and175 or189 Hamiltonian paths,
   or has at least2163 paths. Equality2163 is attained on9 vertices.
   The proof combines an orbit-profile split, a module in a strong quotient,
   and the exact two-path-cover cost. Thus every strong realization of203 or217
   is asymmetric, using the inherited Busch order bound as well. Neither
   assertion claims existence. [THM-5099](../../../01-canon/theorems/THM-5099-value-refined-scc-recurrence-and-sharp-seven-symmetry-gap.md),
   [independent gap audit](audit/strong_gap_audit.md),
   [asymmetry audit](audit/palette_and_asymmetry_audit.md).

## Use the functions

From repository root, standard Python3 only for the delivered functions:

    python3 05-knowledge/results/session239_sequence_engines_20261110/generate.py --sequence A000568 --through 100
    python3 05-knowledge/results/session239_sequence_engines_20261110/generate.py --sequence A000568 --through 200 --binary-precision 16
    python3 05-knowledge/results/session239_sequence_engines_20261110/generate.py --sequence A051337 --through 30
    python3 05-knowledge/results/session239_sequence_engines_20261110/generate.py --sequence A000088 --through 20

The same CLI accepts A000273, A001174 and A054946. JSON terms are strings to
preserve large integers in consumers. Offsets match OEIS: A001174/A054946
start at1; the others start at0. The internal palette engine also supplies
the natural empty-object term at0. A051337's extra a(0)=1 is an OEIS
convention and is kept separate from the nonempty strong-component series.

The full-count algorithm remains finite-horizon and potentially expensive.
The residue flag returns residues, not compressed full integers. The
normalized few-cycle function is `precision/precision.py:normalized_low_bits`.

## Verification and inheritance

- Four counting engines agree through80; the external A000568 reference
  agrees through76; historical C/GMP independently agrees at80,100,120.
- The optimized and inherited compressed engines agree through200.
- Independent divisor-square Davis evaluation checks all0..45; residue/horizon
  controls extend to120, including composite moduli and zero residues.
- The normalized congruences are independently checked through119; the proof
  covers all stated odd n. The low-bit program checks n1001 without computing
  its full A000568 value.
- Every nonempty labeled tournament through6 is checked for the SCC/value
  recurrence, with independent orbit BFS and319 target-divisor controls.
- The sharp gap audit checks all9! vertex orders of each of the eight
  seven-circulant triangle inflations: six values2163 and two2205.
- Palette code is checked against direct partition sums through12; a separate
  audit exhausts every color assignment and permutation through4 for all
  nine realizable palettes with q<=4. CLI prefixes and offsets are tested.

[Reproduction commands](REPRODUCE.md), [independent audit](audit/README.md),
[original plan](PLAN.md). Core/tail factorization, Gaussian identities,
non-P-recursiveness, SCC scalar transforms, valuation oddness, and substitution
kernels are inherited and explicitly credited in the proof files. No external
novelty claim is made. The strong gap is an all-order consequence with a finite
exhaustive certificate, not merely a small-order census.

## What remains open

A genuinely faster general Hamiltonian-path spectrum engine, larger-prime
sharp symmetry gaps, and unrestricted Chan–Pak coverage remain open. The
scalar count improvements do not establish which individual path values
occur. Strongness cannot be dropped from the2163 bound: singleton padding
around Q7 keeps175 paths at arbitrarily large orders. LRC(14) is untouched.

## Sources checked for definitions and indexing

[A000568](https://oeis.org/A000568), [A051337](https://oeis.org/A051337),
[A054946](https://oeis.org/A054946), [A000088](https://oeis.org/A000088),
[A000273](https://oeis.org/A000273), [A001174](https://oeis.org/A001174).
Primary mathematical dependencies already in canon are listed by exact file
in the proof packets; timings and new finite certificates are in this packet.
