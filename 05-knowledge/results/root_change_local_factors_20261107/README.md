# Local factors of the two-root-change multiplicity; small-number coincidences; the cutoff barrier

kind-pasteur-2026-11-07-S15. **THM-4916** (PROVED: the closed-form local factors for every odd prime; identities,
mostly credited) and **HYP-9506** (census-scale factors and the cutoff). The independent audit is PASS WITH
CORRECTIONS, applied ([audit/AUDIT.md](audit/AUDIT.md); MISTAKE-734). Pre-registration
[PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 0cda244632. Chan-Pak and the two-root-change cutoff remain OPEN.

Owner prompt: compare small-number coincidences whose structure we may be able to explain (Chan-Pak, Fermat primes,
the square-sum conjecture); derive the 3-adic local factor of the two-root-change multiplicity; creatively construct a
proven cutoff beyond which two root changes always suffice.

**Concurrent same-prompt work, credited.** Session223 (`../session223_root_changes_and_local_factors_20261107/`) has:
- the two-run product rectangle F(t.Z^a O^b) = (2^a + 1)(2^b + 1);
- the square-sum graph on {2^a + 1}, which has only the edge 3-33;
- in checkpoint 2 (1c524e45a9, pushed 32 s before this packet's checkpoint), a PROOF of the fixed-length 3-adic law
  at every precision: mod 3, 34/27 and 47/54, and with one J, 41/36, 8/9 and 35/36. That proof uses the
  Cesaro-projector method;
- the exact two-J membership equation, whose eq. (2) is THM-4916 (1b).

## Results

1. **The local factor, derived for every odd prime (THM-4916; [PROOF.md](PROOF.md), section 2.1).** In the letter-count
   limit:
   - with two root changes, L_p(0) = 1 + (3p-2)(p-1)/(2p^3) and every other class has 1 - (3p-2)/(2p^3);
   - with one root change, L_p(0) = 1 + (3p-4)/(4p^2) and L_p(-1) = 1 - 1/p^2 + 1/(4p).

   The proof is a confinement cascade at the root changes: a block with d = U = 0 (mod p) freezes the value class.
   The audit proved the block-mixing lemma. The method and p = 3 are Session223's; the general-p formulas are new.
   The local factors are not exactly multiplicative across primes (mod 15: ratios 0.994 to 1.009).
2. **Census scale (section 2).**
   - The census is more contrasted than the letter-count model.
   - Pre-registered: R3 (ordering 5 > 17 > 11 > the rest) and R4 (uniform nonzero classes) CONFIRMED. R5 (classes mod
     27) CONFIRMED only through its 5% widening. R2 (per-prime bands) FAILED at 5, 11, 13 and 17.
   - The audit splits the zero-class excess at 3, 5 and 17: about a quarter is Fermat-type multiplier words, a fifth is
     root-rectangle words Z^a O^b J (divisible by 2^a + 1), and about half is seed confinement.
   - Near-misses avoid 3 and are depleted in 5, 7, 11, 17 and 19.
3. **Identities (section 1), with priority.**
   - (1a) is KNOWN (Session206).
   - (1b) is concurrent with Session223's eq. (2). New: the m, Q(m) form and tail confinement.
   - (1c) Fermat-type multipliers and (1e) the J-budget product family V V in F_<=2 are one-line corollaries of known
     identities.
   - (1d) is Session199's K-form with a = -f.
4. **The cutoff (section 3): no proof.**
   - Barrier lemma: every product or multiplier family is composite, so a prime needs the non-product part of the
     K-form.
   - A single seed misses a positive proportion (Session206); finite banks are OPEN.
   - The routes are a Gerbicz-type state-cone blow-up and a second-moment bound.
   - The census to 2^32 (R1) is CONFIRMED: the 140 misses below 3,780,967 are all the misses below 2^32
     ([R1_OUTCOME.md](R1_OUTCOME.md)).
5. **[ATLAS.md](ATLAS.md): the small-number coincidence atlas.** Chan-Pak {7, 21}, the C = 1 and root-change
   exceptions, the Fermat primes, square-sum paths and cycles, nine cubes, untouchable primes, the three-core holes,
   and the prime square-sum delay. Each is classified as structural, constructive induction (Gerbicz), analytic plus
   computation, or heuristic.
   - The Fermat primes are the prime multipliers 2^i + 1 of THM-4916 (1c).
   - At census scale the multiplier words and the root-rectangle words together carry under half of the boost at 3, 5
     and 17 (audit M3).

## Files

| File | Role |
| --- | --- |
| `PREREGISTRATION.md` | R1-R5 |
| `resdp.c`, `dp_bands.txt`, `dp27_classes.txt` | letter-count residue automaton and the registered bands. The L(nz)avg and ratio columns of `dp_bands.txt` come from 4-decimal per-class values; this changes no outcome |
| `dp27_long.txt`, `dp_primes_n1000.txt`, `q3_long.txt` | truncated excerpts of finite-n automaton output (not limits) |
| `exact_limit.py`, `exact_limit.out`, `exact_limit_q15.out`, `closed_form_check.out` | n -> infinity limits by the projector Pi J Pi J Pi |
| `decomp.c`, `decomp_2p24.txt` | split of M_2 by seed confinement level, parity and class |
| `rc_lf.c`, `lf_2p24.txt`, `eval_s15.py`, `PREREG_OUTCOMES_R2_R5.txt` | registered census and its evaluation |
| `R1_OUTCOME.md`, `census_2p32.txt` | R1 CONFIRMED: no miss of two root changes in (10^9, 2^32] |
| `audit/` | independent audit: AUDIT.md, its own code (`identities_audit.py` reproduces the identity checks of section 1, `cesaro_exact.py` and `cesaro_lazy.c` the closed form, `census_mechanism.c` the M3 split) |

## Reproduce

    gcc -O2 -o resdp resdp.c && ./resdp 27 400 2 && ./resdp 5 1000 2
    python3 exact_limit.py 3 2000 && python3 exact_limit.py 23 2000
    gcc -O2 -fopenmp -o rc_lf rc_lf.c && ./rc_lf 16777216 2 8388608 16777216 > lf_2p24.txt
    python3 eval_s15.py lf_2p24.txt dp_bands.txt dp27_classes.txt
    gcc -O2 -fopenmp -o decomp decomp.c && ./decomp 16777216 2 8388608 16777216
