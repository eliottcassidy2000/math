# Local factors of the two-root-change multiplicity; small-number coincidences; the cutoff barrier

kind-pasteur-2026-11-07-S15. **THM-4916** (PROVED identities) and **HYP-9506** (local-factor conjectures and the
cutoff). Pre-registration [PREREGISTRATION.md](PREREGISTRATION.md) was pushed in 0cda244632 before any registered
census value was computed. Chan-Pak and the two-root-change cutoff remain OPEN.

Owner prompt: compare small-number coincidences whose structure we may be able to explain (Chan-Pak, Fermat primes,
the square-sum conjecture); derive the 3-adic local factor of the two-root-change multiplicity; creatively construct a
proven cutoff beyond which two root changes always suffice.

**Concurrent same-prompt work.** Session223 (`../session223_root_changes_and_local_factors_20261107/`) has:
- the two-run product rectangle F(t.Z^a O^b) = (2^a + 1)(2^b + 1);
- the square-sum graph on {2^a + 1}, which has only the edge 3-33;
- a height-weighted 3-adic factor and the cutoff, both under investigation.

It is credited and not duplicated.

## Results

1. **THM-4916, exact identities ([PROOF.md](PROOF.md), section 1).**
   - Seed confinement: tail values are h - e (m mod 2) mod gcd(R,d).
   - The root-change tail in pre-change coordinates. For odd tails, (m - Q(m))D + (m+1)U, so gcd(Q(m)+1, m+1) divides
     the value.
   - Fermat-type multipliers: value(p J Z O^(j-1) Z) = (2^(j-1)+1)(D + 2U).
   - The defect a = D - E - U + d and the expansion K(v,w) = U_v U_w + E_v E_w - d_v d_w + a_v d_w + a_w d_v of
     Session199's form.
   - The two-state product family V V inside F_<=2.
2. **The 3-adic local factor, derived (section 2).** N mod q is the image of the word under a finite residue
   automaton, and both kinds of confinement make the zero classes heavy.
   - The letter-count automaton (`resdp.c`) gives divisor-type factors: L_p(0) - 1 is about 1.5/p asymptotically.
     For p = 3 the limit is about 34/27 against 47/54.
   - At census scale the Fermat-type multiplier primes 3, 5, 17 (and 11) dominate, and the effect decays like
     1/log N.
   - Pre-registered census tests: R3 (ordering 5 > 17 > 11 > rest), R4 (uniform nonzero classes) and R5 (27/27
     classes mod 27) CONFIRMED. R2 (quantitative band per prime) FAILED: the census boost is larger than the
     letter-count automaton at 5, 11, 13 and 17, which points to height weighting.
   - Near-misses are rough: never divisible by 3, mostly 1 (mod 3), often prime.
3. **The cutoff (section 3).** No proof.
   - Proved: explicit families inside F_<=2, and a barrier lemma (product and multiplier families miss every prime, so
     a proof must use the non-product part of the K-form).
   - The two missing routes are a Gerbicz-type state-cone blow-up and a second-moment bound for the bilinear form.
   - R1 extends the census to 2^32 (section 4).
4. **[ATLAS.md](ATLAS.md): the small-number coincidence atlas.** Chan-Pak {7, 21}, the C = 1 and root-change
   exceptions, the Fermat primes, square-sum paths and cycles, nine cubes, untouchable primes, the three-core holes,
   and the prime square-sum delay. Each is classified as structural, constructive-induction (Gerbicz), analytic plus
   computation, or heuristic, with the mechanism of each. Fermat primes enter the root-change grammar as exactly the
   cheap multipliers of THM-4916 (1c).

## Files

| File | Role |
| --- | --- |
| `PREREGISTRATION.md` | R1-R5 |
| `resdp.c`, `dp_bands.txt`, `dp27_classes.txt`, `dp27_long.txt`, `dp_primes_n1000.txt`, `q3_long.txt` | residue automaton and its predictions |
| `decomp.c`, `decomp_2p24.txt` | exact split of M_2 by seed confinement level, parity and class |
| `rc_lf.c`, `lf_2p24.txt`, `eval_s15.py`, `PREREG_OUTCOMES_R2_R5.txt` | registered census and its evaluation |
| `census_2p32.txt` (R1) | the S14 auditor's independent segmented census, to 2^32 |

## Reproduce

    gcc -O2 -o resdp resdp.c && ./resdp 27 400 2 && ./resdp 5 1000 2
    gcc -O2 -fopenmp -o rc_lf rc_lf.c && ./rc_lf 16777216 2 8388608 16777216 > lf_2p24.txt
    python3 eval_s15.py lf_2p24.txt dp_bands.txt dp27_classes.txt
    gcc -O2 -fopenmp -o decomp decomp.c && ./decomp 16777216 2 8388608 16777216
