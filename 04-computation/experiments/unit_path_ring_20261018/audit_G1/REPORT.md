# Audit G1: THM-4700, THM-4701, the unit-path results note, MISTAKE-635

Independent referee report, 2026-10-18. The subagent could not write files, so the session copied this report from its
final message. It is reproduced in substance; the scripts and their outputs are the auditor's own and sit in this
directory.

**Scope.**
- THM-4700, THM-4701, `unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md` (the "results note"), MISTAKE-635,
  and section 2.1 of the merge note, as of commit b9a2cbf698.
- None of these changed in the later session commits up to 1eef3cc9a7, except merge-note section 3.1b, which is outside
  scope.

**Method.**
- The auditor's own PARI/GP 2.17.4 and Python code; no session module was imported.
- No third-party code was run and nothing was downloaded.
- Four web searches, on ordinary result pages only.

## Summary verdict

**THM-4700: passes, with minor corrections.** Every mathematical claim in statements 1–5 is correct.
- **B2 in F holds.** The auditor built the integral basis `Z[zeta15] ⊗ Z[(1+j)/2] ⊗ Z[1, w, w^2/2]`. Its discriminant
  is `2^32 · 3^24 · 5^36 · 11^24 · 37^16`, and PARI's p-maximal discriminants at 2, 3, 5, 11, 37 match, so it is O_F.
- **The fixed-prime census holds.** `idealprimedec` in the degree-48 field for all 29 primes of S below 3200:
  - the count is `8·#roots`, or 8 when the cubic is irreducible mod p;
  - the residue fields are `F_(p^2)` or `F_(p^6)`;
  - it is 24 at 2549, 2879 and 3119;
  - there are no fixed primes above 2, 3, 5, 11, 37 or above 24 further control primes.
- **Everything else reproduced exactly:** the walks, the Pell list, the two-step criterion with `D ∈ {3, 11, 15, 55}`,
  the Heron criterion and its 6/14 census, the recipe denominators, the `delta_19` pole/zero pattern, the reflection
  links, the mixed quadrilaterals, the house numbers 3.70409 and 4.11794, and the per-pair house table.
- The 14 exclusions are also decided completely by `nfroots`, not only by non-residue certificates.

**THM-4701: passes with corrections.**
- Statements 1–5 reproduced, including the triangle and K_4 laws for every prime power `q <= 32`, and `G_4` = Clebsch
  with its Petersen layer equal to Kneser K(5,2).
- Petersen census: exhaustive recount for all primes `q <= 113` and for `q = 9, 25, 27, 49`, using a different field
  model and search order.
- Statement 4(b) is true for every `q = 19 (mod 20)`, but was proved only at four primes (G1-4).
- The "pentagon becomes mu_5" remark is false (G1-2).

**Results note: corrections required.**

**MISTAKE-635 and merge-note section 2.1: pass.**

## MAJOR

**G1-1. Results note, sections 5(a), 6, 7, 8: the persistent-twin heuristic is refuted.**

The d = 8 row is now in `clock_graded_ladder_d8.out`: twin part -0.08773, lag 1 +0.03827. The C port gives twin parts
-0.0786 and -0.0715 at d = 9 and 10.

So all three pre-registered predictions failed:
- P1 (twin part in [-0.106, -0.098]) failed.
- P2 (lag 1 above +0.040) failed.
- P3 failed, but it was implied by P1, since cross_8 was known in advance.

The session recorded this in `ladder_prereg_outcome_20261018.md` and in the merge note, but not in this note. The
auditor checked rows d = 1..8 against the S41 a4 json; they agree to within 3.9e-12.

Replacement text:
- **Section 5(a), last bullet:** the twin part is negative at every computed depth but decays from d = 8; lag 1 peaks at
  d = 7; the cross term is a difference of two slowly decaying parts; the heuristic failed its pre-registered test.
- **Section 6:** add the rows d = 8, 9, 10 and record that P1–P3 failed.
- **Section 7 ledger:** "ladder rows d = 6..10"; mark the heuristic REFUTED.
- **Section 8:** "determine the decay rates of the twin and lag-1 parts".

**G1-2. THM-4701 remark "Six-colour programme" and results note section 3: "the pentagon of ports becomes mu_5,
independent in the Clebsch graph" is false.**

- The port pentagon `{0, 1, 1 + zeta, 1 + zeta + zeta^2, -zeta^4}` has unit sides, so modulo 2 it becomes an induced
  5-cycle of `G_4` (`g1_port_pentagon_mod2.out`).
- No unit-side pentagon can reduce to `mu_5`.
- What reduces to `mu_5` is the pentagon's set of side vectors. Its diagonal vectors are the consecutive two-sums.

Replacement text: "The side vectors of the port pentagon are mu_5 and its diagonal vectors are the consecutive two-sums.
Modulo 2 the side vectors become the independent neighbourhood mu_5 of 0, the ten two-sums become the Petersen layer, and
phi becomes a cube root of unity. The unit-side pentagon itself becomes an induced 5-cycle. So the K_5 and the Petersen
graph are two shadows of the same five directions."

This has no downstream impact.

## MINOR

**G1-3. THM-4700 statement 3 corollary and its proof.** "L(x) = 2 iff N(4 - N)/3 ∈ (F+)^{*2}" is false for x in U,
where L(x) = 1 but the condition holds. Say instead: "x is a sum of two unit vectors iff ...". The uses in statements 2
and 5(e) are unaffected.

**G1-4. THM-4701 statement 4(b).** Injectivity is claimed for every q but proved only at q = 19, 59, 79, 139.
- The absolute norms of all 45 pairwise differences of the ten points are powers of 5, so the reduction is injective at
  every prime not above 5.
- The auditor checked directly at q = 19, 59, 79, 139, 179, 199, 239, 359.

**G1-5. THM-4701 "New here" overstates priority.**
- The golden realization is the 5-fold rotationally symmetric unit-distance drawing of GP(5,2).
  - Žitnik–Horvat–Pisanski (J. Korean Math. Soc. 49, 2012) give rotationally symmetric unit-distance drawings of
    I-graphs, including GP(5,2).
  - Kaenthong et al. (2025) study the rotation angle of the inner star.
  - The spoke condition forces a 90-degree turn, which is exactly R^2 + r^2 = 1.
  - Retype the realization as KNOWN; keep the 2-adic reduction as new.
- The two-sum/Kneser picture is the classical folded-5-cube model of the Clebsch graph: KNOWN.
- Le Anh Vinh (arXiv math/0510092) studies the unit-quadrance graph on F_q^2, which is G_q for q = 3 (mod 4). For the
  triangle and K_4 claims, write "not found in the literature checked" rather than "new".

**G1-6. Results note, sections 1 and 5(a): ladder claims overstated.**
- "Exact cross terms for d <= 12": only d <= 6 come from an exact integer push. The values for d = 7..12 are float64
  with a rigorous tail bound.
- "No gap; the decay is polynomial" rests on 7 points. On d = 8..12 the two fits are close: power exponent 1.594
  (rms 0.021) against geometric ratio 0.852 (rms 0.032). Type it HEURISTIC.
- "Every cross term is negative, so the tail after d = 12 totals at most E_12 (PROVED)" assumes a sign at all depths,
  which THM-4667 says is not known. Replace with: since E_d >= 0, the cross terms after d = 12 sum to at least
  -E_12 = -1.479, so they cannot be eventually <= -C d^(-alpha) with C > 0 and alpha <= 1 (PROVED, elementary).
  The same wording occurs in `ladder_decay_exponents.py`.
- The numbers themselves reproduce: exponent 1.914, ratio 0.8051, local exponents 2.4, 2.4, 1.9, 1.8, 1.1, 1.7.

**G1-7. Results note section 3.** "The gaps must come from Frobenius permuting components" is too strong. Lang–Weil
error terms for a high-degree variety, and the removal of non-injective points, can also produce zero counts. Change
"must" to "could" and type it HEURISTIC.

**G1-8. THM-4700 remark "Eisenstein triangles" and results note section 4.** Only the sqrt(-3) walks are Eisenstein:
1/7, 1/19 and the two-step 1/13. The walks for 1/3, 1/4, 1/17 and 1/23 come from Q(sqrt-11) and Q(sqrt-15).

**G1-9. No script behind `two_step_reciprocals.out`.** Add one, or cite `g1_reciprocals_pythagoras.py`, which reproduces
the list exactly.

**G1-10. Results note section 2, Heron bullet.** The "iff" omits the hypothesis "some z has N(z) = R_1". Without it the
statement is false: an equilateral triangle with a non-norm squared side is a counterexample. THM-4700 has the
hypothesis; copy it.

**G1-11. THM-4701 proof 5.** "Each new vertex has at least two placed neighbours" is false for search vertices 2, 3 and
5, which have one each. Correctness is unaffected; fix the wording.

## Notes

**N1. Reflection links (THM-4700 5(c)).** The construction degenerates (q' = q) only if N(k) is a square in F+. N(b)
and N(c) are not squares in F, so q' ≠ q in all 120 + 120 cases.

**N2. Mixed quadrilaterals (THM-4700 5(e)).**
- The proof should state distinctness: when N(k) ≠ 1, |u_1| = 1 differs from |e_2/k| = d_k.
- The identity (N, N) = (N, -1) is not needed; (N, N - 1) = (N, -1)(N, 1 - N) suffices.

**N3. Typing.**
- The 14 exclusions and the two-step part of 5(d) are rigorous certificate proofs and can be typed PROVED.
- A simpler complete argument for 5(d): only 2 (for b) + 6 (for c) cross-copy displacements have house <= 2, and none
  of them is a sum of two unit vectors.

**N4. Walk-length quartic.** A_x(lambda) and the Cayley parametrization of U were verified exactly. The restriction
"y ≠ 0" drops the degenerate cases x in U and x in U + 2U, which also give L(x) <= 3.

**N5. MISTAKE-635.**
- Its "fix check" is tautological: the corrected script prints 8·#roots by formula.
- The independent confirmation is `g1_fixed_primes.out`.
- The "Error" bullet should also mention the irreducible case.
- THM-4700's verification line says "p <= 2099", but the `.out` now runs to 2999.

**N6. Prior art.**
- Cite Chilakamarri (Discrete Math. 69, 1988) for the components of the rational unit-distance graph, i.e. B2 for Q(i).
- The overring lemma is classical.
- No literature was found on the near-60-degree triples or the C, A branch.
- The Yang–Mills README was not findable by search and was not fetched.

**N7. Pre-registration.** P3 was implied by P1.

**N8. Rigidity singular value.** The quoted 1.19 depends on how the rows are normalized; with rows p_u - p_v it is
0.595. Rank 15 is confirmed. The realization is also induced: exactly 15 of the 45 pairs are at unit distance.

**N9. Walk lengths in statement 2.** These are lengths of the given walks, i.e. upper bounds for L(1/n). They are
minimal only in the two-step rows; for example L(1/3) >= 3.

**N10. Exhaustive Petersen counts beyond the theorem.** These match the session's later C counter exactly.

| q | 53 | 59 | 61 | 67 | 71 | 73 | 79 | 83 | 89 | 97 | 101 | 103 | 107 | 109 | 113 | 49 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| copies | 480 | 5568 | 7680 | 3840 | 6720 | 0 | 768 | 240 | 240 | 3840 | 9600 | 19200 | 5040 | 7680 | 10080 | 0 |

**Not audited:** the session's later additions (Petersen census to q <= 397, the three-step connector screen, and the
d = 9, 10 rows beyond quoting them).

## Reproduced independently

PARI scripts are run as `gp -q -s 800000000 <file>` from inside `audit_G1/`. `g1_F.bin` (2 MB) was deleted; running
`g1_buildF.gp` first recreates it in under 2 s.

| Script | Output | What it checks |
|---|---|---|
| `g1_probe.gp`, `g1_buildF.gp`, `g1_common.gp` | `g1_probe.out`, `g1_buildF.out` | Builds F and certifies the integral basis |
| `g1_fixed_primes.gp` | `g1_fixed_primes.out` | Fixed-prime census |
| `g1_identities.gp` | `g1_identities.out` | Generators, walks, recipe, reflection links, mixed quadrilaterals |
| `g1_imag_quadratic.gp` | `g1_imag_quadratic.out` | sqrt(-D) lies in F exactly for D ∈ {3, 11, 15, 55} |
| `g1_heron.gp` | `g1_heron.out` | Heron census and house numbers |
| `g1_connectors.gp` | `g1_connectors.out` | House table and the exact 1- and 2-step connector census |
| `g1_nineteen.gp` | `g1_nineteen.out` | Primes above 19 in F and F', delta_19 |
| `g1_quartic.gp` | `g1_quartic.out` | Walk-length quartic |
| `g1_spindle.gp` | `g1_spindle.out` | M has 11 unit edges, C has 15 |
| `g1_petersen_zeta20.gp` | `g1_petersen_zeta20.out` | Golden Petersen realization and its reductions |
| `g1_port_pentagon_mod2.gp` | `g1_port_pentagon_mod2.out` | Counterexample for G1-2 |
| `g1_reciprocals_pythagoras.py` | `g1_reciprocals_pythagoras.out` | S, two-step reciprocals, Pell, Berggren tree to depth 9 |
| `g1_gq.py` (modes `local`, `census`, `count113`) | `g1_gq_local.out`, `g1_gq_census.out`, `g1_gq_count113.out` | Finite unit-circle graphs and the Petersen census |
| `g1_gq_small.py` | `g1_gq_small.out` | G_2 = K_4, G_8 = H(3,4) |
| `g1_petersen_rigidity_unitsum.py` | `g1_petersen_rigidity_unitsum.out` | Rigidity rank, unit-sum statement |
| `g1_ladder_check.py` | `g1_ladder_check.out` | Ladder exponents and fits, rows d = 1..8 |

Prior-art pointers found by the auditor:
- Žitnik–Horvat–Pisanski, unit-distance I-graphs (arXiv 1207.2752);
- Kaenthong et al., a property of a particular generalized Petersen unit-distance graph (2025);
- Le Anh Vinh, chromatic number of unit-quadrance graphs (arXiv math/0510092);
- Chilakamarri 1988, components of the rational unit-distance graph;
- the finite Euclidean graph literature (Medrano–Myers–Stark–Terras and follow-ups, e.g. arXiv 0807.2689).
