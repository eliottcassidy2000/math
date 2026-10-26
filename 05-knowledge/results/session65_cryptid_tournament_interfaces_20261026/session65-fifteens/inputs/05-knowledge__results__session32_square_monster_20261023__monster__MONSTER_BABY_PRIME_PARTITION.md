# The Monster/Baby Monster prime partition and its level-two geometric defect

Session 32, 2026-10-23. Status: CITED group/modular facts; PROVED elementary consequences; FINITE-EXACT polynomial and independent finite-field certificates. No claim of a tournament-spectrum explanation or a new moonshine theorem.

## 1. The integer and all three sets have an exact common source

Write M for the Monster and B for the Baby Monster. The ATLAS gives

|M| = 2^46 3^20 5^9 7^6 11^2 13^3 17 19 23 29 31 41 47 59 71,

|B| = 2^41 3^13 5^6 7^2 11 13 17 19 23 31 47.

Sources: [ATLAS Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) and [ATLAS Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/). Only the displayed orders and the 2.B subgroup are used; the old Monster page's historical maximal-subgroup completeness statements are not inherited.

The user's partition is therefore exactly the following canonical partition relative to this pair of orders:

| Set | Primes | Valuation predicate |
|---|---|---|
| A, unchanged common primes | 17,19,23,31,47 | v_p(|M|)=v_p(|B|)=1 |
| D, new primes relative to B | 29,41,59,71 | v_p(|M|)=1, v_p(|B|)=0 |
| C, enlarged common primes | 2,3,5,7,11,13 | v_p(|M|)>v_p(|B|)>0 |

In this example C also equals the primes p for which p^2 divides |M|. This is a second, Monster-only way to recognize that particular stratum, not a universal theorem about arbitrary pairs of finite groups. No established name for this three-way partition was found in the sources searched. Its definition is canonical without needing a special name.

For

N = 2^5 3^7 5^3 7^4 11 13^2 = 39,046,339,332,000,

we obtain the exact identity

|M|/|B| = N × 29 × 41 × 59 × 71 = 194,478,922,284,018,372,000.

The product of the four new primes is 4,980,721. The five A-primes cancel completely; the six C-primes contribute precisely N; the four D-primes contribute one copy each. This explains every supplied set and every exponent of N together. It is not merely a shared collection of small primes.

The centralizer normalization is important. For a Monster involution t of class 2A,

C_M(t) ≅ 2.B, and B ≅ C_M(t)/<t>.

Thus the actual involution conjugacy-class size is

|t^M| = |M|/(2|B|) = 97,239,461,142,009,186,000
       = (N/2) × 29 × 41 × 59 × 71.

A primary account is Cummins–Matias, [(2+)-Replication and the Baby Monster](https://arxiv.org/pdf/1710.01071), SIGMA 14 (2018), 060, introduction p.2 and §5. There is no claim here that B itself is the chosen subgroup, that N is a subgroup index, or that the numerical ratio |M|/|B| alone is a quotient group. A full Sylow p-subgroup fits inside the centralizer for p in A, whereas the index has a positive p-valuation for p in C∪D; this is an index statement, not a new centralizer classification.

## 2. The four missing Baby primes have a genuine modular meaning

Aricheta, [Supersingular Elliptic Curves and Moonshine](https://arxiv.org/html/1809.07421), SIGMA 15 (2019), 007, provides the precise bridge. His Theorem 1.1 gives Q_p(X)/2 = g(X^p)−g(X) for p not dividing the level, where Q counts supersingular geometric points not rational over F_p. Theorem 1.3 identifies rationality with the presence of a prime-order Fricke element in the corresponding Monster centralizer. For X=X_0^+(2), Table 1 gives exactly the odd prime divisors of |B|. Here X^p=X_0(2p)/<w_2,w_p>=X_0^*(2p). Appendix B supplies the eta coordinate below and Appendix E its j-map.

The usual Ogg list concerns X_0(1), equivalently all supersingular j-invariants being F_p-rational, or g(X_0^+(p))=0. Passing to X_0^+(2) adds the coarse Fricke orbit of a cyclic 2-isogeny. Rationality of this orbit does not require a representative oriented isogeny defined over F_p: Frobenius may exchange its two orientations. The four primes in D pass the first test but fail the second. In group terms they occur in M but not in the centralizer 2.B. This is a source-backed arithmetic/group bridge, not a digit analogy.

Two cautions:

- X_0^+(p) is not X_0(p). The latter has genus zero only for prime levels 2,3,5,7,13. In particular the C-prime 11 is a hostile to identifying C with prime genus-zero levels of X_0(p).
- “New relative to B” does not mean “occurs in no other sporadic group.” The prime 29, for example, also occurs in other sporadic orders. We do not need a classification across all sporadic groups.

## 3. Stronger finite-exact result: exactly one Frobenius pair at each excluded prime

The computation below proves, for the four primes in D, that the rationality failure after adding level two is minimal: exactly two nonrational points, constituting one Frobenius orbit. By Aricheta's formula, all four curves X_0^*(2p) have genus one.

Use the rational coordinate

u=(eta(tau)/eta(2tau))^24, j=(u+256)^3/u^2,
w_2(u)=4096/u, s=u+4096/u.

The eta transformation gives w_2 directly; s is a degree-two quotient coordinate. The normalized 2A Hauptmodul is s+24, but the constant shift does not change splitting degrees. These prime-to-characteristic level-two moduli coordinates are used here away from characteristic 2. We count distinct coarse-moduli points, not stack masses or roots with ramification multiplicity.

| p | Total supersingular j | Total supersingular s | Rational s | Unique irreducible quadratic factor in s |
|---|---:|---:|---:|---|
| 29 | 3 | 4 | 2 | s^2−9s−13 |
| 41 | 4 | 5 | 3 | s^2+8s−8 |
| 59 | 6 | 8 | 6 | s^2−10s−6 |
| 71 | 7 | 10 | 8 | s^2−22s−18 |

All coefficients in a row are interpreted modulo p. The discriminants are respectively 17 mod 29, 14 mod 41, 6 mod 59, and 59 mod 71, all nonsquares. Thus the two roots in each row are genuinely nonrational and exchanged by p-Frobenius.

The independent field audit also identifies the endpoints of these isogenies:

| Characteristic p | Unordered pair of supersingular j-endpoints |
|---|---|
| 29 | {2,25} |
| 41 | {3,32} |
| 59 | {28,48} |
| 71 | {41,66} |

In each row the two distinct nonrational s-points join the same two rational j-points. They are a conjugate pair of parallel edges in the coarse unordered supersingular 2-isogeny multigraph. All eight displayed endpoint j-values avoid 0 and 1728 in their respective characteristics, so exceptional endpoint automorphisms do not complicate this particular pair. Frobenius fixes the vertices and exchanges these two edges. Forgetting the level structure therefore erases the entire rationality obstruction in an explicit, visible way. This is a graph with multiplicities, not a tournament, and no orientation is imposed.

### Certificate A: polynomial elimination

`audit_level_two.py` performs these exact operations:

1. Construct the Legendre Hasse polynomial H_p(lambda)=sum binom((p−1)/2,k)^2 lambda^k modulo p.
2. Eliminate lambda using j=256(lambda^2−lambda+1)^3/[lambda^2(lambda−1)^2], then squarefree and make monic. This gives the full supersingular j-polynomial J_p.
3. Independently count points on one elliptic curve at every j in F_p and verify that precisely the trace-zero j-values are the F_p-roots of J_p. Verify the classical total-degree formula.
4. Pull J_p back along j=(u+256)^3/u^2, then squarefree. Verify that this root set is invariant under u↦4096/u.
5. Eliminate u with u^2−s u+4096=0 and squarefree again. Factor the resulting s-polynomial exactly over F_p.

Universe: every Monster prime p>=5, together with hostile non-Monster primes 37 and 43. There are no probabilistic primality tests or floating-point steps. The p=2 case is excluded by the p∤N hypothesis and the prime-to-characteristic level-two interpretation (not a claim that the abstract coarse genus-zero curve has no good model at 2); p=3 rationality is cited from Aricheta rather than silently counted by a p>3 trace test.

For the nine Baby primes in the computed universe (5,7,11,13,17,19,23,31,47), every s-factor is linear. At each of 29,41,59,71 there is exactly one irreducible quadratic and all other factors are linear. Hostiles 37 and 43 already have nonrational j-values and also have a nonrational s-pair. Hence genus one at level 2p alone does not characterize D among all primes; the Monster/Ogg condition is load-bearing.

### Certificate B: direct finite-field enumeration

`independent_field_audit.py` does not use symbolic resultants. For each of the thirteen Monster primes >=5, it constructs F_{p^2}=F_p[z]/(z^2−d) using a verified quadratic nonresidue d. It enumerates all p^2−1 nonzero u-values, evaluates j(u) and s(u), and retains u precisely when j(u) is a rational supersingular j identified by direct elliptic-curve point counting. The full s-set is then deduplicated. Its size, nonrational count, and Frobenius invariance agree with Certificate A in every case.

This path uses the cited Ogg fact that all supersingular j are rational at the selected primes, also checked by complete factorization in Certificate A. It is an independent implementation of the level-two computation, not a separate proof of Ogg's global classification.

Reproduction:

```
python3 audit_level_two.py > level_two_exact.json
python3 independent_field_audit.py > independent_field_exact.json
```

Independent adversarial audit: a separate implementation enumerated Hasse roots over F_{p^2}, then formed j-, u- and s-sets without importing these scripts or initially assuming rational j. It reproduced all fifteen rows and checked p=3 separately, with no mathematical corrections. Its sole literal-code comment correction (integer resultant followed by reduction, rather than arithmetic directly over F_p(j)) has been applied. The four parallel-edge endpoint pairs were also independently verified, including their generic-automorphism condition. The audit clarified that rationality concerns coarse Fricke orbits, not necessarily an oriented isogeny representative over F_p, and that p=2 is excluded by the theorem’s prime-to-level hypothesis rather than by a claim about every model of the coarse curve; both clarifications are incorporated.

The scripts require Python 3 and SymPy and perform no network requests. The output records all exact factorization data, point counts, group-order arithmetic and nonrational extension-field coordinates.

## 4. Inheritance and discarded shortcuts

Closest proved/cited mechanism: Ogg/Aricheta's Frobenius-rationality versus quotient-genus formula, and the 2A centralizer 2.B. Canonical hostiles: p=11 distinguishes X_0(p) from X_0^+(p); p=37,43 distinguish a genus-one level-two defect from the Monster-only complement. Corrected near miss: the factor 2 between |B| and |2.B| changes exponent dynamics without changing the three support strata. Least-used sidecar: explicit supersingular polynomials on the level-two quotient, now supplied rather than inferred from group orders.

Recovered repository routes:

- `07-reflections/supersingularity-triples.md` is historical inspiration, not a theorem source. Its sentence that supersingularity annihilates all torsion is false. For an elliptic curve the formal-group height is 2, not infinity; geometric p-torsion vanishes, but the p-torsion group scheme is not zero and prime-to-p torsion persists. See the standard height discussion in [Chao Li, Supersingular elliptic curves](https://www.math.columbia.edu/~chaoli/docs/MinorThesis2.html). No proof here uses the historical rhetoric.
- `05-knowledge/results/special_primes_curling_genus_zero_langlands_20261015.md`, §4, already records the Ogg list and a Paley half-imbalance/class-number reformulation. That route detects ordinary level-one rationality, not the difference between M and B.
- `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md`, §§2,6, gives audited Hasse/Legendre/point-count methods and emphasizes the extra level-structure coordinate. MISTAKE-666 in `01-canon/MISTAKES.md` limits each computational range and prevents a false “for every p” extrapolation. Its experiment directory is absent from this checkout; no missing historical script was fetched or implicitly rerun.
- `01-canon/theorems/THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md`, `THM-4759-deep-hole-twining-of-the-leech-voa-is-umbral.md`, and `THM-4760-cusp-criterion-for-fricke-twinings-of-the-leech-voa-and-the-two-descent-family.md` carry genuine lattice/VOA-to-Hauptmodul maps. The refined rationality and cusp conditions are crucial. These theorems do not identify our order-ratio N with a lattice shell, twining coefficient, or tournament invariant.

Live concept board and outcome:

| Concept | Preserved information | Missing coordinate / outcome |
|---|---|---|
| Group-order valuations | Exact A/D/C partition and N | Does not determine multiplication or representation |
| Centralizer and conjugacy class | Actual subgroup 2.B and orbit size | The central factor changes the valuation cycle |
| Ogg/Aricheta rationality | Frobenius action on supersingular points | Level-two isogeny structure is essential |
| Leech/VOA principal moduli | Specific traces, cusps and rationality | No map from N or Hamiltonian gaps supplied |
| Tournament gaps {7,21} | Existing combinatorial objects must be retained | Shared integer labels alone give no transport |

## 5. Relation to the proposed cycle and to {7,21}

The exponent graph of N is exactly 13→2→5→3→7→4 plus 11→1. Identifying the output vertex 4 with the prime vertex 13 gives a legitimate five-state graph quotient, but not integer equality. If the identification is instead read modulo 9, the full exponent map does not descend because 11≡2 while their exponents 1 and 5 disagree. The separate cycle note gives the complete quotient audit and shows that replacing N by N/2 changes this five-cycle to a two-cycle.

There is no source-backed or constructed map here from the forbidden Hamiltonian counts {7,21} to B/M. The viable route is to demand more than scalar agreement: an action, a twined trace, a modular transformation law, a representation-valued statistic, or a level-structure/Frobenius observable. The partition and the four single Frobenius-pair defects are exact positive structure; a claim that they explain the tournament gaps remains OPEN / UNSUPPLIED, not proved by this computation.
