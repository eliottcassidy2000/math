# Second-checkpoint scope and provenance

Published commit: `a60f529a50744a6a95a3d327983f9542780894d2`.
Tree: `4f42622a7f7f82a5642c98fe2b1e26b448c67756`.
All 156 cumulative intended Session138 published blobs were verified.
This checkpoint proved sharp order eight, the uniform optimal asymptotic
coefficient `8/n+O(n^(-2))`, and the all-parameter mixed-profile sign chambers.
It did not include the later exact theorem for n>=1500.

Its curated lane bytes, selectors, source pins, [source-integrity record](SOURCE_INTEGRITY_CHECKPOINT2.json),
and [portable replay](REPLAY_CHECKPOINT2.json) remain unchanged.
The original README follows verbatim; its SHA256 is
`e6e0179f5df1fe2d4d7812ff6937d963e5bc0e76ec199852e478f4719ad13845`.

---

# Session138: stability and actual contexts

**The optimal weighted determinant stability coefficient is `8/n+O(n^(-2))`; exact sharp coefficients are proved through order eight.** The continuation also closes an explicit unbounded mixed-profile construction family without solving general zero-rigidity.

[Results and scope](SYNTHESIS.md) · [Verification and checkpoint history](VERIFICATION.md) · [Continuation replay: PASS](REPLAY_CHECKPOINT2.json)

## 1. The optimal asymptotic stability constant

**PROVED, independently audited.** For weighted tournament matrices of order n, let `Delta=H-det(I+A)-|det A|`, `V=sum_(i<j) A_ij A_ji`, and

```text
K_n = inf_(V>0) Delta/V.
```

Write `r=floor(n/2)` and define the fair-pair coefficients

```text
kappa_(2r)   = (4/r)[1-(3/4)^r-(1/4)^r],
kappa_(2r+1) = (4/r)[1-(3/4)^r].
```

The theorem gives, for `n>=4` and `c0=(1-exp(-4))/4`,

```text
kappa_n/[1+n*kappa_n^2/(2*c0)] <= K_n <= kappa_n,
K_n = 8/n + O(n^(-2)).
```

Its near-candidate consequence is also quantitative: if `Delta/V<=kappa_n`, independent binary rounding T has expected gap at most `2n*kappa_n^2`; a greedy matching leaves at most `2n*kappa_n` variance outside it; and `Pr(T nontransitive)<=2n*kappa_n^2`. These are probability and weighted-mass statements, not exact fair-pair rigidity.

[Asymptotic theorem](asymptotic_stability/PROOF.md) · [Independent audit](independent/asymptotic_audit/AUDIT.md) · [Separate resolvent mechanism](stability/RESOLVENT_ROBUSTNESS_CROSSCHECK.md)

## 2. Exact small-order constants and parity

The largest universal b in `Delta>=bV` is:

| Order n | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| Sharp b | 1 | 3/4 | 7/8 | 3/4 | 37/48 | 87/128 |

At orders four through eight, positive-variance equality is exactly an ordinal sum of the maximum number of fair pairs and at most one singleton. The parity effect is structural: a singleton makes the full determinant vanish; an even perfect pair cover contributes the extra subtraction `(1/4)^r`. Exact equality `K_n=kappa_n` beyond order eight remains **OPEN**.

The order-eight theorem retains positive triangle and cycle reserves. Its complete graph certificate is checked by independent orbit construction and subset dynamic programming, then extended to the full real variance cube by the reviewed multiaffine reduction.

[Sharp eight](stability/SHARP_EIGHT_GRAPH_CERTIFICATE.md) · [Sharp low orders and joint regions](stability/SHARP_LOW_ORDER_STABILITY.md) · [Sharp seven](stability/SHARP_SEVEN_GRAPH_CERTIFICATE.md)

## 3. All-order mechanisms retained from checkpoint one

- Multistar rounding gives `Delta(A)>=E Delta(T)+gamma_(n,k)V`, with optimized coefficient asymptotic to `4/n`. Its original scope remains recorded separately from the later asymptotic continuation.
- The exact matching-uncertainty coefficients hold with arbitrary binary cross-edge orientations. Uncertainty-graph Pfaffian parity gives edge-order certificates and a vertex-cover-parameterized evaluator.
- H=1 is exactly the ordinal-pair/singleton locus. Low-order triangle/variance coefficient regions and stronger weak-gap reserves are exact.

[Multistar and matching support](stability/MULTISTAR_STABILITY.md) · [Uncertainty-graph framework](rounding_graph/PROOF.md) · [H=1 classification](stability/H1_STABILITY.md)

## 4. Actual contexts: constructions and precise barriers

The literal one-edge reflector realizes strong tournament counts `6E(U)-H(U)` and `6E(U)+H(U)` for every nonempty U. Its full-profile operators, source-certified inverse, minimal exterior budget, and unbounded fixed-gap families are proved. Bounded actual height/interior source banks fail even profinite universality; this is not an unrestricted spectral omission.

The ordinal-profile theorem proves internal star-edit zero-rigidity in `C3[TT_k[C3],A,B]` when outside attachments stay fixed. The later sign theorem handles three explicit actual strong-Q4 mixed-profile seeds for every positive side multiplicity: two are positive on the thin boundary and negative throughout the thick quadrant, never zero; the third is positive for arbitrary nonempty side tournaments, with sharp lower gap `135618 H(A)H(B)`.

[Reflector](exact_resources/PROOF.md) · [Resource bounds](exact_resources/RESOURCE_BUDGET.md) · [Restricted profile rigidity](zero_rigidity/PROFILE_DICHOTOMY.md) · [All-parameter sign classification](zero_balance_later/NO_ZERO_THEOREM.md)

## Boundaries

The unrestricted positive-odd H-spectrum conjecture, global boundedness of exact spectral width, arbitrary strong-quotient scalar zero-rigidity, and LRC(14) remain **OPEN**. None of the family-specific no-zero results is a theorem for arbitrary strong quotients. Exact optimal stability coefficients beyond eight are not claimed.

[Checkpoint-one scope](CHECKPOINT1.md), manuscripts, selectors, audit pins, and [portable replay](REPLAY_CHECKPOINT1.json) retain their original through-seven scope. Later results have separate source identities and audit scopes. No finite scan is extrapolated to an all-parameter theorem; no external-priority or formal proof-assistant claim is made.
