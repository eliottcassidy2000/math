# First-checkpoint scope and provenance

Published commit: `6fff918603cec8a2b9e233bbb627d703214d6953`.
Tree: `26a5b3de8341668797234e051d6016ea388e0b79`.
All 118 intended published blobs were verified. The first checkpoint proved
sharp unrestricted stability through order seven and an all-order coefficient
asymptotic to 4/n, alongside the context/resource/restricted-rigidity results.
It did not include the later sharp-eight, optimal-asymptotic, or all-parameter
sign-chamber continuations.

The original lane selectors, reviewed source pins, [source-integrity record](SOURCE_INTEGRITY.json),
and [portable replay](REPLAY_CHECKPOINT1.json) remain unchanged. Current combined
status is in [README](README.md), [synthesis](SYNTHESIS.md), and
[verification](VERIFICATION.md).

Below is the exact first-checkpoint README, retained as a historical scope
snapshot. Its original SHA256 is `ff2e96bb7a69cc0cb13522fba2df0306c64a9f654847cbc19a1b37c6ecd4f939`.

---

# Session138: stability and actual contexts

**Weighted determinant stability has the correct order in dimension, and the sharp unrestricted coefficients are established through order seven.** Actual graph contexts also give infinite strong fixed-gap pairs, while restricted profile rigidity and explicit resource bounds explain why these constructions do not solve exact-spectrum completeness.

[Results and scope](SYNTHESIS.md) · [Verification](VERIFICATION.md)

## 1. All-order weighted stability

For a weighted tournament matrix A of order `n>=3`, write
`Delta=H-det(I+A)-|det A|`, `V=sum_(i<j) A_ij A_ji`, and let T be its independent binary rounding. For every `1<=k<=n`,

```text
Delta(A) >= E Delta(T) + gamma_(n,k) V,
gamma_(n,k) = 4[n-k+2-(n+k+2)/2^k]/[n(n-1)].
```

Optimizing k gives a lower coefficient asymptotic to `4/n`, improving the inherited `24/n^2` bound in the large-order regime. Ordinal fair-pair examples give an asymptotic upper bound `8/n` on the optimal coefficient. Thus the order is settled, not the all-order optimal constant; dimension-free positive stability is impossible.

The sharper conjectured coefficients are proved on the entire **matching-uncertainty class**, even with arbitrary binary cross-edge orientations. The H=1 locus has a complete ordinal-pair classification and exact ratio minimum.

[Multistar theorem and matching corollary](stability/MULTISTAR_STABILITY.md) · [H=1 classification](stability/H1_STABILITY.md) · [Uncertainty-graph rounding](rounding_graph/PROOF.md)

## 2. Sharp small-order coefficients

The largest universal coefficient b in `Delta>=bV` is:

| Order n | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|
| Sharp b | 1 | 3/4 | 7/8 | 3/4 | 37/48 |

Orders three and four are inherited from Session137; this checkpoint adds sharp orders five through seven, stronger triangle/cycle reserves, and exact joint coefficient regions. At orders four through seven, positive-variance equality is precisely an ordinal sum of the maximum possible number of fair pairs and at most one singleton.

For `tau`, the total directed-triangle weight, the exact order-five region for `Delta>=a*tau+b*V` is
`a<=1`, `b<=7/8`, `a/4+b<=17/16`; at order six it is the rectangle `a<=1`, `b<=3/4`.

[Sharp low orders and joint regions](stability/SHARP_LOW_ORDER_STABILITY.md) · [Sharp seven](stability/SHARP_SEVEN_GRAPH_CERTIFICATE.md) · [Sharp weak six](stability/SHARP_WEAK_SIX.md)

## 3. Actual contexts and restricted rigidity

| Result | Scope | Read |
|---|---|---|
| One exterior arc reversal realizes `6E(U)-H(U)` and `6E(U)+H(U)` | Actual strong contexts for every nonempty U; unbounded fixed gaps, source-certified inverse, and complete minimal-context classification | [Reflector theorem](exact_resources/PROOF.md) |
| Every bounded actual height/interior bank misses a residue at some `2^r-1` | Resource obstruction, not an unrestricted spectral omission; bounded substitution width does not bound these resources | [Budget theorem](exact_resources/RESOURCE_BUDGET.md) |
| Homogeneous-pair strict comparison and exact quotient-derivative kernel | All-order realizability identities; a strengthened formal-network countercontrol still does not represent a tournament | [Actual-kernel theorem](zero_rigidity/README.md) |
| Ordinal triangle completions are fully profile-comparable; cyclic embeddings are zero-rigid under internal star edits | All orders in the stated module class; changes to outside attachments are excluded | [Profile dichotomy](zero_rigidity/PROFILE_DICHOTOMY.md) |

## Boundaries

The unrestricted positive-odd H-spectrum conjecture, global boundedness of exact spectral width, arbitrary strong-quotient scalar zero-rigidity, and LRC(14) remain **OPEN**. The all-order optimal weighted coefficient is not determined. Order-eight sharp stability is a separate continuation, not included in this checkpoint.

Finite graph certificates prove continuous-cube statements only through their explicit polynomial reductions. Same-code replay is not an independent implementation. Exploratory searches and executables are excluded from the release selectors; any retained hostile control is labeled as such. No external-priority or formal proof-assistant claim is made.
