# ABC and inter-universal Teichmuller theory: primary-source and transfer ledger

> **Freshness:** primary pages and PDFs checked 2026-08-23.  Recheck the
> preliminary formalization and discussion pages before making a later status
> claim.  This file records what sources state and how the repo may consume
> them; it does not adjudicate the IUT controversy.

## Status gate

- **ABC: OPEN.**  In this repository, ABC is a conjectural antecedent, never a
  proved dependency.
- **IUT I--IV: PUBLISHED.**  The four papers appeared in the 2021 PRIMS
  special issue.  Publication is a bibliographic fact, not by itself a repo
  proof status.
- **IUT-to-ABC: CITED CLAIM / OPEN DISPUTE.**  IUT IV claims consequences
  including Vojta, ABC, and Szpiro.  Scholze--Stix give a direct mathematical
  objection; Mochizuki rejects that objection.  The repo does not collapse
  this disagreement into either `PROVED` or `REFUTED`.
- **2026 formalization: PRELIMINARY RADAR.**  The author's April 2026 report
  describes skeletal early-stage work, not a completed public machine-checked
  proof.
- **Repo consequence:** [THM-3833](../../01-canon/theorems/THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md)
  proves implications *from an explicitly assumed ABC schema*.  It does not
  use IUT as an antecedent.

## 1. The exact ABC interface

For pairwise-coprime positive integers `A+B=C`, the form used in the repo is

```text
for every epsilon>0 there is K_epsilon>0 such that
C <= K_epsilon rad(ABC)^(1+epsilon).                  (1)
```

In logarithmic language this compares an Archimedean/projective height with
the sum of logarithms of primes in the conductor support.  It is the rational
`P^1` minus `{0,1,infinity}` instance of the wider Vojta framework, up to the
usual normalizations and bounded terms.

The type of every coordinate matters:

| coordinate | retained information | information deliberately absent |
|---|---|---|
| `rad(ABC)` | the set of finite primes occurring | valuation depths, term slots, signs |
| `C` or `max(|A|,|B|,|C|)` | primitive projective height | ambient common scale after normalization |
| `A+B=C` | additive certificate and slot assignment | no automatic multi-term generalization |
| `K_epsilon` | conjectural uniformity at fixed `epsilon` | no known numerical pruning threshold |

The smallest fibre hostile is

```text
1+2=3,          1+8=9,          rad(ABC)=6 in both.   (2)
```

The radical address is equal, while heights and ABC qualities differ.  For
two cubes the entire common-scale family

```text
(2^k)^3+(3*2^k)^3=28*8^k,          rad(value)=14      (3)
```

collapses to the single primitive equation `1+27=28`.

Classical `S`-unit theory already proves finiteness when the allowed prime set
is fixed.  ABC would supply the particularly simple conditional height shape
`C<=K_epsilon R_S^(1+epsilon)`; it is not needed merely to know that a fixed
primitive `S`-fibre is finite, and the unknown `K_epsilon` is not an effective
numerical cutoff here.

## 2. What the IUT papers present

The following is a typed reading of the papers' own constructions and claims.
It is safe as a source map; it is not an endorsement of the disputed proof.

### Hodge theaters and separated arithmetic structures

IUT I begins from elliptic/valuation/theta data and constructs Hodge theaters,
described as miniature models of conventional scheme theory in which additive
and multiplicative/value-group aspects have been disentangled.  The theaters
are isomorphic at a conventional level, but the horizontal theta-link is
explicitly not compatible with their native ring/scheme structures.  The ring
structure viewed across such a link is therefore called *alien*.

Typed transfer:

```text
source:       one fully realized Hodge theater
target:       linked Frobenioid / prime-strip data in another theater
preserved:    specified multiplicative monoids, divisors and Galois actions
not carried:  native addition, field embedding and full scheme structure
sidecar:      theater/species identity and the link path.                 (4)
```

### Frobenius-like, etale-like, and Kummer interfaces

The papers keep multiplicative/divisor data and fundamental/Galois data in
different symmetry regimes, then connect selected pieces through Kummer
theory.  This is deliberate forgetting and reconstruction, not permission to
identify the two regimes wholesale.  An `N`-Kummer class forgets `N`th powers;
local logarithms forget roots of unity.  Any downstream decoder therefore
needs a torsion/power-kernel and basepoint/provenance sidecar.

### Theta values and their evaluator

IUT II treats values of a reciprocal `ell`th-root theta function at torsion
points.  In the paper's normalization they have the shape

```text
qbar^(j^2),             j=0,...,(ell-1)/2,             (5)
```

up to a specified root-of-unity ambiguity.  The label `j` is an address; the
quadratic evaluator `j -> j^2` is load-bearing numerical data.

### Log-links, multiradiality, and indeterminacy

IUT III arranges horizontal theta/LGP links and vertical logarithmic links in
a noncommutative lattice.  It presents *coric* data readable from a common
core, *uniradial* data native to one spoke, and *multiradial* algorithms that
represent spoke-native data in forms meaningful from several alien spokes,
normally modulo stated indeterminacies.

The paper's Theorem 3.11 packages theta/LGP pilot data in common log-shell
containers.  Its indeterminacies include label/procession automorphisms,
automorphisms acting on mono-analytic/unit-mod-torsion factors, and an upper
semi-compatible log-Kummer containment rather than equality.  Passing to a
holomorphic hull and log-volume retains a one-sided containment/volume budget
while destroying exact-image and inverse-recovery data.

### Claimed height chain

IUT III Corollary 3.12 claims a comparison between a theta-pilot log-volume
and a `q`-pilot log-volume.  IUT IV then combines that comparison with local
estimates and a choice of auxiliary prime to claim a Vojta-type height versus
log-different/log-conductor inequality, from which it claims ABC and Szpiro.

For repo purposes, the arrows are typed separately:

```text
IUT III multiradial representation
  -> IUT III Corollary 3.12 log-volume comparison       CITED CLAIM / OPEN DISPUTE
  -> IUT IV height-conductor estimate                  CITED CLAIM / OPEN DISPUTE
  -> Vojta / ABC / Szpiro                              CITED CLAIM / OPEN DISPUTE. (6)
```

No arrow in `(6)` is a dependency of proved repo canon.

## 3. The public disagreement, stated at its mathematical locus

Scholze--Stix argue that after consistently identifying the relevant
one-dimensional real spaces, the abstract theta pilot does not retain the
arithmetic-degree weights `j^2`; omitting those weights makes the target
inequality essentially empty, while inserting them creates an inconsistent
transport/monodromy and an error too large for the height estimate.

Mochizuki replies that such a loop uses identifications the framework
deliberately forbids: distinct labels and alien ring structures must not be
silently identified.  His account instead places the pilots, modulo the stated
indeterminacies, in a common container with a single log-volume comparison.

The precise unresolved question to preserve is therefore:

> Does the common-container/multiradial construction license the membership
> step used in IUT III Corollary 3.12 while retaining the `j^2`-weighted
> arithmetic-degree comparison?

Calling this merely “different universes” loses the actual objection.  Calling
publication a resolution loses the truth-status distinction.  Calling the
critique a formal refutation would also exceed the repo's evidence.

## 4. Triangular ranks: address, fibre, and evaluator

The user's triangular convention and the IUT theta labels meet at an exact
but sharply typed interface.  For `T(z)=z(z+1)/2`,

```text
j^2=T(j)+T(j-1),
(j+1)^2-j^2=2j+1,
T(z+2)-T(z-2)=2(2z+1).                               (7)
```

Thus the quadratic evaluator can be stored absolutely as `j^2` or
procedurally as an ordered word of odd increments.  The last identity is
twice one increment.  If label order is quotiented by permutations, those
increments cannot be integrated without a starting value and synchronization
sidecar.

For a primitive Pythagorean triple in the convention of
[THM-3756](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),

```text
B+C=(2r-1)^2=8T(r-1)+1.                               (8)
```

The outer ordinal `r` is a valid shell address, but its fibre has
`phi(2r-1)/2` nodes.  The complementary odd-root coordinate
`d=sqrt(C-B)=2s-1` restores the node.  For example, shell `r=3` contains both
`(5,12,13)` and `(15,8,17)`.  Conditional on ABC, every primitive
Pythagorean triple would satisfy

```text
rad(ABC) >= K_epsilon^(-1/(1+epsilon))C^(2/(1+epsilon)), (8a)
```

but the two shell-mates have radicals `390` and `510`.  The bound consumes the
whole node, not the bare outer ordinal; it gives no finiteness of the known
infinite tree.

The analogous safe carrier for `(5)` is

```text
(theater identity, j, evaluator j^2,
 valuation/log normalization of qbar, link word, torsion ambiguity).       (9)
```

Treating the `j`th selected value as task `j` is harmless.  Replacing `j^2`
by `j` or by `1` inside a height/log-volume calculation is not.  If
`m=(ell-1)/2`, then

```text
sum_(j=1)^m j^2=ell(ell^2-1)/24,
average_(j=1)^m j^2=ell(ell+1)/12.                    (10)
```

The evaluator changes the scale of the aggregate; it is not cosmetic syntax.

## 5. Connections to proved repo mechanisms

Every row below is a procedural analogy or a one-way consumer contract unless
an explicit theorem link says otherwise.

| repo result | exact common mechanism | decisive loss / hostile | status of connection |
|---|---|---|---|
| THM-3756 odd-square chart | dense shell ordinal plus complementary fibre coordinate | equal `B+C=25` does not identify the node | **PROVED in repo**; IUT comparison only ANALOGY |
| THM-3825 prime-colour decoder | valuation quotient/remainder separates shell and scale | radical sends every state in `(3)` to one address | **PROVED in repo**; ABC budget is CONDITIONAL |
| THM-3656 Frobenian sieve | prime colours and correlated support lines | radical forgets colours/incidence; “Frobenian” is not “Frobenioid” | ANALOGY / NO TRANSFER |
| THM-3824 Rule 30 tariff | quotient is valid for a fixed operation; an adaptive stratum needs extra state | fixed Smith class does not recover nonlinear exact gap | ANALOGY / OPERATION AUDIT |
| THM-893 eroded-diamond twins | identical intrinsic shadow, distinct ambient frames | no intrinsic invariant separates literally equal shadows | ANALOGY / FRAME SIDECAR |
| THM-2022 whole-face Frobenius | a whole packet plus carry sidecar survives a Frobenius operation | atomwise or vocabulary-only Kummer matching is invalid | ANALOGY / NO OBJECT MAP |
| THM-3832 triangular root-ratio chart | a birational coordinate simplifies one operation while exposing denominator divisors | chart scalarization loses global polynomial regularity at its boundary | **PROVED in repo**; NO ABC/IUT TRANSFER |

The most faithful shared design rule is:

```text
ordinal address
  + fibre/provenance sidecar
  + evaluation morphism
  + ordered operation word
  + explicit quotient-loss ledger.                                  (11)
```

This rule helps audit IUT-shaped transfers; it does not validate IUT.

## 6. Connection contracts and decisive tests

### Radical support to a natural task

```text
source:       primitive additive packet (A,B,C; A+B=C)
target:       prime-support set or its natural-number rank
preserves:    which primes occur
destroys:     exponents, ordered term slots, signs, height
sidecars:     valuation vectors, additive witness, Archimedean height
test:         compare the two packets in (2).                          (12)
```

### Theta rank projection

```text
source:       weighted selected family {(j,j^2)}
target:       rank j
preserves:    selection, order, equality of labels
destroys:     numerical weight if evaluator semantics is erased
sidecars:     w(j)=j^2, normalization, theater/path
test:         replace j^2 by j or 1 and compare (10).                  (13)
```

### Alien-copy comparison

```text
source:       two typed native structures R_dagger and R_double-dagger
target:       a common weaker container
preserves:    only data proved invariant under the chosen transporter
destroys:     native cross-copy ring identity and exact-image data
sidecars:     copy ID, transporter, direction, indeterminacy orbit/hull
test:         demand cycle/monodromy coherence before comparing scalars. (14)
```

### One-sided hull tariff

```text
source:       exact local images
target:       hull plus coarse log-volume
preserves:    containment and a correctly oriented inequality
destroys:     equality and inverse recovery
sidecars:     inequality direction, enlargement, normalization
test:         reverse the link or use equality after hull formation.    (15)
```

## 7. Procedural task queue generated by the bridge

1. **Typed two-copy toy model.**  Build two isomorphic finite rings with
   distinct type tags, glue only their multiplicative monoids into a common
   measured container, and attach a `j^2`-weighted pilot.  Test horizontal--
   vertical order, the transporter of `j^2`, and whether hull containment alone
   proves membership.  This isolates a logical pattern; it is not an IUT
   formalization.
2. **ABC radical fibres.**  Extend THM-3833's `19,314`-row atlas by grouping
   primitive triples by radical support and measuring height/quality spread.
   Determine the minimal slotwise valuation sidecar for round-trip decoding.
3. **Pell-ray quality.**  Audit normalized gcd and ABC quality along the four
   positive two-cube Pell rays of THM-3547/3375.  Record which prime-colour
   coordinates, rather than the scalar radical, drive each trend.
4. **Colour plus valuation excess.**  Combine THM-3656's support-colour sieve
   with THM-3833's conditional valuation-excess budget only on actual compiled
   cube pairs; do not apply ABC to local survivors lacking an additive packet.
5. **Operation-lattice scheduler.**  Rank typed states `(copy, structure
   level, pilot, operation word)` by shortlex while retaining the noncommuting
   word.  Compare `HV` with `VH`; erasing the word is the hostile control.
6. **Transfer firewall.**  Do not seek an LRC or Jacobian consequence until a
   map preserves owner/phase/arrival or Keller-equivalence data, respectively.

## 8. Primary sources

- Shinichi Mochizuki,
  [*Inter-universal Teichmuller Theory I*](https://www.kurims.kyoto-u.ac.jp/~motizuki/Inter-universal%20Teichmuller%20Theory%20I.pdf),
  official author PDF; Hodge theaters, theta-link, alien arithmetic structure.
- Mochizuki,
  [*Inter-universal Teichmuller Theory II*](https://www.kurims.kyoto-u.ac.jp/~motizuki/Inter-universal%20Teichmuller%20Theory%20II.pdf),
  official author PDF; theta values, Kummer and multiradial language.
- Mochizuki,
  [*Inter-universal Teichmuller Theory III*](https://www.kurims.kyoto-u.ac.jp/~motizuki/Inter-universal%20Teichmuller%20Theory%20III.pdf),
  official author PDF; log-theta lattice, Theorem 3.11 and Corollary 3.12.
- Mochizuki,
  [*Inter-universal Teichmuller Theory IV*](https://www.kurims.kyoto-u.ac.jp/~motizuki/Inter-universal%20Teichmuller%20Theory%20IV.pdf),
  official author PDF; claimed height and Vojta/ABC/Szpiro consequences.
- [PRIMS volume 57 special issue](https://ems.press/journals/prims/issues/1507),
  EMS Press; publication record for IUT I--IV, 2021.
- Mochizuki,
  [*Alien Copies, Gaussians, and Inter-universal Teichmuller Theory*](https://www.kurims.kyoto-u.ac.jp/preprint/file/RIMS1854.pdf),
  official RIMS preprint; alien-copy and weaker-structure transport account.
- Peter Scholze and Jakob Stix,
  [*Why abc is still a conjecture*](https://www.math.uni-bonn.de/people/scholze/WhyABCisStillaConjecture.pdf),
  authors' PDF; ABC formulation and direct critique of the IUT III step.
- Jan-Hendrik Evertse,
  [*On equations in S-units and the Thue--Mahler equation*](https://eudml.org/doc/143111),
  *Inventiones Mathematicae* 75 (1984); classical fixed-`S` finiteness.
- Mochizuki,
  [May 2018 comments on the discussion](https://www.kurims.kyoto-u.ac.jp/~motizuki/Cmt2018-05.pdf)
  and the [official discussion portal](https://www.kurims.kyoto-u.ac.jp/~motizuki/IUTch-discussions-2018-03.html);
  author's rejection of the critique.
- Mochizuki,
  [*Formalization of IUT* (April 2026 preliminary report)](https://www.kurims.kyoto-u.ac.jp/~motizuki/Formalization%20of%20IUT%20%282026-04%29.pdf);
  preliminary radar, not a completed public formal proof.

## 9. Does not prove

This intake does not prove ABC, any IUT theorem, the contested
`IUT III 3.11 -> 3.12` implication, Fermat--Catalan unconditionally, LRC(14),
`JC(2)`, a two-cube support asymptotic, or a common theorem merely from shared
words such as *Kummer*, *Frobenius*, *height*, *radial*, *colour*, or *copy*.
