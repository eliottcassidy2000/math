# Three source tests force a rational Hamiltonian back into the polynomial carrier

2026-10-11. **PROVED:** the squarefree-pole lemma and the three-test
criterion below. **CITED / INHERITED:** the earlier completed-flow and
nonrational-time theorems. **FINITE-EXACT:** the companion controls.
**OPEN:** fixed-source rational repairs, compositions with changing
invariants, their finite compatibility, and planar JC. This is a scoped
self-audited research note, not an independent audit or a novelty claim.

## Inheritance and the corrected target

The [September7 universal-genus theorem](planar_jc48_sep07_universal_genus.md)
and [September8 cusp-ideal theorem](planar_jc48_sep08_cusp_ideal.md) already
exclude every nonzero rational scalar time for the polynomial carriers
they treat. Their proofs and the universal theorem's independent audit
were read before this step. The present question is whether **rational
generators** enlarge the universally source-preserving infinitesimal class.

The closest mechanism is the exact polynomial classification in Section4
of [Hamiltonian source preservation](planar_jc48_sep06_hamiltonian.md).
The hostile is the genuine rational flow of x^2 t, which preserves one
source but leaves the polynomial p,y algebra. The corrected near miss is
assuming that changing the word "polynomial" to "rational" automatically
escapes the old carrier while preserving its consumer. The least-used
sidecar is the order of a pole before and after differentiation.

Board: **source quantifier / coordinate algebra / pole order / squarefree
divisor / invariant / terminal compatibility**. The connection sends three
regular source responses to regularity of two coordinate derivatives.
It preserves the full universal source requirement, loses no poles, and
needs characteristic zero and a squarefree Poisson numerator. The cheapest
hostile replaces that numerator by a square.

## P1. A squarefree numerator cannot hide a rational Hamiltonian pole

Let A=K[p,y], char K=0, and let a Poisson bracket on its fraction field
have {p,y}=F/G, with coprime nonzero F,G in A. Suppose F is squarefree.
Then for S in K(p,y),

\[
             \{p,S\},\{y,S\}\in A\quad\Longrightarrow\quad S\in A.
                                                               \tag{1}
\]

**Proof.** If S has a pole of order m>=1 at an irreducible h, write
S=h^{-m}a, with a a unit in the local ring A_(h). At least one coordinate
partial derivative h_z is nonzero modulo h: in characteristic zero a
nonconstant polynomial has a nonzero derivative of smaller degree.
Differentiating gives

\[
  \partial_z S=-m a h_z h^{-m-1}+h^{-m}\partial_z a.
\]

The first term has exact h-order -m-1, and the second has order at least
-m. Hence this derivative has a pole of order m+1. Multiplication by F/G
can cancel at most one order because F is squarefree. The corresponding
Hamiltonian coordinate derivative still has a pole, contradicting (1).
Thus S has no prime denominator in the UFD A and is polynomial. QED.

Both hypotheses matter. For {p,y}=p^2 and S=1/p one gets
{p,S}=0, {y,S}=1, both polynomial. In characteristic ell>0, the
nonpolynomial S=p^{-ell} has zero coordinate derivatives even for the
ordinary bracket. The lemma is about rational functions, not multivalued
logarithmic Hamiltonians, whose derivatives may have simple poles.

## P2. Three observables characterize the rational universal carrier

Keep the actual source convention

\[
 p=t(1+x^2t),\quad y=xtp,\quad D=p^3-y^2,\quad
 u=y^2/D,\qquad\{p,y\}=-D/p.
\]

Put A=K[p,y] and G_H=-u/2+H for H in A. For **any rational** S in K(p,y),
the following are equivalent:

1. {G_H,S} is in A for every H in A.
2. {G_0,S}, {G_p,S}, and {G_y,S} are in A.
3. S is in K+p^2 D A.

**Proof.** The first implication is immediate. Subtract the H=0 condition
from those at p and y. Both {p,S} and {y,S} are polynomial. Since D is
squarefree and coprime to p, P1 gives S in A. The polynomial calculation
can now be used legitimately:

\[
 \{p,S\}=-(D/p)S_y,\qquad \{y,S\}=(D/p)S_p.
\]

Their regularity implies p divides S_p and S_y. Expanding S in powers of
p, these conditions say that the p^0 coefficient is constant and the
p^1 coefficient vanishes. Thus S=c+p^2R.

For E=2p partial_p+3y partial_y, the remaining source test is

\[
                  \{G_0,S\}=-\frac{py}{2D}E(S).         \tag{2}
\]

Because D is coprime to py, regularity forces D|E(S). The Euler derivation
on A/(D)=K[z^2,z^3] has kernel exactly K: every positive homogeneous
degree is nonzero in characteristic zero. Thus S=c'+DR'. Evaluation at
the origin identifies c'=c. The relatively prime divisors p^2 and D
then give S-c in p^2DA.

Conversely, for S=c+p^2DR, both S_p and S_y are divisible by p, while
E(S)=p^2D(10R+2pR_p+3yR_y). All three responses are polynomial.
The Leibniz rule makes {H,S} polynomial for every polynomial H, proving
the full universal condition. QED.

The **new domain extension** here is from S in A in the inherited
classification to S in Frac(A). The same three tests were already used
in that polynomial proof. No additional coordinate-degree or denominator
bound is assumed.

Every nonconstant S passing these tests is consequently in the already
settled completed carrier. Its nonzero scalar-time flow has no pair of
rational source images by the September7 theorem. Merely permitting
rational S therefore supplies no new universally preserving single flow.

## P3. The sharp surviving choices are narrower requirements

Coordinate preservation alone gives the slightly larger exact class

\[
       \{p,S\},\{y,S\}\in A\quad\Longleftrightarrow\quad S\in K+p^2A.
                                                               \tag{3}
\]

For example S=p^2 passes (3), but {G_0,S}=-2p^3y/D has a genuine cusp
pole. Conversely S=D passes the G_0 and G_y tests but fails G_p, and
S=pD passes G_0 and G_p but fails G_y. These are exact witnesses to the
distinct regularity obligations; they are not a minimality assertion
about every possible three-observable presentation.

The inherited outside-carrier example S=u=y^2/D=x^2t has the rational
flow

\[
 x\longmapsto x/(1-\lambda x),\qquad
 t\longmapsto t(1-\lambda x)^2.
\]

It fixes G_0=-u/2. But {p,u}=-2yp^2/D and {y,u}=-3py^2/D are not in A.
Thus a **particular-source** condition can allow rational motion that
the universal condition prohibits. For the actual nonconstant H of the
finite-response problem, one must solve its own {G_H,S} regularity and
source/depth equations; the H=0 example does not do so.

This is the next productive branch: preserve the particular required
source with its labelled companion, or change the coordinate algebra and
explicitly pay descent back to the original one. A finite composition may
also restore regularity after intermediate losses, but its individual
three-test failures and the exact final return must be recorded. The
single-flow genus theorem does not exclude such a composition.

There is a useful comparison to the
[Collatz rank interface](collatz_atomic_rank_interface_20261011.md).
A common rank controls every permitted controller action; a stopped state
still owes a ROOT bridge. Here a common source-algebra condition controls
every H; weakening it leaves a specific regularity bridge unpaid. This
is a comparison of quantified interfaces, not a map from curve genus to
Collatz descent. In both cases the proposed weaker construction must
retain the exact object on which its final consumer operates.

## Controls and reproduction

[Program](../../04-computation/experiments/planar_rational_carrier_saturation_20261011.py)
and [output](planar_rational_carrier_saturation_20261011.json) use exact
symbolic arithmetic, pole-bearing hostile denominators, both source
charts, and positive polynomial carriers. The general denominator
argument is the proof of P1, not an inference from the finite sample.
Normal and optimized runs agree on 239 always-active exact gates, including
18 independent literal source-chart comparisons. No independent audit is
claimed for this new packet.

```text
python3 -B 04-computation/experiments/planar_rational_carrier_saturation_20261011.py
python3 -B -O 04-computation/experiments/planar_rational_carrier_saturation_20261011.py
```
