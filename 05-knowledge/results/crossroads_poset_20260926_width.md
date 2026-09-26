# Width, insertion laws, and the arithmetic height cut

**Status:** the attached Kahn–Saks manuscript remains **CLAIM UNDER REVIEW**.
Its main finite proof has been read in detail; this audit found no concrete
internal gap. This is not certification of its cited dependencies or its
asserted quantitative rate. Sections 4–6 below are elementary **PROVED**
distributional obstructions/recovery statements, with **FINITE-EXACT** controls.
No Collatz convergence result is claimed.

Source reviewed: Max Aires, *Proof of the Kahn–Saks Conjecture*, September 25,
2026, 25 pages, user-supplied `C:/Users/Eliott/Downloads/Kahn_Saks_Proof.pdf`.
SHA256: `83bc04b19ca1790913f16183322764a5537ea1aa1101d1e4b8a0cbbd43681897`.
Critical formulas on PDF pages 19 and 22–24 were rendered and visually checked.

## 1. Inheritance and the live board

The closest repository mechanism is the exact finite parity-cylinder law and
the growing-word count in
[THM-4495, no-descent count](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md).
The canonical hostile is the all-odd infinite address: every finite prefix is
realized by positive integers, while its compatible dyadic limit is -1.
[THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md)
is the corrected near miss: arbitrarily long strict growth and retained prime
motifs are actual positive-integer statements, but the source changes with the
requested length. The little-used sidecar here is the *probability law after
the source-height cut*, not merely its set of finite orders.

The working board has five objects: uniform linear extensions; the order
polytope; uniform insertion into a chain; positive-slope parity cylinders; and
height-selected positive starts. Their exact intersections, rather than a
forced tournament, produce the useful obstruction below. The anchor is the
attached proof audit; the niche is conditioning; the wildcard is a width-two
insertion gadget inside the six-step Collatz language.

## 2. What the manuscript actually claims

For a finite poset P, let L(P;Z) be the restriction of a **uniform random
linear extension of the entire P** to a selected set Z. This need not be the
uniform linear-extension law of the induced subposet on Z. Already a<c with
b isolated has P(a before b)=2/3, although {a,b} is an antichain.

Theorem 2 claims that large enough width, with k,t,epsilon fixed, produces
either a k-set whose restriction law is within epsilon of all k! orders
uniformly, or a (t+1)-set whose law is within epsilon of uniform insertion of
one vertex x into one fixed order y1,...,yt. These are genuinely different
alternatives. An isolated vertex and a t-chain have exactly the latter law;
its total variation distance from all (t+1)! orders uniformly is

    1 - 1/t!.

Its comparisons are P(x before yi)=i/(t+1). For k=2,t=1 the alternatives
coincide, giving the claimed width-to-balance theorem. For larger parameters,
nearly uniform insertion does not supply a nearly uniform permutation law.
The manuscript's recursive counterexamples in Section 5 address precisely
this distinction. No theorem about width can be inferred just from an
exponentially large number of linear extensions.

## 3. Main proof audit: finite coordinates and conditioning

The following checks concern the submitted argument, conditional on its cited
inputs. The statements of those inputs were inspected in the manuscript;
their full original proofs were not independently re-established here.

* **Order-polytope coordinates (pages 4–5):** F=U_(f(x)), where f is a uniform
  extension and the uniform order statistics U are independent of f. The rank
  scale Z=(n+1)F is essential. The covariance, window, and height formulas do
  not concern an arbitrary weighted distribution on extensions.
* **Selection and doubling (pages 5–6):** the finite partition proof of
  Theorem 8 first controls collisions, then invokes the fixed-k finite-valued
  Komlos selection result. Lemma 9 uses pairwise small-ball estimates to
  obtain a uniformly bounded covering number. Subtracting a common random
  coordinate preserves all orders; no independence of coordinates is used.
* **Antichain/gap reduction (pages 12–14):** Lemma 11 applies the window bound
  to B minus its large-window subset, not just to B. This justifies the
  asserted at-least-half selection. Lemma 12's summed beta densities have
  bound n/(n+1) in rank units on each half-line. Lemma 13 adds comparisons on
  the same ground set, conditions on their joint event E of probability
  at least 1-M^-4, and obtains exactly the uniform order-polytope law of the
  enlarged poset. Its conditional mean estimate retains the same n+1 scale.
  The selected antichain of the enlarged poset is also an antichain of P.
* **Padding and local counts (pages 15–16):** adding m^2-chain padding leaves
  the old relative-order law unchanged. The conditional beta variance
  increases with n. The maximizing pair has a vertex in the old component,
  while R<=m+1; this places the chosen anchor at least m^2 ranks from the
  endpoints. Exchangeability of uniform spacings then gives the truncated
  binomial count used in (35). The finite grid handles arbitrary intervals
  and endpoint conventions. This is a finite estimate before taking R large.
* **Weighted means and the small-ball step (pages 17–20):** the fourth-moment
  tail bound and the doubling net give a bounded-size probabilistic cover,
  uniform in P. Prekopa is applied on the actual affine support of W, so a
  singular ambient distribution is allowed. Jensen is applied to log Phi
  under the truncated g-weighted probability law. The shift from its mean v
  to u_g has norm less than r/4. In (39), a witness w may depend on g: this is
  legitimate because the conclusion is the deterministic count estimate for
  u_g, and the constants are uniform. No one witness for all g is claimed.
* **Transport and finite coin detection (pages 20–23):** the endpoint strips
  in (44) account for entering and leaving trajectories and atoms at interval
  endpoints. The backward error recursion fixes the windows and tolerances
  before increasing R. The number q of tests may depend on P, but q<=Q, so
  the finitely many q-dependent thresholds can be bounded uniformly. The
  auxiliary coin flips detect deterministic measurement vectors; they do not
  replace the original extension by independent extensions for comparisons.
* **Final insertion coupling (pages 23–24):** all errors e_i and the insertion
  variable V are taken from the same W. Their dependence is explicitly
  allowed. A density bound separates V from its fixed quantiles, Chebyshev
  controls the errors, and their union bound gives TV<=2t eta+2Kt sqrt(eta).
  Distinct quantiles give distinct selected vertices. Universal padding
  vertices are excluded by epsilon<1/[2(t+1)], as in the earlier reduction.

The externally cited load-bearing inputs include Aires–Kahn, *Balancing
extensions in posets of large width*, arXiv:2509.11549, Theorem 8.1; Haqi,
*On the Gap of Finite Posets*, arXiv:2608.12678, supplying (11); continuous
Shepp XYZ; and logconcavity/Prekopa facts. Root handles the independent public
source audit. The quintuple-log quantitative bound displayed as (49) is only
asserted in a remark in this manuscript. It was **not validated** here.

The exact finite census checks the stated XYZ covariance, window-variance,
ideal-height, and antichain-window inequalities through six vertices. That
does not verify a theorem with unspecified large-width thresholds.

## 4. The height cut can destroy the poset law

Use the shortcut T(n)=n/2 for even n and T(n)=(3n+1)/2 for odd n. A parity
word records odd steps as 1. The root poset-bridge lane gives the general
fixed-endpoint correspondence with two chains; only the following explicit
instance is needed for this obstruction.

The six-letter words with five ones and 3^(number of ones in each prefix)
greater than 2^(prefix length) are exactly:

| word | least start modulo 64 | position of B1 among A1<...<A5 |
|---|---:|---|
| 110111 | 27 | after A2 |
| 111011 | 39 | after A3 |
| 111101 | 47 | after A4 |
| 111110 | 31 | after A5 |

They are the four linear extensions of the poset A1<...<A5, A2<B1. Uniform
sampling of the four residue cylinders is exactly its uniform extension law.
Uniform sampling of positive starts n<=31, conditioned on this endpoint and
prefix property, selects only starts 27 and 31, hence only the extreme two
insertion positions, each with probability 1/2.

**Proved obstruction.** No poset on these same six event vertices has exactly
those two linear extensions. Every relation of any such poset must be true
in both orders. Their common comparisons generate precisely A1<...<A5 and
A2<B1, which also permit both missing intermediate orders. This does not
exclude representations with extra auxiliary vertices or nonuniform weights.

There is an analytic obstruction as well. Independently of the selected rank
order, draw six independent Uniform(0,1) variables, sort them into the
dependent order statistics U_(1),...,U_(6), and put Z=7U.
For fixed ranks i,j,

    E[Z_i Z_j] = 7 min(i,j)(max(i,j)+1)/8.

Averaging over the two selected orders gives exactly

    Cov(Z_A2-Z_A3, Z_B1-Z_A2) = 1/4 > 0.

In original unit coordinates the covariance is 1/196. Thus the submitted
paper's nonpositive XYZ inequality fails for this height-selected law. The
full four-order positive control has the same covariance -5/32 in rank units.
This is a failure of a proposed transfer, not a counterexample to Shepp or
to the manuscript: the selected law is outside the stated hypotheses.

The support is also nonconvex. In coordinate order (A1,A2,A3,A4,A5,B1), the
two points

    p=(.1,.2,.4,.5,.6,.3), q=(.1,.2,.3,.4,.5,.7)

lie strictly in the selected two chambers, but their midpoint has order
A1,A2,A3,A4,B1,A5, a missing chamber. Hence uniform measure on the selected
union is not logconcave. Both the combinatorial and analytic structures have
been lost by the actual height cut.

Source/target map: parity letters label successive A and B events; the target
is their total order. The map preserves the fixed-length parity word, its
endpoint odd count, and its prefix slope signs. If the word and arithmetic
carry are retained, its least residue can be reconstructed. Forgetting that
carry and height selection loses the actual-source law. The needed sidecar
is a proved estimate on the selected extension weights, not an added arrow
chosen merely to make a tournament.

## 5. Finite uniform laws need not form an infinite uniform law

Let W_L be all L-letter words with strictly positive logarithmic slope at
every nonempty prefix. Under the uniform law on W_3, prefix 110 has probability
1/2 because W_3={110,111}. Under the uniform law on W_4, its probability is
1/3 because W_4={1101,1110,1111}. Their projected laws have total variation
distance 1/6.

Thus the finite laws conditioned on survival to L are not projectively
consistent. Conditioning a longer prefix reweights each shorter word by its
number of allowed continuations. This does not prevent subsequential weak
limits or Doob-type reweighting constructions; it prevents identifying these
particular finite uniform laws as the marginals of one process without a new
argument. Even a correctly constructed measure on infinite dyadic addresses
would still need a separate positive-integer realization argument.

The manuscript avoids this mistake: its crucial comparisons live inside one
finite poset law, its changed-poset conditioning is explicit, and its final
coupling uses one joint random vector. Applying only its endpoint balance
conclusion to a sequence of newly sampled Collatz prefixes loses that feature.

## 6. Strongest elementary survivor: many complete periods restore uniformity

Let a nonempty set C consist of m distinct residue classes modulo N. Sample
uniformly from positive integers n<=X whose residues lie in C, assuming this
set is nonempty. Write X=qN+s, 0<=s<N, and let t be the number of selected
residues represented in {1,...,s}; residue zero has q positive representatives.
Each selected class has q or q+1 representatives, so its total variation
distance from the uniform law on C is exactly

    TV = t(m-t) / [m(qm+t)].

Indeed each of the t larger weights differs from 1/m by
(m-t)/[m(qm+t)], and each other weight differs by t/[m(qm+t)]. For q>=1,

    TV <= 1/(4q).

This applies unchanged to any fixed-endpoint parity-word class, since the
finite parity-to-residue map is bijective. Thus 2^L=o(X) suffices for uniform
order-law recovery, uniformly over the chosen nonempty word class. At X a
multiple of 2^L the law is exactly uniform. For the hostile above q=0,m=4,t=2,
the distance is 1/2. When L is comparable to log_2 X, q need not grow and the
estimate does not remove the obstruction. It gives no control for one fixed
integer as L tends to infinity.

For a given pair, a balance lower bound b under the uniform law therefore
transfers as b-TV under the height-selected law, since each comparison is an
event. This retains an actual quantitative use of a uniform-order theorem
when its other hypotheses hold; balance still does not itself force descent
of a specified integer.

This is the precise stopping boundary: a finite uniform-extension theorem may
be transferred in a regime with many complete residue periods and a proved
TV error. The same transfer at source-sensitive horizons needs a new weighted
or arithmetic theorem. The explicit positive XYZ covariance rules out simply
reusing the entire order-polytope package unchanged.

## 7. Reproduction and scope

Run from the session worktree:

    python -X utf8 04-computation/experiments/crossroads_poset_20260926_width.py
    python -O -X utf8 04-computation/experiments/crossroads_poset_20260926_width.py

Script: [exact probe](../../04-computation/experiments/crossroads_poset_20260926_width.py).
Output: [exact stdout](crossroads_poset_20260926_width.out).

The census contains every distinct transitive relation on n<=6 vertices whose
comparisons agree with the fixed natural topological order: 5,231 posets and
252,266 listed extensions. Every finite poset type of those sizes has such a
labelling, but this is not a census of every labelled poset. It checks 601,302
ordered XYZ triples, 76,193 pair-window bounds, and 96,427 each of the ideal
and antichain bounds. Positive controls include exact insertion gadgets;
the nonuniform antichain restriction and the height cut are hostile controls.

All 8,190 words through length 12 have their affine-carry residue checked
against a direct integer trajectory. The height-selection obstruction is
first in the explicitly bounded universe L<=6, all endpoint counts, and
initial cuts of least residues (eight nontrivial cuts); no broader minimality
claim is made. The many-period formula is independently checked against direct
integer sampling in 151 cases. Ordinary and optimized runs have identical
stdout; every check uses an explicit exception rather than `assert`.
