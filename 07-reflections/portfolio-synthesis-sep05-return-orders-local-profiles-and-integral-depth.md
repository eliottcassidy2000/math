# Cross-thread synthesis: return order, cycle profiles, and integral depth

**Research date: 2026-09-05. Current outcomes are routed to canonical
THM-4416--4419 below.** This is a synthesis and provenance document; the
linked audited theorems control mathematical status. External priority is
not claimed. This session mapped the maintained portfolio and followed
selected proof dependencies, rather than independently verifying all 3,695
theorem files present at startup.

## Results and their actual scope

| Lane | Audited outcome | Remaining boundary |
|---|---|---|
| Signed-cycle spectrum | [THM-4416, D5/D6 gap](../01-canon/theorems/THM-4416-even-graph-cumulative-d5-d6-spectral-gap.md): the single-edge minimum and equality classification hold at every admissible order for both D=5 and D=6 | D>=7 and unrelated tournament H>=disc |
| Laurent first return | [THM-4417, width-two parabolic bound](../01-canon/theorems/THM-4417-width-two-laurent-first-return-parabolic-critical-bound.md): first nonzero CT(f^m) occurs by M+N if min(M,N)=2; a distinct-critical-point count is sharper | General min(M,N)>=3; no uniform final Gaussian moment bound |
| LRC local triple | [THM-4418, sharp pair arithmetic and 44/13 tail](../01-canon/theorems/THM-4418-lrc14-sharp-pair-arithmetic-and-forty-four-thirteenths-tail.md): T_ab<=12/77 and N_ab<=24b/13 sharply; every odd ternary-unit a<b<c with c/b>=44/13 has strict certificate below 6/77 | Comparable-speed triples at unbounded height, physical entry, synchronization, LRC(14) |
| Integral observations | [THM-4419, precision and dyadic triples](../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md): exact precision/kernel/tensor consequences and all-depth three-node dyadic Smith formula | General multiscale clusters, higher jets, moving source modules |

The complete-residue two-jet formula and the all-height separated-contact
identity were independently obtained by a concurrent session before this
session's checkpoint. They are one shared closure each, credited in the
theorems. The width-two arbitrary-support theorem, D5/D6 closure, sharp LRC
pair bounds/tail, and dyadic-triple/precision additions are distinct from
that incoming work. The smaller trinomial census is corroboration only.

## 1. What the repository's strongest mechanisms retain

The maintained canon records NC2/GMC(2) through
[THM-2022, lowest balanced face](../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md).
Its mechanism is a nonzero seed, a selected lowest face, and good-prime
preservation of that whole face. Scalar noncancellation at a later power is
not obtained by separately declaring all underlying atoms positive. The new
return-order theorem strengthens its seed where the Laurent face has width
two on one side; it does not replace the subsequent coefficient-dependent
prime choice.

For LRC(14),
[THM-4409, third-sheet network](../01-canon/theorems/THM-4409-lrc14-third-sheet-component-network-certificate.md)
restores physical contact incidence lost by a pair-only quotient.
The [concurrent raw-carrier synthesis](../05-knowledge/results/synthesis_20260905_lrc_sparse_transport.md)
then identifies the remaining difference as pointwise pair minimization
versus selecting one pair after summation. Our uniform tail proves one
fixed pair suffices in an infinite region. It does not erase the remaining
selector issue on comparable speeds.

The signed-cycle problem has a different native object: a switching class
of a signed complete graph. Its cycle parity genuinely defines the
observable; there is no reason to replace it by a tournament. The old
[THM-4083, D3/D4 deletion proof](../01-canon/theorems/THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md)
already had the right induction, but the needed joint cycle profile was
available in a separate finite Fourier atlas. Reconnecting those objects
closed the next two entire diameter layers.

For integral observations,
[THM-4080, one-scale two-jet partition](../01-canon/theorems/THM-4080-confluent-two-jet-single-scale-smith-partition.md)
uses weighted minors. The missing prime-boundary direction was invisible
modulo p, yet its integral derivative trace divided by p remains a unit
direction. Both concurrent derivations found this one-step saturation.
The deeper dyadic triple shows that adding a node can change existing Smith
factors, so the larger-cluster problem cannot be handled by appending an
unchanged prefix.

Other active routes fit the portfolio without becoming consequences of
these results. The [incoming formal triple continuation](../05-knowledge/results/synthesis_20260905_transgression.md)
distinguishes tangent compatibility, formal integration, and polynomial
termination in a planar-JC chart. The Sun binomial work distinguishes
[universal modular support, THM-4027](../01-canon/theorems/THM-4027-sun-two-four-six-eight-universal-modular-solubility.md)
from [an exact global hole, THM-4026](../01-canon/theorems/THM-4026-sun-two-four-six-eight-binomial-counterexample.md).
Neither an integral precision formula nor an average-count argument supplies
those missing global predicates.

## 2. The six-concept board and inheritance pass

| Object/representation | Preserved question and operation | Hostile or lost coordinate | Cheapest decisive test |
|---|---|---|---|
| Physical comb contacts / raw roof | Bound triple mass using fixed-pair capacities | (1,19,79): exact flow still exceeds physical mass | Retain each carrier's preferred pair before summation |
| Switching class / cycle-count vector | Minimize cumulative eigenvalue after deletion | Antibalanced signing has no negative even cycles | Exhaust small joint profiles and restore subset multiplicities |
| Laurent polynomial / small-root product | Find first nonzero coefficient after cancellation | z^-2+z^-1+z-z^2 waits until m=4 | Compare moment order with exact rational iteration |
| Integral jet matrix / Smith filtration | Recover coefficients at finite precision | p=2,e=1 gives (0,0,2,2), not (0,0,1,3) | Divide the missing derivative trace and test a unit minor |
| Trinomial support / numerical semigroup | Identify channels under addition | A first slice can omit later carry channels | Test both canonical residue carries before applying a two-ray theorem |
| Labelled triple / formal coefficient locus | Preserve a collision under deformation | A tangent direction can fail at second order | Solve the next coefficient equation in the actual source module |

Each of the first four proof notes records its closest proved mechanism,
canonical hostile, corrected near miss, and least-used sidecar. In particular,
MISTAKE-496 supplies the trivial-character exclusion; THM-2070 blocks
support-to-coefficient inference; THM-4010 blocks determinant-to-Smith
inference; THM-4396 blocks pair information from certifying equality alone.

The incoming trinomial classification changed the fifth concept materially:
our height-35 probe was already contained in its height-60 evidence, while
its exact two-carry law explains why arbitrary two-channel slice arguments
can fail. The root-swap proof avoids that slice assumption and handles any
number of terms. The incoming formal continuation changes the sixth concept
from a tangent equation to a whole formal locus; it supplies no automatic
polynomial termination or observer-lattice inclusion.

## 3. Connections that became proofs

| Source -> target | Map and preserved predicate | Destroyed information / needed sidecar |
|---|---|---|
| Local cycle profile -> all-order spectrum | Restrict to six/seven vertex sets and sum with exact binomial occurrence counts; preserves parity contributions to deletion | The scalar cumulative count forgets lengths; retain c3,...,c7 until weighting |
| Width-two moment -> rational parabolic map | R=z²f, T=-R(0)z/R; ord(T²-id)=2m_*+1 | Bare contact order forgets phases; retain the exact leading coefficient and full T |
| Parabolic iterate -> first-return bound | Petal cycles require distinct eligible critical points | Total degree includes critical poles outside the basin; remove gcd(R,R') and the backward orbit |
| Sheet-phase pair -> omitted-comb tail | Exact periodized roof and component count bound every possible third-tooth contact | Mass alone forgets how many pieces can be touched; retain N_ab as well as T_ab |
| Residual rank -> integral observation precision | One divided trace saturates the missing row; Smith factors measure all p-adic losses | Rank alone loses depth; retain every determinantal divisor |
| Independent one-variable grids -> rectangular jets | Tensor their invertible Smith transformations | Coupled geometric or arithmetic source modules are lost; keep the actual product degree box |

The strongest new connection is the Laurent-to-dynamics map: it identifies
the exact target integer with the number of petal cycles and then counts
the critical points they require. The other advances similarly retain data
until the operation that consumes it. This is a research pattern across
subjects, not a theorem identifying their objects or conjectures.

## 4. Sharp boundaries and the next useful tests

1. **LRC anchor.** The next local target is c/b<44/13, with full raw-carrier
   addresses. Use the new T_ab/N_ab pair bounds together with incoming
   THM-4413 owner transversality; test whether the same pair minimizes each
   active roof. A lower bound on a live component and an upper bound on
   total mass have different directions and need not contradict each other.
2. **Moment niche.** The remaining general width starts at three. A product
   of three small roots divided by one root is not a local root swap.
   Construct the actual correspondence and identify a justified critical
   budget before extrapolating the width-two proof. Odd binomials and
   coefficient-cancellation examples remain mandatory controls.
3. **Spectrum continuation.** For D>=7, probe mixed seven/eight-vertex
   profiles against the exact deletion coefficients. The antibalanced
   example still prevents paying an odd-layer surplus with an even layer
   alone. A sampled larger-order minimum would not replace a base proof.
4. **Integral wildcard.** The dyadic triple solves the first repeated-residue
   corner. By residue-cluster CRT it also gives consecutive-node dyadic
   partitions at n=5 and n=6: respectively (0,0,0,0,2,2,2,2,5,7) and
   (0,0,0,0,2,2,2,2,5,5,7,7), independently checked with integer Smith forms.
   General overfull p-clusters require their actual deeper collision tree.
5. **Formal continuation.** Compare the incoming triple-locus compensator
   with the permitted source-normal coefficient space. Finite-precision
   invertibility on our rectangular grid does not establish that inclusion.

These are bounded next experiments, not advertised one-lemma prize reductions.

## 5. Verification, provenance, and maintenance

THM-4416 uses independent Walsh and direct-parity implementations across
2,130,944 switching classes, plus a separately reviewed all-order induction.
THM-4417 has two algebra reviews, primary Milnor source verification, and
exact convolution/composition controls. THM-4418 has 216,039 exact checks,
independent physical/arithmetic constructions, and a separate proof review.
THM-4419 has 441,971 primary gates, independent integer Smith computations,
and all 69 symbolic minors for the dyadic triple. Saved normal and optimized
outputs retain explicit failure gates. No new Lean formalization is claimed.

Source and output pairs are preserved in each theorem's evidence manifest.
The original shared checkout's five untracked outputs were left untouched.
The isolated worktree began at 3eb2b8a66e56, read incoming work through
05d61bbedf, and pushed evidence/reservations at 09c725138c before promotion.

Used existing method cards: search the statement before the method; correct
the object before sharpening the technique; expose the obstruction before
choosing the scale; compute the repair quotient. The new local-profile card
in META-PATTERNS records distinct signed-cycle and LRC evidence, together
with counterindications. Detailed moment-source sidecars were moved behind
links to keep the mandatory startup surface within its existing byte budget.

LRC(14), the wider Laurent first-return problem, D>=7 cumulative gaps,
general higher/multiscale interpolation, and the unrelated flagship
conjectures remain open in the scopes stated above.
