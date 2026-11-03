# What survives composition, averaging and a changed clock

Status: derived statements are labelled below; literature statements are
CITED with their scope. No claim of literature-wide originality is made.

## Inheritance and the live concept board

Anchor: the positive-source route and minimum mean carrier for tournament
attachments. Niche: mixed-weight Keller composition. Wildcard: prime first
arrival and MRB endpoint approximation. Closest proved mechanisms are the
prior endpoint-block transfer, THM-1330's composition monoid, and session31's
Mellin reconstruction. The81/27 source pair is the canonical hostile.
Corrected near misses are signed observations treated as positive gadgets,
linear orbit scaling treated as composition-closed, and density treated as
a first-arrival law. The least-used relevant sidecars are the Paley7
quadratic norm1451 and the MRB negative-moment boundary.

Five live concepts were compared after each computation:

| Concept | New result | What it changes elsewhere |
|---|---|---|
| Allowed future operations | parity carrier571/2668 | internal allocation can be discarded only in the specified fixed-attachment language |
| Actual positive source laws | uniform-birth variance lower bound21/4 | a closed linear mean recurrence is insufficient for the target certificate |
| Composition and retained coordinates | orbit cocycle and degree recurrence5delta+2 | linear source families must be closed under composition before calling them a classification |
| Sampling clock and phase | prime onset, finite cyclic waiting identity, CRT prefixes | density and first arrival need different data and quantifiers |
| Moving analytic endpoint | MRB error transport below | fixed moment accuracy does not control a functional approaching a singular boundary |

These are instances of the existing META-PATTERNS cards “Controlled
forgetting and unlabeled quotients require a sidecar” and “Re-evaluate a
certificate after every fibre-changing operation”. No new universal
meta-pattern is promoted on the strength of numerical coincidences.

## A useful MRB criterion recovered from the old moment bridge

Session31's audited
[tree moment measure](../session31_mirror_mrb_20261022/tree-moments/TREE_MOMENT_BRIDGE.md)
has

    M(p-1)=p^(p-2)/Gamma(p),   Re p>0,
    g(k)=k^(1/k)=k^2/[Gamma(1/k) M(1/k-1)].              (1)

This is a full-kernel re-expression, not a new evaluation of MRB. Its
endpoint expansion is M(1/k-1)=k-log k+gamma+O(log^2(k)/k). A fixed finite
positive quadrature with nodes away from zero has a finite limiting
negative moment B; replacing M by it instead gives g_N(k)~k/B. The MRB
summands fail even to tend to zero. The old continuous endpoint model
repairs this and yields an absolutely summable O(log^2(k)/k^2) correction.

Here is a precise criterion for further approximations (**PROVED**).
Suppose M_tilde(1/k-1)=M(1/k-1)(1+epsilon_k), with real |epsilon_k|<=1/2
eventually, and define g_tilde by (1) with M_tilde. Then

    g_tilde(k)-g(k)=-g(k) epsilon_k/(1+epsilon_k).       (2)

Since g(k)->1, absolute summability of the kernel correction is equivalent
to sum_k |epsilon_k|<infinity. In particular a certified bound
|epsilon_k|<=C/k^(1+eta), eta>0, suffices. Mere relative convergence to zero
does not suffice: epsilon_k=(-1)^k/k for k>=3 makes the alternating
correction asymptotic to-1/k, so that corrected series diverges. These
values are positive moment approximations individually; no claim is made
that this hostile sequence comes from a single positive Mellin measure.

The strongest useful next algorithmic target is a joint schedule N=N(k)
whose certified relative moment errors obey a summable bound, or a better
continuous endpoint model. Fixed-order quadrature success at a regular
mirror point does not establish that bound. The MRB irrationality and
CM-defect questions remain OPEN.

[Crandall's primary paper, section7.5](https://www.marvinrayburns.com/UniversalTOC25.pdf)
already supplies the eta-derivative representation. The repository's
[entire deformation and audited zero](../session31_mirror_mrb_20261022/ENTIRE_MRB_DEFORMATION.md)
give a separate sensitivity test; neither is recast here as a new modular
identity.

The earlier annotated12/13 prompt has a concrete bridge already proved in
[CLOCKS, U1/U2](../weighted_crt_clocks_20261028/CLOCKS.md): for
u=(5+sqrt29)/2, the level-two eta/Fricke expression samples g58^24=u^12.
This is+1 modulo13 and-1 modulo11, with different lift behaviour. The
primitive multiplication clock recovers the trace information lost by the
twelfth-power clock. The Zagier formula supplies the exponent; the split
and inert quadratic rings explain its residue effects. This is the useful
analogy with the new carrier: an observer may forget data only after its
allowed next operations have been fixed. No identification of these rings
with the tournament carrier is asserted. :codex-annotation{index="1"}

## The ETH result: a new algorithmic landscape with explicit input size

The primary
[formalization-scope document](https://raw.githubusercontent.com/openai/math/main/lean/docs/378.md)
states exact deterministic fixed-k CNF counting in time
(L+1)^O_k(1) 2^O_k(n(log log n)^2/sqrt(log n)), polynomial space; its
linked October25 manuscript states a sharper bound. This is CITED here:
we read the stated machine/input scope but did not build its Lean project
or independently audit the complete algorithm. The older repository
[import packet](../external_preprints_mub6_local_smoothing_eth_merge_20261027.md)
is provenance, not a substitute for that scope.

Consequences and limits below are direct deductions from the stated
algorithm, not new claims about its proof. A subexponential SAT algorithm
invalidates ETH-based lower bounds when their reductions and parameters
apply. It does not imply polynomial-time SAT, or automatically supply a
short positive construction for every desired Hamiltonian count.

There are at least2^binom(n,2)/n! tournament isomorphism classes, since one
class contains at most n! labelled orientations. Thus explicitly listing
them still requires2^Omega(n^2) output objects. Computing the scalar
A000568 value by Burnside, listing its classes, and constructing one source
with H=h are different tasks. The new carrier recurrence concerns the
third task's observations; it is not an enumeration of classes.

A proposed SAT compilation must account for its auxiliary variables and
bit lengths. The straightforward circuit for exact Hamiltonian counting
has2^n subset-DP states. Encoding that circuit as CNF retains this cost;
subexponential runtime in the much larger CNF variable count need not beat
direct tournament search. The promising concrete use is exact bounded
search after symmetry/port reduction, with the source orientation and an
independently verified count emitted as the certificate. A uniform proof
still needs a recurrence that closes across all scales.

## Typed connections and stopped implications

| Source -> target / map | Preserved predicate | Lost data / sidecar | Cheapest decisive test |
|---|---|---|---|
| Old forest -> endpoint/parity type | every fixed-attachment mean | internal allocation lost; named endpoints and parity retained | all64 birth cuts at order6; passes170752 pullbacks |
| Actual tournament law -> uniform cover profile | conditional first moment | selected port and cross moments lost | order4 prefix variance propagated; uniform route fails at every n>=5 |
| Keller map -> quotient plus orbit factor | composition, determinant, generic degree | quotient alone loses x scaling | square K; degree7 factor refutes linear closure |
| Prime ordering -> residue density | asymptotic class proportions | phase, first-hit index and dependence lost | factor all nine pre-1451 candidates |
| Paley7 norm -> signed prime candidate | 1451=38^2+7=5*290+1 | primeness and earlier exclusions separate | trial division and explicit preceding factors |
| MRB moment approximation -> kernel | exact formula(2) | endpoint relative error retained | summability rather than fixed-moment agreement |

The productive common principle is operation-sensitive information, with
the exact retained coordinate proved separately in each setting. There is
no single claimed bijection between Keller maps, primes, Collatz or MRB.
The counterexamples narrow the possible proof strategies rather than
establishing global conclusions by analogy.
