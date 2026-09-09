# Tournament continuation: homogeneous contacts, output sheets and weighted cycles

**PROVED scoped results + FINITE-EXACT controls + independent audits,
2026-09-08.** Four theorem files turn the prior cross-concept bridges into
tournament statements. General LRC(14), the Hamiltonian/Pfaffian inequality
`H >= disc`, and smooth-fluid singularity/regularity questions remain
**OPEN**. A universal weighted cubic bound is unproved here; no external
priority claim is made.

The central change is to ask **which operation a tournament representation
supports, and what information that operation still needs**. Strong
components organize directed cycles, but homogeneous pairs control cheap
transport. A projective determinant tournament records relative signs,
but misses a physical output sheet. Weighted triangles measure an exact
instantaneous shear production, but scalar block statistics lose arbitrary
boundary contacts.

This continues [the recovered-concepts pass](cross_concepts_20260908.md)
and its [Euler bridge](euler_bridge_20260908.md). The
[external-source audit](../reference/OPENAI-EULER-AUDIT-2026-09-08.md)
remains the scope authority for the supplied OpenAI repository and paper;
this tournament pass does not complete their full PDE verification.
The session inherited `c53b407b72bc` and reserved THM-4463--4466 in
`ffa8059da` before promotion. Rebased commit hashes can change; the exact
source and output bytes are pinned in the theorem frontmatter and manifest.

## 1. Portfolio and recovered concepts

The assigned anchor is tournament continuation, with LRC as an eventual
consumer requiring its actual legal operations. The niche is checksum
orientation transport. The wildcard is weighted Euler shear production.
The live board has six concepts: strong components, homogeneous modules,
context cut metrics, projective output sheets, cyclic triples and boundary
response kernels.

| Prior mechanism or obstruction | What it changes in this pass |
|---|---|
| THM-2256 / automorphism-contact-dichotomy-for-quotient-frustration | Its bounded-versus-linear alternative now has an intrinsic homogeneous-pair criterion |
| THM-4145 / rooted-homogeneous-pair-expansion-two-defect-formula | The same pair that creates a free transport contact needs mixed two-ear data for Hamilton continuation |
| THM-2225 / checksum and THM-2294 / anchored Plucker tournament | Rotation descends to a signed action on sections with exact periods; quotienting loses the output bit |
| THM-1862 / order-join and THM-1926 / strong-core cycle localization | Weighted cyclic mass localizes, but transitive interface products still affect production |
| THM-2013/2016 / unit triangle counts and reducibility ceiling | The unit formula is inherited; a weighted reducible positive witness identifies the failed extension |
| THM-4459 / shear moments and THM-4460 / continuation states | Preserve the actual spatial and continuation data instead of promoting a scalar summary to a complete state |

The new theorem files below provide the exact full paths for these older
results and their corrections. In particular, the old open wording after
THM-4145 is superseded by THM-4162/4163; only its proved expansion mechanism
is inherited. The new work does not advertise the old score formula,
Plucker identity or strong-component facts as new discoveries.

## 2. Homogeneous pairs classify cheap transport

[THM-4463 / homogeneous-tournament-contacts-and-pinned-scale-floor](../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
defines the external row disagreement

```
h(u,v) = number of other vertices seeing u and v differently.
```

For the forced-pair contact functional,
`w(id,tau)=sum_v h(v,tau(v))`. The zero-disagreement graph is a disjoint union
of consistently directed paths: their vertex sets are the maximal
transitive modules. Relative to an automorphism layer sigma, every zero
contact has `sigma^(-1) tau` equal to a matching of adjacent swaps along
these paths, and its diagonal cost is the number of swaps.

Consequently, outside automorphism axes, the integral forced Hall envelope
at block scale N obeys an exact structural dichotomy:

```
some homogeneous pair exists:  phi_R(N) = 1 for every N >= 1;
no homogeneous pair exists:    phi_R(N) >= N, and phi_R(N) = Theta_R(N).
```

Strong connectivity does not determine this alternative. A strong cyclic
three-block substitution with block sizes 2,1,1 has a homogeneous pair;
a directed join of a 3-cycle and a singleton has none. Distinct
automorphisms have contact at least three, which supplies the sharpened
linear lower bound in the Hall-layer proof.

Now retain an integral pseudometric D from actual pinned exterior cuts.
Let `c_D=min_(homogeneous pairs u,v) (1+2D(u,v))`. When a pair exists,

```
min(N,c_D) <= phi_R,D(N) <= c_D;
phi_R,D(N) = c_D once N >= c_D.
```

One exterior pin can distinguish every homogeneous pair and still give
the constant floor three. Thus separation and scale coercivity are
different properties. A sufficient repair is to make every homogeneous
pair's metric separation grow linearly with N.

**Remaining gap:** this classifies the forced envelope, not the complete
transport cost G. Free-axis internal costs and the residual `G-F_R`
still require a separate consumer. The rooted pair-expansion mechanism
suggests a precise next test: evaluate that residual on one adjacent swap
while retaining the mixed two-ear response. No claim about `H >= disc`
or canonical LRC configurations follows from the envelope alone.

## 3. A tournament can return while the output flips

[THM-4464 / checksum-projective-tournament-and-output-sheet-cocycle](../../01-canon/theorems/THM-4464-checksum-projective-tournament-and-output-sheet-cocycle.md)
takes the compatible cyclic checksum response of length 2q and identifies
antipodal phases only after choosing one representative from each pair.
Those q distinct projective directions give a tie-free determinant
tournament. The antipodal ties in the original 2q directions are retained
as the reason a section is necessary.

After q physical steps, the section has been negated. The address and
the full valued determinant tournament return, but the tracked checksum
verdict flips. One retained global sign recovers the output when the
transported section contains the tracked input; a fixed arbitrary section
also needs that input's relative sheet.

For the **labelled section action**, the least possible period is
`2^v2(q)`. If q is a power of two, every section tournament has period q.
For odd q there is exactly one fixed section tournament, the alternating
regular cyclic one. Thus odd and dyadic behavior differ by an exact signed
rotation law, not just a resemblance between their pictures.

There is a second loss: at four vertices, positive vertex weights cannot
generally restore rank-two determinant amplitudes from their signs. The
small witness has determinant magnitudes `(1,1,1,1,2,1)` and satisfies
the Plucker cancellation `1-2+1=0`; its signs alone give `1-1+1 != 0`.
The preserved data for a weighted downstream consumer therefore include
edge amplitudes as well as the output sheet.

**Remaining gap:** the signed action is exact, but it does not itself
improve a causal coin decoder or provide LRC owner transport. Literal word
orbits can be longer than their 2q checksum response; no shorter word
period is inferred.

## 4. Weighted triangles give the exact shear cubic

[THM-4465 / weighted-tournament-shear-production-and-contact-kernel](../../01-canon/theorems/THM-4465-weighted-tournament-shear-production-and-contact-kernel.md)
uses an actual pairwise relation: coordinate i points to coordinate j
when the nonnegative matrix entry M_ij is positive, with no opposing
positive entry. Put S=sym(M), K=skew(M). Then

```
F = 4 tr(S K^2)
  = 3 * (weighted cyclic-triple sum)
    - (weighted transitive-triple sum).
```

In the Euler velocity-gradient equation this is the instantaneous
production of squared skew norm, while the gradient lies in this cone.
In three dimensions it is also the vorticity contraction `omega^T M omega`.
The identity does not assert that an Euler trajectory stays in the cone.

Unit weights recover the old score-variance formula and reducibility
ceiling. Arbitrary weights change the answer: a reducible seven-vertex
example with all weights in `[1/2,1]` has positive F. Even a strong
tournament can have negative F. Strongness is therefore not the weighted
production statistic.

Constant cross-block contacts close exactly on `(size, edge mass, F)`.
A general source or sink attached with vector b instead charges

```
Q_B(b) = sum_(u<v) w_uv b_u b_v.
```

Two oriented three-vertex blocks have the same `(size,edge mass,F)=(3,4,6)`
and the same internal weight multiset, yet the identical exterior vector
`b=(2/3,4/3,2)` changes their full productions to `+2/9` and `-2/9`.
The weighted boundary quadratic form is therefore necessary even to
predict the next sign. The symmetric magnitude matrix is not assumed
positive semidefinite.

## 5. A sharp norm bound survives a large substitution class

[THM-4466 / sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes](../../01-canon/theorems/THM-4466-sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes.md)
puts `E=sum M_ij^2` and lets C be the weighted cyclic-triple sum. It proves

```
F <= 3C <= E^(3/2)/sqrt(3)
```

whenever all positive cyclic triangles share an edge, and throughout the
recursive constant-contact substitution grammar with transitive or
common-edge quotients. This includes all weighted tournaments on at most
four vertices, but also infinitely many larger strong examples.

The mechanism is a common-edge Cauchy bound followed by an exact cyclic
composition law. Quotient edge weights rescale as
`lambda_ij sqrt(n_i n_j)`, so both cyclic mass and squared energy split
into child and quotient sectors. Convexity closes the norm bound.

Equality is not confined to one triangle. An edge of weight two with
four two-edge returns of weights one has `E=12`, `C=8`, `F=24`, attaining
equality on six vertices. Unused edges have zero weight in this witness.

**Remaining gap:** the bound for arbitrary cyclic supports and arbitrary
weights is unproved here. A separate 6,000-draw exact scout at orders
5--10 found no counterexample, which does not establish the general
inequality or its status in the external literature. The next proof
obligation is control of overlapping cycles outside the common-edge and
constant-contact classes, preserving the multi-return equality family.

## 6. Typed connections and decisive tests

| Source -> target map | Preserved predicate | Destroyed information / needed sidecar | Cheapest decisive test |
|---|---|---|---|
| Incidence rows -> homogeneous-pair graph | Zero forced contact | Full reversal cost / internal G and mixed two-ear response | One adjacent swap in a rooted strong expansion |
| Checksum phases -> projective section determinant tournament | Relative signs and valued determinants under common rotation | Absolute output / section root bit and tracked input sheet | q steps: same tournament, opposite verdict |
| Oriented shear matrix -> weighted triangle sums | Exact instantaneous F | PDE realization / pressure, spatial compatibility and time evolution | Direct matrix contraction against cyclic/transitive products |
| Blocks -> `(n,e,F)` | Constant-contact substitution law | Nonuniform contacts / boundary quadratic form | Same scalar state, exterior signs `+2/9` and `-2/9` |
| Substitution tree -> energy sectors | Cyclic mass and squared energy | Arbitrary overlapping-cycle geometry / actual edge amplitudes | Common-edge equality and a nondecomposable weighted support |

The concurrent `1d099f2c6` checkpoint supplies a compatible but separate
warning: its [complete LRC zero-clock classification](continuing14_20260908_lrc_zero_classification.md)
rejects connected possible graphs that force distinct speeds to coincide.
Its [line-state repair](continuing14_20260908_no3_line_states.md) needs exact
continuation costs after gluing. These are real instances of restoring
coordinates lost by a quotient; they do not manufacture a tournament
orientation on a symmetric zero relation. Any LRC tournament consumer
must first declare its intrinsic observable, ties and legal operation.

## 7. Reproduction and stopping state

The four theorem files contain full proofs, quantifiers, equality/failure
boundaries, dependency slugs, independent audits and frozen source/output
hashes. The [manifest](tournament_continuation_20260908_manifest.json)
pins all four implementations and outputs. Run each listed script normally
and with `python3 -O`; every stored output agrees under both modes.

- Contact lane: 124,468 permutation probes and 329,308 Hall multisets.
- Section lane: 223,007 exact gates, including 36,863 sections.
- Production lane: 115,991 exact gates with independent matrix/triple paths.
- Norm lane: 188,936 exact proved-scope gates, plus the separately labelled
  finite scout. The scout is not a proof dependency.

The useful stopping objects are now explicit: the residual transport cost
on homogeneous swaps, the output-sheet cocycle, the nonuniform boundary
quadratic form, and the cyclic norm inequality outside the proved grammar.
They are narrower, testable obligations arising from recovered repo
mechanisms. None is presented as a solution of a major open problem.
