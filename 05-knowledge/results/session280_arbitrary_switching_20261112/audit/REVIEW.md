# Independent review of Session280

**Final verdict: PASS WITH MINOR SCOPE/EXPOSITION REPAIRS for ancestry,
wildcard, coverage, switching and escape.** No substantive mathematical
failure was found in the retained theorems. The result is not a Chan–Pak proof.

## 1. Claims read and independent reasoning

### Hidden active-word ancestry

The new unknown-label theorem is valid on its stated INPUT strip
`0<=d<=A<=W<2A`, with positive odd A and physical parity. The active block
coefficient bounds `2h>=a+1` and `2b>=t+2` give
`2A'-W'>=A>0`, establishing strip invariance. The output has `D=A` odd.
Then `V/D=a+t(W/A)` lies in `(3t/2,3t)`. These are disjoint OPEN intervals
for successive powers of two. Thus length uniqueness is not a numerical
approximation, nor an assumed fixed schedule. Odd D gives the unique canonical
residue `a=V/D mod t`; recomposition and the strip guard decide the formal
predecessor. Strict scalar growth excludes cycles and proves termination.

The physical seed gate remains necessary. The proof correctly distinguishes
formal strip packets, actual root packets, graph isomorphism classes and scalar
Hamiltonian counts. The claim does not cover a first input outside the strip;
its paid initialization paragraph handles that case without silently changing
the source. Different seeds in a bank may be ancestors of each other, and
only fixed-seed or antichain-bank uniqueness is asserted.

One small reproduction typo was noticed: active words through length10 number
511, not1023. It was already corrected by the author before checkpoint2. No
mathematical statement depended on the typo.

### Observation chart and common dyadic lattice

The probe matrix determinant is exactly1 after the two parity normalizations.
The unnormalized count chart has determinant4, as retained in the proof.
The chronological three-observation determinant is exactly `b_2*t_1`; the
nonautonomous scalar recurrence follows by eliminating W and retains its
initial owner. A literal physical pair of equal first two observations is a
proper hostile control rather than a claim that all possible pairs fail.

For the Smith theorem, every exterior matrix maps the ambient integer lattice
into L0 and maps L0 into2L0. The stated parity calculation proves both facts
without a spectral assumption. Nonvanishing products modulo2 force the first
Smith factor1. The minor gcd and determinant then give the claimed factors;
the block100 calculation proves sharpness for every depth, not only the finite
controls. The theorem is explicitly for ambient `(A,W,d)` coordinates, not an
unstated integral change to the physical parity sublattice.

### Switching capacity

The last-letter envelope is a weighted-average bound on actual cone rays.
Independent rational reconstruction verifies every row of the63-letter upper
certificate and the83/98-letter lower certificates. The399/400 exponent is
rigorously obtained using the integer-power inequality, not float evaluation.
The geometric series constant164 is conservative and correct. Counting words
rather than scalar outputs is safe in the scarcity upper bound. Conversely,
reciprocal presentation divergence is not promoted to scalar support divergence.
The countable-bank and sublogarithmic bounded-repair arguments retain the
scale comparison and the fixed alphabet hypothesis. Empty-word seed values
are correctly excluded from the countable-bank nonempty-closure statement.

The ancestry theorem permits an additional exact statement: for a fixed thin
seed, these presentations correspond to DISTINCT ROOT PACKETS. It does not
remove scalar collisions or make their projections cover an interval.

## 2. Independent checks

Run `python3 audit/independent.py` or its absolute path. It imports none of the
candidate verifiers and uses independently reconstructed block coefficients,
integer minors and exact rational growth matrices. Normal and optimized output
are identical.

-47,244 cross-word/small-source image cases, all distinct, with exact length and
 modular-word recovery;
-300 entire hidden variable-length itineraries, up to100 blocks and lengths2..20,
 decoded without supplied labels or lengths;
-16,129 chronological observation determinants;
-sharp repeated100 Smith factors at every depth1..40;
-all63 upper,83 lower and98 lower exact pressure inequalities.

All three candidate verifiers also replayed byte-identically against committed
`results.json` files. These finite checks support, rather than replace, the
all-order algebraic proofs above.

## 3. Remaining global obligation

The rigorous successes concern state observation, ancestry, integer lattice
structure and branch supply. None proves that each prescribed positive odd
integer other than7,21 has a physical root packet in the reachable image. This
remaining target-fiber existence problem is not discharged by a decoder,
packet injectivity, or supercritical presentation pressure.

## 4. Switched contraction and the first explicit bad triple

The switching proof was subsequently read in full. The quotient formula is an
EXACT residual identity, not a derivative estimate. Its reference-dependent K
is affine in the new normalized ray. The positive-semidefinite convexity
identity therefore transports the three extreme-ray inequalities to the whole
reachable cone. Tracking the last letter selects the correct changing metric;
no common fixed Perron vector is assumed.

An independent expansion of all9,408 two-dimensional quadratic inequalities
passes. Each metric is strictly above(1/2)I. The factor2 comparison of native
positive coordinates after the first block is checked separately, and positive
subsequent matrices preserve it. This supplies the otherwise necessary absolute
normalization factor: the squared coordinate widths are at most
`8 B^2 E0 (249/250)^(k-1)`. Independent rational arithmetic confirms the bound
first drops below1 at depth3292 forB96 and this certificate. It stays below1 at
all later depths. The proof does not claim this threshold is optimal.

The hostile three-block product has the stated characteristic polynomial and
opposite signs at-2 and-1. Its expanding transverse root refutes blanket
cone-slice contraction. It does not establish a scalar collision between
physical integer seeds, nor invalidate the unrelated three-count decoder.
The reported minimality is expressly for total length of ACTIVE TRIPLES.

## 5. New multiplicity corollary reviewed

Coverage Section5b correctly combines distinct packet ancestry with the lower
weighted-word certificate. The exact inequalities `512^delta<gamma` give the
divergent packet series at exponent1+delta. The upper assumption
`m(N)=O(N^epsilon)` for epsilon<delta contradicts that series by an ordinary
convergent p-series comparison.

The finite-depth reciprocal sum is at least `gamma^(m-1)/512`, while every
count at that depth is at most`512^m`. The harmonic-sum estimate therefore
forces a large fiber among those counts. Choosing `m=floor(log_512 X)` yields
the stated cumulative-in-X lower bound. The same argument may use the maximum
fiber RESTRICTED TO THAT DEPTH, so dividing by the maximum order`1+9m` really
does bound the number of unrooted classes. One unrooted class has at most one
intrinsic packet per choice of root. No fixed-order theorem or specified-target
existence is claimed. This is compatible with an infinitely sparse scalar
support and must not be marketed as coverage.

## 6. Joint151-language certificate and exact packet generating functions

The151-letter metric extension was independently checked on all68,403
transition/ray pairs. Every rational metric is above(1/4)I; the exact initial
energy is`12225700863/1690000000`. With squared rate`9991/10000`, the certified
B96 bound drops below1 at depth15417. The same proof as the56-letter result
applies without a new limiting or spectral assumption.

The pressure witness was independently reconstructed on exactly this151-word
list. All151 row ratios exceed1069/1000. The integer inequality
`1024*1000^128<1069^128` validates delta1/128. Thus stability and supercritical
supply genuinely refer to the same physical language. Known-schedule scalar
recovery remains distinct from unknown-schedule recovery with three probes.
The maximum-fiber and unrooted-class bounds transfer with growth1024 and
order bound1+10m; they do not become pointwise target-existence bounds.

The order polynomials were checked from the actual alphabets. Ancestry gives
unique packet histories, including their order, so the rational word generating
functions count packets exactly. The multivariate substitution sends the
exponents to precisely `(hA+bW-d,aA+tW,A)` and adds the block's vertex count.
It is well-defined on the stated cone-supported, positive-A-degree completion,
not unrestricted formal series with negative powers. The growthA'>=5A gives
coefficientwise uniqueness by iteration. The finite packet set at each fixedA
makes specialization to scalar multiplicities legitimate. This is explicitly
not an A000568 enumeration or closed scalar recurrence.

## 7. Core surgery and paid width doubling

The positive cover injection gives E_i>=h and its endpoint sum yields both
receiver-gap inequalities, including the empty-deletion singleton. Substituting
the old gate differences into the fusion formula gives the claimed sandwich
defect. Its equality classification uses the strict boundQ_L-V_L<Q_L; the
extremal source/sink contexts only ordinal-multiply the count gap. Neutral
padding leaves the unequal deletion coordinate, so it is not a full rearming.

Independent backward incoming-edge subset dynamic programming reconstructs
the three9-vertex seeds with counts269,271,273 and common(h,V,e)=(15,441,15).
It also independently counts the six15-vertex constructions as24875 through
24885 in steps2, with receiving values28119 for the first triple and27639
for the second. The candidate verifier separately replays the universal-law
controls, common frames, unit-context transport and long exact recurrences.

The arbitrary left unit-context update closes on(H,h,V,e), so the repaired
triple really is renewable. The even-k congruence and the even-to-odd binary
carry give the displayed six-count family at every stated scale. The common
receiver mismatch of the two translated triples is nonzero and is multiplied,
not canceled, by a common sandwich. The general width-doubling lemma correctly
states the needed actual shared-interface source and the resulting repair bill.
No repeated generic doubling or eventual interval follows.

## 8. Final small repairs and summary scope

The transferable radix lemma now explicitly assumes integer radixb>=2 and
positive integral source coordinates. Its two displayed equations recover
(t,a,A,W); any additional coordinates need their own exact inverse. This
repairs an omitted abstract-domain/invertibility qualification; the concrete
tournament proof was already valid.

The synthesis was reviewed with requested clarifications that the broad-cone
odd/parity conclusions concern physical packets, that packet injectivity is
not an UNROOTED isomorphism-class injection, and that neutral padding does not
repair the old deletion defect. These are scope clarifications, not changed
numerical or algebraic conclusions.

The optional C++ four-edge search was recompiled and replayed: all58,905
candidates and the three reported matches reproduce byte-identically. Both
independent Python audit programs also agree in normal and optimized mode.
