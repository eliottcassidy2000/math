# Exact audit packet: unknotting-number-one challenge and Keller 3+1+3 fragment

Date: 2026-07-28  
Scope: read/audit only on canonical files; all new artifacts live under
`.scratch/unknotting_mu3_20260728/`.

Audited `origin/main` at
`95aebffc9c165f0807209e6572a97e76cefac88a`; the four relevant blobs were:

```text
a718dbf609bdde937fa223457f4c18338fbc379d  04-computation/unknot1_decider.py
2759b51836b4752848e7f89b6e4f08ef81fcdc85  UNKNOT1-DECIDER-STATUS-AND-DESIGN...
a6f607e4677ca4d05d5bb7ddd4717cfeb3dc2a89  FRAGMENT-DECODE-mu3-depth2-tree...
520670a415838574e92ccc2ffb9bca5f73530877  AK-FORCING-WORKBENCH...
```

During this audit, a concurrent session placed an uncommitted repair of
`04-computation/unknot1_decider.py` in the shared checkout (blob
`f6649f88a64f48c1cc93fcbb20bffe77874843f9`).  I did not edit it.  I
independently verified that incoming repair: all 13 old and 3 hostile checks
pass, the six-crossing witness returns `UNKNOWN`, no executable assertions
remain, and normal / `python3 -O` outputs are identical.  The isolated
proposal below was developed against the audited `a718...` origin blob and
adds a structured upper-bound field plus one extra hostile check; those are
useful ergonomics, not prerequisites for soundness.

The same concurrent session then added a `MISTAKE-318` correction banner
and repaired true-certificate wording in the status note.  Those changes
correct the demonstrated soundness error.  The status note's separate
Alexander/SNF/gate/genus-one claims listed in Section 9 still require
selective correction; consequently the all-document patch below should be
treated as a content proposal against the audited pre-correction blobs, not
blindly applied over the live concurrent edits.

## 1. Bottom line

1. **OPEN:** an always-terminating exact decider for `u(K)=1` is not known.
   The Epoch problem says the general unknotting number is not even known
   decidable.  A worst-case one-hour Boolean algorithm for arbitrary
   100-crossing inputs would be a major theorem, not a routine program.
2. **DEMONSTRATED ERROR:** the current engine has a six-crossing false
   `TRUE_CERTIFIED`.  A crossing change to the unknot proves `u(K)<=1`;
   the program failed to certify that the input was nontrivial.
3. **REPAIR PREPARED:** the minimal sound repair returns true only when a
   one-change unknot certificate is paired with a nontriviality certificate
   for the input.  Otherwise it retains an `u<=1` certificate and returns
   `UNKNOWN`.  All 13 old checks and 4 hostile checks pass; normal and `-O`
   outputs are identical.
4. **FINITE-EXACT:** the owner's example PD is exactly `K11n3`, not merely
   a knot with matching determinant/signature.  Knot Atlas currently gives
   `u(K11n3) in {1,2}`.  The current `UNKNOWN` is the correct frontier.
5. **PROVED/FINITE-EXACT:** the truncated `3+3+1=7` fragment is THM-2473's
   sporadic Keller-map inverse tree.  In positional order it is `3+1+3`.
   The middle fiber has one finite survivor and two sheets escaping to
   infinity; no points collide.  The fixing symmetry is order two and the
   generic monodromy is `S3`, not `mu3`.
6. **IRRELEVANT TO THIS DECIDER:** arXiv:2607.24528 proves Feige's
   probability conjecture.  It can at most inspire randomized-runtime
   heuristics and supplies no knot-theoretic completeness or time bound.
   The actually relevant new paper is Lackenby's arXiv:2607.23350, which
   improves the exact unknot oracle but does not decide `u=1`.

## 2. Exact false-positive witness

The current program returns `TRUE_CERTIFIED` on

```text
D5 =
[[1,11,2,10],
 [6,10,7,9],
 [3,8,4,9],
 [11,5,12,4],
 [7,2,8,3],
 [5,1,6,12]].
```

Yet `D5` represents the unknot.  The following path starts with a
one-crossing unknot and uses only legal reverse Reidemeister moves:

```text
D0 = [[1,1,2,2]]

R1+:
D1 = [[1,1,2,4],[2,4,3,3]]

R2+:
D2 = [[1,3,2,2],[6,4,7,3],[7,4,8,5],[8,6,1,5]]

R2+:
D3 = [[1,1,2,12],[6,12,7,11],[5,10,6,11],
      [2,8,3,7],[9,4,10,5],[8,4,9,3]]

R3:
D4 = [[1,5,2,4],[8,4,9,3],[10,5,11,6],
      [2,10,3,9],[6,11,7,12],[7,1,8,12]]

R3:
D5 = [[1,11,2,10],[6,10,7,9],[3,8,4,9],
      [11,5,12,4],[7,2,8,3],[5,1,6,12]]
```

Thus every `Di` is an unknot diagram and `u(K)=0`.  On `D5`, the input
greedy reducer has an empty move log and stalls at six crossings.  Changing
crossing 4 then reduces by

```text
R2@(3,5), R1@0, R2@(1,4), R1@2,
```

so the old engine emits:

```text
TRUE_CERTIFIED
change crossing #4 (1-indexed), then reduce:
R2@(3,5),R1@0,R2@(1,4),R1@2
```

This certificate is a correct proof of `u(K)<=1`, but the exact value is
zero.  Reproduce the construction and verdict:

```bash
python3 .scratch/unknotting_mu3_20260728/unknot1_true_soundness_witness.py
```

The script replays Spherogram 2.4.1's legal reverse-move construction and
then invokes the current repo engine.  Its complete output is
`unknot1_true_soundness_witness.out`.

### Failure mechanism

The first false implication is

```text
greedy R1/R2 failed  =>  the represented knot is nontrivial.
```

Failure of one simplification strategy is not a topological certificate.
The missing coordinate is the distinction between:

```text
the changed diagram is an unknot  =>  u(K)<=1
the input knot is nontrivial       =>  u(K)>=1.
```

Only their conjunction proves `u(K)=1`.

## 3. Minimal safe repair and regression

The proposed patch is:

```text
.scratch/unknotting_mu3_20260728/proposed/unknot1_decider-soundness.patch
```

It is a patch against audited origin blob
`a718dbf609bdde937fa223457f4c18338fbc379d`, not against the concurrently
repaired shared working copy.

It makes four changes:

1. adds `nontriviality_certificates` and an `upper_bound_certificate`;
2. accepts `det(K)!=1` or `sigma(K)!=0` as sound, conservative proofs that
   the input is not the unknot;
3. when a visible change unknots but input nontriviality is uncertified,
   returns `UNKNOWN` and retains only the `u<=1` certificate;
4. adds the six-crossing hostile regression and replaces four executable
   optimization-sensitive assertions with `require`.

The determinant/signature rule is deliberately incomplete.  A nontrivial
Alexander-polynomial-one knot may have `det=1` and `sigma=0`; the repaired
engine returns `UNKNOWN` even if it visibly finds an unknotting change.
That is a safe loss.  The preferred future repair is an exact unknot oracle
on the input.

Verification:

```bash
python3 .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
python3 -O .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
```

Results:

```text
13/13 original checks pass
4/4 hostile checks pass
17/17 total
suite failures: 0
normal output == python -O output byte-for-byte
```

On the hostile input the repaired result is:

```text
VERDICT: UNKNOWN
upper-bound certificate only: change crossing #4 ... R2,R1,R2,R1
note: exact u(K)=1 is withheld because input nontriviality is uncertified
```

## 4. Audit of `FALSE_CERTIFIED`

The mathematical implications used by the negative lane are correct:

1. Murasugi's bound `|sigma(K)| <= 2u(K)` makes `|sigma|>=4` a proof of
   `u(K)>=2`.
2. The Montesinos trick makes the double branched cover of a `u=1` knot a
   half-integer surgery.  In particular its first homology is cyclic.
3. For order `d`, Lickorish's necessary linking-form condition requires a
   generator whose self-linking is `+2/d` or `-2/d`.

The exact code algebra was independently probed in
`false_certificate_audit.py`:

```text
400  PD relabel/crossing-order invariance checks
1000 symmetric integer presentations compared with independent SymPy SNF
 940 cyclic and 60 noncyclic cases
zero cyclicity disagreements
zero square-orbit decision disagreements
normal and -O output identical
```

For `7_4`, the selected generator has `q=11 mod 15`; its unit-square orbit
is exactly `{11,14}`, disjoint from the required `{2,13}`.  Hence the
reported Lickorish obstruction is exact.

The generator search is incomplete but fail-safe.  The explicit symmetric
presentation

```text
[[ 998137,-666540,-27741],
 [-666540, 445105, 18525],
 [ -27741,  18525,   771]]
```

has Smith invariants `[1,1,105]`, so its cokernel is cyclic.  Standard
coordinate images are `15,21,35 mod 105`; no basis vector or supported
pair with coefficients `+/-1` generates, while `(1,1,1)` does.  The code
finds no generator and safely declines to obstruct.

Two documentation qualifications remain necessary:

- the status note says “Smith normal form,” but the implementation uses
  adjugate entries/determinant divisors;
- passing three named regression cases does not itself give universal
  “decision rights.”  The justification is the theorem plus a correct
  implementation; the gate is only a regression control.

This audit is strong finite evidence for the exact implementation, not a
formal proof of every PD-rotation and Gordon--Litherland convention.

## 5. The prompt example is exactly K11n3

Starting at incoming arc 1 and traversing the supplied PD gives the
over/under encounter sequence recorded in `k11n3_dt_audit.out`.  Pairing
odd and even visits gives the signed DT code

```text
[4,8,10,-14,2,-16,-20,-6,-22,-12,-18].
```

The official KnotTheory character encoding is:

```text
bdeGaHJCKFI.
```

Parsing `DTCode4KnotsTo11.m`, this string has zero-based index 618.  The
package's own loop over crossing number, alternating/nonalternating type,
and table index assigns it uniquely to `K11n3`.

This is independently confirmed directly by Knot Atlas: its K11n3 page
lists the exact same eleven PD tuples (in a different order), the same
signed DT code, determinant/signature `{43,-2}`, and unknotting-number
range `{1,2}`.

Primary data:

- https://katlas.org/wiki/K11n3
- https://drorbn.net/AcademicPensieve/Projects/KnotTheory/KnotTheory/DTCode4KnotsTo11.m

Reproduction:

```bash
python3 .scratch/unknotting_mu3_20260728/k11n3_dt_audit.py
python3 04-computation/unknot1_decider.py --example
```

Every one of the eleven visible flips has determinant unequal to one:

```text
31,59,3,3,53,35,35,9,7,7,9.
```

Therefore none is an unknot and this particular diagram has `u(D)>1`.
It does **not** follow that `u(K)>1`; K11n3 remains exactly on the current
`1` versus `2` frontier.

## 6. Strongest honest general algorithmic statement

### Positive semidecision theorem

Given any input PD `D`, the following procedure is sound and halts on every
`u(K)=1` knot:

```text
1. Run an exact unknot recognizer on D.
   If D is an unknot, return FALSE (u=0).
2. Breadth-first enumerate all finite combinatorial Reidemeister sequences
   starting at D.
3. At each reached diagram E, change each crossing c of E.
4. Run the exact unknot recognizer on flip(E,c).
5. If one is an unknot, return TRUE with the full certificate.
```

Soundness follows from the input nontriviality check and the exhibited
one-change path.  If `u(K)=1`, some diagram `E` of `K` has an unknotting
crossing; Reidemeister's theorem gives a finite path from `D` to `E`, so
fair breadth-first enumeration reaches it.  If `u(K)>=2`, the process may
run forever.  Thus `u=1` is positively semidecidable, not thereby decidable.

### Why visible crossing search is not complete

Taniyama proves:

```text
for every nontrivial K and every natural n,
there is a diagram D of K with u(D)>=n.
```

Taking an `u(K)=1` knot and `n=2` proves that no algorithm restricted to
the input diagram's crossings can be complete.

Primary source: https://arxiv.org/abs/0805.3174

McCoy gives the sharp special-class exception:

```text
if K is alternating and u(K)=1,
every alternating diagram has an unknotting crossing.
```

Hence an input promised to be an alternating diagram admits a finite exact
procedure: certify input nontriviality, flip every input crossing, and use
an exact unknot oracle.

Primary source: https://arxiv.org/abs/1312.1278

### Normal surfaces and surgery do not close the outer search

Lackenby's new hierarchy paper supplies a new exact algorithm for
incompressibility and hence unknot recognition.  This materially improves
the **inner oracle** used on the input and every changed diagram.

It does not supply a finite list of all possible unknotting crossing arcs
or a bound on the size of a diagram in which a `u=1` change appears.
Section 9 explicitly says several steps are not described precisely enough
to estimate running time; its initial iteration bound depends on quantities
`L` and `g` that are not then bounded.  The refined hierarchy later bounds
length by `4 c_q(H)`, but no implementation or worst-case laptop guarantee
for 100-crossing diagrams is proved.

Primary source: https://arxiv.org/abs/2607.23350

The Montesinos/Lickorish route likewise gives necessary conditions on the
double cover.  Outside special classes, no cited converse turns these
conditions into a finite complete enumeration of crossing circles.  McCoy's
alternating converse is notable precisely because the general converse is
unavailable.

The reinforcement-learning system often finds short unknotting sequences
for diagrams up to 200 crossings, but its output is an upper bound.  The
paper explicitly distinguishes `u(D)` from `u(K)` and says no algorithm is
known to compute `u(K)`.

Primary source: https://arxiv.org/abs/2409.09032

Therefore no rigorous under-one-hour guarantee follows from present theory.

## 7. Exact Keller-map fragment resolution

For the sporadic Keller map in THM-2473,

```text
u=1+xy
F=(u^3z+y^2u(4+3xy),
   y+3xu^2z+3xy^2(4+3xy),
   2x-3x^2y-x^3z)
det J_F=-2.
```

At `v*=(-1/4,0,0)`, exact Groebner elimination gives

```text
x^3-x=0,
F^-1(v*)={P-,P0,P+}

P-=(-1,3/2,13/2)
P0=(0,0,-1/4)
P+=(1,-3/2,13/2).
```

The next fibers are:

```text
P-: 21119x^3-404x-208       (three distinct roots)
P0: z=0, y=0, x=-1/8        (one point)
P+: 20929x^3+532x-208       (three distinct roots).
```

Thus the positional count is `3+1+3=7`.

The core eliminant is

```text
E(x)=Lx^3+(4-3bc)x-2c.
```

At `P0`, `L=0` and `E=4x+1/2`, so the cubic drops directly to a linear
polynomial.  Two sheets escape to infinity on the Jelonek nonproperness
surface.  Since `det J_F=-2`, finite ramification/collision is impossible.

The involution `sigma=diag(-1,-1,1)` fixes `P0` and swaps `P-` and `P+`.
This is `Z/2`, while THM-2473 proves generic monodromy `S3`.  Calling the
branch “mu3-fixed” is not the exact symmetry statement.

Reproduction:

```bash
python3 .scratch/unknotting_mu3_20260728/keller_depth2_exact.py
```

This exact match supersedes the fragment-decode note's Arithmetic-Kakeya
preference and its “two preimages coincide” reading.  It also removes the
fragment as evidence for the AK workbench's identification-gluing
hypothesis; that hypothesis may be studied independently.

There is no mathematical bridge from this Keller-map fiber computation to
the knot `u=1` decision problem beyond a generic warning that an apparent
finite branching count can lose sheets at infinity.  No transfer should be
claimed without a source/target map and preserved predicate.

## 8. The linked 2026 probability paper

arXiv:2607.24528 is Nie--Wei, *On Feige's conjecture*.  It proves that for
independent nonnegative mean-one random variables,

```text
P(X1+...+Xn < n+1) >= (n/(n+1))^n >= 1/e.
```

It contains no knot diagrams, crossing changes, branched covers,
Reidemeister moves, surgery, Keller maps, or `mu3` argument.  A
distribution-free lower tail can inform a heuristic randomized search
budget only after a valid runtime random-variable model is supplied; it
cannot prove exactness or a worst-case one-hour bound.

Primary source: https://arxiv.org/abs/2607.24528

## 9. Correction-ready documentation packet

The exact proposed documentation diff is:

```text
.scratch/unknotting_mu3_20260728/proposed/documentation-corrections.patch
```

It was generated against the four audited origin blobs at the top of this
report.  Since `MISTAKE-318` and the soundness repair arrived concurrently,
its status-note hunk now overlaps those good edits; transplant the remaining
corrections rather than applying that hunk wholesale.

It corrects all current truth-surface overclaims found by exact term search:

1. `04-computation/unknot1_decider.py`
   - false `TRUE_CERTIFIED` semantics and branch;
   - missing hostile regression;
   - optimization-sensitive executable assertions.
2. `05-knowledge/results/UNKNOT1-DECIDER-STATUS-AND-DESIGN-klein-S691.md`
   - `u<=1` promoted to `u=1`;
   - Alexander polynomial claimed in the pipeline but not implemented;
   - Smith normal form claimed but not implemented;
   - a three-example gate said to confer “decision rights”;
   - Coward--Lackenby cited for a stronger practical decidability claim than
     that paper itself states;
   - K11n3 identification and exact `u(D)`/`u(K)` boundary omitted;
   - July 2026 Lackenby unknot-oracle paper omitted.
3. `05-knowledge/results/FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md`
   - wrong preferred host;
   - collision rather than escape;
   - `mu3` rather than the exact order-two fixing symmetry / `S3`
     monodromy.
4. `06-writeups/AK-FORCING-WORKBENCH-klein-S691.md`
   - transfers the misdecoded fragment into an AK identification hypothesis.

No other current nonhistorical Markdown/Python file matched the exact
unknot-decider or fragment-overclaim terms.  Historical commits/messages
should retain provenance but must not be used as the current truth source.

A concise mistake-ledger entry, with mechanism, witness, strongest survivor,
repair, and regression, is ready at:

```text
.scratch/unknotting_mu3_20260728/PROPOSED-MISTAKE-ENTRY.md
```

## 10. Primary-source ledger

- Epoch open problem and exact status:
  https://epoch.ai/frontiermath/open-problems/unknotting-number
- Taniyama, unbounded diagram unknotting number:
  https://arxiv.org/abs/0805.3174
- McCoy, alternating `u=1` diagrams:
  https://arxiv.org/abs/1312.1278
- Coward--Lackenby, genus-one `u=1` anatomy:
  https://arxiv.org/abs/0809.4142
- Applebaum et al., hard unknots/RL upper bounds:
  https://arxiv.org/abs/2409.09032
- Lackenby, new incompressibility/unknot algorithm:
  https://arxiv.org/abs/2607.23350
- Lidman, exact resolution `u(11n102)=2` illustrating the depth of negative
  certificates:
  https://arxiv.org/abs/2606.12431
- Knot Atlas K11n3 data:
  https://katlas.org/wiki/K11n3
- Official KnotTheory DT table:
  https://drorbn.net/AcademicPensieve/Projects/KnotTheory/KnotTheory/DTCode4KnotsTo11.m
- Nie--Wei, Feige conjecture:
  https://arxiv.org/abs/2607.24528

## 11. Reproduction commands

```bash
python3 04-computation/unknot1_decider.py --test
python3 .scratch/unknotting_mu3_20260728/unknot1_true_soundness_witness.py
python3 .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
python3 -O .scratch/unknotting_mu3_20260728/proposed/unknot1_decider.py --test
python3 .scratch/unknotting_mu3_20260728/false_certificate_audit.py
python3 .scratch/unknotting_mu3_20260728/k11n3_dt_audit.py
python3 .scratch/unknotting_mu3_20260728/keller_depth2_exact.py
```

All exact audit scripts use explicit `require` checks.  Their normal and
`python3 -O` outputs were compared and are identical.
