# Independent review of the fixed homogeneous-clone obstruction

**PASS: proof review and independent exact replay, normal and `python -O`.**
Reviewed target:
[`native/clone_escape/PROOF.md`](../native/clone_escape/PROOF.md).
No correction to its mathematical claims was needed. This review does not
modify the producer's code, matrices, receipts, or existing audit files.

## Scope and the all-order reduction

The reviewed statement excludes a **fixed marked strong tournament** (R,x)
whose homogeneous transitive inflation satisfies either

```
H(R[x := TT_t]) = 63*2^t-35,
H(R[x := TT_t]) = 63*2^t-29
```

on an infinite integer tail. The fixed quotient may have arbitrary exterior
arcs. Its shape and incidences cannot depend on t. This is stronger than the
earlier fixed root-attachment exclusion, but it does not exclude changing
quotients, differently incident prefix blocks, other clock substitutions,
or individual tournaments having these counts.

The run formula is correctly typed. Each Hamilton path decomposes uniquely
into r nonempty clone runs, an optional old prefix/suffix, and r-1 nonempty
old middle paths from O to I. The latter require at least two vertices since
I and O are disjoint. The ordered old pieces have count
`A_r=[U] L M^(r-1) Q`, with the constant1 in both L and Q. The t labelled
clones are assigned surjectively to the r ordered runs; their order inside
each run is forced. Thus the factor is `r! S(t,r)`, without a further
factorial division. Permuting the distinct disjoint middle paths explains
the separate divisibility `(r-1)! | A_r`.

For the greatest index m with A_m nonzero, the coefficient of m^t in the
exponential expansion is A_m>0. It cannot be cancelled by a term of smaller
index. Dividing an eventual identity by m^t, or using linear independence
of finitely many exponential sequences on a tail, therefore forces m<=2.
The remaining identity is `A1+A2*(2^t-2)`. Hence A2=63 and A1=91 or97,
even when the original claim was only eventual. There is no unsupported
substitution of t=1 into an assumed eventual identity.

A strong nontrivial quotient stays strong under this inflation. A weak
quotient retains its condensation cut under inflation, so cannot produce
a strong family. The singleton exception to preservation of strongness
is correctly excluded by the proof's word "nontrivial".

**CITED dependencies:** the Moon--Busch strong minimum gives125 from order10
onward. The strong-deletion theorem supplies a strong order8 deletion of
each strong order9 target; the insertion injection gives H(parent)<=H(child).
The corrected bibliographic route is **THM-1370**,
[`h-spectrum-omits-7-21-all-n`](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
§3, citing J. W. Moon, *On Subtournaments of a Tournament*, CMB9(3),
297--301 (1966), Theorem2, DOI10.4153/CMB-1966-038-7. The bounded extension
argument is also explicit in **THM-4744**,
[`interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve`](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md),
Proofs3--4. These are inherited cited results, not new results established
by the finite program.

## Independent computation

`clone_review.cpp` deliberately uses different implementations from the
producer:

| Obligation | Producer | Independent reviewer |
|---|---|---|
|Hamilton path counts|forward path extension DP|last-vertex deletion DP|
|Strongness|Floyd-style bit closure|directed search from every vertex|
|Order9 isomorphism quotient|minimum degree-block bit code|partial adjacency-preserving bijections|
|Run vector A|squarefree polynomial convolution|literal clone counts, then triangular Stirling inversion|

The independent code starts from the complete inherited order8 catalog,
reselects all135 strong parents with H<=97, and recomputes all256 birth cuts
for each. Results:

```
all births:                       34,560
strong births:                    34,290
births with H91 or H97:                39
ordinary target classes:       2 (H91), 7 (H97)
```

It checks a bijection between its nine independently deduplicated classes
and the producer's nine representatives. No producer canonical key is used
in this comparison. All old rows and added rows are checked to be actual
loop-free tournament orientations.

The producer's order9 canonicalizer also passes source review: its accumulator,
minimum, map key, and shifted unit are all explicitly `uint64_t`. Nine
vertices require36 pair bits. Five of the actual nine canonical codes exceed
2^32, so the width issue is exercised by this certificate, rather than being
only a hypothetical concern.

The reviewer checks that every physical mark of every catalogued H91/H97
base appears exactly once in the receipt. It recomputes all run coefficients
for740 marked cases, consisting of404 strong cases and336 supplementary
weak cases. The659 cases through order8 and81 order9 cases agree with the
producer. The strong cases are:

| Order | H91 classes | H97 classes | Marks |
|---|---:|---:|---:|
|7|15|6|147|
|8|12|10|176|
|9|2|7|81|

For each mark, the reviewer constructs the actual inflated tournaments at
t=1,...,floor((|R|-1)/2)+1 and counts their Hamilton paths directly. The
matrix `r! S(t,r)` is triangular with nonzero diagonal t!, so exact integer
back-substitution uniquely recovers the full allowed A vector. The largest
graph counted has13 vertices; 64-bit count arithmetic is sufficient since
13!<2^63. Every inferred vector equals the producer's vector, and no target
case has A2=63. This also independently checks every literal t=2 count.

Two additional controls bring the total run-vector comparison to742:

* The attained22 family has A=(89,63,0,0), so the slope63 filter accepts a
  realizable positive case.
* The inherited higher-run hostile has A=(9,10,2), detecting the genuine
  third exponential mode.

## Context and reproducibility

The root README correctly keeps attachment-chart failures separate from
forbidden Hamiltonian counts and keeps reusable two-branch closure OPEN.
Its next-target discussion can now cite this stronger fixed-quotient
exclusion. No inference of source identity from equal H, or native legality
from a scalar response identity, is used in this proof.

Run from the repository root:

```
python 05-knowledge/results/c4_port_returns_20261103/audit/check_clone_review.py
python -O 05-knowledge/results/c4_port_returns_20261103/audit/check_clone_review.py
```

The script compiles into a temporary directory and writes
`clone_review.json` / `clone_review_optimized.json`. Their contents are
identical. Producer normal/optimized receipt pairs are also checked for
byte equality. The review receipts pin all mathematical producer inputs
and review code by SHA256. No larger ordinary order9 catalog is assumed
or generated.
