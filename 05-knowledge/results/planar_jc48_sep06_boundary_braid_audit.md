# Independent audit: the equality sextic's rational braid certificate

**Status: PASS — analytic transport, actual consumer, exact source and frozen
replay independently audited.** This audit concerns
[the braid proof](planar_jc48_sep06_boundary_braid.md) and its stated
six-sheet, single-three-cycle obstruction. It does not identify the full
affine complement group or infer a global Jacobian-conjecture result.

## 1. The polynomial is the actual curve

Independently of the producer, I extracted its literal coefficient table
and formed

\[
 \operatorname{Res}_t(t^4+t^3+t^2-u,
 2t^6+3t^5+2t^3+2t^2-v).
\]

The result equals the declared monic quartic F(u,v) exactly. Direct
substitution of the two parametrizing polynomials into F gives zero.
An independent discriminant calculation reproduces

\[
 -u^5(256u^2+11u+12)
 (1225u^4+4639u^3+14945u^2+26313u+15435)^2.
\]

Thus the certificate is attached to the actual equality sextic in the
[resolution supplier](planar_jc48_sep06_resolution_budget.md), not to an
abstract polynomial with a similar local passport. I inherit the audited
whole-curve singularity classification and the
[odd-cusp necessary passport](planar_jc48_sep06_odd_cusp.md); I do not
reinterpret the finite paths as a proof of those prior results.

## 2. Exact root isolation and consistent path transport

I checked `base_taylor`, `taylor`, `rouche`, `radii` and `tube` line by
line. Their polynomial coefficient expansion is exact in Q(i). The upper
norm |Re z|+|Im z| bounds the usual absolute value from above; the lower
norm max(|Re z|,|Im z|) bounds it from below. The strict Rouché inequality
retains the constant term, every higher vertical term, and every nonzero
base-variation term. It applies uniformly to the whole straight base
segment, not just its endpoints.

Each accepted large disc contains exactly one root. The factor 1/32 in
the radii makes all four discs disjoint, and monic vertical degree four
then accounts for every root. In particular the actual base path avoids
the discriminant. The auxiliary endpoint disc has radius at most both
adjacent large-disc radii divided by 64. Combined with the displacement
bound strictly below half the old radius, this identifies the same
labelled root across successive segments.

Both the actual root path and the polygonal centre path lie in the same
four convex disjoint discs on every segment. Straight interpolation
therefore gives a collision-free configuration homotopy. At segment
joints the homotopies agree: they interpolate the same actual endpoint
root and the same centre. Each loop has the identical labelled initial
centre row, and the final unordered centre row equals the initial one
exactly. The common basepoint displacement is consequently one coherent
identification of the actual fibre for all four loops. It is not an
independently chosen conjugacy gauge for each loop.

The witness really has four centres in every row. Its four rational base
paths are closed at 2+i and have 15,987, 1,330, 1,824 and 2,362 accepted
segments. No claim about which critical values their diamonds enclose
is needed: a certified closed path already supplies a necessary relation.
The verifier's default mode uses exact fractions only; the optional
NumPy producer proposes data but supplies no mathematical authority.

## 3. The marked braid convention and affine necessity

The rational multiplication by 1+i/4 preserves orientation. Pairwise
real-coordinate crossings of a polygon segment have exactly computed
rational times. The verifier checks generic endpoints, distinct times,
adjacency in the current order and nonzero imaginary separation. Thus
the recorded signed chronological word is the actual polygon braid;
cancelling adjacent inverse letters preserves its action.

The note's convention is coherent. Put the meridian access stems below
the ordered punctures. For a counterclockwise half-turn the old left
meridian moves to the final right meridian; the old right one moves to
the conjugate x_right^-1 x_left x_right. This is the inverse-Artin map
on the based loops. Transporting a representation uses its inverse,
giving precisely

\[
 (a,b)\mapsto(aba^{-1},a),\qquad
 (a,b)\mapsto(b,b^{-1}ab)
\]

for the positive and negative letters. Applying these changes in
chronological order agrees with the source. The distinction between a
geometric map on loops and contravariant representation transport is
essential and is stated in the final proof.

The monic-polynomial section pays the actual consumer. On a compact
region filling the four base loops, all roots have a uniform bound;
choose a single constant vertical basepoint outside this bound and below
the projected strands. The section loop contracts in the actual affine
complement. Consequently the necessary monodromy condition is exact
equality of the marked tuples, not equality only up to simultaneous
conjugacy.

The note also correctly proves that the vertical free group F4 surjects
onto the affine complement group. Over the regular base, fibre generators
and section lifts generate the bundle group. The section extends across
the finite critical-value set, killing those base lifts in the whole
complement. Every loop can be perturbed off the corresponding critical
vertical fibres. Monicity prevents roots escaping over a finite base
value. This argument requires neither projective transversality at
infinity nor generation of the entire base group by the four selected
loops.

## 4. Complete finite consumer, independent path and controls

The producer retains all forty single three-cycles. Fixing the first to
(012) is simultaneous relabeling, preserving every Hurwitz relation and
the generated action's transitivity. All 40³=64,000 remaining assignments
are retained before the four filters.

I separately implemented the same consumer with literal permutation
tuples, independent of the producer's indexed conjugation table. It
reproduces the survivor counts

\[
 30400,\quad760,\quad19,\quad1.
\]

The sole survivor is four equal three-cycles, so it moves only three
labels. The constant tuple is a positive consistency control. Because
the four actual fibre meridians generate the affine group, this cannot
be the tuple of a transitive six-sheet cover.

As a supplemental convention check, global word reversal, global sign
reversal, and their combination all give the same counts and the same
constant survivor. I replayed these controls. They corroborate, but do
not replace, the explicit below-stem convention in Section 3.

This is stronger than merely failing to find a transitive numerical
tuple. The root paths, strand words, finite universe and actual
necessity map are separately justified. Conversely the finite bank does
not prove that the whole affine group is cyclic: it covers only this
meridian cycle type and this six-label target.

## 5. Reproduction and exact pins

I read the entire source and replayed both commands:

```text
python3 -B 04-computation/planar_jc48_sep06_boundary_braid.py
python3 -B -O 04-computation/planar_jc48_sep06_boundary_braid.py
```

Both independently completed all **172,168 always-active gates** and
reproduced the frozen 635-byte output, including all four words and both
witness hashes. The literal file sizes and SHA256 values were checked:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| [Source](../../04-computation/planar_jc48_sep06_boundary_braid.py) | 9,080 | `459fc3a3920f2335331ee69704accb170a3b3c5ca31f6903c009672faf6d6046` |
| [Output](planar_jc48_sep06_boundary_braid.out) | 635 | `295b812bdfb38a6b0b2a80dddf3fab26de80d3ed00a6d30414f1cc05951a3b3d` |
| [Compressed witness](planar_jc48_sep06_boundary_braid_certificate.json.gz) | 1,393,532 | `a06d7434da196de9a0f721fd913362dac7fd678358673c8d929d7f2c2a2c6951` |
| Decompressed JSON | 5,621,433 | `5e4c61c9a64e1a30ab841de524be966583ee2a29929c1544aea6f4e2f4675868` |

**Final text audit: PASS.** No mathematical correction is needed. The
source's retained candidate header and transport warning are honestly
interpreted as a boundary on code-only inference; this audited analytic
proof supplies the separately required transport. The primary proof may
now record the independent audit and freeze.
