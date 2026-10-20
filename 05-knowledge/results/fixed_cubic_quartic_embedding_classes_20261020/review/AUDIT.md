# Independent audit: fixed cubic/quartic embedding classes and cylinders

20 October 2026.

## Verdict and scope

**PASS within the stated hypotheses, including the exact-count refinement.**
The main note proves the at-most-two bound and generation criterion; the
separate `EXACT_EMBEDDING_COUNT_APPENDIX.md` proves the exact count. The
fixed algebra has exactly two ambient plane-embedding classes when
`c != 0` and `2a = 3b`, and exactly one otherwise. Every fixed cylinder
with degree semigroup `<3,4>` has no entire generating pair whose
`(x,y)`-Jacobian is a nonzero element of `C[x]`.

There is no conclusion about arbitrary moving cubic/quartic families,
generic-only generation, globalization of all actual source algebras, the
whole generic conductor-six Keller sector, or the planar Jacobian conjecture.
The inherited fixed-cylinder theorem and the sharp triple-curve witness
were read and checked in their exact roles. No general embedding-equivalence
classification or relative lifting theorem is being smuggled in.

All reviewer writes are confined to this review directory. No producer
file, shared mathematical checkout, or remote repository was modified by
the reviewer. The parent owns publication. Final promotion changed only status lines,
review pointers, and the README; stripping those bookkeeping additions
recovers both prepromotion proof-file hashes exactly. All checkers,
receipts, and source records remain unchanged.

## 1. The source input applies to every embedding

A fresh independent web search returned the indexed original text of
Abhyankar–Sathaye, *Uniqueness of plane embeddings of special curves*,
1996, printed page 1062, Section 2, at:

https://scispace.com/pdf/uniqueness-of-plane-embeddings-of-special-curves-omujaliiec.pdf

The official DOI is https://doi.org/10.1090/S0002-9939-96-03254-6.
The original paragraph identifies the initial characteristic entries with
the two coordinate pole orders and explicitly describes triangular degree
reduction until neither initial entry divides the other. This is the
load-bearing all-embedding input, not merely a numerical classification
of possible abstract semigroups. No newly downloaded byte-identical copy
of that paper is claimed.

The full local author PDF of Sathaye–Stenerson, *On Plane Polynomial
Curves*, dated September 10, 2014, was independently read. Section 3.1(4),
printed page 9, gives the nonprincipal bounds on the first two and later
entries. The definitions, gcd-drop condition, and primitive-generator
property in Sections 3 and 3.1 support the rest of the finite reduction.
The author source is:

https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf

For the present algebra, the function field really is `C(y)`: its index
in `C(y)` divides both polynomial degrees three and four. The normalization
is consequently `C[y]`, and any plane embedding of the same algebra has
that normalization and the same unique place at infinity. The source
hypotheses therefore survive changing the generating pair.

The triangular reduction terminates by the strictly decreasing sum of
degrees. A constant terminal coordinate would make the algebra polynomial
and normal, contradicting its proper inclusion in its normalization.
The remaining coordinate degrees lie in `{3,4,6,7}`. The five unordered
nondividing initial pairs are the ones listed in the note. Every pair
involving seven already has gcd one and cannot acquire a missing primitive
generator later. The other possibilities give `(4,3)` or `(6,4,3)`.
Swapping the first two entries preserves their gcd, lcm, later bounds,
and conductor contribution, so the ordering convention is harmless.

## 2. Free modules, filtrations, and the two candidate pairs

For a monic polynomial `p` of degree `n`, `C[y]` is finite free over
`C[p]`. The other polynomial `q` is integral. If the joint function field
is `C(y)`, its minimal polynomial over `C(p)` has degree `n`; its coefficients
belong to `C[p]` by normality and integrality. Thus the claimed power basis
of `C[p,q]` is valid. The field hypothesis is essential and is not dropped.

Distinct degree residues make leading cancellation impossible. Counting
missing degrees computes the actual vector-space codimension, not merely
an upper bound. Applying this to `1,k,k^2` over `C[h]` proves the complete
semigroup `<3,4>` and delta three. It also gives all low-degree filtration
spaces used in the affine normalizations.

In degrees `(3,4)`, the two coefficients that must be nonzero make the
affine coordinate change invertible. In degrees `(4,6)`, scaling, subtracting
the first coordinate from the second, and deleting constants give exactly
`p=k+beta*h`, `q=h^2+gamma*h` without changing generation.

The basis changes to `1,q,r,qr` have determinant one over `C[p]`, including
when the correction to `r` contains a multiple of `p*q`. Nonzero beta
gives residue-separated degrees `0,6,11,17` and delta seven. With beta
zero, the wrong gamma gives `0,6,9,15` and delta six. Both contradict
generation of the delta-three algebra. This checks the universal pair
classification, not just the sampled parameters.

The extension from the displayed normal form to every fixed subalgebra
with semigroup `<3,4>` is sound. Chosen degree-three and degree-four elements
generate a subalgebra with precisely the same three gaps. Inclusion and
equal finite codimension force equality. All subsequent source changes
are affine over `C`, not unsupported rational changes over `C(x)`.

## 3. Exact generation and exact number of classes

When `c=0`, both coordinates of the second candidate are even, while `h`
is odd and nonconstant. Hence they cannot generate. When `c != 0` and
`2a != 3b`, the further corrected element has degree seven and leading
coefficient `-c(2a-3b)/2`. Coprime degrees four and seven supply the
primitive-field condition rather than assuming it. The new basis degrees
`0,6,7,13` yield delta five and rule out generation. On `b=2a/3`, the
stated universal recovery identity exactly recovers `(c^3/8)h`.
Therefore the generation criterion is an iff.

It remains important not to conflate generation with inequivalence.
The refinement addresses this correctly. Every intrinsic automorphism
extends to an affine change of the normalization variable and preserves
the degree filtration. Its effects on `h` and `k` therefore have the
triangular affine form with nonzero diagonal coefficients stated in the
note. Thus every such automorphism has an ambient affine lift in the
standard embedding. This justifies reducing unparametrized image
equivalence to a prescribed lift of the second candidate.

On the two-generator locus, let `s` be a root of

    s^3+(4a/3)s-c=0.

The endpoint quadratic is `y^2-s*y+s^2+a=0`. Its two roots have equal
`h` and `k` values. Their common `h` value is `-c+(a/3)s` and their squared
separation is `-3s^2-4a`. The resultant with the cubic is `-27c^2`, so
no endpoint pair collapses when `c != 0`.

If `a != 0`, the cubic has at least two distinct roots. Indeed a monic
cubic with just one distinct root is a cube, and its zero quadratic
coefficient would force that root to be zero and then `a=c=0`. Distinct
sums give distinct singular image points because their `h` values differ.
At each singular point, any prescribed ambient lift agrees to first
order with `(V,U^2+(3c/2)U)`: its difference from this map is divisible
by the implicit curve equation, whose value and gradient both vanish.
The forced determinant `c/2-(2a/3)s` cannot be the same nonzero constant
at these different points. Hence no ambient lift exists.

If `a=0`, the scaling by `d^3=c` reduces both parametrizations to the
already audited triple-curve example. Its local determinant four and
period ratio `-28/5` disagree. The intrinsic-automorphism argument closes
the same equivalence gap there. Consequently there are exactly two
classes on the stated locus, and exactly one elsewhere.

## 4. Fixed-cylinder and ramification quantifiers

Entire equality `C[x,P,Q]=C[x] tensor A` specializes to generation on
every complex fiber. The inherited finite-period lemma explicitly allows
a finite number of embedding classes. Its conormal determinant is a unit
of `C[x]`, hence a nonzero constant. A finite set of reference embeddings
and intrinsic normalization automorphisms gives a finite set of possible
period values; the period polynomial is therefore constant. No algebraic
family of ambient automorphisms is required.

The integration-by-parts sign is correct. Equal endpoint values of every
element of the fixed algebra kill the boundary term, giving
`I'(x)=(s-r)j(x)` whenever the Jacobian is `j(x)`.

The new ramification split is valid and strictly broader than the
collision-only statement. A common zero of `h'` and `k'` forces the
y-derivative of every cylinder element to vanish there. Evaluating the
Jacobian on that fixed line forces `j(x)=0`. Generation is not even needed
for this branch of the argument.

The two branches exhaust the proper normalized plane algebras in question.
A unibranch singular plane-curve germ has vanishing coordinate derivatives
on its normalization; if one derivative were nonzero, that coordinate
would be a local formal parameter and the germ would be smooth. A
multibranch point instead supplies a collision. Delta three guarantees
some singularity. In particular the cusp is handled by ramification, not
falsely declared to have distinct colliding normalization points.

## 5. Independent exact route and reproducibility

The independent checker imports no producer code and uses explicit
exceptions, so optimization does not remove its checks. Its main route
is different from both the producer's residue bases and conductor formula.
In the free `C[u]`-module `C[u,y]/(p(y)-u)`, it computes the determinant
of the columns `1,q,q^2,q^3` in the normalization power basis. When nonzero,
this determinant's u-degree is the codimension of `C[p,q]` in `C[y]`, by
Smith normal form over the one-variable PID.

The universal determinant has 2,959 terms. Its u-degree-seven coefficient
is `9*beta^2`; at beta zero the u-degree-six coefficient is
`-(2*gamma-3c)^2`. Both independent index calculations confirm the main
exclusions for every parameter specialization in their stated branches.
For the remaining candidate, dividing by the corresponding determinant
for `1,h,h^2,h^3` gives exactly

    -(c^2/64) * L(u)^2,
    L(u)=8a^3-28a^2b+32ab^2-12b^3+(8a-12b)u-c^2.

The standard determinant has degree three and leading coefficient minus
one. This independently yields delta five off `2a=3b`, delta three on it
when `c != 0`, and nonprimitivity when `c=0`. The determinant's sign has
no effect on codimension.

Further controls include a 27-case rational parameter box, an axiom-aware
recursive enumeration of characteristic sequences, universal recovery
and collision identities, exact sharp-example periods, a composite of two
entire triangular automorphisms, a generic-only generating hostile with
variable period, and a cusp ramification control.

The independent checker passes 70 checks normally and with `-O`, with
byte-identical output. The main producer checker passes 79 checks and the separate exact-count
checker passes 15 checks, each normally and with `-O`, with identical
output matching its frozen receipts. The eleven-file producer manifest
verifies in full. Final hashes are pinned in `REVIEW_SUMMARY.json`. Reproduce from the packet directory:

    python3 review/independent_checks.py > /tmp/fixed-independent.json
    python3 -O review/independent_checks.py > /tmp/fixed-independent-O.json
    cmp /tmp/fixed-independent.json /tmp/fixed-independent-O.json
    python3 verify_finite_embedding_classes.py
    python3 -O verify_finite_embedding_classes.py
    python3 verify_exact_embedding_count.py
    python3 -O verify_exact_embedding_count.py
    sha256sum -c SHA256SUMS.txt

The source theorem and the universal proofs, not these finite controls,
establish the all-embedding and all-fixed-cylinder quantifiers.

## 6. Dependency and publication-layout checks

All four hashes in `SOURCE_HASHES.json` were independently checked against
the staging files. The three referenced mathematical dependency files also
exist at their canonical repository paths and match those hashes exactly.
The cited source PDF is a reading aid and is excluded from publication.

The note's sibling links are intentionally canonical publication-relative
links. They do not resolve from the nested staging directory, which is
not the publication root. Their targets are present in the shared results
parent and verified in `dependency_verification.json`. The publication
owner should retain that sibling layout when copying this packet.

No mathematical flaw or scope contraction was required. During the audit,
the generation iff, exact-count refinement, and ramification split were
independently verified before their inclusion; they do not turn a
fixed-algebra theorem into a moving-family theorem.
