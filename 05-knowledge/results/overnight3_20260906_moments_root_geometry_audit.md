# Independent audit of the complete-row path proof and two-row signs

**Status: INDEPENDENT AUDIT PASS.** The complete-row simple-negative-root
theorem is PROVED using the explicitly cited classical inputs. The
two-row sign bank is FINITE-EXACT. Unbounded cross-row separation remains
OPEN. The conclusion is independently concurrent with the already promoted
THM-4436; the audited report correctly presents another proof of that same
conclusion rather than a second theorem discovery.

Audited [report](overnight3_20260906_moments_root_geometry.md),
[producer](../../04-computation/overnight3_20260906_moments_root_geometry.py), and
[certificate bank](../../04-computation/overnight3_20260906_moments_root_geometry_certificates.json).
The only requested wording change was to say **positive mass** before
deducing `n>0` in the return-fibre consumer. The producer applied it; the
abstract factorial theorem correctly continues to allow `n=0`.

## 1. Mathematical proof audit

The complete binomial-line lemma is valid for arbitrary positive integers
`p,q`, without a coprimality assumption. From `S_i=(-qi,pi)` to
`T_j=(U-qj,V+pj)`, the displacements are
`U-q(j-i)` and `V+p(j-i)`, giving exactly the stated binomial entry.
The level `px+qy` increases along every directed unit edge, with all sources
on one level and all sinks on the other. All paths for any chosen finite
minor lie in a bounded rectangle, so the determinant involution is finite.

Source and sink orders are compatible: an inverted pairing would reverse
the x-order of two paths between the common boundary levels. Parameterizing
the polygonal paths continuously by level forces an intersection. Grid
edges can intersect only at vertices or along shared edges, so this is the
vertex intersection needed for cancellation. Surviving families have the
increasing pairing and positive sign. No missing density of discrete level
values is assumed when `p,q` have a common factor. The zero-level case
`U=V=0` is correctly separated as the constant sequence.

Translating the complete support is exact: the new offsets are nonnegative
and the new second offset lies in `[0,p)`, so all negative indices vanish.
This does not claim that an arbitrary prefix of a PF sequence is PF.

The abstract factorial factorization also passes. The alpha support begins
at zero and extends through at least `h`; the beta support is exactly
`-L,...,h`. Multiplying alpha by `z^L` aligns both complete supports. Their
ordinary coefficientwise product has overlap precisely `0,...,h` before
the shared shift. Unequal polynomial degrees cause no truncation error.
The argument holds at `n=0` and for noncoprime `A,B`, since those restrictions
are absent from this abstract coefficient identity. Constant rows are
handled separately and have no cancellation parameter.

I independently opened the primary
[Brändén paper](https://arxiv.org/pdf/math/0403364). Section 2 supplies the
finite PF/root characterization. Theorem 3.8 uses the minimum of the two
degrees for the factorial-weighted coefficient product. Section 7 identifies
ordinary coefficientwise multiplication after applying `Gamma={1/k!}`,
and states simplicity away from zero. Here the second factor has strictly
negative roots, the product is nonzero, and removing its exact origin
factor leaves positive constant term. All cited hypotheses are satisfied.

For every positive nonempty balanced mass, reducing the third count modulo
`A` gives the canonical `t`, and then `N>=0`. The balance equation gives
`a*n=b*N+c*t>0`, so `n>0`; this rules out a further hidden lower truncation.
It does not assume that the first return equals `gcd(a+b,a+c)`. The supplied
`(-1,1,2)` control is an appropriate boundary: its first return is at mass
two, although that gcd is one.

The scalar invariant and nonzero monomial prefactor preserve exact moment
vanishing. Individual simple negative roots imply the stated phase wall,
but do not couple different masses. The report explicitly preserves that
limitation. Its common-root example using two simple negative-rooted
polynomials is logically valid and is not presented as a trinomial witness.

## 2. Independent finite certificate replay

The independent [audit program](../../04-computation/overnight3_20260906_moments_root_geometry_audit.py)
imports neither the producer nor SymPy. It performs no root approximation
or root isolation. It reconstructs the declared 221-support universe from
its six direction pairs, five degrees, corner residues and first-two-valid
parameter rule, including the two explicit additions.

For each mass, the auditor enumerates the third-monomial count and solves
the literal balance equation

```text
(a+b)*y+(a+c)*z=a*m,       x=m-y-z.
```

Every nonnegative solution contributes its exact integer multinomial
coefficient. This reconstructs **all 442 coefficient rows**, their contents,
their complete supports and their second-row Laurent shifts independently
of the producer's channel formula.

The auditor constructs its own exact rational Sturm chains. Each supplied
first-root interval is strictly negative and disjoint from the others;
it contains exactly one root by Sturm variation. A degenerate interval
must be an exact rational root with nonzero derivative. The interval count
equals the polynomial degree, proving that **all 1,015 first roots** are
covered. Separate Sturm counts and squarefreeness checks verify that
**all 2,242 second-row roots** are also simple and strictly negative.

Independent rational polynomial division computes `Q mod P`. Interval
Horner evaluation reproduces each recorded rational enclosure exactly,
and every enclosure excludes zero. This proves the finite coprimality
claims without trusting a sampled evaluation near a root.

The monomial shift is essential. The exact sign counts are:

| Second Laurent lower exponent | Compressed Q sign | Raw sign | First roots |
|---:|---:|---:|---:|
| 0 | negative | negative | 552 |
| -1 | positive | negative | 463 |

Thus every true Laurent second-row value is negative, while 463 compressed
values have the opposite sign. Positive content normalization does not
change signs; multiplication by `(-1)^low` restores the removed monomial.

The first-degree counts are `(2:44,3:45,4:44,6:44,8:44)`. All four carry
signatures occur with exactly the claimed counts:
`(0,0):60`, `(0,1):60`, `(1,0):50`, `(1,1):51`.
The prefix-truncation hostile also checks arithmetically: the beta-prefix
discriminant is `-20`, whereas the true row discriminant is `10780`.

## 3. Reproduction and evidence

Run from the repository root:

```text
python3 -u -c 'import runpy; runpy.run_path("04-computation/overnight3_20260906_moments_root_geometry_audit.py", run_name="__main__")'
```

The [preserved audit output](overnight3_20260906_moments_root_geometry_audit.out)
records 13,570 explicit gates. Normal and optimized Python runs produce
identical bytes. There are no floating-point constants or assertions that
optimization could disable.

SHA-256 values over LF bytes:

```text
audit source:
b2e515b6010c150bf5f8dccac8a81eccd7650c963599f9701f62e3ead23b8480
audit output:
308f1d48ab3d37ac7a3c0e98967afe4c60dcde1348c5b2c2b960382d49998ce8
audited producer certificate:
3768e534d170bae5dba0ac3d9c817ddbef789b2ac1a43519362495a934b5d9eb
audited producer source:
61d7f81c02d56d0809dbc4081d7961325a43558a6073ebacc702a3ca553f676d
```

No producer or canon file was changed by this auditor. Acceptance of the
analytic proof and finite signs supplies no unbounded cross-row sign law,
cross-mass coprimality theorem, or additional trinomial-width closure.
