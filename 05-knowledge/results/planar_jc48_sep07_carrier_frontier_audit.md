# Independent audit: arbitrary polynomial coefficients in the y-linear carrier

**Status: INDEPENDENT ANALYTIC AND EXACT SOURCE AUDIT PASS.** The theorem in
[the primary note](planar_jc48_sep07_carrier_frontier.md) is accepted, including
its exact genus formula, geometric integrality, and all nonzero scalar-time
nonrationality consequence. The qualifier is the displayed y-linear invariant;
neither higher y-degree multipliers nor distinct-invariant compositions have
been classified here. JC(2) remains OPEN.

I read the full proof, the standalone source and the incoming constant-slope
and fixed-input flow mechanisms. Independently, the cubic in y has coefficient
sequence -p2 B,-p2 A,p5 B,p5 A-c. The usual cubic discriminant gives precisely
p4 N from the note. At a root of B, comparison of the constant, linear and
quadratic powers of c prevents cancellation of the asserted fixed factor.
The degree parity 2deg B+3 versus 2deg A prevents cancellation at infinity.
Thus neither a generic-coefficient assumption nor squarefree B is hidden.

For local orders a,b, the lower hull has vertices selected from
(0,0),(2,a+eta),(3,b+eta), eta=0 or2. The condition that the quadratic point
lies below the cubic edge is exactly 3a+eta<2b. Its length-two slope gives
ramification a mod2 and its remaining length-one slope is integral. The
single cubic edge has defect two precisely when 3 does not divide b+eta.
At equality, 3a+eta=2b forces a even; all slopes are integral. Transcendental c
makes the initial cubic separable, so there is no additional ramification.
The source's independent lower-hull computation covers both locations and
364 order pairs, including A=0. It retains the actual linear coefficient,
which lies above the relevant hull.

I checked the finite-simple-branch argument separately. A polynomial has
only finitely many critical values even if its critical locus has curves:
the restriction of its differential vanishes on each irreducible curve
component, forcing constancy in characteristic zero. The generic value
therefore has no singular point. The triple-root equations reduce to the
nonzero polynomial A2+3B2p3; at an ordinary double root, the first discriminant
variation is -4q3 times the perturbation's constant coefficient. A repeated
remaining discriminant root would therefore be a forbidden critical point.
Every residual zero contributes exactly one to ramification.

At infinity the two branches with y asymptotic to plus/minus p^(3/2) form
one index-two place for every degree comparison. The third root is unramified,
either the small Hensel root or the root asymptotic to -A/B. This index-two
place has a second role: the cubic function field cannot be a degree-three
constant extension of a rational p-field. Since the constants extension has
degree dividing three, it is trivial. This justifies geometric integrality
before applying geometric Riemann--Hurwitz.

The independent ramification ledger gives

    2g-2=-6+(n-sum nu+sum rho+1),

exactly the stated formula. The budget sum nu<=2deg B gives g>=deg B+5 when
deg B>=1. At constant B the mandatory p=0 defect two gives g>=6. As separate
hand controls: (A,B)=(0,1),(0,p),(0,p2) have genera6,6,8; the balanced common
factor pair ((p-1)2,(p-1)3) has n=25, removed multiplicity6 and remaining
special defect2, hence genus9; (p2,p4) has n=29, removed multiplicity8 and
special defect0, again genus9. These check sharpness, the failed constant-B
extrapolation, and both exceptional local geometries.

The flow consequence uses the invariant I, rather than f(I). Its first
logarithmic coefficient is W=s8(A(s2)+s3 B(s2)), nonzero by opposite parity.
For a first outer term f_k I^k, the p-displacement has coefficient
2 lambda k f_k s W^k at tau^k. Subsequent derivation iterates have larger
order. The same coefficient at every nonzero integer multiple of lambda
excludes finite order. The already audited fixed-input comparison transfers
this literal calculation to the actual completed source. If both images
were rational they would induce a nonconstant selfmap of the geometrically
integral genus-at-least-six generic curve. Riemann--Hurwitz makes it an
automorphism, contradicting its finite geometric automorphism group. This
does not assert that each image separately is nonrational.

I inspected the verifier's always-active checks and replayed it normally and
with python -O, independently of the producer's replay. Both1018-gate outputs
are byte-identical to the frozen383-byte output. The24 named polynomial pairs
retain common factors, zero coefficients, irreducible factors, both balanced
thresholds and three outer Hamiltonian orders. Their good specializations
only control named examples; the analytic proof above pays the all-parameter
claim. No research producer is imported by the verifier.

    python3 04-computation/planar_jc48_sep07_carrier_frontier.py
    python3 -O 04-computation/planar_jc48_sep07_carrier_frontier.py

Frozen SHA256:

    source 04f6008ed12dc33dbd2ef2ce7429fab8a97bc28f77a5cb5c5d232f1a6277db68
    output 5b27c5e31ea9bfdd6a4fe3976a00b37effd780f33fa8a4faa0c72e88eac6cef2
    semantic f65798e59de1a2cce5a34b8fe387ee5b2087216146e231b094b27d7fbda5fb56

No proof correction was required. A later genus theorem about the particular
completed-response supplier, if established, will be a separate dependency:
that supplier is not y-linear and is not covered merely by this result.
