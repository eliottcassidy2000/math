# Independent audit: no repeated beta roots and the reduced boundary

**Verdict: PASS for the all-model analytic theorem and exact identities.**
The complete producer proof was read. The repeated-root exclusion is valid
in the full weak two-interlacer model, including canceled poles and roots
at zero. Compactness gives a uniform positive mutual beta-root gap. Neither
a distance of beta roots from zero nor a strict reference in every fibre
follows. The reduced extremum and zero-atom statements have the stated scope.

## Independent proof of the exclusion

For a monic weak interlacer A of the nonnegative-root B, cancel common
factors in A/B. Weak ordered interlacing removes every higher-order pole
and leaves nonnegative simple residues of total mass one. Thus ordinary
and shifted moment matrices are PSD even before beta simplicity is known.
This supplies the necessary conditions without circular approximation by
strictly interlacing data.

I independently reconstructed the moments through order5 by expanding the
finite geometric inverse of the reversed B, rather than the producer's
forward division recurrence. The determinants were then evaluated by
summing signed permutations. The resulting identities are exactly

    det C_shifted2=(x-75)/3,
    det C_shifted3=-x(x-75)^2/27+15(x-75)y/7
                   -16y^2/49+(x-75)z/3,
    det D_ordinary3=-49x^2/144+269x/4-18y/7-3132.

These are moments of C/B and D/B; they are not Newton sums of the roots of
the quartic numerators. At a repeated zero of B, y=z=0. The first PSD
condition gives x>=75, the second forces x=75, and the third is -37/16.
No positive y floor is needed for this contradiction.

For a repeated positive beta root r, weak interlacing forces C(r)=D(r)=0.
The linear coefficient system has determinant r/12, so solving it is
legitimate precisely in this branch. Independent elimination gives

    x=24r^3/7-36r^2+108r,
    y=3r^4-28r^3+63r^2,
    B'(r)=4r^2(2r^2-14r+21)/7.

The two possible positive roots are(7-sqrt7)/2 and(7+sqrt7)/2. Reduction
modulo the quadratic gives x=126-12r, y=735/4-49r and D_ordinary3=5r-75/4.
For the upper root, y=49/4-(49/2)sqrt7<0; for the lower root, the moment
determinant is -5/4-(5/2)sqrt7<0. These are exact algebraic signs. They
exclude every repeated positive root and hence every higher multiplicity.

The argument only used necessary equations at the putative repeated root.
It need not reconstruct z or prove that both algebraic candidates actually
give beta polynomials. Rejecting both necessary candidates is sufficient.

## Quantifiers, compactness and the extremum consequence

The ordered beta vectors with sum13 and square sum59 form a closed bounded
set. Coefficients of both prescribed quartics depend continuously on those
vectors. Requiring their real ordered roots to weakly interlace is closed,
so the full model is compact. It is nonempty by the exact central controls.
The minimum adjacent beta gap is continuous and strictly positive everywhere
by the exclusion. Its minimum on this compact set is positive. This is an
existence statement for a uniform constant, not an explicit lower bound.

The incoming independently audited fixed-phase boundary theorem in
`long_frontier2_sep06_boundary.md` was read in full. Its compactness and
indefinite-Hessian argument puts every response extremum on a zero, repeated
or shared-interlacer boundary. The present theorem removes the repeated
alternative in precisely the full model. Zero and common-interlacer nodes
remain legitimate. The original first-phase condition P(-s)=0 is retained
when applying that result separately to each fixed s>0.

For a simple B node, its A/B residue is A(r)/B'(r), so a shared A-root
means that this one weight disappears. It does not mean that beta nodes
collide. The ordered middle-resultant selector applies only to fibres with
its stated simultaneous positive-definite reference. The proof supplies
no such reference for every fibre, and a determinant zero outside the
joint PSD intervals does not become admissible merely by being a boundary.

## Zero-atom reduction and hostile controls

At z=0 the simple beta zero has residue3/7 for C/B and1/7 for D/B.
Cross-multiplication independently verifies

    C=(3/7)G+(4/7)v Ctilde,
    D=(1/7)G+(6/7)v Dtilde,

with exactly the quartic and two cubics in the producer. Removing these
fixed positive atoms leaves normalized positive measures. Conversely,
positive residue representations for the two cubic/quartic ratios, combined
with those atoms, imply weak interlacing of C,D with vG. Allowing repeated
positive G roots when formulating the converse creates no loophole: the
same full-model theorem excludes them. Thus this is a faithful iff geometry
reduction, with the original anchors retained.

The geometry reduction does not replace the original binomial carriers.
The original cubic phase equation on z=0 and the full carried doubled
response remain the target; no lower-height factorial-row sign is imported.

The independent source checks the native hostile
B=v^2(v-3)(v-5)^2 and C=v(v-2)(v-5)^2. Their literal root lists weakly
interlace, while the D residue determinant is -37/16. This refutes the
B-only and C-only simplicity extensions. At both(84,35,0) and(84,35,1),
independent Sturm counts, square-free gcds and positive full Hankel leading
minors verify the complete model. The zero control confirms why mutual
root separation must not be reported as separation from the origin.

## Reproduction and raw evidence

The referee imports no producer or repository implementation. It uses exact
rational polynomials and algebraic numbers, geometric-series inversion and
permutation determinants. Its49 always-active gates pass in both normal
and optimized Python, with byte-identical raw LF stdout. Computations
confirm complete polynomial identities and native controls; the analytic
arguments above supply the universal model and compactness quantifiers.

    python continuing5_20260906_simple_beta_boundary_audit.py
    python -O continuing5_20260906_simple_beta_boundary_audit.py

- Referee source SHA256: `3ffa1995a3fd5c9df3d11ca447c81ce41ed6647f68674716f315036597daae29`.
- Referee output SHA256: `27981148b4724708886ba05b6dd1ec3baa7c2bd1a6cdd3275f95a6dfcfb1305a`.

The producer proof read has SHA256
`d0890541d81121948ffeefac4861a13283d59b70dd6c00d8c3fed2aeb217fa01`.
Its source and transcript are pinned respectively by
`979a4540c112b049e032bd1320ed78e7991a5996b0a5a4991cd88e185bbcc185`
and `8f65fb22877e43128560c0dcc3f4fb4c346111adb0d2175324a72cff4a73b7d6`.
No mathematical repair is requested. Status promotion and repository
filing are separate from these frozen outside artifacts.

The [raw-byte checkpoint manifest](continuing5_20260906_manifest.json) pins
the filed source, report and identical normal/optimized transcript. Any
candidate-report hashes above identify the pre-promotion bytes.
