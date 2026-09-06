# Independent referee: the entire original zero-beta sign boundary

**Status: PASS, with no mathematical repair requested.** I accept the coefficient-packet theorem and its C-only zero-boundary corollary in [continuing6_20260906_zero_boundary.md](continuing6_20260906_zero_boundary.md). The proof retains the original phase and all response carriers. The statement holds for every phase, including the unbounded third branch, not merely for a finite root bank. The original two-interlacer shared-root residual and general sign question remain open.

The independent verifier imports no producer code. It reconstructs the native Laurent expressions and formal residue moments, compares all 107 certificate coefficients, and additionally proves the nine univariate chart signs by an independent exact Sturm path. Normal and optimized runs pass **251 always-active exact gates**, with byte-identical raw LF output.

## 1. The predicate, normalization and necessary packet

The producer fixes z=0 and keeps the complete binomial14 odd/even polynomials O,E and the prescribed beta,C_raw,D_raw carriers. In particular

    P=O star beta,
    Q=(O^2+t^(-1)E^2) star(beta^2+2t C_raw D_raw).

The second interlacer remains in the response even when its geometric interlacing hypothesis is dropped. My source obtains O and E by dense coefficient extraction from (1+u)^14, expands the two Laurent factors directly, and reconstructs their coefficientwise product. It recovers q_(-1)=28, the degree-eight polynomial T=sQ(-s), and

    P(-s)/2002=-(12/7)y s^3+x s^2-10s+1/11.

There is no change of phase variable or lower-height replacement after removing a beta zero. The inverse and mixed carries survive.

I independently expand the rational function C/B at infinity and obtain

    det H_shifted,2=(x-75)/3,
    det H_ordinary,3=(x-75)(135-x)/9-8y/7.

Weak C-interlacing with nonnegative beta nodes gives a positive residue measure, including cancellations at repeated nodes. Its ordinary and shifted moment forms are positive semidefinite. Together with y>=0 these give exactly the necessary packet

    x>=75, 0<=y<=(7/72)(x-75)(135-x).

This implies 75<=x<=135 and y<=175/2. The converse geometry claim is neither needed nor made. In particular, this larger coefficient domain is lawful for proving a sufficient response sign.

## 2. Exhaustion of original phases

The endpoint inequalities are valid uniformly on the packet: the original normalized phase is positive at1/110, negative at1/90 and63/1000, and positive at1/8 and9/16. The first three follow from the containing rectangle. The latter two follow from substituting the curved upper bound for y and completing the square in x. The resulting minima are strictly positive even if the vertex lies outside a particular admissible slice.

For y>0 the phase has degree three with negative leading coefficient. The three separated sign changes therefore give all three roots, each positive and simple, in

    (1/110,1/90), (63/1000,1/8), (9/16,infinity).

There can be no extra real or complex root because these already exhaust the degree. At y=0 the degree is two with x>0 and the first two sign changes exhaust it. The disappearing third branch is not assigned a finite root. Since P(t)=2002p(-t), this proves that all roots of P are strictly negative and simple throughout the packet.

## 3. Original elimination and the complete sign proof

Solving the original phase for y gives the displayed Y(x,s), valid for s>0. The independent reconstruction agrees coefficient by coefficient with

    sQ(-s)=-F(x,s)/968,
    (72/7)s^3(Y-c(x))=J(x,s)<=0,

where F has exactly the eleven reported terms and its x^2 coefficient is 264385 s^5(19601s+31920)>0. The sign orientation is correct: the target requires F>0. The positive quadratic leading coefficient would make an endpoint-only argument invalid in general; the producer correctly retains its discriminant and derivative.

For the first phase branch, the raw T is separately convex in x,y on the containing rectangle. Its maximum is therefore bounded by its four corner values. Each negative corner polynomial is positive on the whole closed interval, not just at its endpoints.

On [63/1000,11/100], the discriminant of F in x is strictly negative, while its leading coefficient is positive. Thus F>0 for every real x there. On [11/100,1/8], J(105,s)>0 and J is increasing for x>=105, so packet feasibility forces x<105. The certificate F_x(105,s)<0 and increasing derivative of the quadratic imply F decreases throughout x<=105; hence F(x,s)>=F(105,s)>0.

On [9/16,3/5], J(84,s)>0 and J decreases throughout x<=84. Thus feasibility forces x>84. The strictly positive coefficients of F(84+X,9/16+u) then cover that portion of the third branch. The positive coefficients of F(75+X,3/5+u) cover every later s. These two complete quadrants prove the unbounded tail directly; no empirical tail cutoff or asymptotic approximation is used.

## 4. Independent exact certificate path

For each of the nine univariate chart polynomials I reconstruct its native power coefficients from the independently derived T,F,J. I then compute a rational Sturm sequence. Both endpoints are strictly positive, and the endpoint sign-variation counts agree, excluding every interior root. This supplies a separate sign proof from the producer's positive transformed-coefficient method.

The degrees and variation pairs are:

| Chart | Degree | Left/right variations |
|---|---:|---:|
| -T(75,0,s) | 6 | 4/4 |
| -T(75,175/2,s) | 8 | 5/5 |
| -T(135,0,s) | 6 | 4/4 |
| -T(135,175/2,s) | 8 | 6/6 |
| -disc_x(F)/s^4 | 6 | 1/1 |
| J(105,s) | 3 | 1/1 |
| F(105,s) | 6 | 2/2 |
| -F_x(105,s) | 6 | 1/1 |
| J(84,s) | 3 | 1/1 |

I also independently reconstruct and compare every homogeneous coefficient, including the padded degree-eight first charts, and all 21 monomials in each quadrant. Exactly 107 positive coefficients are covered; no monomial or constant term is omitted. The source checks the derivative identities used to orient both feasible-x barriers.

## 5. Nonvacuity, hostiles and final boundary consequence

At (x,y,z)=(84,35,0), the remaining beta quartic has four positive roots, and both independently reconstructed full residue Gram matrices have positive leading principal minors. Thus the full-model zero boundary is nonempty. The repeated configuration (75,0,0) satisfies

    B=v^2(v-3)(v-5)^2,
    C/B=2/(3v)+1/(3(v-3)),

while D(5)!=0. It is a valid C-only point in the packet even though the two-interlacer simplicity theorem does not apply. The new sign proof introduces no hidden simplicity assumption.

The two essential projections are separately challenged. At (84,1050/11,0), s=1/6 is an original phase with positive T, but the ordinary C Gram determinant is -639/11. At the actual full-model point (84,35,0), T is positive at s=1/4, but the normalized original phase equals335/176 rather than zero. The exact positive response fractions match the producer. Thus neither the packet nor the original-zero condition can simply be dropped.

Combining this result with the already proved full-model simplicity and fixed-phase boundary-extremum theorems is sound: a nonnegative response would give a nonnegative maximum, attained by a possibly different admissible shared-root shape at the same original phase s. This is an extremum transfer, not a pointwise claim that the original shape lies on the boundary. The entire z=0 alternative is now strictly negative, and repeated beta roots were already excluded. Such a remaining adverse shared-root configuration must therefore have z>0 and simple positive beta roots. No sign on that residual is supplied here.

## 6. Reproduction and frozen pins

Run `continuing6_20260906_zero_boundary_audit.py` normally and with `python -O`. The captured `.out` and `_optimized.out` are byte-identical raw LF files. Adjacent files and the standard sibling results layout are supported; no producer implementation is imported.

- Audit source SHA256: `9c147bdad50fe34b650bd43faf3ab9217a48af41cbf11f672de052bf2aeea427`.
- Audit output SHA256: `ecd931bf78deb70bdd85b9698304c9fd4372f5aff9f4cc3e03cca8ed5c78d5b1`.
- Producer source: `dd28188852e5fb955d19d886687f2cb791e25cbf1a6f1c3b2bf66158cc467e6b`.
- Producer output: `6fe1ddf5e2182828405968b34c8236bc8ccea8f5b046e843781dcf26362fef0c`.
- Producer certificate: `f6c169d4190f96465f93d195fe4517537a9c40e85370ae6a736e0e2886e9b12e`.
- Audited prepromotion report: `0bd6098695c32374707e6cb24f24669ba6722b01b2d7b1bcfe725f157c1ead3a`.

All referee files were written outside the repository. No producer, maintained file, namespace or Git state was changed. Parent owns promotion and integration.
