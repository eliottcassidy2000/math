# Independent audit: the degree-eight two-measure geometry decoder

**Verdict: PASS, analytic theorem and stated finite hostile independently audited.**
No mathematical repair is required. The decoder is an exact equivalence for
the stated real coefficient pencil with x,y,z>=0, including singular Hankels,
zero roots, repeated roots and common factors. The degree-seven hostile is a
separate finite-exact obstruction. Neither result proves the remaining
original-phase response sign theorem.

The producer under review is
`continuing4_20260906_moments_packet.md/.py/_certificate.json`. This referee
does not import or execute the producer. Its source uses SymPy rational
expansion, symbolic quotient matrices, all principal minors and exact Sturm
counts. The producer uses a separate standard-library implementation. This
report also audits the parent's compact-interval fibre corollary below,
added to the producer report after its original 128-gate source freeze.

## 1. Universal proof: why the eighth moment is enough

For a monic degree-d denominator B and a real numerator A of degree d-1, the
coefficients of A/B at infinity obey the B recurrence. Define L on
E=R[v]/(B) by its first d moments. The recurrence proves that L(v^j mod B)
equals the formal moment m_j for every j. Consequently the d-by-d ordinary
Hankel matrix through moment 2d-2 is exactly the Gram matrix of L(pq) on E.
This identification uses the fixed denominator; it is not a general claim
about arbitrary truncated moment sequences.

Multiplication by v is selfadjoint by associativity in E. If the Gram matrix
is positive semidefinite, its radical is its nullspace and is invariant under
multiplication. The quotient is a finite-dimensional positive Hilbert space
with a real selfadjoint multiplication operator. The spectral measure of 1
therefore represents A/B as a sum of nonnegative-residue simple real poles.
Equality follows by the complete expansion at infinity. A singular Gram
matrix causes no gap: invisible repeated and nonreal denominator factors may
be canceled. No moment-nine positivity or flat-extension hypothesis has been
silently added in the degree-five application.

The independent source explicitly checks both symbolic companion identities
H T=T^T H and their formal-moment expansions. It also retains the hostile
A=(v^2+1)^2, B=(v-1)(v^2+1)^2: its rank-one Hankel is positive although B has
four nonreal roots with multiplicity. Thus the lemma alone concerns the
reduced ratio and cannot decode an arbitrary denominator.

The prescribed two numerators supply the missing cancellation obstruction.
At a common nonzero zero r their coefficients must satisfy

    x=(24/7)r^3-36r^2+108r,
    y=3r^4-28r^3+63r^2.

Putting r=a+ib, b!=0, reality gives b^2=3a^2-21a+63/2 and then
(2a-7)(2a^2-14a+21)=0. The linear branch requires b^2=-21/4; the quadratic
branch requires b^2=0. Both are impossible. This excludes every common
nonreal zero for all real x,y, independently of positivity assumptions.

If both Hankels are positive, any nonreal B root would have to cancel from
both ratios, contradicting this calculation. All B roots are therefore real.
For x,y,z>=0 and u>0, B(-u) is strictly negative, which places all roots in
[0,infinity). A reduced positive-residue ratio has strictly interlacing simple
poles and numerator zeros. Restoring a common real factor adds the same
root-counting increments to numerator and denominator and preserves weak
interlacing, including multiplicities. Conversely, weak interlacing of the
monic real polynomials gives nonnegative residues after cancellation and hence
the two positive-semidefinite Hankels. Both directions and their boundary
types are valid.

## 2. Exact finite hostile and its actual representing measures

The referee independently reconstructs

    x=78071/1000, y=601/50,
    z=127473806477/203203019250, s=57/2, M=707/100.

Literal rational expansions of C/B and D/B agree with every frozen moment.
For each channel every principal minor, not only the leading minors, is
strictly positive in the ordinary H4, shifted K4, upper M H4-K4, and the
degree-six quadratic localizer. Both H5 determinants are negative and agree
exactly with the producer's displayed fractions.

There is a complete degree-seven positive-measure interpretation. Let
T=H4^-1 K4. In the H4 inner product, T is selfadjoint and 0<T<M I. The first
three monomial shifts are exact, so 1 is cyclic. It has four distinct interior
eigenvalues and strictly positive spectral weights summing to one. Splitting
j into two exponents at most three recovers moments j<=6, and
<T^3 1,T T^3 1>=m7 recovers the seventh moment. This establishes genuine
representing probability measures, not merely necessary matrix tests.

Independently of the producer's moment linear system, the referee constructs
the quadrature denominator as det(v I-T) and numerator as the (0,0) entry of
H4 adj(v I-T). These match the frozen rational quartics and cubics coefficient
by coefficient. Exact Sturm counts put all four roots in (0,M), and the two
denominators are coprime. Disjointness refers only to these canonical
four-atom constructions; it says nothing about every possible pair of
degree-seven representing measures.

The exact omitted-square identity is particularly informative:

    native m8 - quadrature m8 = det(H5)/det(H4) < 0.

The truncated positive measures therefore cannot obey the native recurrence
at the first omitted moment. The independent engine verifies this identity
for both channels. The surrogate also satisfies all four strict Newton
inequalities, the stated residue floor and slope, and the support endpoint
consequences. It is not an old endpoint-violating surrogate.

The referee multiplies the literal O,E,beta,C_raw,D_raw Laurent polynomials,
then takes their coefficientwise product. Every Q coefficient, including the
lowest nonzero coefficient q(-1)=28 and the full 2t C_raw D_raw term, matches.
Deleting that cross term changes the response. The exact original phase and
response are

    P(-57/2)=0,
    Q(-57/2)=1473043757735348617617612017/28089600000 > 0.

Exact Sturm counts give four simple positive original phases and exactly
three simple positive real B roots, with one nonreal conjugate pair.
Therefore separate bounded-support moments through degree seven cannot
certify the desired response sign. The claimed sharpness is restricted to
this hierarchy of separate initial moment packets; arbitrary cross-channel
or other nonlinear sidecars have not been ruled out.

## 3. Positive and singular controls

At (x,y,z)=(2127/25,27144/625,3564/625), the literal B roots are
(1,3,9,22,30)/5. Direct residues are positive, both H5 matrices are positive
definite, and all their principal minors pass.

At (648/7,54,27/7), the exact common factor is v-3. Both H5 matrices are
positive semidefinite of rank four; their radicals are the reduced quartic B.
Separate rational root intervals prove both reduced strict interlacings.
Restoring v-3 proves the original weak interlacings, so strict positivity
would indeed discard a valid point.

At the repeated boundary (75,0,0), B=v^2(v-3)(v-5)^2 and
C/B=2/(3v)+1/[3(v-3)] is positive. The second channel has ordinary three-by-three
determinant -37/16. This confirms the need to retain both channels.

## 4. Additional audited compact-interval fibre corollary

For every fixed x,y>=0, the set of z>=0 satisfying the full weak model is a
compact interval, possibly empty or a singleton. Moments 0,...,4 do not depend
on z. The exact z coefficients of moments 5,...,8 are

    C: (1,14,130,904+4x/3),
    D: (1,15,147,1067+19x/12).

All these moments are affine in z, as checked symbolically. Each H5 condition
is therefore a one-parameter positive-semidefinite matrix pencil. Their
intersection with z>=0 is closed and convex. The proved decoder identifies
z with the product of five nonnegative roots whose sum is 13, so
0<=z<=(13/5)^5 by AM-GM. The set is bounded and hence compact. The fixed square
sum 59 prevents equality in AM-GM, but strictness is unnecessary for
compactness. This establishes interval geometry for every fixed x,y without
an endpoint formula, a claim of nonemptiness, or a response-sign consequence.

## 5. Reproduction and frozen identities

Run the independent source from the repository root:

    python 04-computation/continuing4_20260906_moments_packet_audit.py
    python -O 04-computation/continuing4_20260906_moments_packet_audit.py

It reads the producer source, transcript and JSON only for frozen byte pins
and exact certificate comparison. It never imports the producer. Dependency
lookup supports the adjacent outside layout and the filed results layout.
SymPy 1.14.0 was used. All 601 gates remain active under optimization; normal
and optimized runs produce byte-identical LF transcripts.

- Referee source SHA256: `5befcd5655e488e882c346478d3c16fa1664922cfcde9bcf1f566d4f28172b8d`.
- Referee output SHA256: `8d172c7a17b41f1b23d3aeda000cbc4486622943a1d226d281559e40fde88d52`.
- Producer source pin: `1b7d02c20c631f95b8a539af8948bf92a9b6dc1504ddc05d5a00b6b13892a9bc`.
- Producer output pin: `54e693133895b249e31aeeed325d786ee4ed6c46d8343fa082e802adcbf6932e`.
- Producer JSON pin: `e396e221c4cd8e4b668f0f07a572ed5314b25cb7c66a8ac68eb439d84979170f`.

Files: the computation source, the retained identical normal/optimized
transcript, and this report. The independent referee phase made no Git
or maintained-file changes; root owns the canonical integration.

Integration correction: the outside producer freeze normalized CRLF in its
capture wrapper. Root added explicit stdout LF configuration before filing
and reran raw normal/optimized output against the unchanged LF transcript
and JSON. The source pin above is the corrected filed source; mathematical
code and the128-gate universe are unchanged. The final manifest pins this
referee after updating that source pin.
