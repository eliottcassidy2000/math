# Independent audit of the middle-resultant-root selector

**Verdict: PASS.** The exact PSD interval selector, the full two-channel
intersection and all three strict-reference examples are correct. The
channel switches, positive lower endpoint and nonnested one-channel models
are exact. No mathematical repair is requested. The strict-reference
hypothesis remains essential to the proved selector; no reference supplier,
global nonemptiness or carried-response sign theorem is inferred.

## 1. Universal proof audit

The formal series of A/B, for A=C,D, gives affine moments through order8.
Its nonzero z slopes are exactly the four displayed in the producer.
Consequently the slope matrix has a zero first row and column and the
remaining block

    S = [0 0 0 1; 0 0 1 a; 0 1 a b; 1 a b c].

The determinant is1. An independent direct proof of its inertia writes
S=[[0,A],[A^T,D]], with A=[[0,1],[1,a]] invertible. The variable change
u=u'-A^(-T)Dv/2 removes the lower-right block. After replacing v by A v,
the form is 2u' dot v, with two positive and two negative directions.
Thus the full slope has inertia(2,2,1zero) for all parameters. This is an
all-parameter congruence argument; the finite congruence controls merely
check its arithmetic.

At a positive-definite reference H(z_*), congruence gives

    H(z) ~ I+(z-z_*) H(z_*)^(-1/2) H' H(z_*)^(-1/2).

The symmetric coefficient has exactly two positive, two negative and one
zero eigenvalue. Its positive eigenvalues give two determinant roots below
z_*; its negative ones give two above. Repetition within either pair is
allowed. All five normalized eigenvalues are nonnegative precisely on the
closed interval between the largest lower root and the smallest upper root.
This proves the middle-root rule and includes coincident determinant roots.
The leading quartic coefficient is1, because the slope vanishes on the
first row and column and its lower4 block has determinant1, while H_00=1.

The resultant identity is also correctly normalized. At simple B roots,
the Vandermonde factorization gives det H=product A(b_i), since the product
of B'(b_i) has sign(-1)^10=+1 times the squared Vandermonde. Polynomial
continuation gives the identity on repeated-root boundaries. There is no
absolute-value substitution. The rank loss at an endpoint equals the
multiplicity of its determinant zero, directly from the affine normalized
eigenvalue factors.

For x,y>=0 and a simultaneous reference z_*>=0, intersect the two exact
PSD intervals with z>=0 and apply the proved degree-eight iff decoder.
This gives precisely the complete weak-model fibre. A reference at z_*=0
is permitted if both matrices are positive definite; it need not mean all
beta roots are strictly positive. The separate control at(84,35,0) has one
simple zero beta root, four positive roots and two positive-definite Hankels.
The producer's wording keeps this distinction.

## 2. Independent exact reconstruction

The standalone referee uses Python's standard library only and imports no
producer or SymPy. It independently performs literal rational series
division and reconstructs each of the six monic quartics from five exact
Hankel determinant evaluations. A second path reconstructs it from five
literal9-by-9 Sylvester determinants. The universal affine-rank argument
bounds the degree by4, so these are polynomial identities, not sampled
agreement. Every rational coefficient agrees with the producer JSON.

Separate rational Sturm chains verify all24 isolating intervals, all four
real roots of each quartic and the selected middle ordering. All30 strict
reference leading minors agree exactly and are positive. Both exterior
positive-determinant components fail positive definiteness, as predicted.
The strict reference points also have five positive beta roots and four
positive interlacer roots by independent root counts.

The active upper channel changes from C at(x,y)=(155/2,9) to D at
(155/2,37/4). At(311/4,21/2), the D lower endpoint lies in
(163/2101,9/116)>0, the D upper endpoint in(177/530,176/527), and C's
lower endpoint is negative while its upper endpoint is strictly larger.
Therefore these are the actual full-fibre endpoints. The fibre contains
z=1/5 and excludes0; direct D matrix evaluation at0 is indefinite.

Both one-channel hostiles survive independent checks:

- D is positive definite at(155/2,9,1/10), while det H_C<0.
- C is positive definite at(155/2,37/4,13/50), while det H_D<0.

Their surviving full-rank positive moment representations already imply
simple real beta roots, and the coefficient signs put those roots in the
positive half-line. Independent literal root counts confirm this. They
are genuine one-interlacer models, not the earlier nonreal-beta surrogate.
Neither channel contains the other globally.

## 3. Frozen reproduction

The [referee source](../../04-computation/continuing5_20260906_pencil_selector_audit.py) and
[transcript](continuing5_20260906_pencil_selector_audit.out) pass **182
always-active exact gates** in normal and optimized modes. Raw stdout is
explicitly LF; both captured outputs are byte-identical without newline
normalization. The JSON is validated as data, never imported as a program.

From the source directory:

    python continuing5_20260906_pencil_selector_audit.py
    python -O continuing5_20260906_pencil_selector_audit.py

When filed in04-computation, the referee finds its certificate in
../05-knowledge/results if it is not beside the source.

Reviewed/frozen identities:

- Producer source: `e5b7e67ba1fa123c3570d267e6fe924312890b4eb326f1656ad4350206aa0b5d`.
- Producer report read: `b51c8ed90e2fe2a70797700c439b0e36994f7b58aa38244055f2613a2d5838e8`.
- Producer JSON: `4f4ea7bdcb823f05b6d74b818597a8e4a4a007cc87f8bcf0fc5bfa9033761b8d`.
- Audit source: `823fc51f2b58b5ef60c5aceeca12dc8bec383b07fbc45924c101c2439dfbbdb4`.
- Audit output: `a69ae411ad56821b1fad8fe8c347640785ad9cc67a5afb51459f68e59d454ee4`.

No producer, maintained file or Git state was edited by this referee.
The universal proof and finite controls are accepted in precisely their
stated strict-reference and weak-boundary scopes.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
