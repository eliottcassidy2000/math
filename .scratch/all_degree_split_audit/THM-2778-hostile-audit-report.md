# THM-2778 all-degree complete split closure: independent hostile audit

## Verdict

**PASS, subject to the scope and text gates listed below.**

I found no mathematical obstruction to promoting the proposed all-degree
closure.  The finite top atlases used in THM-2759/2767 are not needed
uniformly.  The missing replacement is the following recurrence-tail fact:
on the top boundary, simultaneous vanishing of the two top fluxes and the
top response forces the quartic series itself to be a square, and hence
forces the unique point `P_infty`.  Every other top-boundary point is
response-active and can therefore be pulled back to the unique finite
source pole before applying the exact-prefix identities.

The audited statement is:

> Let `M=4k-2`, `k>=1`, and let `(P,Q)` be a planar Keller pair in the
> complete chosen-sheet split polynomial exact-square-prefix reduced Faber
> chart, with normalized top coefficient of `E_M` and every allowed lower
> seed retained.  Then `(P,Q)` is a polynomial automorphism.

This remains a theorem about the stated terminal chart.  It does not derive
that chart for an arbitrary Keller pair and does not prove `JC(2)` or
`DC(2)`.

The canonical theorem must not depend on the reserved THM-2764.  Its
load-bearing dependencies should be proved results such as THM-2129,
THM-2181, THM-2230, THM-2723, the all-degree face identity in THM-2752, and
the exact-prefix exclusion in THM-2760.  The load-bearing part of THM-2760
was included in this audit; see Section 8.

## 1. Exact setup and the rational-primitive dichotomy

On the chosen split sheet write

```text
P=(w^2+d)^2+(qw-s),
d=C-B^2/(4U^2),  q=A/U,  s=AB/(2U^2)-E.
```

In full Faber gauge the reduced target is

```text
Q=E_M+sum_(1<=j<M, 4 does not divide j) a_j E_j,
```

with constant `a_j`.  THM-2129/2723 give constants
`lambda,W,kappa`, `kappa!=0`, such that

```text
Phi_Q=lambda,  Psi_Q=W,  U R_Q'=kappa.
```

If `U` is constant, `(x,z)->(x,w)` is a polynomial source automorphism and
THM-2181 closes the monic depressed quartic.  Otherwise THM-2723 gives

```text
U=u_0(x-a)^m,  m>=2,
R_Q=r_0+r_1(x-a)^(1-m),  r_1!=0.
```

Thus `R_Q` has exactly one pole, the finite point `a`.  All original
polynomial coefficients `A,B,C,E,U` are regular there.

With `wt(h,d,q,s)=(1,2,3,4)`, homogenize as

```text
F=Phi_M+sum a_j h^(M-j)Phi_j-lambda h^(M+1),
G=Psi_M+sum a_j h^(M-j)Psi_j-W h^(M+2),
N=R_M+sum a_j h^(M-j)R_j.
```

The affine response is `N/h^(M+3)`.

## 2. The recurrence-tail square theorem

Put

```text
f(z)=1+2dz^2+qz^3+(d^2-s)z^4,
f(z)^(M/4)=sum_(n>=0)c_n z^n.
```

The coefficient recurrence is

```text
n c_n =
 sum_(r=2)^4 f_r (((M/4)+1)r-n)c_(n-r).              (1)
```

At a top-boundary point, `F=G=0` gives

```text
c_(M+1)=c_(M+2)=0.
```

If the top response also vanishes, then

```text
0=R_M=4c_(M+3)+2d c_(M+1)
```

gives `c_(M+3)=0`.  At `n=M+4`, the only predecessor outside this
three-term zero tail is `c_M`, but its multiplier is

```text
4((M/4)+1)-(M+4)=0.
```

Hence `c_(M+4)=0`; recurrence `(1)` then inductively gives
`c_n=0` for every `n>M`.  Therefore `F_0=f^(k-1/2)` is a polynomial and

```text
F_0^2=f^(2k-1).
```

The exponent `2k-1` is odd.  Unique factorization in `C[z]` implies that
every irreducible factor of `f` has even multiplicity, so `f=g^2`.
Since `f(0)=1`, `deg(g)<=2`, and the coefficient of `z` in `f` is zero,
one may write

```text
g=epsilon(1+dz^2),  epsilon in {+1,-1}.
```

Coefficient comparison gives

```text
q=0,  s=0.
```

Projectively `d` is then nonzero, so the point is exactly

```text
P_infty=[d:q:s]=[1:0:0].
```

Consequently every top-boundary point other than `P_infty` has
`R_M!=0`.  On every normalization branch above such a point,
`N/h^(M+3)` has a pole.  No transversality, reducedness of the ambient
complete intersection, or finite Groebner atlas is required for this
conclusion.

## 3. Every non-`P_infty` boundary point is excluded

Let `X` be the normalization of the reduced integral closure of the actual
generic coefficient image.  The image is nonconstant because `R_Q` is
nonconstant.  The rational map from `P1_x` extends to a finite surjective
morphism `P1_x -> X`.  Thus every boundary branch of `X` has a physical
source preimage.

At `q=0`, direct perturbation of

```text
((1+dz^2)^2-sz^4)^(k-1/2)
```

shows

```text
Phi_M=0,
Psi_M=4 binom(k-1/2,k)(-s)^k.
```

The coefficient is nonzero, so a top intersection with `q=0` has `s=0`
and is `P_infty`.

Now take a top-boundary point with `q!=0`.  It is response-active unless it
is `P_infty`, so its source preimage is `a`.  On a finite local index cover
put

```text
beta=hB/(2U).
```

The exact polynomial-prefix identities are

```text
beta^2+d=h^2C,
q beta-s=h^4E.                                        (2)
```

Because `C,E` are regular at `a`, `beta` is DVR-regular.  If `omega` is
its residue, `(2)` gives

```text
d=-omega^2,  s=q omega.                               (3)
```

If `omega!=0`, the exact-prefix nonzero-root theorem THM-2760 says the two
top fluxes cannot vanish together.  If `omega=0`, then `d=s=0`.  The two
top fluxes vanish together exactly when

```text
k=2 mod 3, equivalently M=6 mod 12.
```

The other residue classes are therefore excluded already.  The surviving
root-zero lane is removed uniformly in Section 4.

## 4. Uniform root-zero section obstruction

Suppose `k=2 mod 3` and a source branch maps to

```text
P_q=[d:q:s]=[0:1:0].
```

On the `q=1` index cover put `tau=v(h)>0`.  Identity `(2)` gives
`v(beta)>0`; therefore

```text
v(q_aff)=-3tau,
v(beta_aff)>-tau,
q_aff=q/h^3,  beta_aff=beta/h.                        (4)
```

For every Faber seed `E_j`, its value at the original polynomial section
`z=0` is weighted homogeneous of weight `j` in

```text
wt(beta_aff,C,E,q_aff)=(1,2,4,3).                     (5)
```

This follows directly from

```text
E_j=sum_(i=0)^j c_i w^(j-i)
```

because `c_i` has weight `i`, and the substitutions
`d=C-beta_aff^2`, `s=q_aff beta_aff-E` preserve weight.

Since `M` is divisible by three on this lane, specializing
`beta_aff=C=E=0` gives the exact pure term

```text
binom(M/4,M/3) q_aff^(M/3).                           (6)
```

Its coefficient is nonzero: `M/4` is a half-integer, so no factor in the
generalized binomial coefficient vanishes.  It has valuation exactly
`-M tau`.

Every other weight-`M` monomial either contains `C` or `E`, making
`a+3r<M` for its `beta_aff^a q_aff^r` polar part, or contains
`beta_aff`, whose inequality in `(4)` is strict.  Hence every other top
monomial has valuation strictly greater than `-M tau`.  Every monomial in
every lower `E_j`, `j<M`, has valuation at least `-j tau>-M tau`.
Thus the pure term `(6)` is an uncancellable pole of the polynomial
`Q(x,0)`, a contradiction.

This includes `M=6`, where `(6)` is `(3/8)q_aff^2`, and `M=18`, where it
is `(-21/1024)q_aff^6`.

## 5. Uniform full-bank slope-four lemma at `P_infty`

Work on a finite `d=1` cover.  Write

```text
tau=v(h)>0,  u=v(q),  v=v(s),  b=min(u,v),
k=(M+2)/4.
```

Zero coordinates are assigned valuation `infinity`.  The claim is

```text
b<4tau  ==>  v(N)<(M+3)tau.                           (7)
```

Thus every sub-slope-four branch is response-polar.

### 5.1 Exact row data

For an even row `j=4ell+2`, put `k_j=(j+2)/4`.
Expanding on `d=1`,

```text
((1+z^2)^2+(qz^3-sz^4))^(j/4),
```

all perturbation orders below `k_j` have degree at most `j`.  At order
`k_j`, the relevant quotient is

```text
(qz^3-sz^4)^k_j/(1+z^2).
```

It follows that

```text
ord(Phi_j)=ord(Psi_j)=k_j,
ord(R_j+Phi_j/2)>=k_j+1,                              (8)
(M-j)+4k_j=M+2.                                      (9)
```

The zero polynomial in `(8)` has infinite order.  This convention matters
at `M=2`, where `R_2+Phi_2/2` is identically zero on `d=1`.

Up to a common nonzero scalar, the top faces are

```text
Phi_face =
 sum_(a odd) binom(k,a)(-1)^((a-1)/2)q^a(-s)^(k-a),

Psi_face =
 sum_(a even) binom(k,a)(-1)^(a/2)q^a(-s)^(k-a).      (10)
```

They are the odd and even parts of `(-s+iq)^k`.  If both vanished, then
both `(-s+iq)^k` and `(-s-iq)^k` would vanish, forcing `q=s=0`; hence they
are coprime.  `Psi_face` contains pure `s^k`; the pure `q^k` term lies in
`Phi_face` for odd `k` and `Psi_face` for even `k`.  The top response face
is `-Phi_face/2`.

For every odd row `j`,

```text
Phi_j(1,0,0)=c_j=4 binom(j/2,(j+1)/2)!=0,
Psi_j=O(q,s),
R_j(1,0,0)/c_j=(j+1)/(2(j+3)).                       (11)
```

The last equality follows from the adjacent-binomial ratio

```text
binom(j/2,(j+3)/2)/binom(j/2,(j+1)/2)=-1/(j+3).
```

Odd gaps `g=M-j` are distinct.  Treat the `lambda` row as a Phi-only row
of gap `M+1` and response ratio zero.  Equations `(8)`--`(9)` show that
every lower even flux and response is strictly later than the top order
`kb` whenever `b<4tau`.

### 5.2 Exhaustion of valuation cases

If `v<u`, the pure `s^k` term is the unique top `Psi` term of order `kb`.
Any odd `Psi` row early enough to cancel it has order at least
`g tau+b`; then `g tau<kb`, and that row's constant `Phi` term is uniquely
earlier in the first equation.  Otherwise the top `Psi` term is unique.
Thus no branch exists.

If `u<v` and `k` is even, the same argument uses the pure `q^k` term in
`Psi`.  If `u<v` and `k` is odd, the pure `q^k` term lies in `Phi`.
Let `g` be the least active odd/`lambda` gap.  Orders
`g tau<kb` and `g tau>kb` each leave a unique first-flux term.  At equality
the first flux may cancel, but the surviving response multiplier is

```text
-(1/2+(j+1)/(2(j+3)))=-(j+2)/(j+3)!=0,               (12)
```

or `-1/2` for the `lambda` row.  This proves `(7)`.

If `u=v=b`, coprimality in `(10)` leaves at least one top face nonzero.
If `Phi_face=0`, `Psi_face!=0`; an odd `Psi` term early enough to cancel it
again creates a strictly earlier unique constant `Phi` term.  If
`Phi_face!=0`, the first equation forces `g tau=kb`; all odd `Psi` terms
then start after `kb`, so the second equation forces `Psi_face=0`, and the
response survives by `(12)`.

If no odd coefficient and no `lambda` row is active, the same pure-face or
coprime-face argument leaves an uncancelled top flux.  The `W` term has
order `(M+2)tau=4k tau>kb`.  This covers unequal slopes, equal slopes,
one zero coordinate, and the no-active stratum.

### 5.3 The exact prefix upgrades the pole lemma to slope four

If `b<4tau`, `(7)` localizes the source point to `a`, so `(2)` applies.
At `P_infty`,

```text
beta_0^2=-1.
```

The second identity in `(2)` forces

```text
u=v=b,  s_0=beta_0 q_0.
```

For `beta_0=+i` or `-i`, exactly one of `-s+iq` and `-s-iq` vanishes.
Consequently both the odd and even faces in `(10)` are nonzero.  Cancelling
the first flux requires a least odd/`lambda` row with `g tau=kb`, but its
`Psi` term starts at `kb+b`; all lower even rows and the `W` row are also
later.  The nonzero top `Psi` term cannot cancel.  Hence no physical
sub-slope-four branch exists:

```text
min(v(q),v(s))>=4v(h),
v(q/h^3)>=v(h)>0.                                    (13)
```

Finite index covers and source ramification multiply all valuations by a
positive integer, preserving every comparison above.

## 6. Projective regularity and the vertical contradiction

Sections 2--4 exclude every boundary point except `P_infty`; Section 5
makes `q_aff=q/h^3` regular and zero at every normalization branch above
that point.  It is regular on the affine chart as a coordinate function.
A positive-dimensional complete curve cannot be contained in the affine
chart `h!=0`; otherwise its affine coordinate functions would all be
global constants.  Hence the projective image has a boundary point.

Therefore `q_aff` is a global regular function on complete integral `X`,
vanishes at a boundary point, and must be identically zero:

```text
q=A/U=0.
```

Thus `A=0` and `s=-E` is polynomial.  At `q=0`, every odd `Psi_j`
vanishes, while

```text
Psi_(4ell+2)=
4 binom(ell+1/2,ell+1)(-s)^(ell+1).
```

The normalized top row makes the second-flux equation a nonzero polynomial
over `C` in `s`; hence `s` is constant.

For odd `j`,

```text
deg_d Phi_j(d,0,s)=(j+1)/2,
[d^((j+1)/2)]Phi_j=4 binom(j/2,(j+1)/2)!=0.          (14)
```

These degrees are distinct.  If any odd coefficient is active, the
first-flux equation is a nonzero positive-degree polynomial over `C` in
`d`, so `d` is constant.  Then all response inputs are constant,
contradicting `U R_Q'=kappa`.

If no odd coefficient is active, the first-flux equation forces
`lambda=0` (and in particular a lone nonzero `lambda` gives no branch).
Every even `R_j` vanishes at `q=0` by parity, so `R_Q=0`, again
contradicting `U R_Q'=kappa`.  The nonconstant-`U` branch is empty.

## 7. Low-degree, coordinate, and scheme hostiles

### `M=2`

There is no exceptional low degree.  Directly,

```text
Phi_2=2q,  Psi_2=-2s,  R_2=-dq.
```

The tail-square, root classification, full-bank slope argument, and
vertical triangle all remain valid with `k=1`.  The only bookkeeping corner
is

```text
(R_2+Phi_2/2)|_(d=1)=0,
```

which has infinite rather than zero local order.  Independently, the
reduced quadratic-member case is already closed by THM-2202/THM-2071.  A
canonical proof may either retain `k=1` with the infinite-order convention
or cite that prior closure and start the new argument at `k=2`.

### Weighted coordinates and ramification

The `d=1` and `q=1` charts are finite index covers of weighted-projective
charts, not globally faithful affine coordinates.  This causes no loss:
all identities are weighted homogeneous, and every valuation is multiplied
by the same positive ramification index.  The proof never equates index-cover
point counts with coarse orbit counts.

### Reducible or nonreduced ambient response scheme

No integrality or reducedness of the ambient complete intersection is used.
A map from reduced `P1` kills nilpotents.  Its generic image lies in one
reduced irreducible projective curve; the proof normalizes precisely that
closure.  Unrelated components and embedded structure cannot add a boundary
branch to the actual image or alter a pulled-back valuation.

## 8. Audit of the THM-2760 dependency

The nonzero-root exclusion was independently replayed through `k=30` from
the specialized recurrence

```text
B_r(x)=(1-x^2)^2+r x^3(1-x).
```

After dividing the two top coefficients by `r^k`, their residual
polynomials were coprime in every case.  This is finite hostile evidence;
the all-degree quantifier comes from THM-2760's Schur--Szego/Jacobi proof.

I also checked the load-bearing external statement against the TeX source
of Kostov--Shapiro, *On Schur--Szego composition of polynomials*,
arXiv:math/0605377:

1. Proposition `th:mult` gives product-root multiplicity
   `m_P+m_Q-n` when this is positive, and explicitly says equality zero is
   not a root.
2. Proposition `prop:comp` preserves hyperbolicity when the second factor
   has all negative roots.
3. Theorem `th:ordpart(ii)` says every remaining `B`-root is simple.

THM-2760 proves both Jacobi factors have simple negative roots.  For
`n>=2`, input multiplicities `1+1` never exceed `n`, so there are no
`A`-roots and all output roots are simple `B`-roots.  The case `n=1` is
immediate.  Thus the cited multiplicity theorem is used with the correct
strict inequality and supports the required squarefreeness conclusion.

The root-zero classification itself is elementary and was independently
checked: at `d=s=omega=0`, the series is
`(1+qz^3)^(k-1/2)`, so both flux indices miss multiples of three exactly
for `k=2 mod 3`, while the response index is then a multiple of three with
nonzero half-binomial coefficient.

**Dependency verdict for the load-bearing subset: ACCEPT.**

## 9. Exact companion and reproduction

Files:

```text
.scratch/all_degree_split_audit/audit_all_degree_split_closure_thm2778.py
.scratch/all_degree_split_audit/audit_all_degree_split_closure_thm2778.out
```

The companion imports no proposed discovery code and uses no Python
`assert` nodes.  It checks:

- `k=1..30`, i.e. `M=2..118` in steps of four;
- the terminal recurrence multiplier and eight induction steps per degree;
- exact odd/even top faces, their gcd, pure monomials, response syzygy, and
  both exact-prefix residues;
- an independent recurrence construction of the nonzero-root residual gcd;
- every root-zero congruence and section coefficient in the range;
- every odd-row constant/response ratio and every lower-even slope alignment
  in the largest bank;
- full section reconstruction through degree `22`, including weight grading
  and pure-`q` coefficients;
- direct saturation of the common top triple in degrees `2,6,10,14,18`;
- the full `M=2` bank and the zero-polynomial/infinite-order convention.

It passes `2384` explicit gates.  Normal and optimized Python are
byte-identical and match the stored output.

Reproduction:

```bash
python3 .scratch/all_degree_split_audit/audit_all_degree_split_closure_thm2778.py
python3 -O .scratch/all_degree_split_audit/audit_all_degree_split_closure_thm2778.py
```

LF-normalized hashes:

```text
script_sha256 = e640136a4071ca5b75ab84ad02381755bb2360f6c69e8696e5e7050e129a644f
output_sha256 = 2c8e79888ba5de13d64f90f7b201bbba2f0f1d015f4f1e89a4c0930086836f91
```

## 10. Promotion gates and failure boundary

Promotion is sound if the canonical theorem:

1. states the complete chosen-sheet split polynomial exact-square-prefix
   chart and normalized top `E_M` coefficient explicitly;
2. retains every odd row, every lower even row, `lambda`, and `W`;
3. proves the tail-square theorem rather than claiming unverified uniform
   transversality;
4. treats the projective root-zero lane with the section coefficient `(6)`;
5. declares zero-polynomial local order infinite at `M=2`, or dispatches
   `M=2` to the prior quadratic closure;
6. works on the reduced integral closure of the actual generic image;
7. names THM-2760 as load-bearing (or incorporates its proof), and does not
   depend on reserved THM-2764;
8. states clearly that the result closes this terminal chart only.

No fatal counterexample, missing boundary chart, cancellation stratum,
coordinate-gauge exception, ramification exception, or component-scope
failure survived the audit.
