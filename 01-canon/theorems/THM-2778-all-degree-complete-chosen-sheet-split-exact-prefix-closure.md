---
id: THM-2778
title: "All-degree complete chosen-sheet split exact-prefix closure"
status: >
  PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  planar Keller pair in the complete chosen-sheet split polynomial exact-
  square-prefix terminal family of reduced degree M=4k-2, k>=1, is an
  automorphism.  A recurrence-tail square theorem makes P_infty the only
  common top triple zero; the exact-prefix exclusion and a uniform nonzero
  half-binomial section pole remove every other boundary point; and the full
  odd/even bank forces slope four and a vertical contradiction at P_infty.
  Constant U and reducible or nonreduced ambient response schemes are
  included.  This does not derive the chart for an arbitrary Keller pair or
  prove JC(2) or DC(2).
source: root/all-degree-split-exact-prefix-closure-2026-07-28
audit: >
  all-degree-split-hostile-audit-2026-07-28 independently reconstructed the
  recurrence/UFD tail theorem, the root-zero half-binomial obstruction, all
  odd/even slope faces, the M=2 zero-syzygy edge, projective and component
  scope, and the load-bearing Schur--Szego citation in the primary source;
  replayed 2384 exact gates for M=2..118 under normal and optimized Python
  with byte-identical stored output and no assert nodes: PASS.
external_input: >
  Kostov--Shapiro, On Schur--Szego composition of polynomials,
  arXiv:math/0605377, Propositions th:mult and prop:comp and Theorem
  th:ordpart(ii), through THM-2760
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2752-all-even-zero-first-flux-response-regularization-closure
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
related:
  - THM-2759-split-degree-six-componentwise-monicization-closure
  - THM-2764-all-degree-chosen-sheet-even-zero-flux-componentwise-closure-away-from-six-mod-twelve
  - THM-2767-split-degrees-ten-fourteen-eighteen-componentwise-monicization-closure
script: 04-computation/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.py
output: 05-knowledge/results/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.out
script_sha256: b7994f8e17675b8e61ee227cfe8d2bb7fc37503ef5bef4a8d5b35029d51486e8
output_sha256: 17ea0425ccece5411bf792d6c89074cc4c3cd1f77a3ac9d4affd0e8bf5fce680
hash_basis: LF-normalized bytes
---

# THM-2778 -- all-degree complete chosen-sheet split closure

**PROVED + CITED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Every complete chosen-sheet split polynomial exact-square-prefix terminal
family of reduced degree

```text
M=4k-2,  k>=1,
```

closes componentwise.  The proof combines THM-2760's all-degree exact-prefix gcd
with a new recurrence-tail lemma.  The latter removes the finite top-boundary
atlas entirely: the only common top zero of all three Faber observables is
`P_infty`.

The theorem does not derive the split chart for an arbitrary Keller pair,
treat a nonsplit or nonpolynomial exact prefix, or prove `JC(2)` or `DC(2)`.

## 1. Complete family and primitive dichotomy

Let `(P,Q)` be a polynomial Keller pair over `C` in the complete chosen split
sheet of the polynomial exact-square-prefix chart.  In its standard notation,

```text
P=(w^2+d)^2+(qw-s),
d=C-B^2/(4U^2),  q=A/U,  s=AB/(2U^2)-E.              (1)
```

After the exact target shear, write the full Faber bank

```text
Q=E_M+sum_(1<=j<M, 4 does not divide j) a_j E_j,      (2)
```

where all `a_j` are constants.  THM-2129/2230/2723 give constants
`lambda,W,kappa`, with `kappa!=0`, such that

```text
Phi_Q=lambda,       Psi_Q=W,       U R_Q'=kappa.      (3)
```

If `U` is constant, `(x,z)->(x,w)` is a polynomial source automorphism and
THM-2181 proves `(P,Q)` is an automorphism.  For contradiction, assume the
only other THM-2723 case:

```text
U=u_0(x-a)^m,  m>=2,
R_Q=r_0+r_1(x-a)^(1-m),  r_1!=0.                     (4)
```

Thus `R_Q` has exactly one pole, at the finite point `x=a`, where the original
polynomial coefficients `A,B,C,E` are regular.

Put `alpha=M/4=k-1/2` and define

```text
f(z)=1+2dz^2+qz^3+(d^2-s)z^4,
f(z)^alpha=sum_(n>=0)c_n z^n,

Phi_M=4c_(M+1),
Psi_M=4c_(M+2),
R_M=4c_(M+3)+2d c_(M+1).                             (5)
```

Homogenize with `wt(h,d,q,s)=(1,2,3,4)`:

```text
F_(M+1)=Phi_M+sum a_j h^(M-j)Phi_j-lambda h^(M+1),
G_(M+2)=Psi_M+sum a_j h^(M-j)Psi_j-W h^(M+2),
R_(M+3)=R_M+sum a_j h^(M-j)R_j.                      (6)
```

## 2. The recurrence tail classifies the triple top zero

Coefficient comparison in `f(f^alpha)'=alpha f'f^alpha` gives, for every
`n>=1`,

```text
n c_n=sum_(r=2)^4 p_r ((alpha+1)r-n)c_(n-r),          (7)
```

where `(p_2,p_3,p_4)=(2d,q,d^2-s)` and negative-index coefficients are zero.

Suppose

```text
Phi_M=Psi_M=R_M=0.                                   (8)
```

Then `c_(M+1)=c_(M+2)=c_(M+3)=0`.  At `n=M+4`, the only recurrence term
which might see `c_M` has coefficient

```text
4(alpha+1)-(M+4)=0,                                  (9)
```

so `(7)` gives `c_(M+4)=0`.  Thereafter every input index in `(7)` is at
least `M+1`; induction gives

```text
c_n=0 for every n>M.                                 (10)
```

Hence `f^alpha` is a polynomial.  Squaring yields

```text
(f^alpha)^2=f^(2k-1).                                (11)
```

In the UFD `C[z]`, the exponent `2k-1` is odd.  Since the right side of
`(11)` is a square, every irreducible multiplicity of `f` is even; therefore
`f=g^2`.  Because `f(0)=1` and its linear coefficient is zero, after changing
the sign of `g` one has

```text
g=1+bz^2,
f=(1+bz^2)^2.
```

Comparison with `(5)` gives `b=d`, `q=0`, and `s=0`.  In weighted projective
space the only nonzero such point is

```text
P_infty=[d:q:s]=[1:0:0].                             (12)
```

Conversely `(12)` makes all three top observables zero.  Thus, set-theoretically,

```text
V(Phi_M,Psi_M,R_M)={P_infty}                         (13)
```

for every `M=4k-2`.  No smoothness, Jacobian, irreducibility, or finite
Groebner atlas is used.

## 3. Every other boundary point is response-polar and exact-prefix impossible

Let `X` be the normalization of the reduced integral closure of the generic
coefficient image in `(6)`.  The image is nonconstant, because otherwise
`R_Q` would be constant, contradicting `(3)`.  Properness extends the source
map to a nonconstant morphism

```text
gamma:P1_x -> X,
```

which is finite and surjective.  Hence every boundary branch on `X` has a
physical source preimage.  A boundary point satisfies `Phi_M=Psi_M=0`.
By `(13)`, every boundary point other than `P_infty` has `R_M!=0`; therefore
`R_(M+3)/h^(M+3)` is polar on every normalization branch above it.  Equation
`(4)` places every source preimage of such a branch at the finite point `a`.

On a finite weighted-index cover, put

```text
beta=hB/(2U).
```

The polynomial exact prefix gives

```text
beta^2+d=h^2C,       q beta-s=h^4E.                  (14)
```

Since the right sides and `d,q,s` are DVR-regular at the finite point `a`,
`beta` is regular.  If `omega` is its residue, `(14)` gives

```text
d=-omega^2,          s=q omega.                      (15)
```

If `q omega!=0`, THM-2760 proves that `Phi_M` and `Psi_M` cannot both vanish.
If `q=0`, then `(15)` either gives `P_infty` (`omega!=0`) or the excluded
all-zero homogeneous tuple.  The only remaining case is

```text
omega=0, q!=0.                                       (16)
```

THM-2760 proves that `(16)` is a common-flux point exactly when

```text
M=6 mod12.                                           (17)
```

It is response-active, consistently with `(13)`.

## 4. One half-binomial kills every root-zero survivor

Work on the `q=1` cover at `(16)` and put `v(h)=eta>0`.  Equations `(14)` and
regularity of `C` imply `v(beta)>0`.  In affine variables,

```text
beta_aff=beta/h=B/(2U),       q_aff=q/h^3=A/U,

v(beta_aff)>-eta,             v(q_aff)=-3eta.         (18)
```

At the original polynomial section `z_source=0`, every Faber value `E_j(0)`
is weighted-homogeneous of weight `j` in

```text
wt(beta_aff,C,E,q_aff)=(1,2,4,3).                    (19)
```

At `beta_aff=C=E=0`, coefficient extraction from
`(1+q_aff z^3)^(M/4)` gives

```text
[q_aff^(M/3)] E_M(0)=binom(M/4,M/3).                 (20)
```

Under `(17)`, `M/3` is an integer and `M/4` is a half-integer, so `(20)` is
nonzero.  Its pure `q_aff^(M/3)` monomial has valuation exactly `-M eta`.
Every other weight-`M` monomial is strictly shallower: a positive
`beta_aff` exponent is strict by `(18)`, while a positive `C` or `E` exponent
forces three times the `q_aff` exponent to be less than `M`.  Every monomial
of every lower `E_j(0)` has weight `j<M` and hence valuation greater than
`-M eta`.

Thus the normalized top row supplies a unique pole to the target-sheared
polynomial `Q(x,0)`, a contradiction.  All boundary points other than
`P_infty` are excluded in every degree.

## 5. Uniform full-bank slope-four lemma at `P_infty`

Use a finite `d=1` cover at `(12)`.  Write

```text
eta=v(h)>0,  u=v(q),  v_0=v(s),  b=min(u,v_0),
k=(M+2)/4.                                            (21)
```

Allow zero coordinates by assigning valuation `infinity`.  We prove

```text
b<4eta  ==>  R_(M+3)/h^(M+3) is polar.               (22)
```

Four all-degree identities drive the proof.

1. For every even row `j=4ell+2`, set `k_j=(j+2)/4`.  THM-2752 gives

   ```text
   ord(Phi_j)=ord(Psi_j)=ord(R_j)=k_j,
   ord(R_j+Phi_j/2)>=k_j+1,
   (M-j)+4k_j=M+2.                                   (23)
   ```

   Therefore every lower even flux and response row is strictly later than
   top order `kb` whenever `b<4eta`.

2. Up to a common nonzero scalar, the top leading faces are the odd and even
   parts of `(-s+iq)^k`.  They are coprime.  The `Psi` face contains pure
   `s^k`; pure `q^k` lies in `Psi` for even `k` and in `Phi` for odd `k`.

3. For every odd row,

   ```text
   Phi_j(1,0,0)=c_j!=0,       Psi_j=O(q,s),
   R_j(1,0,0)/c_j=(j+1)/(2(j+3)).                    (24)
   ```

   Indeed `f=(1+z^2)^2` at `(1,0,0)`, so these are consecutive
   half-binomial ratios.  The gaps `M-j` are distinct.  Treat the `lambda`
   row as one additional Phi-only gap `M+1` with response ratio zero.

4. THM-2752 also gives

   ```text
   ord(R_M+Phi_M/2)>=k+1,                             (25)
   ```

   so the top response/first-flux ratio on the leading face is `-1/2`.
   At `M=2` the left side is identically zero, which is interpreted as
   infinite local order; exact certificates must not assign the zero
   polynomial the numerical order zero.

If `s` has smaller valuation, pure `s^k` is the unique top `Psi` term at
order `kb`.  Any odd `Psi` row early enough to cancel it has its constant
`Phi` term strictly before `kb`, where the least gap is unique.  The same
argument applies when `q` has smaller valuation and `k` is even.

If `q` has smaller valuation and `k` is odd, pure `q^k` is in `Phi`.  Let `g`
be the least active odd/`lambda` gap.  Unless `g eta=kb`, either that row or
the top `Phi` term is uniquely earliest.  At equality, first-flux
cancellation leaves the response multiplier

```text
-(1/2+(j+1)/(2(j+3)))=-(j+2)/(j+3)!=0,               (26)
```

or `-1/2` for the `lambda` row.  Lower even and later odd response rows
cannot reach that order, proving `(22)`.

Finally suppose `u=v_0=b`.  If the top `Phi` face vanishes, coprimality makes
the top `Psi` face nonzero; any odd `Psi` competitor would put its constant
`Phi` term strictly earlier.  If the top `Phi` face is nonzero, its equation
forces `g eta=kb`; then the top `Psi` face must vanish, and `(26)` again
leaves a response pole.  With no active odd/`lambda` row, coprimality or the
appropriate pure top monomial makes one flux uniquely earliest.  This proves
`(22)` in all coefficient and zero-coordinate strata.

Every pole in `(22)` pulls back to `x=a`, so `(14)` is regular there.  At
`P_infty`, its residues give

```text
omega^2=-1.
```

If `b<4eta`, the second identity forces

```text
u=v_0=b,       s_0=omega q_0.                        (27)
```

For `omega=+i` or `-i`, one of `-s+iq,-s-iq` vanishes and the other does not;
therefore both odd/even top faces are nonzero.  The first flux can cancel only
against a least odd/`lambda` row with `g eta=kb`.  Its `Psi` row is one order
later, all lower even rows are later by `(23)`, and `W h^(M+2)` is later.
The nonzero top `Psi` term is unique, a contradiction.  With no active
odd/`lambda` row, the nonzero top `Phi` term is already unique.  Hence every
physical branch at `P_infty` satisfies

```text
min(v(q),v(s))>=4v(h),
v(q/h^3)>=v(h)>0.                                    (28)
```

Finite index-cover ramification multiplies every valuation by one positive
integer, so `(22),(28)` descend to the coarse normalization.

## 6. Projective regularity and the vertical triangle

All non-`P_infty` boundary points are excluded by Sections 3--4.  Equation
`(28)` makes `q_aff=q/h^3` regular and zero at every remaining boundary
branch; it is regular on the affine chart.  Every nonconstant projective image
curve meets `h=0`.  Thus `q_aff` is a global regular function on a complete
integral curve and vanishes somewhere, so

```text
q_aff=0.                                              (29)
```

Consequently `A=0` and `s=-E in C[x]`.  At `q=0`, every odd `Psi_j` vanishes,
while every even row `j=4ell+2` is a nonzero scalar multiple of
`s^(ell+1)`.  The normalized top coefficient makes `Psi_Q=W` a nonzero
polynomial equation over `C`; hence `s in C`.

For odd `j`, `Phi_j(d,0,s)` has exact `d`-degree `(j+1)/2` with nonzero
leading coefficient.  The complete odd bank gives the distinct degrees
`1,...,M/2`.  Therefore `Phi_Q=lambda` either forces `d in C`, making every
response input constant and contradicting `(3)`, or is the zero polynomial.
The latter occurs only when all odd `a_j` and `lambda` vanish.  Then only even
response rows remain, and THM-2752's `q^3` divisibility gives `R_Q=0`, again
contradicting `(3)`.

The nonconstant-`U` alternative is empty.  The constant-`U` alternative is an
automorphism by THM-2181.  This closes the complete chosen-sheet split
polynomial exact-square-prefix family for every reduced degree `M=4k-2`.

## 7. Scope and hostile boundaries

1. THM-2760 is load-bearing only for the exact-prefix nonzero-root exclusion
   and root-zero congruence.  Its Schur--Szego multiplicity input was
   independently checked against the primary TeX source: product-root
   multiplicity uses the strict positive threshold, hyperbolicity is
   preserved, and every remaining `B`-root is simple.  Thus two simple-
   negative degree-`n` inputs produce no `A`-root for `n>=2` and only simple
   output roots; `n=1` is immediate.
2. The new tail lemma proves only set-theoretic triple-response support.  That
   is enough: a nonzero top response is a DVR unit on every branch, so no
   transversality or reduced top scheme is needed.
3. The root-zero coefficient is the nonzero half-binomial `(20)`, not a finite
   pattern inferred from degrees `6` and `18`.
4. The full odd bank and `lambda` row are essential at `P_infty`; this is
   stronger than the all-even zero-flux reservation THM-2764.
5. The proof works on the normalization of the reduced integral generic
   image.  A reduced source kills nilpotents; embedded or unrelated ambient
   components do not enter.
6. The conclusion is chart-relative.  Nonsplit/nonpolynomial-prefix families,
   chart entry, degree raising/descent, and other Newton/Jelonek/source-fibre
   branches remain open; `JC(2)` and `DC(2)` remain open.

## 8. Independent exact audit and reproduction

The companion reconstructs the recurrence coefficients without importing
the discovery derivation.  It checks `k=1,...,30` (`M=2,...,118`), including
the terminal multiplier and induction tail, coprime odd/even top faces,
response syzygies, both exact-prefix residues, the nonzero-root residual gcd,
root-zero congruence and section coefficient, every odd-row ratio and lower-
even alignment, full section reconstruction through degree `22`, direct
saturation of the top triple through degree `18`, and the `M=2` infinite-
order convention.  These `2384` finite gates are hostile controls; the
all-degree quantifier comes from the recurrence/UFD proof and THM-2760, not
from extrapolation.

Reproduce from the repository root:

```bash
python3 04-computation/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.py
python3 -O 04-computation/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.py
diff -u \
  05-knowledge/results/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.out \
  <(python3 04-computation/jc2_all_degree_complete_split_exact_prefix_closure_thm2778.py)
```
