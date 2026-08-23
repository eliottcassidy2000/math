---
id: THM-3834
title: "One-sided second-row R2Z2 profiles do not enter the cubic pseudo-plane Darboux packet"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  In the first
  canonical r^2z^2 extension of THM-3821, the remaining one-sided top
  orientation X=0,Y!=0 admits no Darboux pair.  If the preceding first-row
  profile S is nonzero, the 10/7 tower and complete T integral leave an
  unavoidable 2d+1 origin coefficient in r4.  If S=0, the fixed-charge
  shifted 5/2 tower forces f=beta w^2 and Y=alpha e w^5; r4 uniquely fixes
  T=alpha w^5, the next two buckets force p=g=0 through half-integral origin
  resonances, r3z2 makes w constant, and r3z leaves -60alpha w^5 e^3.
  Together with THM-3821, THM-3828, and THM-3829, this proves that the entire
  displayed fixed first r^2z^2 ansatz is empty.  Higher canonical slots and
  the general planar Jacobian problem remain OPEN.
source: jc_zero_debt_lift / cubic-pseudoplane second-row profile lane, 2026-08-23
audit: >
  PENDING INDEPENDENT HOSTILE AUDIT.  The deterministic companion has 40
  active gates checking the Poisson Casimir, monic canonical reduction,
  arm and top buckets, the nonzero-S 10/7 tower, complete T integral, full
  r4 source and symbolic odd origin coefficient, six hostile origin-degree
  controls, the corrected S=0 shifted 5/2 law including its fixed-charge
  term, the exact r4 particular and homogeneous equation, both half-integral
  p/g resonances, the r3z2 constant-w equation, the constant arm, and the
  final nonzero r3z monomial.  Normal and optimized runs byte-match the
  frozen transcript.  No finite-field inference is used.
depends_on:
  - THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate
  - THM-3828-proportional-second-row-r2z2-profile-nonentry
  - THM-3829-misaligned-second-row-r2z2-10-7-profile-nonentry
related:
  - THM-3814-nodal-rz-kummer-profile-degree-gate
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
script: 04-computation/jc2_cubic_pseudoplane_one_sided_second_row_thm3834.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_one_sided_second_row_thm3834.out
script_sha256: 81263a14587121c55ba6d64acecbe220a48e2c9c709943fb127f935fc5f5288e
output_sha256: 596bc49f2c32318cfb4ffe903193a8fec1c4351bd1b9b820ffb2fda205bd5302
semantic_sha256: 5d2a2b84131e0bc6e879f13a5b8b2410540aefdafd5d37752d04d480baef94d6
hash_basis: raw LF bytes
---

# THM-3834 -- the one-sided second row is empty

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  Let `k` be
an algebraically closed field of characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles

```text
f,g,h,kappa,p,q,S,T,X,Y in k[e]
```

consider the fixed first second-row extension

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S+r^2z^2 X,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T+r^2z^2 Y.      (2)
```

If

```text
X=0,                         Y!=0,                      (3)
```

then

```text
{A,C}!=1.                                                (4)
```

Consequently, there is no Darboux pair anywhere in (2).  Indeed:

```text
X=Y=0       is THM-3821;
X!=0        is THM-3828 plus THM-3829;
X=0,Y!=0    is this theorem.                            (5)
```

This is a complete obstruction for the displayed ansatz, not for higher
canonical profiles on the surface and not for planar JC in general.

## 1. Exact one-sided buckets and the arm value

Reduce `{A,C}-1` by the monic relation

```text
z^3=r^2e+r.
```

The remainder has `z`-degree at most two and there is no quotient loss.  The
pure `z` and top buckets are

```text
[z]   =(36e^2f-24e kappa-12f+1)/2,
[r^7] =30e^2(XY'-YX').                                  (6)
```

The top bucket vanishes automatically under (3).  The arm equation is

```text
(3e^2-1)f-2e kappa+1/12=0,                              (7)
```

so in every branch

```text
f(0)=1/12.                                               (8)
```

The `r^6` bucket is

```text
[r^6]=-3e(-7eS Y'+10eY S'+2SY).                         (9)
```

We now keep its two branches separate.

## 2. The branch `S!=0`

Divide (9) only in the rational function field and compare prime residues.
At the prime `e`, and at every other irreducible prime `pi`, respectively,

```text
7 ord_e(Y)-10 ord_e(S)=2,
7 ord_pi(Y)-10 ord_pi(S)=0.                             (10)
```

The nonnegative solutions are

```text
Y=alpha e^6v^10,             S=beta e^4v^7,             (11)
```

with `alpha,beta in k*` and nonzero `v in k[e]`.  The full `r^5` bucket then
has the complete integral

```text
T=(10alpha/(7beta))e^2v^3 f
  +(2alpha/7)e^5v^10+gamma e^4v^7,          gamma in k. (12)
```

The exact companion checks (12) in the unreduced bucket.  If two solutions
differ by `W`, the homogeneous equation is `SW'-WS'=0`; hence `(W/S)'=0`
in `k(e)`, so (12) includes every polynomial solution.

Substitute (11)--(12) into `[r^4]` and cancel its nonzero common polynomial
factor `3e^3v^2/(7beta)`.  The full remaining source `P_4` is frozen in the
companion.  Its first possible origin block is

```text
30alpha f(4efv'-evf'+2fv).                              (13)
```

Every other block either contains at least `e^2v^4` or one of the later
factors

```text
e^3fv^8,        e^4fv^7v',        e^4v^8f',
e^6v^15,        e^7v^14v'.                              (14)
```

Write `v=e^d u`, where `d>=0` and `u(0)!=0`.  By (8), `f` is a unit at the
origin.  The leading coefficient of (13), at order `d`, is

```text
60alpha f(0)^2u(0)(2d+1).                               (15)
```

It cannot vanish in characteristic zero.  All terms in (14), and all
`e^2v^4` blocks involving `kappa`, occur later.  Thus `P_4` cannot be zero,
contradicting the `r^4` bucket.  This closes `S!=0`.

## 3. The shifted `5/2` tower when `S=0`

Set `S=0` before cancelling any profile.  The full `r^5` bucket is

```text
-6e(5eYf'-2efY'+2Yf).                                  (16)
```

The last term is the fixed-charge term; omitting it would give the wrong
unshifted Kummer law.  Since `f,Y` are nonzero, (16) gives

```text
5 ord_pi(f)-2 ord_pi(Y)=0               for pi!=e,
5 ord_e(f)-2 ord_e(Y)=-2.                              (17)
```

Equation (8) says `ord_e(f)=0`, so the UFD solutions of (17) are

```text
f=beta w^2,             Y=alpha e w^5,
alpha,beta in k*,       w(0)!=0.                        (18)
```

This is the shifted `5/2` tower.

After (18), the complete `r^4` equation is, up to the nonzero polynomial
factor `-6beta e w`,

```text
7eT w'-2ewT'+Tw+alpha w^5(3ew'-w)=0.                   (19)
```

It has the exact particular solution

```text
T_0=alpha w^5.                                          (20)
```

If `T=T_0+W`, the homogeneous equation is

```text
7eWw'-2ewW'+Ww=0.                                      (21)
```

Suppose `W!=0` and let `t=ord_0(W)`.  Since `w(0)!=0`, the unique first
coefficient in (21) is proportional to

```text
1-2t.                                                    (22)
```

Its only resonance is `t=1/2`, impossible for a polynomial.  Therefore
(20) is the unique polynomial solution.  No homogeneous branch was lost.

## 4. Two half-integral side-profile gates

The next two buckets are

```text
[r^4z^2]=15e(-2Yp'+pY'),
[r^4z]  =3(-10eYg'+3egY'+2Yg).                         (23)
```

If `p!=0`, put `t=ord_0(p)`.  Using `ord_0(Y)=1` from (18), the first local
coefficient in the first equation of (23) is proportional to

```text
1-2t.                                                    (24)
```

Thus `p!=0` would again require `t=1/2`, and hence

```text
p=0.                                                     (25)
```

If `g!=0`, the second equation in (23) has first local coefficient
proportional to

```text
5-10 ord_0(g).                                          (26)
```

Its only resonance is also `1/2`, so

```text
g=0.                                                     (27)
```

These are local order arguments in `k[[e]]`; they do not divide by `p(0)` or
`g(0)`.

## 5. The terminal two buckets

Substitute (18), (20), (25), and (27) into the full `[r^3z^2]` bucket.  It
reduces exactly to

```text
-10alpha e^2w^4w'.                                      (28)
```

The polynomial domain has characteristic zero and `alpha w!=0`, so (28)
forces `w'=0`.  Hence `w` is a nonzero constant.  Equations (7)--(8) then
give

```text
f=1/12,             kappa=e/8,
Y=a e,              T=a,                 a in k*.       (29)
```

Finally the full `[r^3z]` bucket, with no condition on the unused `h,q`, is

```text
-60a e^3.                                                (30)
```

This is nonzero, the final contradiction.  Thus the `S=0` branch is also
empty, proving (4).

## 6. Exact boundary and next design layer

The four cited theorems now close the complete fixed profile grammar (2):

```text
THM-3821       : X=Y=0;
THM-3828/3829  : X!=0 (aligned and misaligned target rows);
THM-3834       : X=0,Y!=0.                               (31)
```

No target-swap symmetry is assumed in the last line; its asymmetric fixed
charges are treated directly in (16)--(30).

The conclusion stops exactly at the displayed `r^2z^2` slot.  A positive
counterexample search must change the grammar by adding a higher canonical
slot, changing the carrier/arm immersion, or leaving this pseudo-plane
profile.  This theorem does not exclude any of those directions.

## 7. Reproducibility

Run

```bash
python3 04-computation/jc2_cubic_pseudoplane_one_sided_second_row_thm3834.py
python3 -O 04-computation/jc2_cubic_pseudoplane_one_sided_second_row_thm3834.py
```

Both executions must byte-match

```text
05-knowledge/results/jc2_cubic_pseudoplane_one_sided_second_row_thm3834.out
```

The companion uses exact SymPy identities.  It has no random or finite-field
step and no disabled Python `assert`.  The six small fixed-degree origin
replays are hostile controls only; the all-degree proof is the symbolic
coefficient (15).
