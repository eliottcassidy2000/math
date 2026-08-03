---
id: THM-3214
title: "Two-jet pseudo-division locality and Catalan sharpness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The fraction-free reciprocal pseudo-remainder operator consumes exactly two
  coefficient jets per step.  Its kth pivot factors through the initial
  2k-jet and its kth connection through the initial (2k+1)-jet.  Neither
  bound can be lowered: a one-coefficient Catalan deformation preserves every
  earlier input jet and every earlier pivot while changing the indicated
  pivot or connection.  Thus floor(s/2) is the sharp universal
  pseudo-division information budget, although simultaneous nonvanishing
  remains open.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  The exact companion pins THM-3192, checks the universal coefficient formula,
  63 deterministic locality pairs, the Catalan closed orbit through eight
  steps, and sharp pivot/connection deformations through ten steps.  It also
  verifies the leading-order perturbation recurrence coefficient by
  coefficient.  An independent audit rederived the ratio formula, locality
  induction, Catalan orbit, and exact terminal perturbations; it also detected
  and repaired the raw-chain scaling omitted from the first version of `(25)`.
  Normal and `-O` replay byte-match the stored output and both declared
  LF-normalized hashes.
depends_on:
  - THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return
script: 04-computation/factorial_two_jet_locality_catalan_sharpness_thm3214.py
output: 05-knowledge/results/factorial_two_jet_locality_catalan_sharpness_thm3214.out
script_sha256: aec62d92ec6e1147994946fcd135b4ecdbdb2f1f2aa264ebd8a62072d5d8229a
output_sha256: 75d9feaa94b35ed7aba624ff0f3aa88583fe2e0b3d6a9e6c57b16476fc4440c8
hash_basis: LF-normalized bytes
---

# THM-3214 -- two-jet pseudo-division locality and Catalan sharpness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3192 identifies the offset-six walls `H,J,K` with reciprocal-jet
pseudo-division coordinates.  The present theorem removes the offset and the
factorial moments from the local question.  It proves that two new input
coefficients per remainder are sufficient over every commutative coefficient
ring and necessary already on one explicit regular rational orbit.

## 1. The universal two-jet operator

Let `A` be a commutative ring and write

```text
f=sum_(i>=0) f_i z^i,                 g=sum_(i>=0) g_i z^i.       (1)
```

Put

```text
P_1(f,g)=f_0g_1-f_1g_0,                                      (2)

E(f,g)=z^(-2){f_0^2g-[f_0g_0+P_1(f,g)z]f}.                   (3)
```

The coefficients of degrees zero and one inside braces vanish identically,
so `(3)` defines an element of `A[[z]]` without any division in `A`.  More
precisely,

```text
[z^j]E(f,g)
 =f_0^2 g_(j+2)-f_0g_0 f_(j+2)-P_1(f,g)f_(j+1)       (j>=0).   (4)
```

Consequently the `h`-jet of `E(f,g)` is a universal polynomial in the
`(h+2)`-jets of `f,g`:

```text
Jet_(h+2)(f,g)  --->  Jet_h(E(f,g)).                          (5)
```

This remains valid on every pivot wall.  No leading coefficient has been
inverted.

When `f_0` is a unit, let

```text
T_2(u)=z^(-2)(u-u_0-u_1z).
```

Then `(3)` also has the ratio form

```text
E(f,g)=f_0^2 f T_2(g/f).                                    (6)
```

Thus one fraction-free Euclidean step removes exactly the constant and linear
parts of the projective ratio `g/f`.

## 2. Arbitrary iterates and the sharp upper budget

Start from any pair `R_(-1)=G`, `R_0=F` in `A[[z]]` and iterate

```text
R_(k+1)=E(R_k,R_(k-1))                         (k>=0).        (7)
```

Define the kth pivot and connection coordinate by

```text
rho_k=[z^0]R_k,                  chi_k=P_1(R_k,R_(k-1)).      (8)
```

Repeated use of `(5)` gives, for every `h,k>=0`, a universal polynomial map

```text
Jet_(h+2k)(F,G)  --->  Jet_h(R_k).                           (9)
```

Indeed, the assertion is immediate for `k=0`.  If it holds through `k`, then
`Jet_h(R_(k+1))` uses `Jet_(h+2)(R_k,R_(k-1))`; the first entry uses initial
jets through `h+2+2k`, and the second uses no more.  This proves `(9)` by
induction.

Taking `h=0`, and then using the one-jets in `(2)`, yields

```text
rho_k depends only on Jet_(2k)(F,G),
chi_k depends only on Jet_(2k+1)(F,G).                       (10)
```

These statements concern information locality, not nonvanishing.  They are
valid even if one or more preceding pivots vanish.

## 3. A regular Catalan orbit

Work now over `Q[[z]]` and let `C` be the Catalan series

```text
C=1+zC^2=1+z+2z^2+5z^3+14z^4+... .                         (11)
```

Use the initial pair

```text
R_(-1)=1,                         R_0=C.                     (12)
```

Set `epsilon_(-1)=epsilon_0=1` and

```text
epsilon_(j+1)=-epsilon_(j-1).                               (13)
```

Then every iterate has the closed form

```text
R_j=epsilon_j C^(2j+1)                         (j>=0).       (14)
```

For the first step, `C^(-1)=1-zC`, so

```text
T_2(C^(-1))=-C^2,                 E(C,1)=-C^3.              (15)
```

For every later step,

```text
C^(-2)=(1-zC)^2,
T_2(C^(-2))=-C^2,
E(C^m,C^(m-2))=-C^(m+2).                                  (16)
```

The scaling law

```text
E(af,bg)=a^2b E(f,g)                                       (17)
```

then proves `(13),(14)`.  In particular every pivot on this orbit is the unit
`epsilon_j`; there is no hidden singular step in the sharpness examples below.

## 4. The 2k pivot bound cannot be lowered

Fix `k>=1` and deform only the coefficient which is invisible to
`Jet_(2k-1)`:

```text
C_t=C+t z^(2k).                                             (18)
```

Let `R_j(t)` be the chain `(7)` starting from `(1,C_t)`.  For `0<=j<=k`,
the first changed term is

```text
R_j(t)-R_j(0)
 =epsilon_j t z^(2(k-j))+O(z^(2(k-j)+1)).                  (19)
```

To prove `(19)`, start with coefficient one at order `2k`.  Before the last
step the perturbations do not change any constant or linear input term.  If
the first changed term of `R_j` is `lambda_j t z^a`, then `(4)` shows that the
unique contribution at order `a-2` in `R_(j+1)` is

```text
lambda_(j+1)=-epsilon_j epsilon_(j-1) lambda_j.             (20)
```

Since `lambda_0=1` and `(13)` gives
`-epsilon_j epsilon_(j-1)epsilon_j=epsilon_(j+1)`, induction
proves `(19)`.  There are no nonlinear corrections at the displayed order:
through the last required step all perturbed constants and linear terms used
in `(4)` are still unchanged.

At `j=k`, `(19)` gives the exact pivot identity

```text
rho_k(C_t)=epsilon_k(1+t).                                  (21)
```

Thus `C` and `C+z^(2k)` have identical `(2k-1)`-jets and identical pivots
through `k-1`, while their kth pivots are `epsilon_k` and `2epsilon_k`.
Therefore no universal map through `Jet_(2k-1)` can determine `rho_k`.

## 5. The 2k+1 connection bound cannot be lowered

For `k>=1`, instead put

```text
C_t=C+t z^(2k+1).                                          (22)
```

The same induction gives

```text
R_j(t)-R_j(0)
 =epsilon_j t z^(2(k-j)+1)+O(z^(2(k-j)+2))       (0<=j<=k). (23)
```

Hence `R_k` changes first in its linear coefficient, `R_(k-1)` changes first
in degree three, and all their constants are fixed.  Formula `(2)` gives

```text
chi_k(C_t)-chi_k(C)
 =-epsilon_k epsilon_(k-1)t !=0.                            (24)
```

For `k=0`, the same conclusion is immediate from
`P_1(C+t z,1)=-1-t`.  Therefore `chi_k` is not determined by the initial
`2k`-jet for any `k>=0`.

The two sharpness families are regular on both sides: taking `t=1` leaves all
displayed pivots nonzero.  The failure is genuine missing input information,
not division by a vanishing earlier pivot.

## 6. Exact explanation of the offset-six budget

Apply `(10)` to THM-3192's reciprocal top jets.  For the raw universal chain
`mathcal R_(-1)=b`, `mathcal R_0=a`, its exact truncated identities and chart
units give

```text
Jet_3 E(a,b)=U_H Jet_3(r),       Jet_1 E(r,a)=U_J Jet_1(s),

rho_1=P_2(a,b)=U_H r_3                         ~_p H,
rho_2=U_H^2 P_2(r,a)=U_H^2 U_J s_2             ~_p J,
chi_2=U_H^3 U_J P_1(s,r)=U_H^3 U_J U_K K       ~_p K.      (25)
```

Here `~_p` denotes equality up to a unit of `Z_(p)` for `p>=197`.  The extra
powers of `U_H,U_J` record the fraction-free scaling of the raw chain:
`E(cf,g)=c^2E(f,g)` and `P_1(cf,dg)=cdP_1(f,g)`.  All displayed `U` factors
are p-units, so they change neither the chart ideals nor the offset-six
heights `2,4,5`.  Those heights are forced by the universal filtered operator;
they are not accidental sizes in the symbolic formulas.  More generally, an
available top jet of order `s` can determine at most the first `floor(s/2)`
pivots, and the Catalan family proves this ceiling is exact for the universal
pseudo-division problem, not specifically inside the factorial-moment family.

What remains open is geometric rather than local: a sufficient jet does not
force its selected pivot to be nonzero, and `(9)` does not choose a surviving
chart when several Pluecker coordinates cancel together.  Hence this theorem
does not by itself prove arbitrary-offset separation, `NC(2)`, or `GMC(2)`.

## 7. Exact evidence

Run

```text
python 04-computation/factorial_two_jet_locality_catalan_sharpness_thm3214.py
python -O 04-computation/factorial_two_jet_locality_catalan_sharpness_thm3214.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer arithmetic only.  It pins the promoted THM-3192 blob; checks `(4)` on
deterministic coefficient tables; checks `(9)` on 63 pairs of inputs with
identical required jets and adversarial distinct tails; verifies `(14)`
through eight steps; and verifies `(19)--(24)` through ten steps, including
unchanged earlier pivots.  No floating point, randomness, imported executable,
or assertion-sensitive test is used.

**QED.**
