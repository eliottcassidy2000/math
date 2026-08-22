---
id: THM-3583
title: "Universal exponent-two two-by-four weight Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  For every squarefree polynomial
  Sigma of degree at least two, a polynomial Darboux pair on
  c^2 e=Sigma(b) cannot have two nonconstant weight pieces in one output and
  at most four in the other.  Thus a two-piece output forces at least five
  pieces in its mate and at least seven pieces in total.  This is not a
  universal seven-piece floor: the general six-piece 3 x 3 cell remains
  open, although THM-3579 closes its equal-step subclass.  No polynomial
  Darboux pair and no counterexample to JC(2) is claimed.
source: root / delegated Darboux hostile lane, 2026-08-21
audit: >
  The support-component graph, three four-node placements, six R=2/T=2
  ladders, weight-zero and arm-order boundaries, square-homogeneous branch,
  hidden gcd-five extraction, and hidden gcd-four square-root normalization
  were derived independently and checked by exact symbolic identities.
  Ordinary and optimized companion runs agree byte for byte.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
related:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3579-equal-step-three-by-three-danielewski-darboux-nonentry
  - THM-3581-critical-value-multiarm-keller-compiler-and-A13-carrier
script: 04-computation/jc2_danielewski_two_by_four_weight_nonentry_thm3583.py
output: 05-knowledge/results/jc2_danielewski_two_by_four_weight_nonentry_thm3583.out
script_sha256: 2bfd82141735ebacf1a3bd6290dc80aba3a3a13e9649a9f7fd4b77655fe05157
output_sha256: adee46c27ff115f51eef502f2ac4e385e80deebd20582c85ee7a7d344a1b7a82
hash_basis: LF-normalized bytes
---

# THM-3583 -- universal exponent-two two-by-four weight Darboux nonentry

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**  This closes the `2 x 4`
cell left open by THM-3569.  The conclusion is universal in the squarefree
target polynomial, but only for boundary exponent two and only in the sector
where one Darboux output has two weights.

All rings are over `C`.  Let `Sigma in C[b]` be squarefree of degree at least
two and put

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)).                         (1)
```

Use the Poisson bracket and grading

```text
{b,c}=c^2,       {c,e}=-Sigma'(b),       {b,e}=-2ce,
wt(b,c,e)=(0,1,-2).                                      (2)
```

Every weight-`r` piece is uniquely `c^r f(b)` on `c!=0`, with

```text
Sigma^ceil(-r/2) divides f             when r<0,          (3)
```

and

```text
{c^r f,c^s g}=c^(r+s+1) W_(r,s)(f,g),
W_(r,s)(f,g)=s f'g-r f g'.                               (4)
```

As in THM-3569, subtract scalar weight-zero pieces before counting support.
The theorem says that

```text
{P,Q}=1,       #supp(P)=2       implies       #supp(Q)>=5. (5)
```

The statement is symmetric under exchanging `P,Q`.

## 1. Two elementary tools

If `r,s` are nonzero and have the same sign, unique factorization gives

```text
W_(r,s)(f,g)=0
  iff
f=A h^(|r|/d),       g=B h^(|s|/d),
d=gcd(|r|,|s|),                                      (6)
```

up to nonzero constants.  If `r,s` have opposite signs, the logarithmic
derivative equation would make `f^|s|g^|r|` constant, so no negative-weight
coefficient satisfying `(3)` can occur.  If exactly one weight is zero, the
weight-zero coefficient is forced to be constant and is removable.  A
both-zero endpoint is tautological; in the four-node scalar component below,
the only such retained endpoint is the upper endpoint of `L` at `R=T=1`.
It is kept through `(11)--(15)` and then killed by the arm gate `(16)`.
We call `(6)` together with these zero-weight boundaries the common-power
lemma.

We will also repeatedly obtain an Euler factor

```text
E_(h,K)=K h'+2h K'.                                     (7)
```

Because `Sigma|h`, `deg h>=2`.  For nonzero `K`, the leading coefficient of
`E_(h,K)` is

```text
(deg h+2deg K) lc(h)lc(K),                              (8)
```

so `E_(h,K)` has positive degree and cannot divide the unit `1`.

THM-3569 supplies the already proved homogeneous, one-by-arbitrary,
two-by-two, and two-by-three exclusions.  We use those smaller cells only
after explicitly deleting a bracket-zero support component.

## 2. Exhaustive support-component classification

Let the two weights of `P` be `r_0<r_1` and put

```text
delta=r_1-r_0>0.                                       (9)
```

The only `Q` weights that can contribute to the scalar row are

```text
alpha=-r_1-1,          beta=-r_0-1=alpha+delta.         (10)
```

Both must occur: if only one occurred, its homogeneous bracket would itself
equal one, contrary to THM-3569.

Join two `Q` weights by an edge when their difference is `delta`.  Bracket
rows belonging to different connected components of this graph never meet.
A component not containing `alpha,beta` contributes no scalar row, hence its
whole bracket with `P` is zero and the component can be deleted.  The
component containing `alpha,beta` cannot have two or three vertices, by the
smaller-cell theorem.  With at most four `Q` pieces it must therefore have
exactly four consecutive vertices.

The lowest singleton row pairs `r_0` with the least `Q` weight, and the
highest singleton row pairs `r_1` with the greatest.  The common-power lemma
forces the lower pair to be negative and the upper pair to be nonnegative.
Thus, for positive integers `R,T`, write

```text
supp(P)={-R,T-1},          delta=R+T-1.                 (11)
```

Put

```text
s_j=-T+j delta.                                        (12)
```

The scalar complements are `s_0=-T` and `s_1=R-1`.  There are exactly three
ways to place them in a four-vertex path:

| placement | indices | `Q` weights |
|---|---|---|
| double lower `L` | `-2,-1,0,1` | `-2R-3T+2, -R-2T+1, -T, R-1` |
| balanced `B` | `-1,0,1,2` | `-R-2T+1, -T, R-1, 2R+T-2` |
| double upper `U` | `0,1,2,3` | `-T, R-1, 2R+T-2, 3R+2T-3` |

Writing `q_j` for the coefficient of the piece of weight `s_j`, every
internal row has the single uniform equation

```text
W_(T-1,s_j)(F,q_j)+W_(-R,s_(j+1))(f,q_(j+1))=delta_(j,0).
                                                               (13)
```

The two path endpoints give

```text
W_(-R,s_m)(f,q_m)=0,
W_(T-1,s_(m+3))(F,q_(m+3))=0,                         (14)
```

where `m=-2,-1,0` in `L,B,U` respectively.  Equations `(11)--(14)` are the
complete support classification; no cancellation row has been dropped.

At a simple root of `Sigma`, the scalar row is

```text
W_(T-1,-T)(F,q_0)+W_(-R,R-1)(f,q_1)=1.               (15)
```

The first summand can be nonzero there only when `T=2`: for `T=1` its
derivative multiplier is zero, and for `T>=3` the coefficient `q_0` has arm
order at least two.  Symmetrically, the second can be nonzero only when
`R=2`.  Therefore every survivor lies on one of the six ladders

```text
R=2 with L,B,U,                 T=2 with L,B,U.         (16)
```

We now close all six, including their common-power side branches.

## 3. The three `R=2` ladders

### 3.1 Balanced placement

Here

```text
supp(P)={-2,T-1},
supp(Q)={-2T-1,-T,1,T+2}.                              (17)
```

The lower endpoint and its adjacent row integrate to

```text
f=A h^2,                 q_(-1)=B h^(2T+1),
q_0=D h^T+((2T+1)B/(2A)) h^(2T-1)F,       Sigma|h.    (18)
```

For `T>=2`, both summands of `q_0` have arm order at least two.  Each term in
the scalar row `(15)` is then divisible by `h`, a contradiction.  For
`T=1`, the upper endpoint is

```text
W_(0,3)(F,q_2)=3F'q_2=0.                              (19)
```

Thus `F` is constant and its weight-zero piece is removable, reducing to an
already excluded smaller cell.  The balanced ladder is empty.

### 3.2 Double-lower placement

The weights are

```text
supp(Q)={-3T-2,-2T-1,-T,1}.                           (20)
```

The lower endpoint has weights `-2,-3T-2`.  If `T` is odd, their gcd is one
and `(6)` gives `f=A h^2`; neither scalar summand can survive an arm, because
`T!=2`.  Hence `T=2k` is even.  Put

```text
r=2k-1,       u=3k+1.
```

The two endpoints and two internal zero rows give, on the main branch,

```text
f=A h,                    q_(-2)=B h^u,
F=L K^r,                  q_1=M K,
C=uB/A,                   mu=3kC/(2A),
q_(-1)=C h^(3k)F,
q_0=D h^k+mu h^(3k-1)F^2,                 Sigma|h.    (21)
```

Direct substitution factors the scalar row exactly as

```text
(K h'+2h K')
[ AM
  -k r L D h^(k-1)K^(r-1)
  -mu(3k-1)rL^3 h^(3k-2)K^(3r-1) ].                  (22)
```

The Euler factor has positive degree by `(8)`.

There is one possible homogeneous addition to `q_(-1)`:

```text
C_0 h^((4k+1)/2).                                     (23)
```

It is polynomial only when `h=H^2`.  For `k>=2`, both scalar coefficients
then have arm order at least two.  At `k=1`, the complete exceptional branch
is

```text
f=A H^2,                 q_(-2)=B H^8,
q_(-1)=(4B/A)H^6F+C_0H^5,
q_0=DH^2+(6B/A^2)H^4F^2+(5C_0/(2A))H^3F,
q_1=MF,                                                     (24)
```

and its scalar row is

```text
[H(HF)'/(2A^2)]
[4A^3M-4A^2D-15AC_0HF-48BH^2F^2].                     (25)
```

Again a nonconstant factor remains.  This closes every odd, even, and
square-homogeneous boundary of the double-lower ladder.

### 3.3 Double-upper placement

Now

```text
supp(Q)={-T,1,T+2,2T+3}.                              (26)
```

The lower endpoint is `W_(-2,-T)(f,q_0)=0`.  Odd `T` makes `f` a square and
again kills both scalar modes.  Thus `T=2k`, and

```text
f=A h,                 q_0=B h^k,       Sigma|h.      (27)
```

Put

```text
r=2k-1,       u=4k+3,       v=2k+2,       w=2k+4.    (28)
```

The upper endpoint has `d=gcd(r,u) in {1,5}`.  On the main branch, or after
the extraction described below, write

```text
F=L K^r,                     q_3=M K^u,
C=AMu/(Lr),                  q_2=D K^v+C hK^w,
lambda=ADv/(Lr),             mu=ACw/(2Lr),
q_1=E K+lambda hK^3+mu h^2K^5.                        (29)
```

The scalar row is exactly

```text
(K h'+2h K')
[ -krLB h^(k-1)K^(r-1)
  +AE+3A lambda hK^2+5A mu h^2K^4 ].                  (30)
```

It cannot equal one.

The hidden boundary is `d=5`, equivalently `k=3 mod 5`.  Before extraction,
the endpoint has

```text
F=L J^(r/5),                q_3=M J^(u/5).             (31)
```

Here the forced parts of the descending recurrence have the forms

```text
q_2=C h J^((2k+4)/5),       q_1=mu h^2J.              (32)
```

The `q_2` homogeneous exponent `(2k+2)/5` and the `q_1` homogeneous exponent
`1/5` are nonintegral.  Since `k>=3`, `(27)` and `(32)` make both scalar
summands vanish at an arm unless the `q_1` homogeneous mode is nonzero.
Unique factorization then forces `J=K^5`, converting `(31)` to `(29)` and
returning to factor `(30)`.  Thus the gcd-five branch supplies no escape.

## 4. The three `T=2` ladders

The overlap `R=T=2` was already closed above.  Assume `R>2`.  At every arm,
the `f`-summand in `(15)` vanishes, so the coefficient `q_0` must be simple
and `F` must be a unit at that arm.

### 4.1 Balanced placement

The weights are

```text
supp(Q)={-R-3,-2,R-1,2R}.                              (33)
```

The lower endpoint gives

```text
f=A h^(R/d),          q_(-1)=B h^((R+3)/d),
d=gcd(R,3).                                             (34)
```

At an arm let `ell=ord(h)`,

```text
m=R ell/d,             n=(R+3)ell/d.                  (35)
```

In the adjacent row, `W_(1,-R-3)(F,q_(-1))` has exact order `n-1`, because
`F` is a unit.  The other term, `W_(-R,-2)(f,q_0)`, has exact order `m` and
leading coefficient `R-2m`.  This coefficient cannot vanish: `d` is odd and
`m=R ell/d`.

Cancellation would require `n-1=m`, hence `3ell/d=1`.  The only integral
possibility is `d=3,ell=1`, but then `m=R/3`, contradicting the membership
bound `m>=ceil(R/2)`.  The balanced ladder is empty.

### 4.2 Double-lower placement

Here

```text
supp(Q)={-2R-4,-R-3,-2,R-1}.                          (36)
```

Set `d=gcd(R,4)`.  The lower endpoint and the two successive first-order
equations have the following common-power skeleton:

```text
f=A z^(R/d),               q_(-2)=B z^((2R+4)/d),
q_(-1)=C z^((R+4)/d)F + homogeneous,
q_0=lambda z^(4/d)F^2+D z^(2/d)+induced homogeneous.  (37)
```

At an arm let `ell=ord(z)`.  Membership gives

```text
(R/d)ell>=ceil(R/2).                                  (38)
```

The displayed particular term of `q_0` has order `4ell/d`; a homogeneous
addition induced from `q_(-1)` has order `3ell/d`.  The former could be
simple only for `d=4,ell=1`, which violates `(38)`, while the latter cannot
be simple because `d` divides four.  Thus the only possible simple term is
`D z^(2/d)`.  Therefore `2ell/d=1`.  It follows that `d` is even and
`R=2k`.  If `d=2`, put `h=z`; if `d=4`, polynomiality forces `z=h^2`.
This is the hidden gcd-four square-root branch.  In both cases `h` is simple
at every arm and

```text
f=A h^k,                   q_(-2)=B h^(2k+2),
q_1=M F^(2k-1),            Sigma|h.                   (39)
```

The possible half-integral homogeneous mode is excluded by the simple arm
of `h`.  Define

```text
C=2(k+1)B/(kA),            lambda=(k+2)C/(2kA).
```

The two internal coefficients are

```text
q_(-1)=C h^(k+2)F,
q_0=D h+lambda h^2F^2.                                  (40)
```

The scalar row factors as

```text
(F h'+2h F')
[ AMk(2k-1)h^(k-1)F^(2k-2)-D-2 lambda hF^2 ].         (41)
```

Its Euler factor has positive degree.

### 4.3 Double-upper placement

Finally,

```text
supp(Q)={-2,R-1,2R,3R+1}.                             (42)
```

The lower endpoint `W_(-R,-2)(f,q_0)=0` is incompatible with a simple
`q_0` when `R` is odd.  Thus `R=2k` and

```text
f=A h^k,                 q_0=B h,          Sigma|h.   (43)
```

The upper endpoint gives `q_3=N F^(6k+1)`.  Put

```text
u=6k+1,          E_0=ANu,
q_2=D F^(4k)+E_0 h^kF^(6k),
lambda=4ADk,     mu=3AE_0k,
q_1=M F^(2k-1)+lambda h^kF^(4k-1)
                  +mu h^(2k)F^(6k-1).                 (44)
```

Substitution gives the last scalar factorization:

```text
(F h'+2h F')
[ -B
  +AkM(2k-1)h^(k-1)F^(2k-2)
  +Ak lambda(4k-1)h^(2k-1)F^(4k-2)
  +Ak mu(6k-1)h^(3k-1)F^(6k-2) ].                    (45)
```

The nonconstant first factor closes the sixth ladder.

## 5. Exact consequence and failure boundary

Combining the six ladders with THM-3569 gives

```text
PROVED:   (#supp P,#supp Q)=(2,m) or (m,2) forces m>=5;
          if one output has two pieces, the total is at least seven;
          every six-piece Darboux pair, if one exists, is exactly 3 x 3.

PROVED:   THM-3579 separately closes the equal-step 3 x 3 subclass.

OPEN:     the general, necessarily unequal-step 3 x 3 cell;
          all 2 x 5 and larger cells;
          existence of any polynomial Darboux pair on Y_Sigma;
          JC(2).                                          (46)
```

In particular, this theorem does **not** prove a universal seven-piece
minimum: unequal-step `3+3=6` remains live.  It also does not extend
automatically to surfaces `c^n e=Sigma(b)` with `n>=3`; the negative-weight
membership and scalar-arm gate change with `n`.

The hostile rational pair

```text
{b,-1/c}=1                                             (47)
```

still explains the obstruction boundary: the one-piece solution exists only
after retaining a pole.  Degree-one `Sigma` remains the sharp polynomial
exception, while a repeated root makes the Poisson tensor degenerate and
prevents every polynomial Darboux pair for a different reason.

## 6. Exact companion

Reproduce with

```bash
python3 04-computation/jc2_danielewski_two_by_four_weight_nonentry_thm3583.py
python3 -O 04-computation/jc2_danielewski_two_by_four_weight_nonentry_thm3583.py
```

The companion enumerates the `L/B/U` support paths for `1<=R,T<=8`, checks
the `[1,1,2,2,2]` row multiplicity profile, verifies all displayed
first-order recurrences and scalar factorizations through `k,T=8`, exercises
the square-homogeneous boundary, the gcd-five rows `k=3,8`, the gcd-four
rows, and the balanced arm-order hostile through `R=64`.  These are exact
identity and hostile controls.  The universal theorem is the support,
valuation, UFD, and Euler-factor proof above, not extrapolation from the
finite ranges.

**QED.**
