# Two first-return channels: all doubling carries strengthen a positive certificate

**Status: PROVED / INDEPENDENTLY AUDITED.** Root and observer independently
audited the complete theorem, factorial-ratio lemma, affine-fibre bookkeeping,
polynomial certificate, exact detection, and nonvacuity. The finite controls
below corroborate the proof; they do not replace it. Canonical promotion and
shared navigation are root-owned.

## 1. Exact statement and inheritance

Let

```text
f(u)=alpha*u^(-a)+beta*u^b+gamma*u^c,
a>=1, 1<=b<c, gcd(a,b,c)=1, alpha*beta*gamma!=0.
```

Suppose the **complete first nonempty support-return fibre has exactly two
channels**. Then, with `g=gcd(a+b,a+c)`,

```text
CT(f^g) and CT(f^(2g)) cannot both vanish.
```

There is no bound on either endpoint, the return mass, or the coefficient
values. No semigroup freeness and no absence of doubling carries is assumed.
The first nonzero moment is exactly `g` or `2g`, according as the first
moment does not or does vanish. Both outcomes occur on the coefficient
torus. Consequently this whole class obeys the sharp width ceiling
`m_*<=a+c`, since the inherited classification gives `2g<=a+c`.
Reflection and common charge scaling preserve the zero/nonzero moments;
the primitive statement is the precise normalization used below.

The closest proved mechanism is
[THM-2639, free equal-mass two-rung persistent-collision certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md).
Its all-level free-semigroup hypothesis implies a three-channel second row,
but does not cover the present nonfree class. The exact canonical residues
and the two floor carries come from
[the prior trinomial classification, Sections 2--3](synthesis_20260905_moments_trinomial.md).
The canonical hostile is `(-13,1,8)`: its first row has two channels and
its second has five, so silently retaining only sums of first-row channels
is invalid. The least-used sidecar is the **parity of each extra affine
index** after tuning the first scalar sum.

The live board is: first-row multiplicity; primitive affine step; canonical
residues; factorial-ratio comparison; torus monomial quotient; exact
polynomial remainder. The source-to-target map takes the complete affine
fibre to its integer index and factors out a nonzero coefficient monomial.
It preserves every multinomial weight and all cancellation information.
Dropping the extra indices destroys completeness; their signs, not their
absence, supply the successful operation here.

### Concurrent-work and duplication audit

Incoming commit `e10dff3181` proves the entire trinomial smaller-endpoint
strip through eight by thirty symbolic resultants and five finite
opposite-endpoint cases; see
[the complete width-eight report](overnight_20260906_moments_width8.md).
That theorem includes higher first-row multiplicities, so this result does
not subsume it. Conversely, the theorem here has no endpoint bound and is
not supplied by its finite family bank. Their common two-channel cases are
scope intersections, not separate new discoveries.

The incoming report already proves the sharp width-three family below and
contains the two-carry hostile. Those are explicitly inherited controls.
The independently active empty-core session's board also names the same
nonfree two-first-channel target; development is therefore treated as
concurrent, with no separate priority credit presumed. Targeted searches
of theorem/result statements, exact control supports, no/later-carry and
two-rung phrases found THM-2639 and the two incoming notes as the closest
mechanisms. No external-priority claim is made.

## 2. The full first and second affine fibres

Put `A=(a+b)/g`, `B=(a+c)/g`. The inherited classification gives

```text
1<=A<B, gcd(A,B)=1, gcd(a,g)=1,
a=A*(B+r)+B*z,        0<=r<B, 0<=z<A,
x=g-B-r-z>0.
```

Here the coefficient of `B` in the canonical middle count is exactly one
because there are exactly two first-row channels. They are

```text
v=(x,B+r,z),          w=(x+B-A,r,z+A).
```

Their difference `d=(B-A,-B,A)` is the **primitive** integer affine-fibre
step: `gcd(A,B)=1`. Thus no hidden intermediate channel is lost. In fact
the original balance equations at mass `kg` are exactly

```text
A*y+B*z'=k*a,          x'=kg-y-z'.
```

Every nonnegative solution in `y,z'` has `x'>0`, because
`a*x'=b*y+c*z'` and the mass is positive. Therefore the second row is
precisely

```text
2v+j*d,       -epsilon_z<=j<=2+epsilon_y,
epsilon_y=floor(2r/B), epsilon_z=floor(2z/A),
epsilon_y,epsilon_z in {0,1}.
```

The lower and upper bounds come from the third and second coordinates.
The positivity observation ensures that the first coordinate removes no
endpoint in this actual trinomial setting. The three generated channels
have indices `0,1,2`; the only extra indices are `-1` and `3`.

## 3. Elementary factorial ratios, with all strict boundaries

For integers `t>=0,L>=1`, define

```text
R_+(t,L)=binom(2t+L,t)/binom(2t,t),
R_-(t,L)=binom(2t+L,t)/binom(2t+2L,t+L).
```

Then

```text
R_+(t,L)<2^L,
R_-(t,L)<=2^(-L),
L>t  =>  R_+(t,L)<=2^(L-1).                         (R)
```

The first inequality follows from the product
`prod_(i=1)^L (2t+i)/(t+i)`, whose factors are strictly less than two.
For the second, use

```text
R_-(t,L)=prod_(i=1)^L (t+i)/(2t+L+i).
```

Every factor is at most one half. Equality holds precisely when `L=1`.
For the last inequality first set `L=t+1` and rewrite

```text
R_+(t,t+1)=prod_(i=1)^t (2t+1+i)/(t+i) <=2^t.
```

Each factor is at most two. Increasing L strictly decreases
`R_+(t,L)/2^L`, since its consecutive ratio is
`(2t+L+1)/(2(t+L+1))<1`. This proves (R). Equality in its last bound
occurs only at `(t,L)=(0,1)` and `(1,2)`. The strict hypothesis `L>t`
cannot be dropped: `R_+(1,1)=3/2>1`.

Now allow **any integers** `0<A<B`, `x,r,z>=0`, `r<B`, `z<A`, and form
the two vectors v,w above, of common mass `g=x+B+r+z`. This purely
factorial lemma does not need coprimality, positive charges or x>0. Define

```text
N_v=prod_i binom(2v_i,v_i),
N_w=prod_i binom(2w_i,w_i),
M=prod_i binom(v_i+w_i,v_i),       delta=B-A>0.
```

The exact coordinate decompositions are

```text
M/N_v = R_+(x,delta) R_-(r,B) R_+(z,A) <1/2,
M/N_w = R_-(x,delta) R_+(r,B) R_-(z,A) <=1/2.        (H)
```

Indeed the powers of two sum to `delta-B+(A-1)=-1` and
`-delta+(B-1)-A=-1`. The first inequality is strict because delta>0.
The second can be equality: exactly `A=1,B=2,r=1,z=0`, for any x>=0.
Thus its weak inequality is intentional; the sum in (H) is still strictly
less than one.

## 4. The generated-three-channel determinant is strictly negative

Let `mult(u)=|u|!/prod_i u_i!`, and write

```text
p=mult(v), q=mult(w),
C=mult(2v), D=mult(v+w), E=mult(2w),
Delta=C*q^2-D*p*q+E*p^2.
```

All five weights are positive integers. Cancelling factorials gives

```text
Delta/(p^2*q^2)
 =binom(2g,g)*(1/N_v-1/M+1/N_w) <0.                 (D)
```

The strict sign follows from `(M/N_v)+(M/N_w)<1` in (H).
This is the extension of THM-2639's coefficient comparison that does not
require the first vectors to be extreme rays of the whole semigroup.

## 5. Both carries strengthen the full polynomial certificate

For the coefficient monomials `X=alpha^v1 beta^v2 gamma^v3` and
`Y=alpha^w1 beta^w2 gamma^w3`, put `T_1=CT(f^g)`, `T_2=CT(f^(2g))`.
Let `L_-` and `L_+` be the multinomial weights at indices -1 and 3,
respectively, or zero when the corresponding carry is absent. Exactness
of the affine fibre gives

```text
T_1=pX+qY,
T_2=L_-*X^3/Y+C*X^2+D*X*Y+E*Y^2+L_+*Y^3/X.       (F)
```

These expressions are identities in the coefficient torus; all displayed
terms correspond to actual nonnegative exponent vectors when present.
At `T_1=0`, the normalized ratio is `t=Y/X=-p/q<0`. Therefore

```text
T_2/X^2 = Delta/q^2 + L_-/t + L_+*t^3 <0.          (S)
```

The inequality is for the normalized real scalar, not for the generally
complex moment T_2 itself. Each extra carry contributes negatively; no
channel has been discarded. This proves noncancellation immediately.

Here is an explicit identity without Laurent denominators. Define the
positive integer

```text
K=-p*q*Delta+L_-*q^4+L_+*p^4 >0
```

and the homogeneous cubic

```text
Q_3(X,Y)=q^3*L_+*Y^3
 +(q^3*E-p*q^2*L_+)*X*Y^2
 +(q^3*D-p*q^2*E+p^2*q*L_+)*X^2*Y
 +(q^3*C-p*q^2*D+p^2*q*E-p^3*L_+)*X^3.
```

Direct polynomial division of `X*Y*T_2` by `pX+qY` proves

```text
q^4*X*Y*T_2-Q_3(X,Y)*T_1 = K*X^4.                 (B)
```

The right side is nonzero on the coefficient torus. This is a literal
polynomial certificate after substituting the coefficient monomials, not
just a claim that an uncomputed resultant is nonzero.

## 6. Exact detection, unbounded nonfree family, and overlap controls

All support-return masses are divisible by g, and the first moment is the
binomial `pX+qY`. Its monomial quotient `Y/X` is nonconstant and takes every
value in `C*`, so it can be tuned to `-p/q`. Equation (B) then gives first
detection exactly at 2g. Otherwise detection is at g. This proves both
nonvacuity and the exact worst-case detection time for each support.

The equality `2g=a+c` is equivalent to B=2, forcing A=1. Since the first
row has two channels, this leaves a=2 or a=3. The a=2 class lies inside
THM-4417. The a=3 class is precisely the incoming report's already proved
sharp family `(-3,g-3,2g-3)`, `g>3`, `gcd(g,3)=1`. Its control
`(-3,1,5)` is retained; no new equality-family discovery is claimed here.

For an unbounded-endpoint genuinely nonfree family, take n>=1 and set

```text
A=n, B=n+1, r=1, z=0, g=(n+1)^2,
a=n(n+2), b=n(n^2+n-1), c=n(n+1)^2+1.
```

These supports are primitive because `a=g-1`, so `gcd(a,g)=1`; b>0 and
c-b=g. Their first row has two channels, while their first extra channel
appears at rung `k=n+1` through `floor(k/(n+1))`. The semigroup is never
generated by the first slice. For n>=3 its smaller endpoint is
`a=n(n+2)>=15` and grows without bound, outside the incoming strip through
eight. The first example there is `(-15,33,49)`, with first mass 16 and
worst-case detection 32. This family has no doubling carry for n>=2,
but the theorem also covers every actual one/two-carry case.

Readable controls, with signed support `(-a,b,c)`, are

| Support | First mass | Doubling carries | Role |
|---|---:|---|---|
| (-3,1,9) | 4 | none | inherited free equal-mass control |
| (-4,1,11) | 5 | none | nonfree; first extra channel at rung three |
| (-3,1,5) | 4 | upper only | inherited width-three equality control |
| (-13,1,8) | 7 | both | canonical five-channel hostile |

For the last row `(p,q)=(42,210)`, `Delta=-423783360`, and the full tuned
normalized second moment is `-1227968/125`, not the value obtained by
dropping the two extra channels. The positive integer in (B) is
`3821063113728`.

## 7. Failure boundary and the next open lane

The weakest theorem hypothesis is the complete first fibre having exactly
two channels, after the stated three-charge normalization. No condition is
placed on later-generation completeness. The weaker factorial lemma has
the explicit integer/residue hypotheses in Section 3 and no charge content
assumption. There is no Delta=0 case within those hypotheses.

For three or more first-row channels, one can no longer replace T_1 by a
binomial or force `Y/X=-p/q`. For example `(-4,1,6)` has first mass five
and three channels with weights `5,30,10`. Keeping two of those channels
would change the first equation itself. This is the first failed implication
of extending the present certificate to arbitrary multiplicity; it is not
a counterexample to two-rung coprimality. General trinomial two-rung
coprimality, the arbitrary-support Laurent width bound at smaller endpoint
at least three, and uniform Gaussian detection remain open.

The initial no-doubling-carry idea survived but was unnecessarily narrow:
the stronger ratio bound `L>t`, rather than `L>2t`, gives (D) for every
canonical residue. The cheap two-carry hostile then reveals that both
omitted terms have the favorable odd index. This is the mechanism, not an
extrapolation of a positive finite coefficient census.

## 8. Exact reproduction and independent paths

Run

```bash
python3 -B 04-computation/nc2_two_rung_overnight_hexagon_sep05.py
python3 -B -O 04-computation/nc2_two_rung_overnight_hexagon_sep05.py
```

The declared parameter universe is every `1<=A<B<=9`, `gcd(A,B)=1`,
`0<=r<B`, `0<=z<A`, `1<=x<=10`, filtered only by b>0 and `gcd(a,g)=1`.
It has 3,665 rows, of which 3,530 are nonfree. The four doubling-carry counts
are 1,187 / 999 / 794 / 685 for `(epsilon_y,epsilon_z)=00/10/01/11`.

For every row, one path scans the original negative multiplicity and solves
the original charge equation to reconstruct both complete weighted fibres.
A second path constructs the primitive affine line and both carry endpoints.
They agree before any polynomial comparison. Exact rational evaluation at
the first-row zero is then compared with the integer remainder, and all five
coefficients of (B) are checked independently. All 76,923 explicit gates
pass, and normal, optimized and stored output are byte-identical.
Scalar-ratio controls include
2,501 `(t,L)` pairs and the strict-residue failure. Wide controls cover the
nonfree families, delayed carries, and the inherited sharp family.

The all-height proof is elementary and independent of the finite replay.
No numerical coefficient roots, literature-priority, or Lean result is
claimed. The full proof remains a candidate until independent audit.

```text
script af6cb010e4a2f72b731a70a5ef0b31cc5472f7e8347f1a018bcc12bcafbd2ab6
output 2e7b9502365bb22dbefcdc7dce2ee09eb7bb74cf7db28ea7296e443e18b996c5
```
