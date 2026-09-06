# Two first channels force a nonzero second trinomial return, with both carries

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED by root.**
The root independently checked all binomial ratios, the determinant direction,
complete carry-index bounds, and the conclusion for complex torus
coefficients on 2026-09-06. No theorem ID is reserved here. The all-exponent proof below
closes the exactly-two-first-channel stratum of the corrected trinomial
problem. Three or more first channels remain outside its scope.

Canonical namespace: [THM-4432, trinomial two-channel two-rung noncancellation
with carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md).
The linked canon file controls namespace promotion; this note supplies the
audited proof and evidence for the shared concurrent result.

## Inheritance and the research move

The closest proved mechanism is [THM-2639,
free equal-mass two-rung certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md).
The corrected near miss is [MISTAKE-544](../../01-canon/MISTAKES.md), repaired
by [the September 5 trinomial classification](synthesis_20260905_moments_trinomial.md).
The canonical hostile is its `(-13,1,8)`: two channels at the first mass seven
gain **two** extra channels at mass fourteen. The least-used relevant sidecar
is the pair of exact residue carries. They must be retained, but at a tuned
first zero their contributions have the favorable sign.

The active concept board was: minimal support channels; numerical-semigroup
residues; multinomial coefficient ideals; first-moment cancellation phase;
sharp Laurent width; horizontal Wick transport. The proof turns the carry
obstruction into a sign certificate. This implements the protocol cards
“expose the obstruction first” and “find the hidden second coordinate.”

The paper starting the session, Heule–Scheucher,
[*Happy Ending: An Empty Hexagon in Every Set of 30 Points*](https://arxiv.org/pdf/2403.00737),
motivates seeking a small certificate while retaining its validity sidecar.
There is no asserted map from point orientations to Laurent moments. The
actual mathematical map here is the channel-to-one-variable map specified
below, independently of that methodological inspiration. No claim of
literature novelty is made.

The prior Sep05 evidence was initially invisible in this sparse checkout;
it was then read directly by `git show HEAD:<path>`. Its Section 6 leaves
general two-rung coprimality open and gives finite evidence for 2,275 collided
supports. Its Section 3 identifies, but does not dispose of, the two carry
terms. This result adds the all-exponent proof when the first row has two
channels; it does not reclassify the already proved support semigroup.

## 1. Statement and exact first-return formula

Let

```text
f(u)=c0*u^(-a)+c1*u^b+c2*u^c,
a>=1, 1<=b<c, gcd(a,b,c)=1, c0*c1*c2!=0.
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g.
```

Thus `1<=A<B`, `gcd(A,B)=gcd(a,g)=1`. Let a **channel** mean a
nonnegative integer multiplicity vector, before its multinomial weight.
Assume the first nonempty balanced channel row has exactly two channels.
The corrected classification gives its mass as `g`. Equivalently, the
canonical representation has the form

```text
a=AB+As+Bt, 0<=s<B, 0<=t<A.
n=g-B-s-t, h=B-A,
v=(n,B+s,t), w=(n+h,s,t+A), d=w-v=(h,-B,A).
```

Here `n>0`, because `a*n=b*(B+s)+c*t>0`. Put

```text
alpha=g! / (n! (B+s)! t!),
beta =g! / ((n+h)! s! (t+A)!),
z=c0^h*c1^(-B)*c2^A.
```

**Theorem.** The first nonzero constant-term moment is exactly

```text
m_*(f) = g    if z != -alpha/beta,
m_*(f) = 2g   if z == -alpha/beta.                 (1)
```

Consequently `m_*(f)<=2g<=a+c`. The result is uniform over all nonzero
complex coefficients and all exponents in this stratum, including nonfree
balanced semigroups. It uses no real-coefficient, genericity, or bounded
relation-height assumption. Reflection handles two negative charges;
common charge dilation preserves all moments and can be divided out first.

The singleton stratum already has noncancellation at its first support
return. Hence the sharp width bound now follows elementarily for all
trinomial supports with at most two first channels. The residual of this
approach begins at three first channels.

## 2. An elementary short-residue binomial inequality

For integers `H>=1,r>=0`, define

```text
R_H(r)=binom(2r+H,r)/binom(2r,r)
      = product_(i=1)^H (2r+i)/(r+i),
S_H(r)=binom(2r+2H,r+H)/binom(2r+H,r)
      = product_(i=1)^H (2r+H+i)/(r+i).
```

Every factor gives

```text
R_H(r)<2^H,                  S_H(r)>=2^H.          (2)
```

There is a crucial improvement when `r<H`:

```text
R_H(r)<=2^(H-1).                                  (3)
```

To prove it, `R_H(r)/2^H` strictly decreases with `H`, since its successive
ratio is `(2r+H+1)/(2(r+H+1))<1`. It therefore suffices to take `H=r+1`.
Set

```text
Q_r=R_(r+1)(r)/2^r
   =(3r+1)! r! / (2^r (2r+1)! (2r)!).
```

Then `Q_0=Q_1=1` and

```text
Q_(r+1)/Q_r=3(3r+2)(3r+4)/(8(2r+3)(2r+1))<=1,
```

because denominator minus numerator is `5r(r+2)>=0`. This proves (3),
with equality only `(r,H)=(0,1)` or `(1,2)`. The strict first inequality
in (2) holds even at `r=0`.

## 3. The three internal coefficients have a strictly negative determinant

Write

```text
C=(2g)! / product_i (2v_i)!,
D=(2g)! / product_i (v_i+w_i)!,
E=(2g)! / product_i (2w_i)!,
Delta=C*beta^2-D*alpha*beta+E*alpha^2.
```

Set

```text
Nv=product_i binom(2v_i,v_i),
Nw=product_i binom(2w_i,w_i),
M =product_i binom(v_i+w_i,v_i), K=binom(2g,g).
```

Factorial cancellation gives

```text
Delta/(alpha^2*beta^2)=K*(1/Nv-1/M+1/Nw).          (4)
```

The channel coordinates turn the needed inequality into two half budgets:

```text
M/Nv=R_h(n) * S_B(s)^(-1) * R_A(t)
    <2^h * 2^(-B) * 2^(A-1)=1/2,

M/Nw=S_h(n)^(-1) * R_B(s) * S_A(t)^(-1)
    <=2^(-h) * 2^(B-1) * 2^(-A)=1/2.              (5)
```

Here `h=B-A`, and the two applications of (3) use exactly `t<A` and
`s<B`. Thus `M/Nv+M/Nw<1`, and (4) proves

```text
Delta<0.                                         (6)
```

This extends the free-ray determinant from THM-2639. Freeness is replaced
by both endpoint residue bounds, which are available because the complete
first row contains exactly two channels.

## 4. Both extra second-row channels reinforce the sign

For a coefficient vector `c=(c0,c1,c2)`, use `c^r=product_i c_i^(r_i)`.
The first moment is

```text
CT(f^g)=c^v*(alpha+beta*z).                        (7)
```

Every balanced mass is a multiple of `g`, since the balance equation is
`g(A*y+B*z_count)=a*m` and `gcd(a,g)=1`. The complete second channel row is

```text
2v+j*d, -epsilon_t<=j<=2+epsilon_s,
epsilon_t=floor(2t/A), epsilon_s=floor(2s/B),
epsilon_t,epsilon_s in {0,1}.                     (8)
```

To check completeness directly, its positive counts are
`2B+2s-Bj` and `2t+Aj`; their nonnegativity gives precisely the bounds
in (8). The negative count is then positive by the original charge
equation, so it introduces no extra restriction.

Let `K_j=(2g)!/product_i(2v_i+j*d_i)!`, which is strictly positive for
every index in (8). Then

```text
CT(f^(2g))/c^(2v)=sum_(j=-epsilon_t)^(2+epsilon_s) K_j*z^j.  (9)
```

At the tuned value `z=-alpha/beta<0`, the internal terms `j=0,1,2` sum
to `Delta/beta^2<0`. The only possible extra indices are `j=-1` and
`j=3`; each contributes a strictly negative real number. Thus the
**complete** right side of (9) is negative. Since `c^(2v)!=0`, the actual
moment is nonzero. Together with (7) and mass divisibility, this proves (1).

This sign argument does not discard the carries or replace the second row
by a square of the first. The offending channels from MISTAKE-544 are the
two extra odd powers that make the certificate stronger.

## 5. Sharpness, an explicit width-three family, and the method boundary

The delayed value `2g` is attained for every support in the theorem's
stratum by choosing nonzero coefficients with `z=-alpha/beta`.
The span inequality is equality exactly when `B=2`, which forces `A=1`.
In the two-channel stratum this leaves `a=2` or `a=3`.

For a concrete family beyond the width-two theorem, take any integer
`g>=4` with `gcd(g,3)=1` and

```text
f_g(u)=u^(-3)+eta*u^(g-3)+u^(2g-3),
eta^2=-6/(g-2).
```

Here

```text
CT(f_g^g)=binom(g,3)*eta^3+g*(g-1)*eta=0,

CT(f_g^(2g))=
  2g(g-1)(2g-1)(23g^2-47g+20)/(15(g-2)^2)>0.       (10)
```

Equation (10) follows by substituting `eta^2` into the four channels
`(2g-6+j,6-2j,j)`, `j=0,1,2,3`. Its quadratic factor is positive for
`g>=4` (it is positive at four and has positive derivative there onward).
Thus `m_*(f_g)=2g=a+c`. For `g>=7` this is genuinely width three;
the smaller allowed `g=4,5` have width one and two after reflection.
The extra `j=3` term is present, so the free-semigroup theorem alone does
not supply this family.

The exactly-two-channel hypothesis cannot be dropped from the determinant
argument. For the primitive support `(-4,3,10)`, the first mass is seven
with channels `(3,4,0),(4,2,1),(5,0,2)`. Keeping only the first two gives
`Delta=0`. For `(-4,5,14)`, the first mass is nine with channels
`(5,4,0),(6,2,1),(7,0,2)`; the first two give `Delta=126309456>0`.
The first failed implication is applying both short-residue bounds to an
adjacent pair that does not span the complete first slice. These examples
do not refute general two-rung coprimality; they refute this unqualified
pair-compression argument. The missing coordinate is the remaining first
channels, whose polynomial zero is no longer a single negative ratio.

## 6. Connection contract and limitations

**Exact source and target.** A complete two-channel first return and its
second return map to `alpha+beta*z` and the Laurent polynomial in (9),
where `z=c^d`. Dividing by nonzero coefficient monomials preserves their
zero sets. The quotient forgets the overall moment magnitude and phase;
`c^v`, `c^(2v)` restore these. Multinomial weights and both carries are
necessary sidecars. The decisive hostile is `(-13,1,8)`, whose second
row has five channels and whose normalized tuned second moment equals
`-1227968/125`.

**Wick transport.** For integer `H>=c`, set

```text
P(Z,W)=c0*Z^H*W^(H+a)+c1*Z^H*W^(H-b)+c2*Z^H*W^(H-c).
```

With `L(Z^i W^j)=delta_(i,j)*i!`, one has
`L(P^m)=(Hm)!*CT(f^m)`. Hence the same exact first-detection formula
holds on this horizontal Wick face. The common factorial is the sidecar
preserving cancellation. This does not control arbitrary radial data,
general Gaussian detection order, higher-dimensional moment conjectures,
LRC(14), or the sharp Laurent bound for wider general supports.

## 7. Reproduction and exact scope

```text
python3 -B 04-computation/trinomial_two_channel_empty_core_returns_sep06.py
python3 -B -O 04-computation/trinomial_two_channel_empty_core_returns_sep06.py
```

The source is
[the exact verifier](../../04-computation/trinomial_two_channel_empty_core_returns_sep06.py)
and its saved output is
[the matching transcript](trinomial_two_channel_empty_core_returns_sep06.out).
All checks use explicit failures that remain active under `-O`.
Normal and optimized output bytes agree. SHA-256 of raw LF bytes:

```text
source cbc09597b4e6b8ec4cddccbd46e47f20d4c36dd48cba4c818829d2a01be9ae53
output 90a973c99fa3eb8b83d5a424245fed19a7b6a1ca25167785ad4b36c408cfe3c9
```

The support census is every primitive `(-a,b,c)` with `1<=a,b<=40`,
`b<c<=60`, filtered to exactly two channels at mass `g`. It has 559
supports: carry profiles `(0,0),(0,1),(1,0),(1,1)` occur 258, 136, 112,
53 times. Original charge-equation enumeration checks the complete second
row, determinant sign, both half budgets, and the actual nonzero tuned
moment in exact rational arithmetic. The separate binomial lemma is checked
at all `1<=H<=80,0<=r<H` (3,240 cases).

An independent implementation repeatedly multiplies rational Laurent
polynomials on six controls, including all carry profiles and both named
correction witnesses. It checks **all** moments through `2g`, proving
the returned first-nonzero index is `g` at positive coefficients and `2g`
at the tuned rational coefficient point. It also verifies the exact
normalized moment against (9). The width-three expression is checked at
all 65 allowed `g` between four and one hundred. The analytic proof,
rather than any of these finite universes, supplies the all-exponent claim.

## Concurrent integration

After the independent proof and root audit, incoming commit `4d1ad2a390`
reserved [THM-4432, trinomial two-channel two-rung noncancellation with
carries](../../01-canon/theorems/THM-4432-trinomial-two-channel-two-rung-noncancellation-with-carries.md).
That reservation explicitly credits concurrence with this empty-core
overnight session. It arrived as an **UNPROVED EMPTY STUB**, whose reserved
state was not used as a proof dependency; the linked canon owns subsequent
promotion status.
This is one shared mathematical result with concurrent ownership, not
separate novelty for each artifact. The present note remains an explicit
audited proof with independently replayed evidence.
