# Independent audit of the leading-coefficient exactness bridge

**Status: INDEPENDENT ANALYTIC / SOURCE / NORMAL-AND-OPTIMIZED REPLAY
AUDIT PASS. Frozen producer hashes match. No mathematical correction
remains; root owns primary promotion.**
September 8, 2026. Referee: `three_ray_geometry`.

Producer: root. Accepted primary/source stem:
`planar_jc48_sep08_leading_exactness`. Only this audit sidecar is owned by
the referee; no producer, older theorem, source or output was edited.

## 1. The independently checked mathematical statement

Let `F in C(x)[t]` have t-degree n>=1 and nonzero leading coefficient
`A in C(x)`. A rational `G in C(x,t)` with `J_(x,t)(F,G)=1` forces

    dx/a to be exact in K=C(x)(a),   a^n=A.

This is a necessary condition for arbitrary lower coefficients. It gives a
nonzero-residue obstruction at any finite zero of A of order exactly n.
In the specified DG quartic square-prefix setting, the leading coefficient
is `A=N0(x)^2` and one may take `a=sqrt(N0)`. Thus a finite double root of
N0 excludes a rational mate, and therefore any global mate, regardless of
the lower coefficients or M-value there.

The word finite is essential: for example `N0=(x-r)^6`, the binary-octic
profile is six at r and two at infinity, but `dx/sqrt(N0)` is exact.
The proposed leading gate alone therefore does not close that orientation
of the profile. Any binary `(4,2,2)` profile has two distinct double points;
at least one lies in the finite x-chart, so this caveat does not leave an
unpaid `(4,2,2)` orientation. No choice of a different chart is assumed.

## 2. Field, derivation and formal Hensel construction

All calculations are over characteristic zero. Choosing a root a does not
require `Y^n-A` to be irreducible. The field K is the simple field generated
by that chosen root, and may have degree strictly smaller than n.
Because C contains all nth roots of unity, different choices of a generate
the same field and differ by a complex scalar root of unity. Exactness of
dx/a is independent of that choice.

The derivation `d/dx` extends uniquely through this separable algebraic
extension, with

    a' = A'/(n a^(n-1)).

Extend it coefficientwise to the Laurent-series field K((v)), holding v
fixed. Separately `d/dv` holds K fixed. These derivations commute. Every
coefficient operation is formal and finite at each power; no analytic
convergence or specialization away from poles in x is required.

Write `F=A t^n+sum_(i<n) f_i t^i`. The formal equation

    A U^n + sum_(i<n) f_i U^i v^(n-i) = 1

has the constant solution `U(0)=a^(-1)`. Its derivative in U at v=0 is
`n A(a^(-1))^(n-1)=n a`, a nonzero element, hence a unit of K. Formal
Hensel recursion uniquely determines `U in K[[v]]` with that constant
term. Set `T=U/v`. Then

    T=a^(-1)v^(-1)+O(1),       F(x,T)=v^(-n).

The embedding of rational functions is fully paid. Since T has v-order
minus one, distinct powers of T have distinct orders, and the highest
power of any nonzero polynomial in K[t] cannot cancel. Thus substitution
`t -> T` injects K(t) into K((v)), in particular the required C(x,t).
Every nonzero rational denominator remains nonzero and is invertible in
K((v)). Equivalently an algebraic T over K would contradict `F(T)=v^-n`.

## 3. Exact coefficient, sign and index

Let `Gtilde=G(x,T)`. Differentiating `F(x,T)=v^-n` with v fixed gives
`T_x=-F_x/F_t`. The sign of the Jacobian therefore yields

    Gtilde_x = G_x+G_t T_x
             = (G_x F_t-G_t F_x)/F_t
             = -1/F_t.

Differentiating instead with respect to v gives
`F_t T_v=-n v^(-n-1)`, hence

    Gtilde_x = (1/n)v^(n+1) T_v.

The derivative of the leading term of T is `-a^-1 v^-2`. Therefore the
coefficient of v^(n-1) on the right is precisely `-1/(n a)`.
The derivative of the v^0 term of T contributes nothing; all positive
powers contribute strictly later. If b is the v^(n-1) coefficient of
Gtilde, then b lies in K and

    b'=-1/(n a),       d(-n b)=dx/a.

This argument includes n=1. It imposes no bound on the pole order or
numerator/denominator degrees of G. A nonzero constant Jacobian lambda
is reduced to 1 by dividing G by lambda.

## 4. Residue consumer and independently reconstructed controls

At a finite point r where `A=(x-r)^n H(x)` and H(r)!=0, every local branch
of the normalized cyclic field has

    a=(x-r)h(x),       h(r)^n=H(r),       h(r)!=0.

The nth roots of the unit H are formal holomorphic units in
`C[[x-r]]`; there is no local ramification. Hence

    dx/a = (1/h(r)) dx/(x-r) + regular terms.

Its residue is nonzero. The differential of a rational function on a
normal curve has zero residue at every point, as is immediate by
termwise differentiation of its Laurent series. This contradicts the
necessary exactness. Global reducibility of `Y^n-A` or a rational K does
not remove the local pole: it merely selects fewer of the branch choices.

For n=4 and `A=N0^2`, use `a=sqrt(N0)`. If
`N0=(x-r)^2 H(x)` with H(r)!=0, the residue is one of the nonzero values
`+/-1/sqrt(H(r))`, with the corresponding single choice if the field is
already rational. No assertion about M, lower F coefficients, generic
fibre smoothness or compact-fibre pole cancellation is needed.

The referee independently checked the following model controls, all now
included in the accepted producer:

1. For `N0=x^5(x-1)^3`, set `z=sqrt((x-1)/x)`. Then
   `x=1/(1-z^2)` and `a=z^3/(1-z^2)^4`, so
   `dx/a=2(z^-2-2+z^2)dz=d(-2/z-4z+(2/3)z^3)`. This leading gate passes.
2. There is an actual rational pure-monomial mate for that example:
   `A=x^10(x-1)^6`,
   `B=(8x^2-4x-1)/(6x^9(x-1)^5)`,
   `F=A t^4`, `G=B t^-3`. A separate SymPy expansion gave
   `-3A'B-4AB'=1`. This is not asserted to be a global DG section.
3. Passing the leading gate is not sufficient for general lower terms:
   `F=x^2+t^2` has A=1. On its generic fibre F=c, put `h=x+i t`; then
   `-dx/F_t=-i dh/(2h)`, with a nonzero residue, so no rational G can
   have J(F,G)=1. This retains the necessary-only direction.

A supplementary formal observation was sent to the producer and is now
proved explicitly in its Section 5: for a pure
monomial `F=A t^n`, exactness is sufficient for a rational mate. If
`H'=1/a`, average H over the cyclic Galois group to transform under its
character as a^-1. Then `B=-H/(n a^(n-1))` lies in C(x), and
`G=B t^(1-n)` has Jacobian one. This observation is a model theorem only;
it does not infer sufficiency for arbitrary lower coefficients. The final text includes the correct character-weighted average and the
explicit descended Jacobian identity. Both directions are accepted.

## 5. Full primary read: additional residue, descent and scope checks

The final primary was read completely, after the analytic bridge had been
checked independently. Its added two-finite-pole statement is valid for
arbitrary positive integers alpha,beta and distinct finite p,q. The
coefficient of `(x-p)^(alpha-1)` in `(x-q)^(-beta)` is

    (-1)^(alpha-1) binom(alpha+beta-2,alpha-1)
       * (p-q)^(-alpha-beta+1),

which is nonzero. This independently derives the producer's residue
formula and recovers its two-distinct-finite `(4,4)` octic case. A nonzero
constant nth-root factor scales the residue and does not cancel it.
Neither a source-coordinate translation of infinity nor a statement about
all placements of a two-point profile is implicit in this calculation.

For pure-monomial sufficiency the radical extension is indeed normal:
all roots of `Z^n-A` are already in C(x)(a). It is separable and Galois,
even when its degree is a proper divisor of n. The character on a is
faithful on this group and has nth-power one. If `H'=1/a`, then

    Hbar = (1/|Gal|) sum_gamma chi(gamma) gamma(H)

still has derivative 1/a, and a change of index gives
`gamma(Hbar)=chi(gamma)^(-1)Hbar`. Thus
`B=-Hbar/(n a^(n-1))` is invariant and belongs to C(x). Direct
calculation gives `-(n-1)A'B-nAB'=1`. This is a sufficient condition for
`F=A t^n` and not for a general lower-coefficient deformation.

The genuine lower-row control `F=x^3t^2+x`, `G=xt/F` has Jacobian one.
Its inverse branch is
`T=sqrt(1-xv^2)/(x sqrt(x) v)`, so
`Gtilde=sqrt(1-xv^2) v/sqrt(x)` and
`[v]Gtilde=1/sqrt(x)`. Its derivative is
`-1/(2x sqrt(x))`, exactly the coefficient demanded by the theorem.
The chosen square-root branch is the formal one with constant term one.
There is no analytic branch continuation assumption.

The cubic-field example `F=x^2t^3`, `G=-1/(xt^2)` also has Jacobian
one. With x=q^3, the leading root is a=q^2, and `dx/a=d(3q)`.
This is an explicit hostile to demanding the primitive in C(x) instead
of in the correctly retained radical field. It also checks the character
choice in the descent: H=3q transforms as a^-1, and the descended B is
exactly -1/x.

The conic necessary-only hostile is rigorous on the generic fibre over
C(c). Restriction of a rational mate is defined in that generic function
field; an isolated special-fibre denominator problem cannot evade the
nonzero residue in `-i dh/(2h)`. The rational-field pure-monomial passing
control and the conic failing control therefore retain different, clearly
stated scopes.

### Antecedent recovery

I read the actual file
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md).
Its current theorem is stronger than the historical square/parity name,
but its top-coefficient and parity descent are accurately described as
antecedents in the producer. They are in the polynomial quadratic-fibre
setting and are not treated as this arbitrary-rational-mate theorem.

After checking the core-paper route, I independently opened the official
MPIM PDF of Makar-Limanov and Trakhtenberg,
[*Properties of a Jacobian mate*, 2024-33](https://archive.mpim-bonn.mpg.de/5148/1/mpim-preprint_2024-33.pdf).
Its printed pp.4–5 contain the radical-expansion lemma and successive
subtraction of fractional powers from a transformed mate. The producer
credits these as related Newton/edge mechanisms; it does not transport
the present theorem from that paper or assert external priority. This
limited provenance statement is accepted. No downloaded literature file
was placed in the repository.

## 6. Source universe and exact independent replay

I read the complete 5,993-byte source. Its `need` function raises on a
false value and its algebraic `zero` check uses exact rational cancellation.
No Python assertion disappears in optimized mode. The source neither
changes producer files nor relies on a numerical tolerance or hidden
sample filter.

The declared 241 gates decompose as follows:

| Complete control family | Gates |
|---|---:|
| n=1,...,6, e in {0,2,3,4,5,6}: rational monomial mate, primitive, coefficient |108|
| Nontrivial cubic radical field |4|
| Genuine nonconstant lower-row rational mate and formal coefficient |5|
| n=1,...,6 explicit inverses with nonconstant lower terms |24|
| n=1,...,8 multiplicity-n logarithmic root controls |16|
| All ordered finite pole orders alpha,beta=1,...,6 |72|
| All four finite/infinity placement types of `(4,2,2)` |4|
| Passing odd `(5,3)` octic and its actual rational model mate |3|
| Leading-pass / generic-fibre-fail conic hostile |4|
| The infinity-double passing primitive |1|

The omitted e=1 in the positive monomial bank is explicitly the logarithmic
case tested separately. No vanishing control was removed after observing
its result. The source's `cancel` before coefficient extraction pays the
reported symbolic-expression false negative: the final check is an exact
identity of Laurent coefficients, not a numerical or branch-dependent
simplification.

I independently ran both commands from the actual worktree:

```bash
python3 04-computation/planar_jc48_sep08_leading_exactness.py
python3 -O 04-computation/planar_jc48_sep08_leading_exactness.py
```

Both completed with all241 gates. Their outputs are byte-for-byte equal
to each other and the saved frozen output: 459 bytes. The source and
output pins were re-read after both runs; they are unchanged.

| Frozen artifact | Bytes | SHA256 |
|---|---:|---|
| Source |5993|`8613ca7d44fc85c9619b02b8a67401c66aed745c5e7977a949bc51b7cb332446`|
| Output |459|`2dfa3a229dc8436be302c06173181aa89f8549f9b4e547fcd80ddd431d35f418`|
| Primary before root's status-only promotion |10585|`3e485175123fc0a5d850460f7ad025d797e50eac515921b8db84c3224be880e6`|

The gate-label SHA256 is
`151ca181cf03cd30da473d1d548e9e5778e897c79795cea8576a30d00aefcf9c`.
This is a reproducibility pin, not a replacement for the checked values
or the unbounded analytic proof.

## 7. Final finding

Full analytic proof, field types, actual source-coordinate consumer,
finite/infinity scope, pure-monomial equivalence, necessary-only hostile,
complete source universe and both independent frozen replays **PASS**.
No remaining mathematical or source correction is requested. This
excludes the declared finite-root classes, including the entire binary
`(4,2,2)` placement class, and does not close the full quartic problem or
JC(2). Root owns promotion, maintained routing and Git integration.
