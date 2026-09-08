# Independent audit of the all-m polynomial-mate weight gate

**Status: INDEPENDENT ANALYTIC + SOURCE + EXACT-REPLAY AUDIT PASS.**
No correction is required. This audit accepts the complete statement in
[the frozen primary](planar_jc48_sep08_polynomial_weight_gate.md), uniformly
for every integer `m>=1`, every characteristic-zero coefficient field,
and every polynomial mate degree. Its conclusion is necessary, not
sufficient, and its polynomial-source hypothesis is essential.

I read the full primary and source, reconstructed all leading-row cases,
and replayed the source normally and with `-O`. Both outputs agree with
all 553 bytes of the frozen 1,947-gate output. Three additional direct
original-source identities independently recovered the complete `(6,2)`
constant row, the tuned first row, and its actual derivative on `u=0`.
No producer functions were imported and no producer artifact was edited.

## 1. The complete source image is retained

The theorem's source and target rings are typed correctly:

    w=u^m t,
    K[u,t] -> K[u,u^-1,w],
    u^i t^k -> u^(i-mk)w^k.

The inverse exponent for a monomial `u^L w^k` is `L+mk`. Thus its
membership in the original polynomial ring is equivalent to
`k>=max(0,ceil(-L/m))`. Different pairs `(i,k)` give different pairs
`(i-mk,k)`, so cancellation cannot erase this exact image condition.
A finite sum belongs to the source image if and only if each of its
Laurent rows satisfies the stated divisibility. In particular every
nonzero negative row has positive `w` degree and is divisible by `w`.
This is an equivalence for the complete subring, not a necessary-only
heuristic inferred from sample monomials.

The coordinate Jacobian is exactly `u^m`. For rows `u^i f_i(w)` and
`u^L g_L(w)`, their contribution to the original bracket is

    u^(i+L+m-1) [i f_i g_L'-L f_i' g_L].

I checked the sign and the exponent directly from the original source
monomial bracket. If all rows of a proposed mate have exponent `L>=0`,
then every nonzero contribution has exponent at least `m`. The only
possible pair with `i+L=0` is `i=L=0`, whose coefficient is zero.
A mate must therefore have a negative least exponent. This preliminary
step correctly handles a leading constant, a constant proposed mate,
and the corner `m=1`; it does not divide by a possibly zero leading
exponent.

## 2. The constant-row proof has no hidden cancellation

If the constant row `f(w)` of `Fhat` is nonconstant and `L<0` is the
least mate exponent, the unique first bracket term is

    -L f'(w)g_L(w) u^(L+m-1).

Its coefficient is nonzero in characteristic zero. Every higher row of
either polynomial has strictly greater exponent, so the term cannot be
cancelled by a later correction.

Equality to one would force `L=1-m`. For `m=1`, this contradicts the
already established `L<0`. For `m>=2`, the source-image condition makes
`g_L` divisible by `w`; its product with the polynomial `f'` cannot be a
nonzero constant. This argument covers every nonconstant polynomial
`f`, including degree one, without first discarding its constant term.
Thus the constant row really is a scalar in `K`.

## 3. All degrees of the first row and the exact initial mate

After subtracting that scalar, write
`Fhat=u g(w)+O(u^2)`.

If `g=0`, the original polynomial `F-f0` is divisible by `u^2` after
substitution. Evaluation at `u=0` makes both of its original partial
derivatives zero in `K[t]`, so its Jacobian with any polynomial cannot
be one. This is a polynomial-ring argument and works over any
characteristic-zero field; it does not assume algebraic closure or apply
to a rational mate with a pole on that line.

If `g` is nonconstant, let `h=deg g>=1` and `j=deg g_L>=1`. The first
bracket coefficient is

    g g_L'-L g' g_L.

Its degree is exactly `h+j-1`, and its leading coefficient is the
product of the two leading coefficients times `j-Lh`. Since `L<0`,
this integer is strictly positive and therefore nonzero in `K`.
This pays noncancellation for arbitrary polynomial `g`, rather than just
an affine first row. Again all omitted rows occur strictly later.
The exponent would have to be zero, giving `L=-m`, but the coefficient
still has degree at least one. It cannot equal one.

The only remaining case is `g=beta` with `beta!=0`. Its leading
coefficient is `beta g_L'`, which is nonzero because a nonzero negative
row is divisible by `w` and has positive degree. Thus

    L=-m, beta g_L'=1.

Characteristic zero gives `g_L=w/beta+C`; the source-image divisibility
forces `C=0`. The exact least exponent and exact leading mate row are
therefore proved, not merely bounded.

Finally, a source monomial that survives on `u=0` is `t^k`, with weight
`-mk`. The proved minimum weight excludes every `k>=2`. At `k=1` the
coefficient is exactly `1/beta`; the constant is free. Consequently
`G(0,t)=t/beta+gamma`. No higher row can change this conclusion, because
the source monomial map is injective. This checks the final assertion
of the theorem as well as its first-jet condition.

## 4. Square-prefix and actual DG applications

For `w=u^2t`, the displayed hypotheses
`ord N>=4`, `ord P,ord M>=2` ensure that the entire square-prefix
expression belongs to `K[u,w]`. The geometric source of those hypotheses
is explicitly separate from this theorem.

I independently recomputed the first two coefficient rows. At an active
fourfold, the constant row has nonzero `w^4` coefficient `N4^2`, so the
first gate excludes polynomial mates. When `ord N>=5`, its constant row
is

    (P2 w+Q0)^2+M2 w+R0.

Constancy gives `P2=M2=0`. If `ord N=5`, the next row has quadratic
coefficient `2Q0 N5`; hence `Q0=0`, then `M3=0`, and its nonzero
constant requirement is `R1!=0`. If `ord N>=6`, the next row is

    (2Q0 P3+M3)w+2Q0 Q1+R1.

This gives exactly `M3=-2P3 Q0` and `2Q0 Q1+R1!=0`. Higher numerator
coefficients affect only later `u` rows, so these computations retain
all degrees, rather than assuming an exact low-degree truncation.

For the actual all-p sixfold equations in the primary, I also checked
in the original `(u,t)` variables that the tuned derivative is

    F_u(0,t)=-2p(lambda+2dA).

The independently reconstructed chart rows give

    c=nu=0, mu=-2bd,
    beta=-2p(lambda+2dA)!=0.

Thus the declared `p=0` polynomial exclusion follows immediately. The
three degree-four section rows and both degree-two linear-section rows
are verified with `x=u+p`; the argument does not assume translation is
an automorphism of the surface. The application is to the explicit
coefficient class with its stated jets and makes no unproved geometric
entry claim or rational exclusion.

## 5. Positives, hostile controls and scope

The polynomial positives `F=f0+beta u`, `G=t/beta+B(u)` work for every
`m`, every nonzero `beta`, and arbitrary polynomial `B`. They verify the
exact initial mate row, including its unavoidable `t/beta` restriction
on the original line.

Both rational hostile formulas have the stated signs. For all `m>=1`,

    F=u^(m+1)t, G=1/(m u^m), J(F,G)=1.

Its leading negative row has a nonzero constant coefficient, violating
the polynomial source-image condition. For `m>=2`, the second formula
`F=u^m t`, `G=u^(1-m)/(m-1)` similarly has Jacobian one and violates
that condition. The first example covers `m=1`; the second is not used
there. These controls demonstrate why the theorem cannot be promoted to
arbitrary rational or formal mates.

The polynomial `F=u+u^2` satisfies the asserted first-jet shape but has
zero source gradient along `u=-1/2`, so the necessary condition is not
sufficient. Its use is legitimate over every characteristic-zero field.
The source and proof do not assert a stronger classification.

The cited positive-weight canon is only an antecedent. Every argument
needed for the present signed-weight theorem is included explicitly;
there is no reliance on a pending geometric `(6,2)` classification.

## 6. Frozen source audit and exact replay

The source was read in full. It uses explicit exception checks that stay
active under `-O` and imports no inherited mathematical implementation.
Its monomial universe is exactly `m=1..6`, `L=-12..5`, `w` exponent
`0..6`, giving 756 image controls. The separate leading-degree controls
use `h,j=1..6`, `L=-12..-1`, giving 432 parameter rows with both degree
and noncancellation checks. Original-source higher rows are retained in
its bracket controls, and `m=1` is treated separately. These finite
universes corroborate the identities; the all-m and all-degree proof is
Sections 1--3 above and the primary's corresponding argument.

Reproduce using

```sh
python3 -B 04-computation/planar_jc48_sep08_polynomial_weight_gate.py
python3 -B -O 04-computation/planar_jc48_sep08_polynomial_weight_gate.py
```

Both independent runs pass **1,947 gates**, matching every byte of
[the frozen 553-byte output](planar_jc48_sep08_polynomial_weight_gate.out).
All supplied pins were independently recomputed.

| Artifact audited | Bytes | SHA-256 |
| --- | ---: | --- |
| Primary before promotion | 11,099 | `eb14dc46fae493fa76007797fdd725d63ccf98473276fb76e13d86d97917c014` |
| Source | 5,982 | `885b2c01b2498c4dc991b91796ab1ee9365fc190cbd5699bd2ab76676c012509` |
| Frozen output and both audit replays | 553 | `a0a99f9fb55e5e5bfe1075e60d1f800e1a86ef02bb5cfc86aeb052e7f766a6e7` |

Semantic trace:
`7c6566cea7919b18480be705f73791af41e44f876d1d78d66837ab208b277b79`.
The additional three original-coordinate application checks did not
import producer functions. No primary, source or frozen output was
changed by this audit. Parent-owned promotion can use this acceptance.
