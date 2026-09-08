# An exact entire inverse hierarchy without a rational mate

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This theorem concerns the complete all-finite binary `(7,1)` entry on the
fixed DG surface. It proves a formal equivalence and a sharp failure of
rational algebraization, not a Jacobian-conjecture theorem.

## 1. Complete coefficient contract and the all-index statement

Use the original surface and volume

    W=(P1_x x P1_z)\{z=x^2},  t=1/(z-x^2),
    omega=dx wedge dt.

Let `H in L2`, `L in L1`, `F=H^2+L`, with two distinct finite zeros of
multiplicities seven and one in the leading binary octic. Normalize the
nonzero separation and leading scalar by actual source/surface scalings,
and write `u=x-p`. The translation is an ordinary polynomial coordinate;
no automorphism of `W` is asserted. The complete original global rules in
[the DG filtration](planar_jc48_sep08_dg_quadratic.md) give

    N=u^7(u-1),
    P=2u^6-(4p+2)u^5+a u^4+b u^3+c u^2+d u+e,
    Q=u^4-(4p+1)u^3+(a+4p^2)u^2+(b-2ap+4p^2)u+f,
    M=m4 u^4+m3 u^3+m2 u^2+m1 u+m0,
    R=m4 u^2+(m3-2pm4)u+h,
    H=N t^2+P t+Q,  L=M t+R.                        (1)

In particular `P0=e` and `P1=d` are arbitrary and retained. No finite
active-point, first-jet, source-submersion or mate hypothesis is imposed
on (1).

Put `kappa=sqrt(N)`, `K=C(u)(kappa)`, and choose the unique formal inverse

    F(u,T(u,v))=v^-4,
    T=kappa^-1 v^-1+T0+T1 v+T2 v^2+... .             (2)

**Theorem.** The following are equivalent for the entire
coefficient universe (1):

1. `M(1)=0`.
2. Every differential `T_j du`, for `j=-1` and all integers `j>=1`,
   has a primitive in `K`.
3. There is a formal Laurent series
   `Gformal in C(u)((1/t))` with `J(F,Gformal)=1`.

When these hold, a formal mate can be chosen in
`t^-3 C(u)[[t^-1]]`. There is no exactness condition on `T0`.

These equivalent formal properties do **not** imply a rational mate.
The literal global example in Section 6 has `M(1)=0` and no rational
mate. Every formal mate for that example is therefore transcendental
over `C(u,t)`; it is not merely an inconvenient rational expression.

The inherited mechanism is
[the radical leading-coefficient primitive](planar_jc48_sep08_leading_exactness.md),
with its explicit formal substitution and derivation contract. The
[all-finite (5,3) proof](planar_jc48_sep08_five_three.md) supplies the
exact opposite-pair identity, and the
[unit-M boundary analysis](planar_jc48_sep08_quartic_common_root.md)
supplies the actual hostile. The closest corrected near miss is to
assume that increasing inverse order must eventually detect a global
obstruction. The present all-index induction shows exactly why it need
not.

The live concepts are complete global rows, two-point ramified parity,
a closed recurrence ring, coefficientwise trace, and algebraization.
The source is a quartic polynomial with its actual lower rows; the
target is a compatible formal primitive equation. The map retains every
inverse coefficient and the source derivation, but forgets whether the
series is algebraic or rational. Section 6 keeps the compact actual fibre
as the missing sidecar. No external priority claim is made for classical
formal-inverse or residue methods.

## 2. Globality makes the entire recursion bounded at infinity

Center the quadratic by `t=y-P/(2N)` and set

    D0=Q-P^2/(4N), E0=R-MP/(2N), C0=M/kappa,
    B0=M^2/N=C0^2.

The full equations (1), including all low rows, give

    deg(P^2-4NQ)<=8,
    deg(2NR-MP)<=8,
    deg(M^2)<=8.                                    (3)

The first cancellation removes degrees twelve through nine; the second
removes degrees ten and nine. Since `deg N=8`, each of `D0,E0,C0,B0`
is bounded at both points over `u=infinity`. This does not follow from
a free quartic prefix: the induced `Q,R` rows are essential.

Set `W=kappa v y`. The exact centered inverse equation is

    W^4+2D0 v^2 W^2+C0 v^3 W+(D0^2+E0)v^4=1,
    W(0)=1.                                        (4)

Writing `W=1+sum_(n>=1) w_n v^n`, the coefficient at each new order is
`4w_n` plus a polynomial in earlier coefficients and `D0,E0,C0`.
Thus, for **every** `n`,

    w_n in C[D0,E0,C0].

It follows at both infinities that

    T_j=O(u^-4),  j=-1 and every j>=1.               (5)

Indeed `T_j=w_(j+1)/kappa` for these indices and `kappa~±u^4`.
Centering changes only `T0`; its size is irrelevant. Apart from the two
infinities, possible poles of these coefficients lie only over `u=0,1`,
since (4) introduces no other denominator.

## 3. Odd coefficients and the explicit genus-zero field

The field has the rational parametrization

    z^2=(u-1)/u,  u=1/(1-z^2),  kappa=u^4 z.         (6)

The points `u=0,1` are ramified, corresponding to `z=infinity,0`.
There are two unramified infinity points, `z=+1,-1`.
The deck involution is `z->-z`, and sends `kappa` to `-kappa`.
Uniqueness of the formal inverse shows

    tau(T(u,v))=T(u,-v), hence tau(T_j)=(-1)^j T_j.  (7)

For odd `j`, including `j=-1`, the differential `T_j du` is therefore
anti-invariant. At either ramified point choose a local parameter on
which the involution acts by minus one. An anti-invariant differential
has only even powers in its Laurent coefficient multiplying that
parameter's differential. It cannot have a term of exponent `-1`, so
its residue vanishes. At the two infinities (5) makes the differential
regular. There are no other poles.

Since (6) is rational, a meromorphic differential with all residues zero
has a rational primitive: this is ordinary partial-fraction integration
in `z`. Thus every odd-index differential is exact, independently of
`M(1)`. For the leading one a literal primitive is

    du/kappa = d[2z-(4/3)z^3+(2/5)z^5].             (8)

The argument is about the whole normalized field. It does not replace
its two infinity points with one or mistake the field for a disconnected
radical presentation.

## 4. All even coefficients, necessity, and the single surviving condition

The four centered inverse branches are `y(v),y(-v),y(iv),y(-iv)`.
Their sum is zero because the centered quartic has no cubic term.
Consequently every centered coefficient whose index is divisible by four
is zero. In particular

    T_(4ell)=0, ell>=1.                             (9)

For the remaining even indices the exact opposite-pair formula is

    T_(4ell+2)=-M/(4N) [s^ell]V(s), ell>=0,
    (1-E0s)V^2+(D0 B0/4)s^2 V^4
                  +(B0^2/64)s^3 V^6=1,
    V(0)=1.                                        (10)

This follows by factoring the depressed quartic into its opposite-root
quadratics; the formula and its normalizations are also recorded in
Section 5 of [the all-finite (5,3) proof](planar_jc48_sep08_five_three.md).
Its coefficient recursion has leading coefficient two, so at **every**
order

    [s^ell]V in C[E0,D0 B0,B0^2].                   (11)

Now assume `M(1)=0`. At the simple root `u=1`,

    D0 has pole order at most one,
    E0 is regular,
    B0=M^2/N vanishes to order at least one.

Hence all three generators in (11) are regular there. The factor `M/N`
in (10) is regular as well. Every even coefficient is rational in `u`,
has no pole at `u=1`, and by (5) is `O(u^-4)` at infinity. Its only
possible pole on `P1_u` is at zero. The residue theorem makes that residue
zero. Rational partial fractions therefore give an exact differential
already in `C(u)` for every such index. Together with Sections 2--3,
this proves the all-index sufficiency in statement 2.

Conversely `T2=-M/(4N)`. Its rational residue at `u=1` is `-M(1)/4`;
on the ramified point of (6), the residue is twice that value. If `T2 du`
is exact in `K`, this residue is zero and `M(1)=0`. Thus statements 1
and 2 are equivalent. Neither a finite scan of indices nor a bound on
the number of future coefficients enters this proof.

## 5. Formal Laurent mates and descent to the original coefficient field

The formal Hensel and substitution argument used in
[leading exactness](planar_jc48_sep08_leading_exactness.md) applies to
(2). The derivation `partial_u` extends uniquely to `K` in characteristic
zero. If a rational or formal Laurent mate is written in the inverse
coordinate as `Gtilde`, the chain rule gives

    partial_u Gtilde |_v = (1/4)v^5 partial_v T.     (12)

For the forward implication, choose primitives in `K` using statement 2
and define coefficients by

    g_(j+4)'=(j/4)T_j, j=-1 or j>=1,
    g4=0,
    Gtilde=sum g_(j+4) v^(j+4) in v^3 K[[v]].       (13)

The index-zero term in the right-hand side of (12) vanishes, explaining
why `T0` imposes no condition. Equation (13) solves (12) coefficientwise.

Put `q=1/t`. The inverse series to `T` is

    v=F(u,1/q)^(-1/4) in q K[[q]],
    v=kappa^-1 q+O(q^2).

Composition of (13) is well-defined: every coefficient receives only
finitely many terms. The formal chain rule now gives
`J(F,Ghat)=1` in `K((q))`. The two derivations on this field are
`partial_u` holding `q` fixed and `partial_t=-q^2 partial_q`.
Both commute with the deck involution and with coefficientwise trace.
Consequently

    Gformal=(Ghat+tau(Ghat))/2 in q^3 C(u)[[q]]

still has `J(F,Gformal)=1`. This proves statement 3 from statement 2.
It does not assume any convergence.

Conversely, substitute `q=1/T(u,v)` into any formal Laurent mate in
`C(u)((q))`. Its lower bound on exponents makes this composition valid
in `K((v))`. Equation (12) and coefficient extraction then prove all
exactness assertions of statement 2. This establishes the full formal
equivalence.

A useful boundary follows from the same trace principle. If a formal
mate is algebraic over `C(u,t)`, let `E` be the finite field it generates.
The source derivations extend uniquely to `E`, and normalized field
trace commutes with them. The trace of the algebraic mate would be a
rational mate with the same constant Jacobian. Thus a polynomial with no
rational mate has no algebraic formal mate either. This observation will
apply to the actual example below.

## 6. A global hostile that passes the entire hierarchy

At `p=0`, take the literal functions

    H=u^3(u-1)(1+u^2t)^2,
    L=(u-1)t,
    F=H^2+L.                                       (14)

Their coefficients are exactly (1) with the appropriate lower parameters:

    N=u^7(u-1), P=2u^5(u-1), Q=u^3(u-1),
    M=u-1, R=0.

Thus `M(1)=0`, so all inverse coefficients in the theorem are exact and
formal Laurent mates exist. Both functions are genuinely global on `W`.
In its second chart,

    H=(1-r)b_D^2,
    L=-r(1-r)(1+r^2 b_D),
    F|D=b_D^4.                                     (15)

Nevertheless **there is no rational mate**. At the finite sevenfold point,
`M(0)=-1` is a unit and the normal order of `H` is five. It is an
unbalanced case, `7!=3*5`, of the proved unit-M analysis. More explicitly,
the three determinations `s~u^(14/3)` form one actual normalized branch
with `(ord u,ord s)=(3,14)` and regular differential order two.
At the finite simple point the shared-simple-root lemma gives regularity.
There is no zero of the octic at infinity. At every generic point on `D`,
the actual factor `r^2` in the source volume gives a zero rather than a
pole of the relative form.

On the original smooth generic fibre away from `D`, the relative form is
a unit. Thus on every compact normalized generic component the complete
form `eta=omega/dF` is holomorphic and nonzero. A rational mate would
restrict to a meromorphic primitive on such a component. A pole of that
primitive would make its derivative have a pole, impossible here.
The primitive would therefore be constant on the compact component,
contradicting the nonzero form. This argument does not assume generic
connectedness or infer a rational map from formal coefficients.

The first failed implication is now exact: all coefficient primitives,
even with a descended formal Laurent mate, do not imply rational
algebraization. By Section 5, every formal mate for (14) is necessarily
transcendental over the original rational function field.

For a concrete finite check of the construction, in the coordinate (6)
set

    g3=-(2z-(4/3)z^3+(2/5)z^5)/4,
    g6=1/(48u^6),
    g7=(3/4)[z^3/6-3z^5/10+3z^7/14-z^9/18].

For (14), `T1=0`; these give the required derivatives of the `v^3,v^6,v^7`
coefficients in (13). Expanding them back through order `q^7` gives
rational functions of `u` and satisfies `J(F,Gtrunc)=1+O(q^5)` exactly.
This is a control of descent and depth. The all-order construction is
proved in Sections 2--5, not extrapolated from this truncation.

## 7. Exact verification universe and status boundary

The source is
[planar_jc48_sep08_seven_one_formal.py](../../04-computation/planar_jc48_sep08_seven_one_formal.py).
It checks the complete original-coordinate section constraints and all
cancelled leading rows, the retained arbitrary `P0,P1`, the full
`M(1)=0` factor space, regularity of the three universal recurrence-ring
generators, the actual quadratic parametrization, finite abstract
recursion/parity controls, both global charts of (14), and the descended
formal-Jacobian control through `q^7`. The minimal failed-condition
control `M=1` has actual ramified residue `-1/2` at the simple point.
All gates remain active under `-O`.

The unbounded quantifiers are supplied by the polynomial-ring inductions
and the parity/residue proofs, not by the finite controls. The proof does
not claim a rational-mate criterion for all-finite `(7,1)` or extend the
surface-coordinate normalization to arbitrary projective boundary maps.
The [independent full audit](planar_jc48_sep08_seven_one_formal_audit.md)
accepts all-index exactness, the complete regular recurrence ring,
formal Laurent descent and the necessarily transcendental hostile.
Both93-gate replays and independently reconstructed resolvent/truncation
controls agree. This is a formal criterion, not rational algebraization.


Reproduce from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_seven_one_formal.py
python3 -B -O 04-computation/planar_jc48_sep08_seven_one_formal.py
```

Both executions pass **93 always-active exact gates** and reproduce all
407 bytes of [the frozen output](planar_jc48_sep08_seven_one_formal.out).

| Artifact | SHA-256 |
| --- | --- |
| Source, 7,342 bytes | `723dcb0a947295ceec0483755fcc8487476139ac69513063dd81319700fa8453` |
| Frozen output and both replays, 407 bytes | `dec45c98a41a9b4ab1443a16137ec231a9f44722149517e2c22fb76e98d6e522` |
| Semantic record | `af2ddce9f9de8d962728c9bdeb5ab5d4a0a4cad6e5ba21c06f6ea5c23f92ff2e` |
