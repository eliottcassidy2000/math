# The oriented four-wheel, its collision correction, and the Apéry form

Status: **PROVED** elementary identities over every odd finite field;
**CITED** Ahlgren–Ono prime evaluation; **FINITE-EXACT** extension-field
controls. The general prime-power trace identity is **not promoted here**.

## 1. Inheritance and exact object

The concurrent THM-4652 and HYP-9274 files were reserved empty namespaces at
intake. They were not proof dependencies. The closest external mechanism is
Ahlgren–Ono's evaluation of a quadratic-character fourth hypergeometric
moment. The important hostile is forgetting whether a graph sum permits
collisions. A second hostile is replacing a Frobenius trace over a field
extension by the Fourier coefficient with that field order as its index.
The sidecars retained below are orientation, field presentation, hub marking,
and the two opposite-rim collision strata.

Let $q$ be an odd prime power and let $\chi$ be the quadratic character of
$\mathbf F_q$, with $\chi(0)=0$. Our labelled wheel has hub $0$ and ordered
rim $(x_1,x_2,x_3,x_4)$. Define the **anchored injective signed sum**

$$
W(q)=\sum_{0,x_1,x_2,x_3,x_4\ \mathrm{distinct}}
\chi\!\left[x_1x_2x_3x_4
(x_1-x_2)(x_2-x_3)(x_3-x_4)(x_4-x_1)\right].                 \tag{1}
$$

Let $A(q)$ be the same sum over all four rim variables. Neither expression
is an unsigned embedding count. Allowing the hub to vary multiplies each
sum by $q$, by translation. Reversing $e$ of the eight displayed difference
factors multiplies the result by $\chi(-1)^e$. In particular, orienting the
rim as $a\to c,a\to d,b\to c,b\to d$ on its $K_{2,2}$ presentation reverses
two edges relative to a directed rim cycle, so gives the same sum.

For $q\equiv3\pmod4$ the character determines a Paley tournament; for
$q\equiv1\pmod4$ it determines a symmetric residue relation. The algebra
below applies in both cases; the tournament interpretation does not.

The concept board is: oriented graph weight, collision partition, Legendre
trace moment, modular coefficient, field extension, and ordered-run
context. No equality between the resulting statistics identifies the
underlying objects or their witnesses.

## 2. Exact all-map versus injective identity

Define

$$
S(q)=\sum_{y\in\mathbf F_q}\chi(y)
       \left(\sum_{x\in\mathbf F_q}\chi[x(x-1)(x-y)]\right)^2. \tag{2}
$$

Terms with $x_1=0$ vanish. Scaling all four rim variables by $x_1^{-1}$
changes the weight by an eighth power, which has character one. Thus fix
$x_1=1$ and multiply by $q-1$. Writing the other vertices as $(x,y,z)$,
the weight factors as

$$
\chi(y)\chi[x(x-1)(x-y)]\chi[z(z-1)(z-y)].
$$

The two sign changes used in this rearrangement cancel. Summing proves
$A(q)=(q-1)S(q)$.

Every collision with the hub or an adjacent rim vertex has weight zero.
The only remaining noninjective strata are $y=1$ and $z=x$:

* On $y=1$, the surviving weight is $\chi(xz)$, with $x,z\ne0,1$;
  its sum is $(-1)^2=1$.
* On $z=x$, the surviving weight is $\chi(y)$, with $x\ne0,1$ and
  $y\ne0,x$; its sum is $-\sum_{x\ne0,1}\chi(x)=1$.
* Their intersection has $q-2$ terms, all of weight one.

Inclusion–exclusion therefore gives, for every odd finite field,

$$
\boxed{A(q)-W(q)=(q-1)(4-q)},\qquad
\boxed{W(q)=(q-1)(S(q)+q-4)}.                              \tag{3}
$$

The collision contribution can be negative: at $q=7$, $A=-186$ whereas
$W=-168$. Consequently subtraction of collisions in a signed sum is not
the same operation as subtracting a positive count. For example,
$\chi(d)^2=1$ only on the retained guard $d\ne0$; it is zero at $d=0$.

## 3. A classical theorem closes the prime identity

Write

$$
f(z)=\eta(2z)^4\eta(4z)^4
 =q\prod_{m\ge1}(1-q^{2m})^4(1-q^{4m})^4
 =\sum_{n\ge1}a(n)q^n.                                   \tag{4}
$$

Ahlgren–Ono, *A Gaussian hypergeometric series evaluation and Apéry number
congruences*, J. reine angew. Math. 518 (2000), 187–212, define this form in
(1.5). Their Lemma 2.2 identifies (2), for an odd prime $p$, with
$p^3\,{}_4F_3(1)_p$; Theorem 6 evaluates it as $-a(p)-p$.
[Primary author PDF, pp. 5–7](https://uva.theopenscholar.com/files/ken-ono/files/049.pdf).

Combining that cited theorem with the elementary collision computation,

$$
\boxed{W(p)=-(p-1)(a(p)+4)}\quad\text{for every odd prime }p. \tag{5}
$$

Thus the equation is not restricted to a finite list of tested primes.
The new work here is the explicit wheel-to-moment interface and its
retained normalization, not a new modularity theorem.

For $\mathcal A(n)=\sum_{j=0}^n\binom nj^2\binom{n+j}j^2$, the same
paper's Theorems 3 and 5 give the following useful transported checks:

$$
\frac{W(p)}{p-1}\equiv-\mathcal A((p-1)/2)-4\pmod{p^2},
\qquad
\frac{W(p)}{p-1}\equiv-p^3-5\pmod{32}.                    \tag{6}
$$

Its quoted coefficient bound also gives
$|W(p)+4(p-1)|\le2(p-1)p^{3/2}$.

An independent check uses the four-square coefficient formula displayed
on p. 7 of that paper. A second geometric check counts actual nonzero
quadruples with $\sum_i(x_i+x_i^{-1})=0$; Ahlgren–Ono,
*Modularity of a certain Calabi–Yau threefold*, Theorem 1, relates this
prime-field count $N(p)$ by $a(p)=p^3-2p^2-7-N(p)$.
[Primary author PDF](https://uva.theopenscholar.com/files/ken-ono/files/050.pdf).
The checker counts this equation directly, without using $a(p)$.

## 4. Prime powers retain a different coordinate

All of (1)–(3) remain proved over $\mathbf F_{p^r}$. A natural extension
of (5) would use the trace of the modular local factor:

$$
t_0=2,\quad t_1=a(p),\quad
t_r=a(p)t_{r-1}-p^3t_{r-2}.                               \tag{7}
$$

This differs from the Fourier coefficient recurrence, whose initial term
is $a(1)=1$. In particular $t_2=a(p)^2-2p^3$, while
$a(p^2)=a(p)^2-p^3$.

The smallest tested extension is decisive. In $\mathbf F_9$,

$$
S(9)=29,\quad W(9)=272,\quad a(9)=-11,\quad t_2=-38.
$$

Substituting $a(9)$ into (5) gives $56$, not $272$. Substituting $t_2$
gives $272$. The checker independently constructs the field and counts
the normalized wheel triples; it also uses a second irreducible polynomial
for $\mathbf F_9$ to verify presentation independence of the value.

The eight tested extension orders are $9,27,81,25,125,49,121,169$.
All satisfy $S(p^r)=-p^r-t_r$ in these exact computations. This packet
does **not** turn those tests into an all-prime-power theorem: that would
require a cited or derived full Frobenius identity with its boundary
contributions, beyond a theorem stated only for primes. The elementary
wheel evaluator remains complete for every supplied finite field.

## 5. Interface and next obligation

The source is a labelled, oriented wheel over a marked finite field. The
target is a one-parameter character moment, then a prime-index modular
coefficient. Scaling and translation preserve the signed weight, but
forget the individual vertex assignments. The collision partition restores
the injective predicate; the field degree prevents confusing $a(p^r)$
with $t_r$; the orientation count restores the sign.

For an ordered-run or Hamiltonian expansion, a coefficient multiplying
this wheel must still be derived with the labelled trail placements and
injectivity guards. Equation (5) supplies its value once that coefficient
and normalization are proved. It does not itself count Hamiltonian paths.
In particular, a direct quotient that is not a wheel may reduce to one
after eliminating degree-two vertices; graph-isomorphism inspection alone
does not prove a higher run moment is polynomial in the field order.

The script exports `Field`, `moment`, `wheel`, `prime_wheel`, and a
source-bound `WheelReceipt`. `audit_receipt` checks exact integer types,
the complete field presentation, and recomputed signed counts. The
Frobenius routine is deliberately separate from the proved wheel routine.
Malformed field moduli, Boolean/float aliases, changed counts, and a receipt
from a different field presentation are rejected.

## 6. Reproduction and finite universe

Run:

```text
python -B -X utf8 04-computation/experiments/wheel_apery_extension_20261014h.py
python -B -O -X utf8 04-computation/experiments/wheel_apery_extension_20261014h.py
```

The declared controls include every odd prime at most 199; literal
normalized wheel triples for primes at most 31 and for extension orders
at most 49; eight extension fields; coefficient comparison for every
$1\le n\le199$; all translation/scaling assignments for $p=3,5,7$;
the two opposite-rim collision sums; and 18 malformed-interface controls.
Field multiplication and character identities are checked independently
of the modular formula. The extension comparisons are finite controls,
not evidence of a stronger quantified claim.

Normal and optimized runs agree exactly after LF normalization:
**71,191 checks**, JSON SHA256
`0def2094b474edf3e6e065cfe8038165a952ae1c04d42ef36c4a16be33a4f922`.
The saved [JSON](wheel_apery_extension_20261014h.json) is the complete
finite output; the [script](../../04-computation/experiments/wheel_apery_extension_20261014h.py)
contains the source-bound evaluator and independent controls.
