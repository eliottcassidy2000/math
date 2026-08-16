# The full Fourier--Dirichlet residue table recovers every Jordan rung

**Status: elementary proved synthesis, not canon and not a promotion of
RESERVED THM-3486.**  THM-3485 is the proved recurrence input.  THM-3486 is
currently a provisional proof candidate and remains outside the proved graph.
The calculation below independently audits its analytic core and records a
strictly stronger, characterwise residue statement for its owner to assess.

## 1. Inheritance and the missing operation

The closest proved algebraic mechanism is THM-3485: if

```text
a_n=P_r(n),                    n=r mod p,
Q_j(t)=1/p sum_r zeta^(-jr)P_r(t),
```

then the minimal tail shift polynomial is

```text
chi_a(x)=product_(j:Q_j!=0)(x-zeta^j)^(deg Q_j+1).       (1)
```

THM-3359 is the closest proved analytic mechanism: a periodic support has a
Hurwitz-zeta Dirichlet series whose residue at one is its cycle density.
RESERVED THM-3486 combines the top trivial colour of (1) with that harmonic
residue.  Its sharp hostile is

```text
a_n=(-1)^n n^d.                                          (2)
```

The ordinary critical harmonic coefficient is zero, although (2) carries a
full `(x+1)^(d+1)` Jordan block.  The corrected near miss is therefore not to
ask more of the trivial scalar residue.  The least-used sidecar is the
**complete character bank before taking residues**.

## 2. Exact meromorphic transform

Assume the lane coefficients are complex, write

```text
P_r(t)=sum_(m=0)^d c_(r,m)t^m,
```

and fix a primitive `p`th root `zeta`.  For every character `0<=j<p`, define
initially in the half-plane of absolute convergence

```text
D_j(s)=sum_(n>=1) zeta^(-jn) a_n/n^s,
Re(s)>d+1.                                                (3)
```

For `r in {0,...,p-1}`, let `r+` be its representative in `{1,...,p}`; thus
`0+=p`.  Splitting (3) by lanes gives the finite identity

```text
D_j(s)=sum_(r=0)^(p-1) sum_(m=0)^d
 zeta^(-jr)c_(r,m)p^(m-s) HurwitzZeta(s-m,r+/p).          (4)
```

Equation (4) is a meromorphic continuation to the whole `s`-plane.  Its only
possible poles are the simple poles `s=1,...,d+1`.  Since Hurwitz zeta has
residue one at argument one,

```text
boxed:
Res_(s=m+1) D_j(s)
  =1/p sum_(r=0)^(p-1) zeta^(-jr)c_(r,m)
  =[t^m]Q_j(t).                                          (5)
```

This proves the complete tomography statement: the labelled pole-residue
table

```text
R=(Res_(s=m+1)D_j(s))_(0<=j<p,0<=m<=d)                  (6)
```

is exactly the two-dimensional Fourier coefficient array of the word.  No
asymptotic inference or finite recurrence fit is involved.

## 3. Recurrence reconstruction

For a fixed colour `j`, equation (5) determines whether `Q_j` is zero and, if
not, its exact degree:

```text
d_j=max{m:Res_(s=m+1)D_j(s)!=0}.                         (7)
```

Substitution into proved THM-3485 gives

```text
boxed:
chi_a(x)=product_(j:some residue in row j is nonzero)
                 (x-zeta^j)^(1+max{m:Res_(m+1)D_j!=0}).  (8)
```

Thus the complete twisted Dirichlet pole packet and the finite shift-module
Jordan profile determine one another once the period and character gauge are
retained.  Over rational lanes the Galois grouping of the residue rows gives
THM-3485's cyclotomic descent.

The word “complete” in this paragraph is scoped to an exactly
periodic-polynomial tail.  An arbitrary finite prefix contributes a Dirichlet
polynomial, which is entire and hence invisible to (6).  The residues recover
the tail recurrence, not the transient prefix.

## 4. Relation to the critical harmonic limit

At the top degree no meromorphic continuation is needed.  For every `j`,
ordinary Abelian comparison gives

```text
lim_(N->infinity) 1/log N
  sum_(n<=N) zeta^(-jn)a_n/n^(d+1)
 = [t^d]Q_j.                                             (9)
```

The matching colour contributes its coefficient times `H_N`; every other
top colour contributes a bounded root-of-unity harmonic sum, and lower layers
are absolutely summable.  THM-3486 is precisely the case `j=0` of (9).

At degree zero, (5) also sharpens THM-3359 in a different direction.  The
untwisted residue is only the accepted-cycle density.  The complete `p`
character residues are the finite Fourier transform of the eventual Boolean
cycle, so inverse Fourier transform recovers every accepted residue class
relative to the declared cycle origin.  It still cannot recover the finite
transient or whatever graph, automaton, or tournament generated that Boolean
readout.

If one wishes to avoid analytic continuation at the lower poles, (9) can be
iterated.  After recovering all degree-`m` residues, subtract

```text
n^m sum_j Res_(s=m+1)D_j(s) zeta^(jn)                   (10)
```

from the current residual word and repeat at degree `m-1`.  Fourier inversion
makes (10) exactly the degree-`m` lane layer.  This demodulate--measure--deflate
algorithm reaches the same table (6), but (4)--(5) explain why it works in one
line.

## 5. Exact hostiles and boundaries

1. **Trivial-colour blindness.**  For (2), `D_0` is regular at `s=d+1`,
   while `D_1` has residue one there.  Any observer using only the ordinary
   harmonic sum loses the entire nontrivial block.
2. **Composite-period missing orders.**  A period-four word that depends only
   on parity has zero residue rows at the primitive order-four characters.
   The full bank detects colours of orders one and two only, exactly as in
   THM-3485's composite-period hostile.
3. **Period refinement.**  Replacing `p` by a multiple introduces labelled
   zero rows and permutes the surviving roots, but (8) is unchanged after zero
   rows are omitted.
4. **Finite-prefix loss.**  Adding any finitely supported sequence changes
   every `D_j` by an entire Dirichlet polynomial.  This is the sharp destroyed
   coordinate.
5. **Character gauge.**  Replacing `zeta` by another primitive root permutes
   the rows of (6).  The labelled table changes; the product (8) does not.
6. **Zero word.**  Every residue vanishes and the tail recurrence is one, in
   the convention of THM-3485.

For THM-3484's ternary determinant lanes, the degree-seven coefficient vector
is constant `(-16384,-16384,-16384)`.  The top residue row is therefore

```text
(-16384,0,0),                                            (11)
```

while the degree-six nontrivial rows are nonzero.  Equations (7)--(8) recover
the proved exponent packet `(8,7,7)` without a Hankel fit.

## 6. Connection contract and cross-frontier lesson

| field | exact content |
|---|---|
| source | a complex periodic-polynomial word with labelled period `p` |
| operation | twist by each character, form its Dirichlet series, take every integer pole residue |
| target | the complete Fourier coefficient array, then THM-3485's minimal shift polynomial |
| preserved | period lanes, character labels, every polynomial layer, every nonzero Jordan depth |
| destroyed | finite prefix, address provenance, graph/tournament/LRC/JC meaning |
| load-bearing sidecars | the complete character bank, period, primitive-root gauge, meromorphic continuation or iterative deflation |
| cheapest hostile | `(-1)^n n^d`, which kills the trivial residue but not the nontrivial one |

The same operation-level warning appears in the live LRC endpoint experiment:
a quotient or trivial character can identify a forced bridge even when a
refined character separates it.  That is inspiration, not an object-level
identification.  The LRC bank concerns endpoint-response characters on a
finite relation quotient; (3) concerns a Dirichlet transform of an index
word.  No LRC current, factorial moment, Jacobian divisor, Rule-30 complexity,
or LRC(14) conclusion follows from (5).

The reusable move is narrower: when a scalar critical transform sees only one
isotypic component, demodulate before taking the critical residue.  Retain the
full character label until the native target predicate has been reconstructed.
