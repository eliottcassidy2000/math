# Beyond elementary Jacobians: Cohn cores, factorial repair tails, and the odd/even cycle split

**Status: CITED + INDEPENDENTLY REPROVED for Wright's elementary-Jacobian
criterion; PROVED-SYMBOLIC for the Cohn transport calculation, the rectangle
lemma, and the cyclic parity law; FINITE-EXACT for the companion checks; CITED
for Cohn non-elementarity and the low-entry-degree decomposition; CONDITIONAL
for the motivating Weihrauch recovery implication; OPEN for every
counterexample construction proposed below.**  No planar Jacobian
counterexample and no proof of `JC(2)` is claimed.

The exact companion is
[`jc2_cohn_parity_cycle_repair.py`](../04-computation/jc2_cohn_parity_cycle_repair.py),
with frozen
[`output`](../05-knowledge/results/jc2_cohn_parity_cycle_repair.out).

## 1. Inheritance and the change of search language

The closest mechanism is Wright's 1978 weak Jacobian theorem: for a planar
Keller map, elementary factorability of its Jacobian matrix is equivalent to
polynomial invertibility
([Wright, *The amalgamated free product structure of GL_2...*](https://www.sciencedirect.com/science/article/pii/0022404987900041)).
The reduced-word argument developed independently in this session and in
[`jc2-small-rational-four-gauss-search-boxeph-20260821.md`](jc2-small-rational-four-gauss-search-boxeph-20260821.md):
reproves the direction needed here: if a determinant-one planar polynomial
Jacobian matrix lies in the elementary subgroup `E_2(k[x,y])`, then its
integrated map is tame.  Thus increasing the number of small shear parameters
can never find a counterexample.  The criterion is classical; the compact
leading-monomial proof below is the session's independent audit, not a novelty
claim.

The canonical hostile is Cohn's matrix

```text
C=[1+xy   x^2]
  [-y^2  1-xy],                 det C=1.                 (1)
```

Cohn proved that this matrix is not elementary over `k[x,y]`; see
[Cohn, *On the structure of the GL_2 of a ring* (1966)](https://www.numdam.org/item/PMIHES_1966__30__5_0/).
A recent decomposition theorem gives additional evidence that this is the
correct smallest non-elementary seed: over an algebraically closed
characteristic-zero field, a matrix in `SL_2(k[x,y])` with an entry of degree
at most two is elementary or is elementary-decorated Cohn type, after a
polynomial coordinate automorphism
([Chapovskyi--Kozachok--Petravchuk, Theorem 1](https://arxiv.org/abs/2412.03688)).
This literature result is used only as a search router, not as a dependency
of any new proof here.

The corrected near miss is the finite six-word search: its empty rational
boxes were real computations, but the all-length theorem explains them
without a height bound.  The least-used sidecar is the rectangle-coboundary
law from `LTT-096`: a binary `2 by n` table is of the form
`e_ij=x_i XOR y_j` exactly when its `2 by 2` rectangle parities vanish, and
the potentials are unique up to one global bit.

The new connection contract is therefore

```text
finite-choice rectangle rigidity -> monomial repair graph
source: size-n rectangular answer classes in 2 x n
target: coefficient channels of a transport/curl equation
map: row bit -> repair branch; column -> exponent channel
preserved: rectangle parity, global sign gauge, odd/even cycle monodromy
destroyed: Weihrauch uniformity, polynomial integrability, nonproperness
needed sidecar: exact curls plus the non-elementary SL_2/E_2 class
cheapest test: rank(I+cyclic shift), then exact coefficient transport
```

The concept board after this pull is:

| lane | object | invariant | decisive test | status |
|---|---|---|---|---|
| Anchor | elementary `DF` | class in `SL_2/E_2` | reduced-word leading monomial | elementary class closed |
| Niche | Cohn core `(1)` | non-elementary class | repair-curl cokernel | first one-factor repair has an infinite tail |
| Wildcard | recovery rectangles | global-bit quotient | cardinality of `A x B` | odd rigid, even balanced |
| Arithmetic | even exponent cycle | alternating holonomy | exact rational cycle solve | next search generator |
| Smooth bridge | exponential formal repair | all finite truncations | terminal residual | analytic solution, no polynomial solution |

## 2. Wright's elementary-Jacobian criterion, independently replayed

Put

```text
E_+(f)=[1 f;0 1],       E_-(f)=[1 0;f 1].              (2)
```

### Cited theorem and reproduced direction

Let `k` be a characteristic-zero field and let `F in k[x,y]^2` satisfy
`det DF=1`.  If `DF` is a finite product of matrices `(2)`, then `F` is a
tame polynomial automorphism.

### Proof mechanism

Merge adjacent factors of the same sign.  First suppose that the resulting
alternating word has no constant parameter and that its rightmost factor is
`E_+(f_1)`.  Write the continuant rows as

```text
v_0=(0,1),  v_1=(1,f_1),
v_j=v_(j-2)+f_j v_(j-1).                                (3)
```

With lexicographic order `x>y`, nonconstancy makes the product term in `(3)`
strictly dominate the skipped term.  For every `j>=2`,

```text
LM(v_j,1)=product_(i=2)^j LM(f_i),
LM(v_j,2)=product_(i=1)^j LM(f_i).                      (4)
```

One row of the full word is one of these continuant rows.  If `(f_1)_x` is
nonzero, write `LM(f_1)=x^e y^r`, `e>=1`.  In a row `(P,Q)`, equation `(4)`
gives `LM(Q)=LM(P)LM(f_1)`.  The leading monomial of `Q_x` strictly exceeds
every monomial of `P_y`: for `e>1` it has larger `x` exponent, and for `e=1`
it has the same `x` exponent and larger `y` exponent.  This remains true when
the leading term of `P` is independent of `y`, in which case it simply
vanishes from `P_y`.  Hence `P_y-Q_x` cannot be zero, contradicting that the
row is a gradient.

Thus `f_1 in k[y]`.  Choose `A'(y)=f_1` and
`T(x,y)=(x+A(y),y)`.  The chain rule gives

```text
D(F o T^-1)(z)=N(T^-1 z),                               (5)
```

where `DF=N E_+(f_1)`.  Equation `(5)` is an integrable elementary word one
factor shorter; substitution by `T^-1` preserves polynomiality and
nonconstancy.  The opposite rightmost sign is identical after exchanging
`x,y` and using lexicographic order `y>x`.

Constants do not create an escape.  Zeros delete and adjacent signs merge.
Endpoint constants are affine source or target gauges.  For an internal
`c in k*`, the exact Bruhat identity

```text
E_+(a)E_-(c)E_+(b)
 =E_+(a+c^-1) w_c E_+(b+c^-1),
w_c=[0 -c^-1;c 0]                                      (6)
```

reduces the elementary length after moving `w_c` to the source end using

```text
w_c E_+(h) w_c^-1=E_-(-c^2 h),
w_c E_-(h) w_c^-1=E_+(-h/c^2).                         (7)
```

Precomposition by the constant linear map absorbs `w_c`; parameter
substitution again preserves the hypotheses.  Induction ends at a constant
Jacobian matrix.  Every reduction used only affine maps and triangular
shears, so `F` is tame.  QED.

The failure boundary is load-bearing: `SL_2(k[x,y])` is strictly larger than
`E_2(k[x,y])`, and `(1)` is an explicit witness.  Suslin stability in matrix
size at least three does not remove this two-by-two obstruction.

## 3. The first Cohn repair produces factorially small rationals

The rows of `(1)` have curls

```text
(1+xy)_y-(x^2)_x=-x,
(-y^2)_y-(1-xy)_x=-y.                                  (8)
```

Left multiplication by `E_-(Z)` leaves the first row unchanged and changes
the second curl.  Closing that second row would require

```text
L(Z):=(1+xy)Z_y-x^2 Z_x-xZ=y.                          (9)
```

Give `x^i y^j` weight `j-i`.  Every term of `L` lowers this weight by one, so
the right side forces the weight-two ansatz

```text
Z=sum_(i>=0) a_i x^i y^(i+2).                          (10)
```

Coefficient comparison gives

```text
2a_0=1,                 (i+2)a_i+a_(i-1)=0,
a_i=(-1)^i/(i+2)!.                                      (11)
```

The degree-`D` truncation therefore satisfies exactly

```text
L(Z_D)=y+(-1)^D x^(D+1)y^(D+2)/(D+2)!.                 (12)
```

It is never a polynomial solution, however high the cutoff.  Formally the
unique weight-two solution is

```text
Z_hat=(exp(-xy)-1+xy)/x^2.                             (13)
```

This is the sharp correction to the initial “small rational coefficients”
intuition.  The correct structural seed really does generate tiny rational
numbers, but factorial denominators record an infinite transport tail rather
than a finite polynomial map.  A height-only search sees increasingly good
near misses and can mistake the shrinking coefficient in `(12)` for
convergence toward a polynomial solution.  Exact terminal-residue tests must
come before height ranking.

The smooth umbilic construction provides an illuminating but non-equivalent
comparison.  A flat exponential envelope can suppress every finite jet while
retaining a winding index.  Here the exponential occurs as the formal inverse
of an algebraic transport operator.  Smooth completion is available; the
polynomial category forbids the infinite tail.  No analytic completion of
`(9)` is a Jacobian-conjecture candidate.

### 3.1 One shear on each side never closes even one row

The path obstruction survives an arbitrary right upper decoration.

**Theorem.**  For every `W,Z in k[x,y]`, the lower row of

```text
E_-(Z) C E_+(W)                                        (13a)
```

is not closed.  Symmetrically, no `E_+(R) C E_-(U)` has a closed upper row.
Thus a non-elementary-core search needs more than one elementary correction
on at least one side before it can even reach both curl equations.

For `(13a)`, put

```text
A=1+xy,  P=-y^2+AZ,
Q=1-xy+x^2Z+WP.
```

Direct differentiation gives the exact identity

```text
P_y-Q_x=-y+L(Z)-partial_x(WP).                         (13b)
```

Suppose first that `deg W=d>=1` and `deg Z=N>=1`.  The unique term of degree
`d+N+1` in the closure equation is

```text
-partial_x(xy W_d Z_N),                                (13c)
```

which is nonzero in characteristic zero because its argument is nonzero and
divisible by `x`.  Hence one of `W,Z` must be constant.

If `Z=z` and `W` is nonconstant, `(13b)` becomes

```text
(y^2-z-zxy)W_x-zyW=xz+y.                              (13d)
```

For `z=0`, this demands `W_x=1/y`.  For `z!=0`, an `x`-degree `m>=2` has
top coefficient `-zy(m+1)w_m`; `m=1` forces `w_1=-1/(2y)`; and `m=0`
cannot supply `xz`.  All are impossible in `k[x,y]`.

It remains that `W=c` is constant.  If `c=0`, equations `(10)--(12)` apply.
For `c!=0`, the leading homogeneous equation for
`Z_N=x^N p(y/x)` is

```text
t(2+ct)p'=(N+1)(1+ct)p.                               (13e)
```

When `N` is even its nonzero solution has a half-integral exponent; when `N`
is odd it is a scalar multiple of
`[t(2+ct)]^((N+1)/2)`, whose `t`-degree is `N+1>N`.  Neither is the leading
form of a degree-`N` polynomial.  A constant `Z` also fails directly.  This
proves the lower statement; source/target exchange gives the upper one.

The two same-sign layouts die even earlier.  In
`E_-(Z) C E_-(U)`, the untouched top row would require
`x^2 U_y-x=0`, hence `U_y=1/x`.  In `E_+(R) C E_+(W)`, the untouched bottom
row would require `y^2 W_x-y=0`, hence `W_x=1/y`.  Therefore no decoration
with at most one elementary factor on each side of `C` can be a Jacobian
matrix.

## 4. The exact odd/even mechanism

### 4.1 Rectangles

More generally, a size-`n` rectangle in a `g by n` target has height dividing
`gcd(g,n)`.  Coprimality of gauge order and branch count forces a singleton
gauge row.  For `g=2`, let `R=A x B` be a nonempty rectangle in
`{0,1} x {0,...,n-1}` with `|R|=n`.  Since `|A|` is one or two,

```text
n odd:   |A|=1, |B|=n, so R is a full row;
n even:  either a full row or |A|=2, |B|=n/2.           (14)
```

There are exactly two full rows and, when `n` is even,
`binom(n,n/2)` balanced rectangles of the second type.  This proves the
finite geometry used in the motivating recovery draft.  It does **not** by
itself prove the stated Weihrauch implication: that conclusion also needs
the draft's finite-commitment theorem and the following faithful-incidence
hypotheses:

1. the raw oracle answer set, before the backward functional, is a rectangle;
2. it has exactly `n` distinct cells, rather than merely an `n`-element
   decoder image of a larger rectangle;
3. the committed decoder preserves tooth/stripe incidence bijectively; and
4. the commitment has recognizable positive evidence, or an effective
   exact-cardinality lower bound supporting negative enumeration.

Without exact cardinality, odd rigidity already fails for a four-cell
`2 by 2` subrectangle inside `2 by 3`.

There is a sharper orthogonal-partition theorem.  Suppose two rectangular
size-`n` tooth blocks partition `2 by n`, while `n` rectangular size-two
stripe blocks also partition it, and every tooth/stripe intersection is a
singleton.  For odd `n`, the teeth are exactly the two rows and the stripes
the `n` columns, up to row swap and column permutation.  For `n=2m`, there is
one further type:

```text
teeth = 2 x S and 2 x S^c,       |S|=m;
stripes = rowwise perfect matchings S -> S^c.           (14a)
```

Conversely every pair of rowwise perfect matchings realizes `(14a)`.  With
target coordinates fixed and partitions internally unlabeled, there are
`n!/2` such even geometries, plus the standard grid.  Their essential
sidecar is not parity alone but the cycle type of the relative matching
permutation `phi_1^-1 phi_0`.  At `n=2` this is only the transposed degenerate
grid; `n>=4` supports genuine matching holonomy.

A binary table satisfies all `2 by 2` rectangle parities precisely when

```text
e_ij=x_i XOR y_j.                                      (15)
```

The potentials in `(15)` are unique modulo simultaneous toggling of all
`x_i` and all `y_j`.  Thus the information loss is exactly one global bit,
not an arbitrary quotient.

### 4.2 Cycles

Rescale the Cohn coefficients by

```text
b_i=(i+2)! a_i.
```

The homogeneous part of `(11)` becomes

```text
b_i+b_(i-1)=0.                                         (16)
```

If `n` **abstract normalized** channels are cyclically identified, `(16)` is the matrix
`I+S_n`, where `S_n` is the backward cyclic shift.  Exact linear algebra gives

```text
det(I+S_n)=2, rank=n       for odd n;
det(I+S_n)=0, rank=n-1     for even n,                  (17)
```

and the even kernel is spanned by

```text
(1,-1,1,-1,...,-1).                                    (18)
```

Equations `(14)` and `(17)` expose the same unweighted support mechanism in
two representations.  Odd cardinality prevents a balanced two-row split and
odd unit monodromy kills the alternating mode.  Even cardinality permits two
rows over half the columns and even unit monodromy preserves one alternating
mode.  Its overall sign is the analogue of the global bit in `(15)`.

The multiplier sidecar is load-bearing.  A weighted cycle

```text
alpha_i c_i+c_(i-1)=0
```

has determinant

```text
product_i alpha_i-(-1)^n.                              (18a)
```

Thus an even support cycle has a kernel only when its multiplier holonomy is
`+1`; an odd cycle needs holonomy `-1`.  Factorial rescaling flattens the
**interior** of the Cohn path but not a hypothetical wrap seam.  The actual
first-`n` Cohn weights are `2,3,...,n+1`, so their product is `(n+1)!` and
their cyclic determinant is `(n+1)!-(-1)^n`, never zero.  The first small
rational even model is instead the reciprocal two-cycle with weights
`2,1/2`.

This is a precise algebraic bridge, but its scope must not be inflated:
rectangle classes are finite-choice answer sets, while exponent cycles are
coefficient-transport graphs.  The map preserves parity and the one-bit
gauge; it forgets computability uniformity on one side and multiplier
holonomy plus polynomial integrability on the other.

## 5. A new counterexample generator

The elementary theorem and `(12)--(18)` change the enumeration order.

1. **Choose a non-elementary core.**  Start with `(1)` or another certified
   representative of `SL_2(k[x,y])/E_2(k[x,y])`.  Do not enumerate a longer
   elementary word.
2. **Decorate without erasing the class.**  Multiply on the left and right by
   short elementary matrices with coefficients in a small rational bank.
   These operations preserve the non-elementary double coset while changing
   the two row curls.
3. **Build the monomial transport graph.**  Vertices are exponent channels;
   directed edges are nonzero coefficient transfers in the curl equations.
   Record the two row labels and the global sign gauge.
4. **Reject trees and one-ended ladders.**  They reproduce `(11)` and leave a
   terminal residual.  Rank this rejection before coefficient height.
5. **Compute weighted holonomy.**  Reject every cycle whose edge-gain product
   differs from `(-1)^n`.  Support parity alone is not a certificate; in
   particular every cyclic closure of the raw factorial Cohn ladder dies.
6. **Retain reciprocal even cycles first.**  Begin with edge gains such as
   `2,1/2`, then solve the exact rational coefficient schemes while
   saturating away elementary and constant-factor branches.  Matching
   holonomy and multiplier holonomy are both necessary, never sufficient.
7. **Restore the sidecars.**  Impose both curls, exact determinant, the Cohn
   class certificate, source/target gauge normalization, the degree-spectrum
   gates, and finally collision or nonproperness.
8. **Route positive-dimensional coefficient schemes arithmetically.**  Use
   primary decomposition first; genus-zero branches get rational
   parametrization, genus-one branches get exact point/incidence and
   Mordell--Weil searches.  The elliptic machinery belongs here, after the
   parity/cokernel filter, not on the dense `8690`-slot coefficient space.

The cheapest first box is a two- or four-channel reciprocal-gain exponent
cycle.  Its sign pattern is fixed up to global reversal, but its multiplier
product must be pinned exactly before exponent placements are enumerated.  A
useful success is not a tiny residual: it is an exact finite cycle whose two
curls vanish and whose matrix remains outside `E_2`.

## 6. Other bridges and failure controls

### Coherent/dephased networks

The alternating vector `(18)` is a phase sidecar.  Squaring amplitudes or
retaining only conductances erases its global sign and, on more complicated
cycles, its holonomy.  This matches the earlier warning from dephasing:
coefficient magnitudes can pass every polygon inequality while the lost phase
and Segre-cycle equations still obstruct an actual Jacobian pair.  The parity
filter should therefore be applied before, not instead of, exact phase and
rank-one constraints.

### The direct planar rectangle transfer fails

The most tempting current planar cell is an exact hostile, not a positive
example.  In
[`THM-3613-three-by-four-size-seven-ray-parity-gate.md`](../01-canon/theorems/THM-3613-three-by-four-size-seven-ray-parity-gate.md),
the `W003` odd-tail gate leaves seven placements.  Mapping them to
`C_2 x C_7` by scalar-arm orientation and fibre index gives

```text
{+} x {1,2,3,4}  union  {-} x {1,2,3}.                 (19)
```

This set is one cell away from a balanced `2 by 3` rectangle but is not a
rectangle.  Coefficient and regularity coupling has destroyed product
independence, so odd cardinality alone cannot be imported into `JC(2)`.
For any finite planar atlas, the cheapest product test is to compute the two
index sets `I_+,I_-`: the incidence is rectangular exactly when one is empty
or `I_+=I_-`.  Rank near misses by `|I_+ triangle I_-|` and attack the
unmatched bracket first.

Even matching holonomy also needs a symmetry firewall.  If row flip lifts to
affine source/target reflections intertwining the map, the fixed-line
injectivity theorem forces automorphy; a commuting effective even-order
action is likewise in the scope of the cited Miyanishi symmetry obstruction.
The useful even cell is therefore one with nontrivial relative matching
permutation and no polynomial involution realizing its combinatorial row
swap.

### Characteristic two

Odd-cycle determinant `2` is load-bearing.  In characteristic two the odd
operator in `(17)` is also singular and the distinction collapses.  The
present search is over `QQ` or `C`, consistent with the planar conjecture.

### Acyclic truncation hostile

For every `D`, equation `(12)` is an exact near solution with a single failed
monomial.  It is the canonical hostile for any numerical residual norm,
floating-point solver, or bounded coefficient comparison.  The terminal
coefficient tends rapidly to zero while exact solvability never occurs.

### Even kernel hostile

The existence of `(18)` does not solve an inhomogeneous system.  Indeed, a
unit forcing on one vertex of an even cycle violates the alternating
left-kernel compatibility.  Even parity supplies ambiguity or obstruction,
not automatic solvability.  Worse, `(18)` disappears as soon as the weighted
holonomy differs from one; the raw Cohn multipliers do exactly that.  A
construction must both tune the multiplier product and pair terminal residues
so their alternating sum vanishes.

## 7. Reproduction and honest frontier

Run

```bash
python3 04-computation/jc2_cohn_parity_cycle_repair.py
python3 -O 04-computation/jc2_cohn_parity_cycle_repair.py
diff -u \
  05-knowledge/results/jc2_cohn_parity_cycle_repair.out \
  <(python3 04-computation/jc2_cohn_parity_cycle_repair.py)
```

The companion checks the determinant and curls of `(1)`, thirteen exact
factorial truncations, the formal exponential identity, cycle determinants
and ranks for `3<=n<=12`, the nonsingular factorial-weighted Cohn cycles, the
reciprocal `2,1/2` positive kernel, all balanced-rectangle counts in that
range, all orthogonal-grid counts through `n=6`, and all binary
rectangle-coboundary tables for `2<=n<=6`.

The honest frontier is now sharply separated:

- elementary-Jacobian Keller maps are tame: **CITED + INDEPENDENTLY
  REPROVED**;
- Cohn non-elementarity: **CITED**;
- the Cohn one-factor repair no-go and odd/even cycle law:
  **PROVED-SYMBOLIC + FINITE-EXACT**;
- the full odd-`n` Weihrauch recovery implication: **CONDITIONAL** on the
  unreproduced finite-commitment and class-typing steps of the supplied draft;
- existence of a holonomy-tuned exponent-cycle correction with both curls,
  retained non-elementary class, and nonproperness: **OPEN**.

The next calculation should enumerate exact reciprocal-gain two-, four-, and
six-channel repair cycles around `(1)`, not longer elementary words and not a
larger dense coefficient height box.
