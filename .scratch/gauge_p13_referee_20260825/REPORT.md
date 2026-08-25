# Hostile referee: rational-edge gauge, tournament-zeta tangent, and the `p=13` shell

**Verdict: PASS WITH REQUIRED NOTATION AND SCOPE REPAIRS.**  The two proposed
matrix/valuation results are correct, and the `C6`--Berggren table is exact.
The candidate must not be promoted unchanged, however, because it suppresses
a load-bearing distinction between the diagonal and principal Kubota--Leopoldt
branches.  Its continued-fraction cylinder also needs an explicit convention
that permits finite termination.  These are repairable statement defects, not
counterexamples to the underlying elementary results.

No canon file was edited.  An independent standard-library verifier is

```text
.scratch/gauge_p13_referee_20260825/independent_audit.py
```

and runs as

```text
python -B .scratch/gauge_p13_referee_20260825/independent_audit.py
python -B -O .scratch/gauge_p13_referee_20260825/independent_audit.py
```

The transcripts are byte-identical.  Script SHA-256 is
`a5b0861f1056d880dc8b7af341cf3b4493a8370ed67eb6b628c9cf6e73b9ba6e`;
the UTF-8 transcript SHA-256 under the PowerShell `Out-String` replay is
`5b6ceb4b4f241762d344b3463afd55795128063515362471e45f23b7fdef3a0a`.

## 1. Exact verdicts

| Candidate | Verdict | Required boundary |
|---|---|---|
| `W=DAD^{-1}` and graph-zeta invariance | **PASS** | finite matrix; weights are nonzero in a field, or units in a commutative ring |
| `v_p(zeta_T(x)-1)=3v_p(x)` when `p` does not divide `c3` | **PASS; sharp iff boundary** | `x!=0`, `v_p(x)>=1`; if `p|c3`, the valuation is strictly larger (or infinite) |
| `U_13/{+-1}` to the `r=7` Berggren shell | **PASS** | a bijection of labels, not a group-compatible ancestry map |
| Stern/Legendre row is the order-two `C6` character | **PASS** | it is a character/color row, not a tournament orientation |
| trivial orbit sum targets `zeta_13(3)` | **FAIL AS NOTATED** | true for `zeta_13^Delta(3)=L_13(3,omega^{-2})`, false for Beukers' principal `L_13(3,1)` convention |
| finite Stern/continued-fraction cylinders are arithmetic-type blind | **PASS AFTER DEFINITION** | use the real interval/prefix class allowing terminating continuations; a strict infinite symbolic cylinder has no rational members |
| finite `p`-adic balls are arithmetic-type blind | **PASS** | nonempty open balls in `Q_p`; membership alone is the observer |

## 2. Rational-edge gauge

Let `R` be a commutative ring, let `A in M_n(R)`, and let
`d_1,...,d_n in R^x`.  Put

```text
D=diag(d_1,...,d_n),             W_ij=A_ij d_i/d_j.
```

Then, entry by entry,

```text
W=D A D^(-1),
I-uW=D(I-uA)D^(-1)                         in M_n(R[u]).
```

Consequently

```text
det(I-uW)=det(I-uA),
W^k=D A^k D^(-1),
tr(W^k)=tr(A^k).
```

This works for every finite directed multigraph, with loops and digons, and is
not tournament-specific.  In a field, “nonzero” is sufficient; over a general
ring the correct hypothesis is **unit**, not merely nonzero.  Zero weights,
infinite graphs without a separately defined determinant, and noncommutative
coefficient systems are outside this statement.

For the weighted Bowen--Lanford product, each closed-walk weight is one because
the vertex ratios telescope.  Thus the weighted and unweighted graph zeta
functions agree.  The mechanism is literal conjugacy, not the false inference
“an arbitrary exact additive cochain makes an unrelated determinant
singular.”  Hence there is no conflict with MISTAKE-409.

### Candidate-script audit

The candidate script's line 117 checks an expression against the identical
expression used to define it.  The substantive determinant and trace checks
are valid, but use only three Stern tournaments and `D=diag(1,...,n)`.  The
independent verifier instead checks all `512` three-vertex digraphs, including
loops and digons, against three positive/negative/fractional rational weight
rows, four spectral parameters, powers through five, and a multiple-edge
hostile.  It performs `6,148` determinant and `7,680` power/trace gates and
rejects a zero-weight row.

## 3. The `p`-adic tangent

For a finite tournament, THM-1926 gives the first nonconstant determinant
coefficient:

```text
P_T(u)=det(I-uA)=1-c3(T)u^3+u^4 Q(u),       Q in Z[u].
```

Let `x in Q_p` be nonzero with `v_p(x)=m>=1`.  Since `P_T(x)=1 mod p`, it is a
`p`-adic unit, and

```text
zeta_T(x)-1
  =(1-P_T(x))/P_T(x)
  =x^3(c3(T)-xQ(x))/P_T(x).
```

If `p` does not divide `c3(T)`, the parenthesized factor is a unit, proving

```text
v_p(zeta_T(x)-1)=3m.
```

Conversely, if `p|c3(T)`, both terms in `c3-xQ(x)` are divisible by `p`, so
the valuation is strictly greater than `3m`, or is infinite when the
difference vanishes.  Thus for `m>=1` the displayed hypothesis is the exact
iff boundary, not just a convenient sufficient assumption.  If merely
`x in p^m Z_p`, the unconditional conclusion is the weaker inclusion
`zeta_T(x) in 1+p^(3m)Z_p`.

A minimal necessity hostile is the order-four tournament with

```text
A=[[0,0,0,1],
   [1,0,0,0],
   [1,1,0,0],
   [0,1,1,0]],
P_T(u)=1-2u^3-u^4,             c3=2.
```

At `p=2`, `x=2`, one has `P_T(2)=-31` and

```text
zeta_T(2)-1=-32/31,            v_2=5>3.
```

The phrase “the first p-adic digit is the triangle count” must therefore say
“when `p` does not divide `c3`,” or be read only as a formal coefficient
statement.  The independent verifier checks the coefficient `-c3` by a
Leibniz determinant for all `1,099` labelled tournaments through order five
and checks `102,700` valuation gates for every labelled tournament through
order six, including nontrivial unit factors and `m=2` controls.

## 4. The `U_13/{+-1}`--Berggren shell

The finite dictionary is correct.  The group `U_13` is cyclic of order twelve,
and quotienting by `{+-1}` gives `C6`.  Each antipodal pair has exactly one odd
representative

```text
d in {1,3,5,7,9,11}.
```

Under THM-3756, set `r=7` and `s=(d+1)/2`.  Then

```text
(A,B,C)=(13d,(169-d^2)/2,(169+d^2)/2),
B+C=169,                    C-B=d^2.
```

Because `13` is prime and `1<=d<13`, these are precisely the six primitive
ordered triples in the complete outer shell `r=7`.  Independent forward
generation by the affine `L,M,R` child maps gives

| `d` | triple | word | `(d/13)` |
|---:|---|---|---:|
| 1 | `(13,84,85)` | `LLLLL` | `+1` |
| 3 | `(39,80,89)` | `ML` | `+1` |
| 5 | `(65,72,97)` | `RM` | `-1` |
| 7 | `(91,60,109)` | `LLR` | `-1` |
| 9 | `(117,44,125)` | `LRR` | `+1` |
| 11 | `(143,24,145)` | `RRRRR` | `-1` |

The class of `2` generates the quotient and has odd-representative cycle

```text
1,11,9,5,3,7.
```

Since `(-1/13)=+1`, the Legendre symbol descends to the quotient.  Since
`(2/13)=-1`, its values on that generator cycle alternate.  It is exactly the
unique nontrivial order-two character of `C6`, the `k=3` row under the usual
six-point DFT indexing.

The map to triples is a label bijection, not a homomorphism into Berggren
ancestry; the candidate correctly records this loss.  The nonidentity class
`{+-5}` is self-inverse because `5^2=-1 mod 13`, so a translation-invariant
Cayley tournament on `C6` cannot exist.  The Legendre row should be called a
quadratic **character/color** row.  No graph on the six quotient vertices is
defined merely by listing its values, so “Paley graph” is unnecessarily
strong terminology.

### Required `p`-adic-zeta notation repair

This is the only material mathematical presentation defect.  Beukers uses

```text
u_a=omega(a)^(-2) H_p(3,a,p),
T_p(a/p)=2p^3 u_a,
T_p(x)=T_p(1-x).
```

Thus `u_(p-a)=u_a`, so the six reflection **labels** are valid.  These formulas
and Corollary 11.3 are in Beukers,
[Irrationality of some p-adic L-values](https://arxiv.org/html/math/0603277v2).

But the sum of the normalized carriers is

```text
sum_(a in U_13) u_a
  =L_13(3,omega^(-2))
  =zeta_13^Delta(3),
```

under the historical report's explicitly chosen diagonal convention.  It is
twice the **trivial-character** orbit sum.  In Beukers' Appendix notation, the
principal branch is instead

```text
zeta_13^principal(3)=L_13(3,1)
  =sum_a H_p(3,a,13)
  =sum_a omega(a)^2 u_a.
```

After antipodal descent, `omega^2` is an order-six character of `C6`, not the
trivial character.  The candidate report never states the diagonal convention
and repeatedly writes bare `zeta_13(3)`.  Promotion must restore the superscript
or spell out `L_13(3,omega^(-2))`; otherwise the asserted target DFT row changes.

## 5. Arithmetic-type blindness of finite cylinders

The intended real statement is true after defining the object.  Let `w` be an
admissible finite simple-continued-fraction word and let `C(w)` be the set of
real numbers whose canonical finite or infinite expansion begins with `w`,
allowing a finite continuation to terminate.  Then:

- appending a finite admissible tail gives a rational point of `C(w)`;
- appending an eventually periodic infinite tail gives a quadratic irrational
  by Lagrange's theorem;
- the set of infinite tails has cardinality continuum, while the algebraic
  numbers are countable, so `C(w)` contains transcendental points.

The same holds for the corresponding nonempty Stern interval.  If instead
“cylinder” means a subset of the strict infinite symbolic shift
`N_(>0)^N`, rational numbers are not elements at all; they are finite words or
boundary/terminating completions.  The theorem must choose one convention.

For `p`-adics, every nonempty open ball in `Q_p` contains a rational point
because `Q` is dense.  Every such ball has cardinality continuum, whereas the
elements algebraic over `Q` form a countable set.  Hence every ball contains
uncountably many elements transcendental over `Q`.  A bounded digit/residue
prefix therefore cannot decide rationality or transcendence.  This is a
statement about the cylinder-membership observer; exact labels, recurrence,
height, and global compatibility may retain additional information.

Neither cardinality argument proves anything about the arithmetic type of a
named `p`-adic zeta value.  It proves only that finite-prefix membership cannot
do so by itself.

## 6. Canon and correction collision audit

1. **Gauge:** THM-4057,
   `01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md`,
   equation `(CW5)`, already proves endpoint-ratio telescoping.  THM-1926,
   `01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md`,
   owns the determinant/zeta dictionary.  Their conjunction makes the new
   gauge theorem an immediate but useful general corollary.  No exact proved
   duplicate of the full diagonal-similarity statement was found.
2. **Tangent:** THM-1926 owns `P_T(u)=1-c3u^3+...`; the valuation equality is a
   new one-line local corollary.  THM-4093,
   `01-canon/theorems/THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent.md`,
   is only a `RESERVED / UNPROVED EMPTY STUB` and supplies no dependency.
3. **Shell:** THM-3756,
   `01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md`,
   Theorems 3.1/4.1 already give the odd-residue outer shell; THM-4057 `(CW29)`--`(CW30)`
   gives the same root/triple chart; and THM-4059,
   `01-canon/theorems/THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution.md`,
   Section 6 already proves that the Stern sign on `U_13` is Legendre.  Section
   5 of the candidate is therefore an exact explicit specialization/composition,
   not a new independent mechanism.
4. **Arithmetic type:** THM-4088,
   `01-canon/theorems/THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md`,
   is also an empty reserved stub.  MISTAKE-219 already states the broad
   firewall that a finite prefix supplies no arithmetic type.  THM-3848 and
   THM-4072 provide stronger nearby finite-prefix versus infinite-predicate
   hostiles, but not this exact rational/quadratic/transcendental cylinder
   statement.
5. **Corrections:** MISTAKE-409 is respected because conjugacy is exhibited;
   MISTAKE-418 is respected because the shell uses the primitive odd-root
   chart; MISTAKE-209/219 forbid the named-irrationality promotion that the
   candidate correctly declines.

## 7. Promotion wording

The safe theorem statement is:

> Let `A` be the adjacency matrix of a finite digraph over a commutative field
> and let every vertex weight `d_i` be nonzero.  With
> `W_ij=A_ij d_i/d_j`, one has `W=DAD^{-1}` and therefore
> `det(I-uW)=det(I-uA)` and identical graph zeta.  For an unweighted finite
> tournament with triangle count `c3`, every nonzero `x in pZ_p` satisfies
> `v_p(zeta_T(x)-1)>=3v_p(x)`, with equality if and only if `p` does not divide
> `c3`.  Separately, the six labels `U_13/{+-1}` identify with the complete
> primitive Berggren shell `r=7`, and the Stern sign is its order-two `C6`
> character.  The trivial orbit sum targets the diagonal branch
> `L_13(3,omega^{-2})`, not the unlabeled principal branch.  Finite real or
> `p`-adic cylinders do not determine arithmetic type, under the explicit
> cylinder conventions above.

With those repairs, the mathematical core passes hostile review.
