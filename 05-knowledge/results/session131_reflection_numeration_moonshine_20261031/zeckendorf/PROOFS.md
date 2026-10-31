# Session 131: reflection gluing, nonconsecutive digits, and the carry at zero

**Status:** PROVED elementary maps and identities; FINITE-EXACT 17-state
minimality certificates and pinned order-five specialization. No novelty claim
for classical numeration. These are integral **metagraph cycles**, not cycles
inside a tournament or Hamiltonian cycles in a square-sum graph.

## 0. Inheritance and connection contract

Inherited at source commit `41a9d96041ca6f970158450972c5295695387e53`:

* `session130_blue_black_theorems_20261031/gain_cohomology/THEORY.md`, Sections
  1–3: vertexwise-fixing reflection; parallel black twins; integral cycle
  lattice; eigenspace splitting index `2^10`.
* `session125_origin32_numeration_20261031/numeration/PROOFS.md`, Sections 2–5:
  Chung–Graham (CG) normal form and terminating positive carries.
* Its `BEATTY_BRIDGE.md`, Sections 1–3: the canonical shift is
  `P(n)=ceil(beta*n)`, `beta=(3+sqrt(5))/2`.
* `session125_rooted_forest_numeration_20261031/PROOFS.md`, Sections 2–4:
  explicit CG–tiling–rooted-path-forest bijection and the actual operation
  `w -> w0 <-> T -> TMM <-> append one rooted singleton`.

Closest proved mechanism: the index-two sum/difference lattice in each twin
pair. Hostile inherited example: one black-loop twin with coefficient one,
whose reflection projections are half-integral. Corrected near miss:
nonconsecutive numerals are not independent reflection eigencoordinates.
Least-used sidecar: the lower-boundary carry of CG normalization.

**Source -> target:** an integral reflection cycle -> canonical signed digit
words for its integer coordinates, plus a finite binary gluing code. The map
preserves exact coefficients, cycle addition, and reflection; it forgets
nothing when signs, positions, and the graph boundary are retained. Independent
numeric shifts destroy gluing. The repair is a precisely bounded carry, not a
new square-sum adjacency law.

## 1. Exact integral gluing coordinates and their carry law

Take any finite multigraph whose involution fixes vertices and fixed edges
pointwise and swaps the other edges in equally oriented pairs. Write `B` for
the fixed-edge graph, `K` for one edge of each twin orbit, and `Q=B union K`.
Let an integral chain have fixed-edge coefficients `b` and twin coefficients
`x_j,y_j`. Define

```
s_j = x_j+y_j,             d_j = x_j−y_j.
```

The chain is an integral cycle exactly when

```
∂B b + ∂K s = 0,          s_j ≡ d_j (mod 2) for all j,
```

with integral `b,s,d`; its inverse is `x=(s+d)/2, y=(s−d)/2`.
The boundary equation is separate and must not be dropped.

Choose the unique `epsilon_j in {0,1}` congruent to `s_j,d_j`, and write

```
s=2u+epsilon,             d=2v+epsilon.
```

Then the exact inverse and cycle condition become

```
x=u+v+epsilon,            y=u−v,
∂B b + 2∂K u + ∂K epsilon = 0.                  (1)
```

Thus `(b,u,v,epsilon)` with (1) is a bijective integer normal form. Negative
coordinates cause no problem: the remainder convention is always 0 or 1.

For two such normal forms the transported group operation is

```
b'' = b+b',
u'' = u+u'+epsilon*eta,
v'' = v+v'+epsilon*eta,
epsilon'' = epsilon XOR eta,                   (2)
```

where products are coordinatewise. This follows from
`epsilon+eta=(epsilon XOR eta)+2epsilon*eta`. It proves closure of (1), since
(2) reconstructs addition of the original cycles. Reflection acts by

```
R(b,u,v,epsilon) = (b,u,−v−epsilon,epsilon).     (3)
```

The simultaneous `epsilon*eta` terms are genuine integral carries. On a
black-loop pair, adding the state of one loop to itself gives a nonzero carry:
`(u,v,epsilon)=(0,0,1)` added to itself gives `(1,1,0)`. Deleting these carries
would give zero instead of twice that loop.

The sum of the integral plus and minus eigensublattices is exactly the
subgroup `epsilon=0`: invariant chains have `v=epsilon=0`, while
anti-invariant chains have `b=u=epsilon=0`. Conversely a cycle with
`epsilon=0` splits by (1). Hence epsilon is the actual gluing obstruction,
not an analogy with digit parity.

### The allowable epsilon code

Contract each connected component of `B` in `Q`, retaining black loops and
multiplicities, to obtain `H=Q/B`. The realizable epsilon vectors are exactly
`Z_1(H;F_2)`.

Necessity: sum (1) over every blue component modulo two. Sufficiency: a binary
cycle in `H` lifts to a binary cycle of `Q` by routing its even demands inside
each connected blue component. Every binary graph cycle lifts to an integral
cycle (use a spanning-forest fundamental basis). Its black coordinates give
the required epsilon, and then (1) defines `u`.

Thus the image has dimension `|E(K)|+c(Q)−c(B)`. In the pinned order-five atlas,
black edges are ordered

```
(1,4), (2,4), (2,5), (2,6), (2,7), (2,8),
(4,5), (4,8), (5,7), (6,6), (6,8), (8,8).
```

The blue components are `{1}`, `{2}`, and `{0,3,4,5,6,7,8,9}`. Consequently

```
epsilon_0=0,             epsilon_1+...+epsilon_5=0 mod 2,             (4)
```

with the last six bits free. This is a 10-dimensional binary code. Every one
of its 1024 words is attained, not merely counted. The checker enumerates the
actual image of integral quotient fundamental cycles and matches (4).

## 2. Two numeration languages, one exact parity gate

Use `F_0=0,F_1=1`. CG uses weights `A_i=F_(2i+2)=1,3,8,...` and digits 0,1,2,
with no factor `2 1* 2`. Zeckendorf uses `W_i=F_(i+2)=1,2,3,...`, binary digits,
and no consecutive ones. All strings in this section are MSB-first; leading
zero padding is allowed for synchronization. A signed integer is represented
by a separate sign and the canonical magnitude word; zero has sign zero.
Signs do not affect parity.

The parities are

```
A_i mod 2 = 1,1,0,1,1,0,...,
W_i mod 2 = 1,0,1,1,0,1,... .                  (5)
```

For either system, a pair of individually legal words for `s,d` lifts to an
integral twin pair **iff the represented magnitudes have equal parity**.
Neither nonconsecutiveness nor last-digit agreement replaces this gate.
For example, `s=1,d=2` have legal Zeckendorf words `1,10` but cannot lift;
`s=3,d=1` do lift although their CG last digits are `0,1`.

Here is a complete finite recognizer for the gate. State is `(l,r,n,p)`:
`l,r` are the two individual legality states, and `n,p` are modulo-two values
of the difference of the prefixes and their one-place shifts. Start at
`(0,0,0,0)`. Read a digit pair `(a,b)` and put `t=(a+b) mod 2`.

For CG, `l=1` means Locked. A 2 from a locked word goes to a rejecting sink;
otherwise update each legality state by `0 -> Free`, `2 -> Locked`, with
`1` preserving the old state, and update

```
(n,p) -> (p+t, n+p+t) mod 2.                   (6)
```

For Zeckendorf, `l,r` are the preceding digits. Consecutive ones go to the
sink; otherwise new legality states are `a,b`, and

```
(n,p) -> (p+t, n+p) mod 2.                     (7)
```

Accept exactly nonsink states with `n=0`. The recurrence for CG is
`A_(i+2)=3A_(i+1)−A_i`; that for Zeckendorf is `W_(i+2)=W_(i+1)+W_i`.
They give (6)–(7) directly, including their different digit injection terms.
Thus these machines recognize exactly the canonical equal-parity relation
for words of arbitrary length, not just a finite tested range.

Each machine has 16 live states and one sink, and **17 is minimal** for this
padded pair language. The portable checker supplies a reachability word for
every state and a distinguishing suffix for every one of the 136 unordered
state pairs. Checking those finite witnesses proves the lower bound by
Myhill–Nerode; the displayed machine proves the upper bound. The certificate
is regenerated without external packages in normal and optimized Python.
This is an exact finite-state theorem, not a guess from state counts.

Apply the gate to every twin pair and retain (1)/(4), or use the parent
packet's boundary-preserving integral exchange-block basis. One must not
replace the rank-23 lattice by one scalar numeral.

## 3. What 'even and odd interwoven, nonconsecutive' actually means

For a Zeckendorf word define its even- and odd-indexed Fibonacci strands

```
e_i = z_(2i+2),           o_i = z_(2i+3),       i>=0.
```

They are binary. The exact nonconsecutive conditions are

```
e_i o_i=0,               o_i e_(i+1)=0.        (8)
```

Each strand by itself can contain arbitrary adjacent ones: its weights were
already separated by one Fibonacci index. The forbidden adjacency lives
**between** the two strands. Naively interleaving two independently legal
binary strands loses precisely (8); `e_0=o_0=1` gives adjacent Fibonacci
weights 1 and 2, whose sum normalizes to 3.

There is a reversible, value-preserving map to signed even-Fibonacci digits:

```
c_i = e_i + o_(i−1) − o_i,       o_(-1)=0.      (9)
```

Indeed `F_(2i+3)=A_(i+1)−A_i`, so the sum telescopes to the same integer.
Conditions (8) ensure `c_i in {−1,0,1}`. Unlike a scalar rank identification,
this map has the following explicit two-state inverse, reading LSD-first:

| previous `o` | `c` | output `(e,o)` | next state |
|---|---:|---|---:|
| 0 | 0 | (0,0) | 0 |
| 0 | 1 | (1,0) | 0 |
| 0 | −1 | (0,1) | 1 |
| 1 | 0 | (0,1) | 1 |
| 1 | 1 | (0,0) | 0 |
| 1 | −1 | forbidden | — |

Start and finish in state 0, adding the final carry digit if needed. The
transition equations are exactly (8)–(9); they prove mutual inverseness.
Equivalently every −1 opens a block `−1 0* 1`, which must close before another
−1. For example `2=F_3` becomes `−A_0+A_1`; its LSD digits are `(-1,1)`.

This gives a genuine nonconsecutive/interwoven representation theorem.
The target is a **signed even-weight normal form**, not the positive CG
normal form: the positive-only rewrites of Section 4 do not apply to its
negative digits. Converting onward to positive CG can always evaluate the
integer and use the inherited greedy encoder. Likewise these local signed
digits are distinct from the separate sign-magnitude encoding used for
arbitrary negative cycle coefficients.

But Fibonacci index parity is not the same as reflection eigenvalue, nor
integer value parity. To use (9) for a reflection cycle, encode each of its
integer coordinates and retain the parity gate of Section 2. These distinct
constraints are not silently identified.

## 4. A carry theorem detecting exactly why independent shifts fail

For a nonnegative digit array `a`, let its raw one-place-shift value be
`T(a)=sum_i a_i A_(i+1)`. For canonical digits of `n`, inherited arithmetic
says `T(a)=P(n)=ceil(beta*n)`.

Positive CG normalization uses these inherited value-preserving rewrites:

```
3 A_0 -> A_1,
3 A_i -> A_(i+1)+A_(i−1)                         (i>=1),
2 A_i + sum_(l<j<i) A_j + 2 A_l
      -> A_(i+1)+A_(l−1)                       (i>l>=0),
```

with `A_(-1)=0`. The final rule removes a forbidden `2 1* 2` block. Every
rewrite strictly decreases the nonnegative coefficient sum, so arbitrary
choices terminate in the unique canonical word.

**Boundary-carry theorem.** Under any such normalization, every rewrite with
lowest index at least one preserves `T`; a rewrite whose lowest index is
zero decreases `T` by exactly one. Therefore

```
number of zero-boundary rewrites = T(input) − P(value(input)).       (10)
```

Proof: shift each rewrite by one index. The higher identities remain the
same recurrence/block identity. For `3A_0 -> A_1`, the shifted difference is
`3A_1−A_2=9−8=1`. For a block with `l=0`, its shifted identity has an extra
`A_0=1` on the right, which the shift of the unshifted output omits. Hence
again the decrease is exactly one. Sum over the terminating rewrite path.

For a sum of `r>=1` canonical nonnegative numerals, (10) becomes

```
C(n_1,...,n_r)=sum_j P(n_j)−P(sum_j n_j),
0 <= C <= r−1.                                (11)
```

The bound is the ceiling-sum inequality (and remains valid when some inputs
are zero). In particular **adding two canonical CG words has at most one
zero-boundary rewrite, regardless of rewrite strategy**. Its presence is
exactly `kappa(a,b)=P(a)+P(b)−P(a+b) in {0,1}`.

For positive `a,b`, let `rho(n)=P(n)−beta*n in (0,1)`. Then
`kappa(a,b)=1` iff `rho(a)+rho(b)>1`; equality cannot occur because
`beta*(a+b)` is irrational. This gives the exact phase criterion, not merely
a bound. No transcendental approximation is needed in the checker:
`P(n)=(3n+isqrt(5n*n))//2+1` for `n>0`, and `P(0)=0`.

### Reflection-compatible growth and its sharp correction

Use any integral equivariant normal form with exchange-block coefficients
`x,y` (constructed in the parent lane by a blue-first spanning forest).
Define the odd extension `K(n)=sign(n)P(abs(n))`. Coefficientwise `K` on that
basis is a well-defined boundary-preserving, reflection-equivariant map.
It is not additive and is basis-dependent; no canonical graph growth is
asserted.

In eigen coordinates `s=x+y,d=x−y`, its exact formula is

```
s' = K(x)+K(y)  = K(s)+kappa_Z(x,y),
d' = K(x)−K(y)  = K(d)+kappa_Z(x,−y),          (12)
```

where `kappa_Z(a,b)=K(a)+K(b)−K(a+b)` always lies in `{−1,0,1}`. For equal
signs this follows from (11) and oddness. For opposite signs, subtract the
smaller magnitude first and reduce to the equal-sign case. Thus (12) is
an exact bounded-carry repair of naive independent eigenstrand shifts.
Both new sums have equal parity because their sum and difference are twice
integers. Reflection changes `d'` to `−d'` and preserves `s'`.

A minimal positive-twin failure (by `x+y`) is `x=2,y=1`, hence `(s,d)=(3,1)`.
Independent CG shifts give `(8,3)`, which has no integral lift. The repaired
shift is `(9,3)`, from twins `(6,3)`. Its correction `+1` is the unique bottom
carry `3A_0 -> A_1` when adding the canonical words of 2 and 1. The same
failure is a valid cycle witness on one of the actual black-loop pairs.

### The smallest universal shift block is exactly three digits

There is also an exact repair by grouping shifts. Let `K^k` mean iteration,
not a numerical power. In CG parity registers, appending zero applies

```
M = [[0,1],[1,1]] over F_2,       M^3=I.
```

For a liftable eigenpair the difference register starts `(n,p)=(0,p)`.
If `p=0`, all iterated simultaneous shifts remain liftable. If `p=1`, the
states cycle `(0,1) -> (1,1) -> (1,0) -> (0,1)`, so the liftable shifts are
**exactly** `k=0 mod 3`. Thus three is the smallest positive universal block
length that preserves integral gluing. The sharp witness `(s,d)=(3,1)` gives

```
(3,1) -> (8,3) -> (21,8) -> (55,21).
  lifts     fails     fails       lifts
```

The identities `A_(i+2)=3A_(i+1)−A_i` and
`A_(i+3)=8A_(i+1)−3A_i` give, also for the odd extension,

```
K^2(n)=3K(n)−n,          K^3(n)=8K(n)−3n ≡ n (mod 2).                (13)
```

These hold because zero-appending preserves canonical legality at every
stage. Applying independent `K^3` to the two eigen coordinates of every
exchange block is consequently an integral, reflection-equivariant map
which preserves its gluing residue. It still differs from twinwise `K^3`:
the two outputs differ by `8*kappa_Z(x,y)` and `8*kappa_Z(x,−y)` in the
respective eigen coordinates. Restoration of integrality does not make
the operations identical or additive. The fixed/trivial and sign blocks
can be treated by the same odd scalar map in the parent's integral basis.

**Zeckendorf corollary.** The same exact return-period classification holds
for appending zeros to genuine nonconsecutive Zeckendorf words, since (7)
has the same zero-input matrix `M`. If `J(n)` is this shift with odd extension
to negative integers, the Fibonacci recurrence gives

```
J^2(n)=J(n)+n,           J^3(n)=2J(n)+n ≡ n (mod 2).                (14)
```

The sharp witness is `(3,1) -> (5,2) -> (8,3) -> (13,5)`. Thus the common
three-step compatibility is a concrete transported modular operation, not
just the observation that both systems contain Fibonacci numbers. Formula
(14) does not identify the Zeckendorf shift `J` with the CG Beatty shift `K`.

**Scope boundary.** Applying `K` independently to arbitrary edge
coefficients generally destroys the cycle boundary; only applying it in a
cycle basis preserves that boundary. On CG words for rooted path forests,
`P` really does append a singleton. Transporting that operation to chosen
cycle-basis coefficients gives a chosen arithmetic dynamics, not an
isomorphism to square-sum vertex arrival or Hamiltonian gluing. The preceding
positive result and obstruction explain both the relation and its limit.

## 5. Reproduction, hostile controls, and exact scope

Run, using Python's standard library only:

```
python3 check_bridge.py > normal.json
python3 -O check_bridge.py > optimized.json
cmp normal.json optimized.json
python3 check_bridge.py --test-failing-gate       # must exit nonzero
python3 -O check_bridge.py --test-failing-gate    # must exit nonzero
```

The checker regenerates both full DFA minimality certificates, compares
integer evaluation with the modular registers, checks reversible interwoven
digits, tests both high-first and low-first normalization, compares bottom
carry counts with exact Beatty ceilings, verifies lattice reconstruction,
reflection and addition, and enumerates all 1024 realized parity words.
Explicit exceptions retain all gates under `python -O`. Universal claims
rest on the proofs above; finite tests verify implementations and specialize
the pinned finite atlas. Receipts state the precise tested universes.
