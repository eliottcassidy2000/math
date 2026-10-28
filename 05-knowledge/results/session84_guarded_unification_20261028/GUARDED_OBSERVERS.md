# Guarded observers: the exact quotient and the arithmetic rank are different minima

**Session 84, 2026-10-28. Status:** elementary statements below **PROVED**;
verification windows **FINITE-EXACT**. The quotient and linear-rank theorems
are standard in substance, not claimed as new. The explicit arithmetic
applications, exact minors, and separation of the two notions of state are
this note's contribution to the session; no external priority claim is made.

## Results first

1. The coarsest exact deterministic observer of a guarded operation system
   identifies precisely the objects with identical output **and definedness**
   in every continuation/context. This is the standard contextual-congruence /
   Myhill–Nerode construction with an explicit failure symbol.
2. On rectangles where every multiplication is coprime, the aliquot kernel
   `K_s(a,b)=s(ab)` has rational separation rank at most **2**, attained already
   at rows `(1,2)` and columns `(1,3)`, with determinant `−1`.
3. For nonsquarefree products, the FUS defect kernel `K_δ(a,b)=δ(ab)` has rank
   at most **4**, attained on two disjoint three-prime blocks. A concrete
   `4×4` minor has determinant **4**.
4. These ranks do **not** measure bits or complete guarded state. Under named
   coprime multiplication, the exact guard-aware quotient for scalar `δ`
   on all `n≥2` is `(rad(n),τ(n))`. For `s`, contextual equivalence is equality of `n`:
   one common fresh-prime continuation, together with the present output,
   reconstructs the integer. Adding named prime-exponent promotions to the
   defect language likewise refines its quotient all the way to equality.
5. The Berggren residue observer is already a continuation congruence; its
   failure is that the branch sign does not descend through that quotient.
   The sign itself is a two-state word observer. This is not a proof of
   unbounded branch-history memory.

**Notation.** `σ(n)=sum_{d|n} d`, `s(n)=σ(n)−n`, `τ(n)=#Div(n)`,
`ω(n)=# {p:p|n}`, `rad(n)=product_{p|n}p`. FUS uses **counts**, never `σ`.
All prime letters in a shape are distinct unless explicitly stated otherwise.

## 1. The standard guarded-context theorem

Let `X` be a set, `O:X→Y` an observation, and `T_a:D_a→X` partial
operations indexed by named labels `a∈A`. A word `w∈A*` acts successively
when every step is defined. Set

```
R_x(w) = O(T_w x)  if the whole continuation is defined,
         ⊥        otherwise,
```

where `⊥` is a new symbol not in `Y`, and the empty continuation is defined
with `R_x(ε)=O(x)`. Define

```
x ≡_O y  iff  R_x(w)=R_y(w) for every w∈A*.
```

### Theorem 1 (guarded observational quotient; standard)

`≡_O` is the coarsest equivalence relation such that

* equivalent inputs have equal present output;
* `x≡_O y` implies `x∈D_a ⇔ y∈D_a` for every named operation;
* if that operation is defined, `T_a x≡_O T_a y`.

Consequently `Q=X/≡_O` carries well-defined partial updates and an output
map. If another sufficient deterministic carrier `h:X→Z` has output,
guard decisions, and updates factoring through `h`, then

```
h(x)=h(y) ⇒ x≡_O y.
```

There is a unique map `h(X)→Q` sending `h(x)` to `[x]`, preserving the
reachable updates, guards, and output. Thus `Q` is the canonical coarsest
sufficient quotient, not just one possible list of statistics. A finite
exact carrier exists iff `≡_O` has finitely many classes; in that case the
number of classes is its minimum possible number of states.

**Proof.** Equality of residual functions is an equivalence. The empty word
gives equal outputs. The word `a` tests whether `T_a` is defined, because
an actual output cannot equal `⊥`. If defined, every continuation `w` from
`T_a x` is the continuation `aw` from `x`, so residual equality gives the
update congruence. Conversely, any relation satisfying the three bullets
preserves definedness and output along every word, by induction on length;
it therefore refines `≡_O`. Apply this to `ker h`. The factor map is
well-defined, and every quotient class is hit. ∎

For an arbitrary partial algebra with binary or higher-arity compositions,
replace words by **one-hole contexts with fixed parameters**, including the
identity context. The same statement and proof hold. To show compatibility
with a `k`-ary operation, replace its arguments one at a time using those
contexts; each replacement preserves definedness and the resulting class.
The family of contexts must be closed under the substitutions being used.

### What this inherits, rather than discovers

The right-congruence construction is Myhill–Nerode in its standard setting;
see [Cornell's theorem and coarsest-congruence proof](https://www.cs.cornell.edu/courses/cs682/2008sp/Handouts/MN.pdf).
The explicit failure symbol merely incorporates partial-operation guards.
The repository already has the single-step kernel criterion and linear
sidecar theorem in **THM-840**, §§1–1a. Theorem 1 organizes their all-futures,
guard-aware form; it is not a claim to have invented contextual equivalence.

Three questions must not be conflated:

* **Static descent:** does a target `P` factor through a proposed observation
  `O`? Exactly when `ker O⊆ker P`.
* **Continuation closure:** is `ker O` preserved by the named updates and
  their guards? Exactly when `O` itself carries those updates.
* **Minimal complete quotient:** which inputs agree in every future
  experiment? This is `≡_O`, possibly much finer than current-output equality.

A target may fail static descent through a quotient that is perfectly
continuation-closed. Conversely, a current output may be easy to compute
and still fail continuation closure.

## 2. Linear arithmetic observers: a different minimization problem

Let `A×B` be an **admissible rectangle**, meaning every pair is in the
operation's domain, and let `K:A×B→Q` be its scalar observation. Its
separation rank is the least `d` allowing

```
K(a,b) = Σ_{i=1}^d u_i(a) v_i(b).
```

### Theorem 2 (rank factorization; standard)

The separation rank is taken in `{0,1,2,…,∞}`: `∞` means that no finite
factorization exists. It equals the row-span dimension when that dimension
is finite, and is `∞` otherwise; equivalently, it is the supremum of finite
submatrix ranks in this extended set.
A nonzero `d×d` minor rules out every factorization with fewer than `d`
scalar channels; an explicit `d`-term factorization proves the matching
upper bound.

**Proof.** A `d`-term expression puts all rows in the span of its `v_i`.
Conversely a basis of the row span expresses every row with coefficient
functions `u_i`. For infinite index sets, any finite independent family of
rows is detected by finitely many coordinate evaluations, yielding a
finite minor of the same rank. ∎

This is ordinary matrix factorization. On full word monoids with
`K(u,v)=f(uv)`, its transition-compatible version is the classical
Carlyle–Paz/Fliess Hankel-rank minimal weighted-automaton theorem; see
[Kiefer's self-contained notes](https://arxiv.org/abs/2009.01217), which
explicitly claim no new results. We do **not** silently apply that theorem
to a partial matrix of legal and illegal arithmetic contexts: here rank is
computed on stated fully admissible rectangles.

These are linear scalar channels with rational values of unbounded size,
not bits. A countably infinite set can be encoded nonlinearly in one
integer, and polynomially many arithmetic registers need not determine a
guard. Any claim about a minimum must name this representation class.

## 3. Aliquot multiplication: rank exactly two, full contextual separation

For coprime `a,b≥1`, divisor-sum multiplicativity gives

```
s(ab) = σ(a)σ(b) − ab.
```

Thus the kernel has separation rank at most 2, with row feature
`(σ(a),a)` and bilinear coefficient matrix `diag(1,−1)`.
The coprime rectangle with rows `(1,2)` and columns `(1,3)` is

```
       1  3
   1 [ 0  1 ]
   2 [ 1  6 ],      determinant −1.
```

Therefore the universal bound 2 is sharp. Some individual rectangles
(singletons, for example) have smaller rank: the statement is **not**
that every coprime rectangle has rank exactly 2. The feature updates under
a known legal multiplier are linear:

```
(σ(n),n)ᵀ ↦ diag(σ(b),b) (σ(n),n)ᵀ.
```

The one-coordinate current observation does not update deterministically:
`s(2)=s(3)=1`, but multiplication by the common legal factor 5 gives
`s(10)=8` and `s(15)=9`.

### Proposition 3 (exact guarded aliquot quotient)

On positive integers with `T_b(n)=nb` defined iff `gcd(n,b)=1`, two inputs
are observationally equivalent for `O=s` iff they are equal.

**Proof.** Suppose `n≡_s m`. The empty continuation gives `s(n)=s(m)`.
Choose any prime `r` dividing neither input; such a prime exists. Then

```
s(nr) = (r+1)σ(n)−rn = (r+1)s(n)+n,
s(mr) = (r+1)s(m)+m.
```

Equality of continuation outputs implies `n=m`. The converse is immediate.
In fact `n=s(nr)−(r+1)s(n)` reconstructs the input. ∎

Rank 2 therefore coexists with **no nontrivial set-theoretic observational
quotient**. There is no contradiction: linear factorization dimension is
not the cardinality or bit size of the exact deterministic state.

## 4. Correct FUS counts and the four-dimensional defect observer

For `N≥2`, let `D(N)={d:d|N,1<d<N}`. Let `F` count all its elements,
`U` its prime elements, and `S` its squarefree elements. Then

```
F(N) = τ(N)−2,
U(N) = ω(N)−[N prime],
S(N) = 2^ω(N)−1−[N squarefree].
```

Primes belong to the squarefree set: `U⊂S`, so `F=U+S` is not a partition
identity. These definitions, the correction terms, and the balanced shapes
are inherited from the recovered divisor-shape source (§7 below).

Define `δ=F−U−S`. For **nonsquarefree** `N`,

```
δ(N)=τ(N)−2^ω(N)−ω(N)−1.
```

Do not apply this reduced expression to primes or squarefree composites.
The full expression has the additional terms `[N prime]+[N squarefree]`.
The balanced shapes are exactly `p`, `p³`, and `p²qr`.

For completeness, the short inherited classification mechanism is:
a squarefree composite has `δ=−ω`; for nonsquarefree support size `r`,
balance requires `τ=2^r+r+1`. At `r=1` this gives exponent 3; at `r=2`
it would require a product of two integers at least 2 to equal 7, impossible;
at `r=3` the minimum nonsquarefree divisor count is 12 and forces shape
`(2,1,1)`; at `r≥4`, even the minimum `3·2^(r−1)` exceeds `2^r+r+1`.
Primes have `(F,U,S)=(0,0,0)` and are balanced separately.

For coprime multiplication, write `q(n)=2^ω(n)` and

```
v(n)=(τ(n),q(n),ω(n),1)ᵀ.
```

The first two coordinates are multiplicative characters on legal products;
`ω` is additive. A legal multiplier `b` acts by

```
       [ τ(b)   0     0    0   ]
T_b =  [  0    q(b)   0    0   ],
       [  0     0     1   ω(b) ]
       [  0     0     0    1   ]

v(nb)=T_b v(n),       δ(n)=(1,−1,−1,−1)v(n)  for nonsquarefree n.
```

### Proposition 4 (exact defect separation rank)

On any coprime rectangle all of whose products are nonsquarefree,

```
K_δ(a,b)=δ(ab)
       =τ(a)τ(b)−q(a)q(b)−ω(a)−ω(b)−1
       =v(a)ᵀ C v(b),

C = [1  0  0  0]
    [0 −1  0  0]
    [0  0  0 −1]
    [0  0 −1 −1],        det C=1.
```

The separation rank is at most 4, and this bound is attained on disjoint
prime blocks. Take rows with patterns `(p²,p³,p²q,p²qr)` and columns with
the same patterns on three completely different primes. Their feature
matrices both equal

```
V = [ 3 2 1 1 ]
    [ 4 2 1 1 ]
    [ 6 4 2 1 ]
    [12 8 3 1 ],      det V=2.
```

Consequently

```
K = V C Vᵀ
  = [ 2  5  6 15 ]
    [ 5  9 12 27 ]
    [ 6 12 15 34 ]
    [15 27 34 73 ],      det K=4.
```

Concrete rows are `(4,8,12,60)` and columns `(49,343,539,7007)`.
All 16 products are legal and nonsquarefree. Hence no three-channel
rational separation can reproduce this defect kernel. ∎

Again, `δ` alone is not a continuation congruence. The smaller hostile is

```
δ(4)=δ(12)=−1,       δ(4·5)=−1 ≠ 0=δ(12·5).
```

Even the balanced subset is not preserved uniformly by a fresh prime:

```
δ(8)=δ(60)=0,        δ(8·7)=1 ≠ 3=δ(60·7).
```

The rank computation measures precisely an arithmetic response, not the
memory needed to determine whether an arbitrary named prime is fresh.

## 5. Exact guard memory and the minimal FUS quotient

### Proposition 5 (support is exactly the coprimality guard observer)

For positive integers `n,m`,

```
(gcd(n,b)=1 ⇔ gcd(m,b)=1 for every b≥1)
    iff rad(n)=rad(m).
```

**Proof.** Equal supports give equal coprimality tests. If the supports
differ, any prime in their symmetric difference is a distinguishing test.
∎

This cannot be replaced by a fixed finite state space: fixing a prime
`q` and `k` further distinct primes `p_1,…,p_k`, the nonsquarefree numbers
`q² product_{i∈S}p_i`, one for every subset `S`, have `2^k` distinct guard
behaviours. A fixed-length binary encoding for these behaviours needs at
least `k` bits. This is a support-specific statement; it does not follow
from the rank-4 lower bound.

In particular, even all four arithmetic coordinates can miss a guard:
`v(4)=v(9)=(3,2,1,1)`, but multiplication by 2 is legal only for 9.

### Proposition 6 (exact guarded defect quotient)

Take `X={n≥2:n is nonsquarefree}`, observation `δ`, and all named
multipliers `T_b(n)=nb` with guard `gcd(n,b)=1`. Then

```
n ≡_δ m   iff   (rad(n),τ(n))=(rad(m),τ(m)).
```

**Proof.** Contextual equality gives equality of every one-step guard,
hence equal radicals by Proposition 5. Thus `ω` and `2^ω` agree.
Equality of the present `δ` then gives equal `τ`.

Conversely, equal radicals and divisor counts give equal present `δ`
and equal guards. If a common multiplier is legal, then

```
rad(nb)=rad(n)rad(b),       τ(nb)=τ(n)τ(b).
```

The resulting pairs are again equal. Induction proves equality along
every legal continuation, including agreement on failures. ∎

This carrier does not retain the full exponent profile: `12=2²·3` and
`18=2·3²` both have `(rad,τ)=(6,6)` and are equivalent here. It is the
canonical coarsest quotient for the **named operation family**, not a
universal arithmetic state. If overlapping multiplication were permitted,
the common multiplier 2 would produce `δ(24)=1` and `δ(36)=2`, breaking
this equivalence. Operation scope is part of the theorem.

The same `(rad,τ)` quotient is exact for the vector observation `(F,U,S)`
on all `N≥2` under coprime multiplication: its radical supplies `ω`, while
`τ=2^ω` detects squarefreeness and `(ω,τ)=(1,2)` detects primality. Conversely,
guards determine the radical and `F` determines `τ`.

**Full-domain scalar corollary (boundary corrections retained).** In fact,
the scalar `δ` also has exactly this quotient on all `N≥2`. To prove the
extra boundary cases rather than assume them, let `r` be a fresh prime.
The corrected formulas give

```
c(N):=2δ(N)+ω(N)−δ(Nr)
     =[N squarefree]+2[N prime] ∈ {0,1,3}.
```

Indeed `Nr` is not prime, is squarefree iff `N` is, and has doubled
`τ` and `2^ω` and incremented `ω`. Substitution proves the displayed
identity. The three values distinguish nonsquarefree integers, squarefree
composites, and primes. Equal guarded residuals recover equal radicals,
then the same `c(N)` by a common fresh prime, hence the same two indicators
and the same `τ` from the corrected current-defect formula. Conversely,
`(rad,τ)` determines the indicators, output, guards, and legal updates.
Thus the full-domain contextual equivalence is again equality of
`(rad,τ)`. For example `p` and `p³` have equal current defect and support,
but a fresh-prime continuation distinguishes their boundary types.

### Proposition 7 (prime promotions force full exponent memory)

Keep the nonsquarefree domain and the coprime-multiplication language of
Proposition 6, and add, for each named prime `p`, the partial promotion

```
P_p(n)=pn,    defined iff p|n.
```

Then contextual equivalence for `δ` in the enlarged language is equality
of integers.

**Proof.** Old contexts already recover `(rad(n),τ(n))`. For each prime
`p|n`, write `e=v_p(n)`. Promotion changes the local divisor-count factor
from `e+1` to `e+2` and leaves `ω` unchanged. Therefore

```
δ(pn)−δ(n) = τ(pn)−τ(n) = τ(n)/(e+1),
e = τ(n)/(δ(pn)−δ(n))−1.
```

The denominator is a positive integer. Thus the common old state and all
named promotion responses recover every exponent on the recovered support,
hence `n` itself. Equality clearly suffices. ∎

An explicit hostile is `972=2²3⁵` versus `288=2⁵3²`. Both have
`(rad,τ,δ)=(6,18,11)` and are indistinguishable in the old language, but
promotion at the named prime 2 produces increments 6 and 3, respectively:
`δ(1944)=17` versus `δ(576)=14`. The smaller `12,18` hostile above is
another instance. This is a precise operation-language refinement, not an
analogy: adding experiments refines the contextual quotient, and here the
refinement is strict all the way to the identity relation. It does not mean
that every larger operation family necessarily increases memory.

## 6. Berggren: the contrasting failure is target descent

The recovered canon **THM-4808(c)** proves that for every fixed `M≥1`,
the root Euclid pair `(2,1)` and the pair

```
(16M+2,4M+1) ≡ (2,1) mod M
```

have opposite branch characters `χ=(-1)^(number of A letters)`. The
second is reached by the root-first word `A D^(M−1) U D`.

The observation consisting of the Euclid vector modulo `M` already
updates under each branch matrix modulo `M`: it is a genuine forward
congruence. What fails is **χ's descent through that congruence**. Adding
`χ` gives a finite forward observer, and `χ` alone needs exactly two
reachable states, since `A` flips it and `U,D` preserve it.

Thus it is false to infer “the sign requires storing the whole word”
from residue non-descent. The integer/tree position may determine its
word, and a sign can depend on all letters while admitting a two-state
streaming update. This is distinct from the unbounded support-guard
requirement proved for arbitrary named prime multipliers above.

## 7. Provenance and verification

### Recovered repository evidence

Sources below are pinned at commit
`ae7c073e864d05be74fed85ba0a2bc164e026a6d`.

* `07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md`:
  motivation only. It explicitly calls itself **REFLECTION, not canon**.
* `01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md`:
  proved operation-kernel criterion and minimum linear sidecar rank
  `dim N(ker R)`. This is direct repository prior art, not a new theorem here.
* `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md`, §3:
  the actual divisor-count definitions, correction terms, balanced shapes,
  and prime-adjunction recursion. The unrelated claims elsewhere in that
  source are not imported here.
* `01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md`,
  conventions and §(c): the explicit every-modulus hostile and branch matrices.
* `05-knowledge/results/session83_aliquot_divisor_grammar_20261028/FIBER_SHAPE_CENSUS.md`:
  arithmetic/source conventions and finite-census context; none of its census
  entries is promoted into an infinite-family statement here.

### Executable controls

Run `python3 check_observers.py`. It uses Python's standard library and exact
integers/rationals, not floating point. Explicit `require` checks remain
active with `python3 -O` as well as normal execution. The checker prints
**PASS** for:

* literal proper-divisor counts and balanced-shape classification for
  `2≤n≤5000`;
* the multiplication identities on all 8,771 coprime pairs `1≤a,b≤120`,
  including 2,464 with nonsquarefree left factor;
* both displayed exact minors and their determinants;
* both defect hostiles, the aliquot hostile, the guard hostile, and the
  operation-extension hostile;
* 1,868 fresh-prime reconstruction checks;
* the full-domain fresh-prime boundary classifier on `2≤N≤500`;
* prime-promotion reconstruction of every exponent for nonsquarefree
  `2≤n≤500`, including the `972,288` hostile.

These finite checks are reproducible controls. The universal claims are
proved above; they are not extrapolations from these windows.

### What the synthesis does and does not establish

The common mechanism is **contextual distinguishability**. A hostile pair
proves non-descent or failed update closure; a separating minor proves a
linear-channel lower bound; a guard test exposes omitted admissibility
memory. These are three related experiments, not interchangeable theorems.
No cross-domain isomorphism, uniform finite-state theorem, optimal bit
encoding, or new general automata theorem is claimed.
