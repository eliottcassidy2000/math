# Session259: diagonal strata, exact birth clocks, and the joint sum–product lift

Status: **PROVED in this note** for the elementary propositions with proofs below;
**INHERITED** for the stated THM-2422 definitions and laws; **FINITE-EXACT** only
for the accompanying small experiment; **OPEN** for the proposed future uses.
No claim of external novelty. No external literature theorem is needed.

## 0. Inheritance and the live concept board

The authoritative starting point is
`01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`.
Its labelled binary operation fibre, not its existential simple shadow, retains
parents, multiplicities, diagonal pairs, and birth times. Its swap-fixed loci
are `z=2a` and `z=a²`. Its hostile examples are prime squares 4 and 9, whose
strict proper-factor fibres cannot distinguish them from primes. Its corrected
near miss is the Fibonacci selected spine versus the full dyadic synchronous
closure. None of those inherited results is re-audited here.

The underused sidecar is **the pair itself, observed simultaneously through
sum and product**, suggested by
`05-knowledge/hypotheses/HYP-3003-summand-multiplicand-farey-basis-merge.md`.
That document and HYP-3004's dichotomy recursion atlas are synthesis, not
proofs. Their useful prescription is to specify what each quotient forgets.

`07-reflections/the-summand-alternation-is-an-involution-parity-law-kps-S131.md`
suggests the fixed-point viewpoint and the blue/black analogy. Here we use the
abstract involution mechanism; no historical SC census or rigidity claim is
needed. In particular `product-graph-sc-spine-fractal-dimensions.md` supplies no
numerical premise.

Live board and session portfolio:

1. Anchor: one-integer insertion into the labelled parent set; exact new-edge law.
2. Anchor: input birth versus output birth; keep both clocks.
3. Niche: the discriminant chart joining the two operation graphs losslessly.
4. Niche: involution parity as the first member of a prime-arity hierarchy.
5. Wildcard: the noncommuting doubling/squaring action on that chart.
6. Wildcard: blue/black fixed/free decompositions, with an explicit no-transfer boundary.

The primary gain is not another binary count. It is a reusable calculus in which
inserting one integer, changing arity, and moving between doubles and squares
are exact operations on the same retained parent data.

## 1. One new integer creates a star before either projection

Let `P_N={(a,b):1<=a<=b<=N}`. Define the sum/product chart

    Φ(a,b)=(s,p)=(a+b,ab),      Δ=s²−4p.

**Proposition 1 (lossless joint chart and birth).** Φ is a bijection from
positive unordered pairs onto

    J={(s,p) in Z²: p>=1, Δ=d² for an integer d>=0,
                    s>d, s≡d (mod 2)}.

The inverse is `a=(s−d)/2, b=(s+d)/2`. The least parent-prefix containing the
pair is `B(s,p)=(s+d)/2`. Thus

    Φ(P_N)={j in J:B(j)<=N}.

Proof: the displayed inverse consists of positive integers and has sum s and
product p; conversely the discriminant of a pair is `(b−a)²`. Uniqueness follows
from the ordered representative `a<=b`. This also proves the birth formula.

**Proposition 2 (the newly born strip).** Adjoining N to `[N−1]` adds exactly

    {(s,p)=(N+a,Na):1<=a<=N},
    equivalently p=N(s−N), N+1<=s<=2N.

There are N weak new pairs and N−1 strict new pairs. Their last point is
`(2N,N²)`, the one diagonal event. Proof: every new pair contains N, and the
other coordinate ranges independently through `[N]`.

The nested joint clouds therefore grow by straight finite strips. Their two
projections have very different collisions, although they arise from precisely
the same N new objects. A sum target receives a new pair only from
`N+1,...,2N`; a product target receives one only from `N,2N,...,N²`.
The exceptional unit pair `(1,N)` lands at product N and is not an ancestry
edge from two strictly earlier nodes. For proper-factor ancestry delete a=1;
for strictly distinct ancestry also delete a=N. These are different deletions.

There is a useful third clock. If a graph is exposed only on integer vertices
`[T]`, a pair and its operation target are all present first at

    T_+(a,b)=a+b,          T_×(a,b)=ab.

For addition the target is later than each parent. For multiplication this is
also true when `a>=2`, but equality holds at `a=1`. Consequently parent-prefix
insertion and integer-target exposure must never be silently identified.
The chart retains all three quantities `B,s,p` without ambiguity.

Example: (3,5) enters the parent cloud at 5, its sum edge at 8, and its product
edge at 15. The same edge has three honest ages depending on the question.

## 2. Doubles and squares are branch coordinates, not equivalent subsets

The diagonal is the discriminant-zero locus

    D={(2a,a²):a>=1}.

Projection to s produces the even numbers; projection to p produces the
squares. This identifies their **parent witness** a, not the two target sets
as subsets with the same natural-number filtration. At integer target cutoff T
there are `floor(T/2)` additive diagonal events but `floor(sqrt(T))`
multiplicative diagonal events. At parent cutoff N there are exactly N of each.
Thus “even versus square” is an exact branch-locus equivalence only after the
clock is transported. The discrepancy is linear versus quadratic target time.

For a restricted set S of available parents, the diagonal loci are `2S` and
`S²`; if targets must also lie in S, they become `S∩2S` and `S∩S²`.
This gives the precise thin-set version of the earlier twin-rank doubling
sidecar, without asserting any density or infinitude for such intersections.

## 3. A noncommuting dynamics on the growing joint cloud

Define, for positive integers c,

    D_c(s,p)=(cs,c²p),
    Q(s,p)=(s²−2p,p²).

**Proposition 3.** These maps preserve J and correspond respectively to
`(a,b)->(ca,cb)` and `(a,b)->(a²,b²)`. They obey

    Δ(D_c j)=c²Δ(j),        Δ(Qj)=s²Δ(j),
    B(D_c j)=cB(j),         B(Qj)=B(j)²,
    Q D_c = D_(c²) Q.

Proof: expand `(ca)+(cb)`, `(ca)(cb)`, `a²+b²=(a+b)²−2ab`, and
`(b²−a²)²=(a+b)²(b−a)²`. The larger positive parent stays the larger.
The composition law follows by squaring ca and cb.

Both operations preserve the diagonal/strict dichotomy; neither collapses a
strict pair onto the diagonal. But doubling and squaring do not commute:
`QD_2=D_4Q`, not `D_2Q`. On positive parent pairs the missing factor of two
cannot vanish. Example `(1,2)`: square after doubling gives `(4,16)`, while
double after squaring gives `(2,8)`.

This is a precise replacement for treating `*2` and “square” as parallel
recursions in HYP-3004. At the parent birth clock they act as `B->2B` and
`B->B²`; at the discriminant they act by different multipliers. A scale label
is necessary to compare different words. No statement about Collatz
termination follows, and no Collatz carrier is introduced just for resemblance.

More generally, define `L_0=2, L_1=s`, and

    L_k=s L_(k−1)−p L_(k−2).

Induction gives `L_k=a^k+b^k`. Parent kth powers act exactly by

For every integer k>=1,

    Q_k(s,p)=(L_k(s,p),p^k),
    Δ(Q_k j)=Δ(j) [sum_(i=0)^(k−1) b^(k−1−i)a^i]²,
    B(Q_k j)=B(j)^k,       Q_k D_c=D_(c^k) Q_k.

The divided-difference factor is positive, including when a=b by the displayed
sum, so every positive power preserves strictness. This is a polynomial lift
of power recursion on the merged graph, not a graph isomorphism of all J:
only pairs of kth powers lie in its image.

## 4. The same insertion calculus in two convolution algebras

For a finite parent set S, use the formal basis `[a]`. In the additive algebra
put `[a]*[b]=[a+b]`; in the multiplicative algebra put `[a]*[b]=[ab]`.
Let `F_S=sum_(a in S)[a]`. Define the operation-specific Adams map

    ψ_j([a])=[ja] (additive),       ψ_j([a])=[a^j] (multiplicative).

These are multiplicative endomorphisms of their respective algebras. Then

    O_2(S)=F_S*F_S,
    H_2(S)=(F_S*F_S+ψ_2(F_S))/2,
    E_2(S)=(F_S*F_S−ψ_2(F_S))/2.

O counts ordered pairs, H unordered weak pairs, E unordered distinct pairs,
with their operation targets retained as coefficient addresses. The fractions
are integral because swap acts freely off the diagonal. This is THM-2422's
binary mechanism in a common algebraic language, not an isomorphism between
integer addition and multiplication.

If n is not in S, then

    H_2(S∪{n})−H_2(S)=[n]*F_S+ψ_2([n]),
    E_2(S∪{n})−E_2(S)=[n]*F_S.

These identities prove Proposition 2 in both projections, for any restricted
parent set, not just a full prefix. They also show exactly where restoring
diagonals changes the growth process. In ordinary additive generating
functions `ψ_j` is substitution `x->x^j`; in a multiplicative Dirichlet
polynomial it is `s->js`. The coefficient addresses, and not just the function
values, must be retained.

### Modulo two and modulo four

Reducing the ordered count in characteristic two gives

    F_S*F_S = ψ_2(F_S) (mod 2).

At target z this says odd ordered count iff z=2a with a∈S, or z=a² with a∈S.
It gives no lower bound on the number of off-diagonal witnesses: they disappear
in pairs. The next bit restores their parity:

    E_2(S)(z) mod 2 = ((O_2(S)(z)−ψ_2(F_S)(z))/2) mod 2.

Thus an ordered count modulo four, together with the diagonal bit, determines
strict unordered parity. Example S=[3], additive target 4: ordered count 3,
diagonal count 1, strict count 1. An unordered parity claim cannot be deduced
from ordered parity alone. This is an exact hierarchy of retained bits.

## 5. Higher arity: parity/squares become divisibility/perfect powers

**Proposition 4 (prime cyclic fixed-point law).** For any prime ℓ and finite
S, let O_ℓ^+(z) and O_ℓ^×(z) count ordered ℓ-tuples in S with sum or product z.
Then

    O_ℓ^+(z) ≡ 1_(ℓ|z and z/ℓ∈S)                    (mod ℓ),
    O_ℓ^×(z) ≡ 1_(z=a^ℓ for some a∈S)               (mod ℓ).

Proof: rotate tuple coordinates by one step. Since ℓ is prime, each orbit has
size 1 or ℓ. The size-one tuples have all coordinates equal. Their sum/product
targets are ℓa and a^ℓ, and positivity makes a unique when it exists.
Equivalently this is the Frobenius identity `F_S^{*ℓ}=ψ_ℓ(F_S) mod ℓ`.

This is the exact arity generalization of even/square: at arity three it is
multiple-of-three/cube, and similarly for every prime. It does not say the
unordered ℓ-parent count has the same congruence.

The first composite hostile is arity four, S={1,2}, additive target 6.
There are six ordered tuples, so the count is 2 modulo 4 although no all-equal
tuple sums to 6. A two-periodic orbit survives. Thus replacing “prime ℓ” by
arbitrary arity is false. The missing coordinate is the stabilizer period.

### The full collision strata (not just the all-equal one)

Let H_k and E_k count weak and strictly distinct unordered k-parent tuples.
Products below use the chosen convolution algebra; t records arity. Then

    sum_(k>=0) H_k t^k = product_(a∈S) (1−t[a])^(−1),
    sum_(k>=0) E_k t^k = product_(a∈S) (1+t[a]).

These are formal series; H_0=E_0 is the algebra identity (target 0 for addition,
target 1 for multiplication). Expanding each factor proves the identities.
Inserting n gives, without any asymptotic assumption,

    H_k(S∪{n})−H_k(S)=sum_(j=1)^k [n]^{*j}*H_(k−j)(S),
    E_k(S∪{n})−E_k(S)=[n]*E_(k−1)(S).

The j summand records **how often the newest integer is repeated**. This is a
complete higher-arity new-node growth law. The binary diagonal is just j=2.
The ordered insertion law similarly has binomial factors `binom(k,j)`.

Cycle decomposition supplies the equivalent explicit formula

    H_k=sum_(λ⊢k) (1/z_λ) product_j ψ_j(F_S)^{*m_j},
    z_λ=product_j j^(m_j) m_j!,

where λ has m_j parts equal to j. Proof: quotient ordered tuples by coordinate
permutations; a tuple fixed by a j-cycle is constant along that cycle, giving
ψ_j(F_S). Counting permutations of each cycle type gives k!/z_λ.
This is the exact symmetric-power/plethystic bridge. It is a counting identity
in characteristic zero; division by k! must not be imported blindly modulo ℓ.

For example at arity three put C=F_S^{*3}, T=ψ_2(F_S)*F_S, D=ψ_3(F_S).
Then

    H_3=(C+3T+2D)/6,          E_3=(C−3T+2D)/6.

In addition T counts `2a+b=z`; in multiplication it counts `a²b=z`.
Thus doubles/squares are precisely the **two-equal collision stratum** inside
ternary growth, while triples/cubes give the all-equal stratum. Higher arity
does not abandon the earlier dichotomy; it nests it inside a partition-indexed
family of collision strata.

## 6. What transfers to blue/black and what does not

For any finite set X with involution σ, there is an exact decomposition into
fixed points and two-element orbits:

    |X|=|Fix σ|+2·(# free orbits).

For parent tuples σ is swap. In the S131 tiling analogy the relevant X is a
line fibre and σ is the specified reflection; “blue” refers to its fixed part
and “black” to its free part. The transferable theorem is this orbit equation,
and its fibrewise version whenever σ preserves the fibre. It is not a map
from individual integer pairs to SC nodes or tilings.

For growth, the honest extra hypothesis is an equivariant nested family
`X_N⊆X_(N+1)`. Only then is the new shell itself σ-stable and satisfies the
same fixed/free count law. The arithmetic parent-prefix family has this
property automatically. An arbitrary list of graph isomorphism classes or a
quotient that changes with N need not. No SC enumeration or interleave law is
inferred. A positive future target would be an actual equivariant insertion
map for the tiling family; until supplied, that transfer is **OPEN**.

## 7. Connection and loss ledger

| Source → target | Exact map/predicate retained | Lost coordinate and repair |
|---|---|---|
| Pair → joint chart | (a,b)↦(a+b,ab), bijection and birth | Only ordering lost; canonical a<=b repairs it |
| Joint chart → summand fibre | project s, preserve sum target | Product/difference and individual pair lost; retain p or d |
| Joint chart → factor fibre | project p, preserve product target | Sum and individual pair lost; retain s |
| Diagonal → doubles/squares | a↦(2a,a²), same fixed witness | Natural target clock differs; retain a or B |
| Pair power → chart polynomial | Q_k, preserve strict/weak and powers | Not surjective; retain kth-power image condition |
| Ordered fibre → mod ℓ | cyclic orbit cancellation | Nonfixed witness mass lost; keep integer count |
| Binary → k-parent growth | symmetric/exterior series, all collision strata | Two coordinates no longer reconstruct tuple for k>=3 |
| Involution pairs → blue/black | fixed/free orbit equation only | No object-level or growth map supplied |

The higher-arity loss is substantive: equal sum and product need not identify
unordered triples. `(1,6,6)` and `(2,2,9)` both have sum 13 and product 36,
but birth clocks 6 and 9. Thus the joint (s,p) chart is lossless in arity two
only. For k parents retain all elementary symmetric coordinates, equivalently
the monic polynomial with those positive integer roots, or retain the tuple.
This is the cheapest hostile against exporting the binary chart unchanged.

## 8. Consequences and next questions

**Direct proved growth information:** one new parent produces a labelled star;
its weak/strict binary difference is one point; its higher-arity difference is
partitioned by the multiplicity of the newly inserted integer. All formulas
hold for arbitrary finite restricted parent sets, so sparse-set experiments
need no new combinatorial framework.

**New exact bridge:** the joint discriminant chart makes summand and factor
fibres two projections of one growing integer object, and supplies compatible
birth, scaling, and power maps. It resolves the dense-even/sparse-square
mismatch by transporting the clock instead of claiming similar densities.

**Open/useful experiment:** for a restricted S, compare the earliest and latest
parent births inside a fixed target fibre before and after Q_k. The chart
transport controls images exactly, but a power target can have other parents
not in the image. Quantifying that extra fibre is a genuine growth question.

**Open/useful higher-arity question:** which collision partition first changes
the birth profile when weak parents are allowed? The exact insertion identities
make this finite experiment unambiguous and distinguish diagonal effects from
ordinary additional ancestry.

**Stopped tangents:** no forced Collatz transfer, square-sum Hamiltonicity,
SC census prediction, or LRC consequence is asserted. They require an actual
predicate-preserving carrier beyond these identities.

## 9. Reproduction

Run from repository root:

    python3 05-knowledge/results/session259_arithmetic_graph_growth_20261111/dichotomies/scout.py
    python3 -O 05-knowledge/results/session259_arithmetic_graph_growth_20261111/dichotomies/scout.py

The script uses direct tuple enumeration separately from algebraic coefficient
formulas. It records bounded universes, insertion checks, prime-arity controls,
and the two named hostile examples. The proofs above, not finite computation,
supply the all-N claims. Receipt and SHA-256 manifest are adjacent.
