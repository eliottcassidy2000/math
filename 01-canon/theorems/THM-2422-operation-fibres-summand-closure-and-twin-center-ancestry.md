---
id: THM-2422
title: "Swap-fixed operation fibres, dyadic summand closure, and twin-center ancestry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The commutative addition and multiplication fibres have exact diagonal
  corrections supported, respectively, on even targets and perfect squares.
  Distinct-summand closure from {2,3} is N minus {1,4,6}; its synchronous
  stages have one last transient hole and then the exact frontier
  M_(t+1)=2M_t-1. The proposed immediate A014574 recurrence fails first after
  startup at 348=312+36, but every normalized twin center k>=3 checked through
  center 100,000,000 has two distinct earlier twin-center parents. That last
  statement is FINITE-EXACT only; its all-n extension remains OPEN.
source: codex-2026-07-26-operation-fibre-ancestry
depends_on:
  - THM-361-product-sum-defect-normal-form
  - THM-362-natural-operation-graph-shadows
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
related:
  - HYP-1994-twin-goldbach-necklace
  - HYP-3003-summand-multiplicand-farey-basis-merge
external: "OEIS A014574, https://oeis.org/A014574"
script: 04-computation/twin_center_additive_ancestry_thm2422.py
output: 05-knowledge/results/twin_center_additive_ancestry_thm2422.out
script_sha256: f8de71d40d646494864dcf178732f199ab4782c88efadb6bde71b37e935c833d
output_sha256: 153af3af5517d78f1532d4d5c8b458fbd03cdd0f198df2cc8612d12a3420cb8f
hash_basis: working-tree bytes (LF)
---

# THM-2422 -- operation fibres, summand closure, and twin-center ancestry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem separates three objects that had repeatedly been blended:

```text
the existential simple operation shadow,
the labelled commutative parent fibre over one target,
the dynamics obtained by repeatedly adjoining such targets.             (1)
```

The simple additive shadow is just the order relation by THM-362. It cannot
see parity, a missing parent, or the difference between a Fibonacci spine and
the whole closure. Those structures live in the labelled fibre.

## 1. The swap-fixed sidecar

For `z>=1`, define the weak and strict unordered additive fibres

```text
A_z^<= = {{a,b}: 1<=a<=b, a+b=z},
A_z^<  = {{a,b}: 1<=a< b, a+b=z}.                            (2)
```

For `z>=2`, define `M_z^<=` and `M_z^<` by replacing `a+b=z` with
`ab=z`. Put

```text
delta_+(z)=1 when z is even and 0 otherwise,
delta_x(z)=1 when z is a perfect square and 0 otherwise.     (3)
```

Then

```text
|A_z^<=| = floor(z/2),
|A_z^< | = floor(z/2)-delta_+(z)
          = floor((z-1)/2)
          = (2z-3-(-1)^z)/4,                                (4)

|M_z^<=| = (tau(z)+delta_x(z))/2,
|M_z^< | = (tau(z)-delta_x(z))/2.                            (5)
```

Here `tau` is the divisor-counting function. The mechanism is the same in
both formulas: quotient ordered parent pairs by the transposition
`(a,b)<->(b,a)` and retain its fixed set. The fixed equation is

```text
2a=z                       for addition,
a^2=z                      for multiplication.               (6)
```

Thus the exact multiplicative analogue of additive parity is the sparse
perfect-square wall, not another period-two sequence.

If the unit pair `{1,z}` is deleted from the multiplicative fibre, then

```text
|{{a,b}:1<a<=b, ab=z}| = (tau(z)-2+delta_x(z))/2,
|{{a,b}:1<a< b, ab=z}| = (tau(z)-2-delta_x(z))/2.             (7)
```

The ordered proper-divisor witness set has size `tau(z)-2`, so

```text
z is prime iff tau(z)-2=0.                                   (8)
```

More sharply,

```text
the weak proper factor fibre is empty iff z is prime;
the strict proper factor fibre is empty iff z is prime or a prime square.
                                                                  (8a)
```

Indeed, if composite `z` has no strict proper pair and `p` is its least prime
divisor, then `p<=z/p`; strict emptiness forces equality, so `z=p^2`.
The diagonal sidecar is therefore essential: `4=2^2` is the first composite
that the strict fibre alone misclassifies as an atom. It is simultaneously
the first absent additive node unlocked by restoring `4=2+2`. The value
`9=3^2` is the next, non-startup multiplicative hostile.

### Proof

For `z=2m`, the strict additive first parent is
`a=1,...,m-1`; for `z=2m+1`, it is `a=1,...,m`. This proves
(4). Divisors of `z` pair under `a<->z/a`; exactly one divisor is fixed
precisely when `z` is a square, proving (5). Removing `{1,z}` gives (7),
and (8) is the definition of multiplicative atomicity.

Put `r(z)=|A_z^<|`. Summing (4) gives

```text
sum_(z<=N) r(z) = floor((N-1)^2/4),                           (9)

r(N+1)-r(N) = 1 if N is even, and 0 if N is odd,             (10)

Delta^2 r(N)=(-1)^(N+1).                                     (11)
```

Equivalently,

```text
sum_(z>=1) r(z)X^z = X^3/((1-X)(1-X^2)).                     (12)
```

This is the precise sense in which the growing summand fibre is alternating:
it is a period-two boundary quasipolynomial. The infinite simple graph itself
does not alternate.

## 2. The startup module and the full distinct-summand closure

Let `Cl^<(S)` be the smallest set containing `S` and closed under

```text
a,b in S, a<b  =>  a+b in S.                                 (13)
```

Then

```text
Cl^<({2,3}) = {2,3,5} union {z:z>=7}
            = N minus {1,4,6}.                               (14)
```

Indeed,

```text
5=2+3,             7=2+5,             8=3+5,                 (15)
```

and for every `z>=9`,

```text
z=2+(z-2)                                                   (16)
```

uses distinct already-live parents. Conversely, `1` has no positive
parents, the only strict parent of `4` is `{1,3}`, and the strict parents
of `6` are `{1,5}` and `{2,4}`. Hence those three nodes never appear.

Restoring equal parents changes the answer sharply:

```text
4=2+2,              6=3+3,
Cl^<=({2,3})=N minus {1}.                                    (17)
```

There is a useful recursive description. For a nonseed target let `h(z)=1`
when it is absent and `h(z)=0` when it is present. Then

```text
h(z) = AND_(a<b, a+b=z) (h(a) OR h(b)).                       (17a)
```

A target is missing exactly when the missing set hits every labelled parent
edge. With `1` absent and `2,3` forced present, the sole strict edge
`1+3=4` forces `h(4)=1`; the two edges `1+5=6` and `2+4=6` then force
`h(6)=1`. At `7`, the edge `2+5` escapes the missing set, and (16) supplies
an escaping edge thereafter. This AND-of-ORs ancestry law is the precise
hypergraph analogue of the user's transitivity intuition; the simple shadow
forgets every clause.

The target packet `{1,4,6}` also has an exact operation-weld reading from
THM-361:

```text
1 = the multiplicative unit and absent additive source,
4 = 2+2 = 2*2, the unique binary product-sum target,
6 = 1+2+3 = 1*2*3, the unique distinct positive product-sum target. (17b)
```

For target `4`, deleting the unit-parent edge `1+3` and the diagonal `2+2`
removes every route. For target `6`, the edges `1+5` and `2+4` inherit the
two earlier holes while `3+3` is the deleted diagonal. This is an exact
shared operation-fibre boundary, not a claim that multiplication causes the
additive closure.

For every `z>=13`, the missing module `{1,4,6}` deletes exactly the three
otherwise-valid parent pairs containing those entries. Therefore the
parent-fibre size inside (14) is

```text
r_Cl(z)=r(z)-3,                  z>=13.                       (18)
```

The irregular beginning propagates as a constant labelled-fibre deficit,
not as missing later vertices. This information disappears completely in
the existential shadow `x<z`.

## 3. Synchronous closure has an exact dyadic self-similarity

The asynchronous closure (14) hides the time at which a parent becomes
available. Retain it by defining

```text
S_0={2,3},
S_(t+1)=S_t union {a+b:a,b in S_t, a<b}.                     (19)
```

Direct calculation gives

```text
S_1={2,3,5},
S_2={2,3,5,7,8},
S_3={2,3,5,7,8,9,10,11,12,13,15},                            (20)

S_4={2,3,5} union [7,28].                                    (21)
```

Thus `14` is the last transient hole. For example, the new upper interval in
`S_4` has the witnesses

```text
14=3+11, 16=3+13, 17=2+15, 18=5+13, 19=7+12,
20=7+13, 21=8+13, 22=7+15, 23=8+15, 24=11+13,
25=10+15, 26=11+15, 27=12+15, 28=13+15.                     (22)
```

Suppose for some `t>=4` that

```text
S_t={2,3,5} union [7,M],              M>=28.                 (23)
```

No new sum exceeds `M+(M-1)=2M-1`. Conversely, for every
`M+1<=z<=2M-1`, choose

```text
a=floor((z-1)/2),          b=ceil((z+1)/2).                  (24)
```

Then `7<=a<b<=M` and `a+b=z`. Hence

```text
S_(t+1)={2,3,5} union [7,2M-1].                              (25)
```

Since `M_4=28`,

```text
M_t=27*2^(t-4)+1,                    t>=4.                    (26)
```

This is the exact toothpick-like self-similarity: after one finite startup
repair, the distance `M_t-1` doubles at every synchronous generation.

The Fibonacci chain

```text
2,3,5,8,13,21,...                                           (27)
```

is a lawful selected ancestry path in (19), but it is not the forced
backbone or the full frontier. Already the Fibonacci target `13` has three
live parent pairs,

```text
13=2+11=3+10=5+8,                                           (28)
```

and its growth rate `phi` differs from the dyadic frontier (26).

## 4. What the A014574 recurrence actually says

Let

```text
C={c>=3:c-1 and c+1 are prime}.
```

By THM-2413 this is OEIS A014574. Apart from the exceptional center `4`,
all centers are divisible by six, so write

```text
C={4} union 6K,
K={k>=1:6k-1 and 6k+1 are prime}
 ={1,2,3,5,7,10,12,17,18,23,25,30,...}.                     (29)
```

Order `K` increasingly as `(k_i)`. The proposed local recurrence has the
exact reformulation

```text
k_i=k_(i-1)+an earlier K-term
iff g_i:=k_i-k_(i-1) belongs to K.                           (30)
```

It is therefore a self-gap incidence condition, not a Fibonacci recurrence.
After the exceptional raw transition `6=4+2`, the A014574 recurrence holds
from `12` through `312` and then fails:

```text
348=312+36,                  36 notin C,                     (31)

58=52+6,                       6 notin K.                    (32)
```

Thus the all-n immediate-predecessor claim is **REFUTED**, with A014574
term 21 as the first post-startup hostile.

There is a stronger repair. For `k in K`, define the distinct parent fibre

```text
P(k)={(a,b):a,b in K, a<b, a+b=k}.                           (33)
```

When this fibre is nonempty, choose

```text
a*(k)=min{a:(a,k-a) in P(k)},       b*(k)=k-a*(k).           (34)
```

Because `K` is ordered by value, `b*(k)` is the latest possible legal
ancestor. This gives a canonical answer to “which previous terms combine.”
At the first local failure, scanning ancestors backward gives gaps

```text
52 -> 6,       47 -> 11,       45 -> 13,       40 -> 18 in K,
```

and hence

```text
58=18+40,                    348=108+240.                    (35)
```

Define the repair depth as the number of necklace vertices strictly between
`b*(k)` and `k`. The immediate recurrence is depth zero; (35) has depth
three.

The exact computation establishes the finite statement

```text
for every k in K with 6k<=100,000,000 and k>=3,
P(k) is nonempty.                                            (36)
```

Its universe contains `440,312` A014574 centers, of which `440,311` belong
to `6K`; all `440,309` checked targets `k>=3` have distinct parents. The
largest observed repair depth is `123`. In contrast, only `185,644` of the
`440,310` normalized nearest-predecessor transitions satisfy (30), so the
finite data distinguish the robust global repair from the failing local
rule.

Statement (36) is **FINITE-EXACT**, not a theorem for all twin centers. The
all-n assertion

```text
K minus {1,2} subset K+_distinct K                            (37)
```

remains **OPEN**. It is a thin restricted-addition problem and is not implied
by the elementary closure (14).

## 5. The prime-sextuple carrier has no local obstruction

A parent triple `a,b,a+b in K` is equivalent to simultaneous primality of

```text
6a-1, 6a+1, 6b-1, 6b+1, 6(a+b)-1, 6(a+b)+1.                 (38)
```

For a prime `p>=5`, scale residues by the invertible factor six and put
`x=6a`, `y=6b`. The forbidden residues are the six lines

```text
x=+-1,             y=+-1,             x+y=+-1.              (39)
```

Each line has `p` points. Parallel pairs are disjoint. Lines in different
directions have twelve pairwise intersections, and there is no triple
intersection: a sum of two signs is `-2,0,2`, never `+-1` modulo `p>=5`.
Inclusion--exclusion therefore leaves exactly

```text
p^2-6p+12                                             (40)
```

allowed ordered residue pairs. This is positive for every `p>=5`; the
forms in (38) automatically avoid `p=2,3`. Hence the six-form family is
admissible.

Admissibility explains why no congruence class can refute (37), but it does
not prove (37). The missing ingredient is a uniform global lower bound for
the representation fibre `P(k)` along the already-prime target set `K`.

There is also a targetwise refinement of (40). For a lawful target residue
`w=6(a+b) notin {+-1}` modulo `p`, the number of ordered parent residues is

```text
p-2,                         w=0,
p-3,                         w=+-2,
p-4,                         otherwise.                     (41)
```

The forbidden first-parent residues are
`{+-1,w-1,w+1}`. Their collisions give exactly the three cases in (41);
summing (41) over the `p-2` lawful target residues recovers (40).

Modulo five this exposes a concrete thinning. Apart from the prime-equals-five
startup `k=1`,

```text
K mod 5 subset {0,2,3}.                                      (42)
```

The ordered parent-channel counts from `{0,2,3}^2` into target residues
`0,1,2,3,4` are

```text
3,1,2,2,1.                                                   (43)
```

Thus a target already known to lie in `K` avoids both one-channel residue
classes, whereas the unrestricted `K+K` problem in HYP-1994 must cover them.
Consistently, ten of the eleven currently recorded HYP-1994 holes are
`1 mod 5`. This explains a real difference between the two problems, but
does not supply the missing global lower bound.

## 6. Scope and connection ledger

The lawful connection is

```text
source:
  the labelled commutative operation cospan;

target:
  the summand closure and the twin-center parent fibre;

map:
  retain the full preimage of each target under (a,b)|->a+b;

preserved:
  parent multiplicity, the swap-fixed diagonal, ancestry time, and missing
  dependencies;

destroyed by the simple shadow:
  the second parent, parity, the module {1,4,6}, and repair depth;

needed sidecars:
  weak-versus-strict diagonal data and the synchronous generation;

cheapest hostile tests:
  4 and 9 for the multiplicative diagonal, 14 for the synchronous startup,
  13 for nonunique Fibonacci ancestry, and 348 for the A014574 local rule.
```

There is no proved causal passage from the missing additive nodes
`{1,4,6}` to the distribution of primes. The weld with A014574 is instead
the exact multiplicative-atom predicate on the six linear forms (38).
Likewise, taking logarithms turns products into sums but sends the integer
lattice to the nonuniform set `{log n}`; it does not identify the summand
and multiplicand fibres.

## 7. Exact verification

Run

```bash
python3 04-computation/twin_center_additive_ancestry_thm2422.py
python3 -O 04-computation/twin_center_additive_ancestry_thm2422.py
```

The companion:

- constructs the exact byte sieve through `100,000,001`;
- checks every A014574 center and canonical parent in (36);
- reproduces all `440,309` canonical pairs by an opposite traversal from the
  latest possible second parent backward;
- records the complete repair-depth histogram;
- verifies (4)--(11), (14), and (18) through target `10,000`;
- independently constructs `S_0,...,S_10` and checks (20)--(26);
- verifies the additive and multiplicative diagonal formulas, with `9` as
  hostile control; and
- enumerates (40) for twelve primes from `5` through `43`.

Every truth-bearing check uses explicit exceptions and runs identically
under optimized Python. The all-n elementary statements are proved above;
the script supplies the separately labelled finite-exact scope (36).

## 8. Independent hostile audit

An independent implementation rebuilt the A014574 centre set and distinct
parent fibres without importing the companion's ancestry routine. It found
exactly `440,312` centres through `100,000,000`, hence `440,311` normalized
values in `K`, and checked all `440,309` targets `k>=3`. Every target had a
distinct earlier-parent representation; the final centre was `99,999,588`.
It also independently located the first post-startup failure of the local
predecessor rule at

```text
348=312+36,
```

and verified that the repaired parent fibre is nonempty there. Normal,
optimized, and stored transcripts byte-match, and both recorded SHA-256
hashes match the theorem header.

The audit does **not** promote (37): distinct-parent ancestry beyond centre
`100,000,000` remains open. This finite/all-n boundary is part of the
audited statement.
