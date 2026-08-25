---
id: THM-4094
title: "Hamiltonian matching deficit and two-prime-lane completeness reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Full deletion
  incidence gives a left-perfect matching and the exact Hamiltonian-path
  deficit as excess compatible
  insertions plus orphan paths; retaining only a selected matching loses that
  deficit, minimally at the transitive-triple/C3 boundary. Using the proved
  omissions 7,21, the proved strong carry atoms 49,63,343, and order-join
  multiplicativity, global H-spectrum completeness is equivalent to strong
  realization of the two infinite lanes p (odd prime p!=7) and 7p (odd prime
  p!=3). Finitary many-sorted
  first-order theories whose incidence fibers are all finite have a uniform
  fiber bound, so the unbounded finite tournament/path universe requires
  external finite/standard semantics. The global H-spectrum conjecture
  remains OPEN.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-001-redei
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-1862-order-join-reduction-principle
related:
  - THM-012b-insertion-decomposition
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-796-three-sorted-recursive-tiling-line-node-incidence
  - MISTAKE-209
  - MISTAKE-219
  - MISTAKE-499
script: 04-computation/hamiltonian_matching_deficit_two_prime_lane_thm4094.py
output: 05-knowledge/results/hamiltonian_matching_deficit_two_prime_lane_thm4094.out
script_sha256: 5b6c779958ad1de276873eabd3585f432bb593e32811fb138dc5d3c19948d616
output_sha256: 5f174f77a873a02bd92bbebe4a5d70b3c329abc8152fb1330d4fe1dc6d98df34
semantic_sha256: a7376b1f1d0bdea71e063caf85415a6555d070130c24304b2b2c6fb24571d76f
hash_basis: raw LF bytes for files; canonical compact JSON for the semantic ledger
audit: >
  ACCEPT after making the finite carry hypotheses explicit, strengthening the
  composable sidecar to full deletion fibers and cut signatures, repairing the
  singleton compactness witness, and scoping deletion notation to |T|>=2.
  The independent hostile audit replayed normal/-O/frozen transcripts and
  checked the deletion identity, minimal hostile, prime extraction, exhaustive
  seven-adic case split, compactness boundary, and exact hashes.
---

# THM-4094 -- Hamiltonian matching deficit and two-prime-lane completeness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem separates three claims compressed by a "matching completeness"
intuition. Full deletion incidence computes the local Hamiltonian-path
response exactly. A selected left-perfect matching retains only monotonicity
and loses the response deficit. Even exact monotonicity, padding, finite
coverage, and multiplicative strong-component closure do not supply the new
strong atoms required for global spectrum completeness.

There is also a precise logical boundary: ordinary finitary many-sorted
first-order semantics cannot require unbounded finite incidence fibers while
excluding infinite fibers in every model. This is a compactness obstruction,
not a disproof of a standard finite existence statement. No literature-
priority claim is made.

## 1. Notation and inherited facts

For a finite tournament `T`, write

```text
Ham(T) = set of directed Hamiltonian paths of T,
H(T)   = |Ham(T)|.                                      (1)
```

Let

```text
S     ={H(T):T a finite tournament},
S_str ={H(C):C a finite strongly connected tournament}. (2)
```

We use three proved facts, cited by ID and full slug because the legacy
`THM-1370` namespace collides.

1. `THM-001-redei` proves that `H(T)` is odd for every tournament `T`.

2. `THM-1862-order-join-reduction-principle` proves that if `T` has strong
   components `C_1,...,C_r` in condensation order, then

   ```text
   H(T)=product_i H(C_i),                               (3)
   ```

   and `H(T_1 ▷ T_2)=H(T_1)H(T_2)`.

3. `THM-1370-h-spectrum-omits-7-21-all-n` proves

   ```text
   7,21 notin S                                         (4)
   ```

   and proves that every other odd value through `609` belongs to `S`. It
   explicitly leaves

   ```text
   S = positive odds minus {7,21}                       (5)
   ```

   as a conjecture.

Source or sink padding is the order join with a singleton. It preserves `H`;
it changes the order of a witness but creates no new value.

## 2. Full deletion incidence and the exact deficit

Fix a tournament `T` with at least two vertices and a vertex `v`. Put

```text
L_v=Ham(T-v),                 R_v=Ham(T).               (6)
```

Define a bipartite relation `D_v` from `L_v` to `R_v` by

```text
P D_v Q
 iff deleting v from the vertex word Q gives exactly P
     and P is still a directed Hamiltonian path of T-v. (7)
```

For `P in L_v`, put

```text
a_v(P)=#{Q in R_v:P D_v Q}.                             (8)
```

Let `O_v` be the number of **orphan** paths `Q in R_v`: paths for which
deleting `v` does not give a Hamiltonian path of `T-v`. An orphan necessarily
has `v` internally between vertices `x,y` with

```text
x -> v -> y -> x,                                       (9)
```

so its failed shortcut marks a directed triangle through `v`.

> **Theorem 2.1 (deletion-incidence deficit).** Every left degree satisfies
> `a_v(P)>=1`; every right degree in `D_v` is zero or one; and
>
> ```text
> H(T)-H(T-v)
>   =sum_(P in Ham(T-v))(a_v(P)-1)+O_v.                 (10)
> ```
>
> Consequently any choice of one neighbor of each left vertex is an
> injection `Ham(T-v)->Ham(T)`, while its unmatched right vertices are
> counted exactly by the right side of `(10)`.

### Proof

Write a path in `T-v` as `P=(x_1,...,x_m)`. Set `s_i=1` when `v->x_i`, and
`s_i=0` when `x_i->v`. There is always a legal insertion position:

- if `s_1=1`, insert `v` before `x_1`;
- if `s_m=0`, insert `v` after `x_m`; and
- otherwise `s_1=0,s_m=1`, so the binary word has a `0 -> 1` transition;
  insert `v` at one such transition.

Thus `a_v(P)>=1`. Deleting `v` from a non-orphan target path produces exactly
one ordered old path, so every right vertex of `(7)` has degree one unless it
is an orphan, when it has degree zero. Therefore

```text
R_v
 = disjoint_union_(P in L_v) D_v(P)  disjoint_union  Orphan_v.       (11)
```

Taking cardinalities gives

```text
H(T)=sum_(P in L_v)a_v(P)+O_v.                          (12)
```

Subtracting `|L_v|=H(T-v)` proves `(10)`.

For completeness, the legal positions are the initial position when
`s_1=1`, the internal `0 -> 1` transitions, and the final position when
`s_m=0`. Binary transition balance gives

```text
a_v(P)=1+#{i:s_i=1,s_(i+1)=0}.                         (13)
```

Thus `a_v(P)-1` is exactly the Type-II count of
`THM-012b-insertion-decomposition`. Finally, choose one neighbor in every
nonempty left fiber. Distinct left fibers are disjoint by uniqueness of
deletion, so the choice is injective and its unmatched targets are exactly
the deficit `(10)`. QED.

## 3. The selected-matching quotient fails minimally

Let the **selected-matching quotient** retain only

```text
(L_v, chosen image M_v subset R_v, the bijection L_v -> M_v),        (14)
```

discarding `R_v\M_v` and the rest of `D_v`. The target cardinality and the
deficit `(10)` do not factor through `(14)`.

The smallest actual hostile uses the common two-vertex base `0->1`.

1. Insert vertex `2` as a sink. The transitive triple has

   ```text
   |L_2|=1, |R_2|=1, left degree (1), O_2=0, deficit 0. (15)
   ```

2. Insert vertex `2` to obtain `0->1->2->0`. The directed triangle has

   ```text
   |L_2|=1, |R_2|=3, left degree (2), O_2=1, deficit 2. (16)
   ```

After selecting one extension, both quotients `(14)` are the same one-edge
matching. Their exact increments are nevertheless zero and two. No hostile
exists below order three because every tournament on at most two vertices is
transitive.

This establishes the precise boundary:

```text
full two-sort deletion incidence: exact;
selected matching quotient:       monotonicity only.                 (17)
```

There is no claim that "two sorts are inherently insufficient." Keeping the
complete right sort and incidence retains the cardinality. The loss is caused
by deleting unmatched witnesses. For one scalar update, the total `(10)` is
sufficient. The count profile `(a_v(P)-1)_(P in L_v)` and the orphan
total/path list are necessary but are not sufficient for composition: a
faithful next-insertion sidecar retains the full fibers `D_v(P)`, the orphan
paths, and their legal and failed cut signatures.

## 4. Support is not object or witness multiplicity

The exact referee records two actual tournaments with `H=9`: the strong
order-five almost-transitive tournament `R_5`, and the reducible tournament
`C_3 ▷ C_3`, whose value is `3*3=9` by `(3)`. They contribute

```text
{9}       to value support,
[9,9]     to an object-indexed term multiset.            (18)
```

In a two-sort model with tournament sort `X`, path sort `Y`, and ownership
incidence `E(x,y)`, the spectrum is

```text
{|E(x,-)|:x in X},                                      (19)
```

not `|Y|` and not the object-indexed multiset of fiber sizes. This is the
MISTAKE-209 support-versus-multiplicity guardrail in the present setting.

## 5. Prime and seven-prime strong-atom extraction

Let `P_odd` denote the odd primes.

> **Lemma 5.1 (ordinary prime lane).** If `p in P_odd` and `p in S`, then
> `p in S_str`.

By `(3)`, a tournament with `H=p` has a product of positive integral strong-
component values equal to the prime `p`. Exactly one component has value `p`;
all other factors are one. This is the required strong witness. QED.

> **Lemma 5.2 (seven-prime lane).** If `p in P_odd` and `7p in S`, then either
> a strong component has value `7p`, or a strong component has value `7`.
> Consequently, by `(4)`, every realized `7p` belongs to `S_str`.

If no component has the whole value `7p`, the nonunit component factors form
a proper factorization of `7p`. One is divisible by seven and is a proper
divisor, hence equals `7` (also when `p=7`). This contradicts `(4)`. QED.

The finite prefix supplies three useful carry atoms. Values `49,63,343` are
attained, and each is strongly attained:

- a proper component factorization of `49` requires a factor `7`;
- a proper component factorization of `63=3^2*7` requires `7` or `21`; and
- a proper component factorization of `343=7^3` requires `7`.

All are impossible by `(4)`, so

```text
49,63,343 in S_str.                                    (20)
```

## 6. Exact two-prime-lane equivalence

> **Theorem 6.1 (global completeness reduction).** Under Rédei parity,
> `(3)`--`(4)`, and the proved strong carry atoms `(20)`, the global conjecture
> `(5)` is equivalent to
>
> ```text
> P_odd \ {7}       subset S_str,                       (21)
> 7(P_odd \ {3})    subset S_str.                       (22)
> ```

The second lane includes `p=7`, hence includes `49`.

### Necessity

If `(5)` holds, every odd prime other than seven is attained, so Lemma 5.1
gives `(21)`. For every odd prime `p!=3`, the value `7p` is odd and is neither
`7` nor `21`; hence it is attained. Lemma 5.2 gives `(22)`.

### Sufficiency

Assume `(21)` and `(22)`. Let `n` be a positive odd integer other than `7`
and `21`. Write

```text
n=7^a m,                    gcd(m,7)=1.                (23)
```

Every prime factor of `m` belongs to `(21)`. We construct `n` by order joins.

1. If `n=1`, use the singleton tournament.  Otherwise, if `a=0`, join the
   strong witnesses for the prime factors of `m`.
2. If `a>0` is even, use `a/2` copies of the strong `49` witness and the
   prime witnesses for `m`.
3. If `a>=3` is odd, use `343`, then `(a-3)/2` copies of `49`, then the prime
   witnesses for `m`.
4. Suppose `a=1`. If `m=3^k`, the exclusions force `k>=2`; use
   `63*3^(k-2)`.
5. In the remaining `a=1` case, choose a prime `p>=5` dividing `m`. It is not
   seven by `(23)`. Use the `7p` witness from `(22)` and the ordinary-prime
   witnesses for `m/p`.

In every case, repeated order join makes the Hamiltonian-path count exactly
the displayed product `n`. Hence every allowed odd `n` lies in `S`; together
with `(4)` this proves `(5)`. QED.

The connection to a prime-place observer is real, not formal. The valuation
`nu_7` is additive under the exact SCC product `(3)`. Ordinary prime atoms
carry `nu_7=0`; the exceptional `7p` lane carries the lone `nu_7=1` that
cannot be separated as `7*p`; and `49,343` transport higher seven-adic
valuation. No analogous exceptional lane is needed elsewhere because the
ordinary prime itself is not forbidden there.

## 7. Current frontier and a construction hostile

The prefix through `609` proves the ordinary-prime lane through `p=607` and
the seven-prime lane through `p=83`, since `7*83=581`. It also supplies the
carry atoms `(20)`. Multiplication advances beyond the literal prefix:

```text
611=13*47 in S.                                         (24)
```

At the scope of this theorem, the first target not forced by the prefix and
multiplication is prime `613`, and the first seven-prime target is
`7*89=623`. These are targets relative to the prefix `(25)`, not asserted
gaps. THM-4097 subsequently gives explicit strong witnesses for both and
extends the current first-unforced lane values to `2,887` and
`7*419=2,933`. THM-4102 then moves them to `14,657` and `14,777`, and
THM-4104 moves the current frontier to `80,407` and `7*11,527=80,689`.

There is a sharp actual-tournament hostile to inference from the finite prefix
and value-controlled operations alone. Put

```text
B={odd m:1<=m<=609, m notin {7,21}}.                    (25)
```

For each `b in B`, choose a tournament `T_b` supplied by the proved prefix,
and let `C_B` be the class generated by order join and source/sink padding.
Then

```text
{H(T):T in C_B}=<B>_times.                              (26)
```

Every generating operation gives the indicated product, and every finite
product is realized by repeated order join. Since `613` is prime and is not a
generator, it is absent from `(26)`. Thus exact finite coverage, padding, and
multiplicative closure do not entail the next prime lane. The class `C_B`
contains genuine tournaments; it is not claimed closed under all arbitrary
vertex insertions.

Equation `(10)` explains what an arbitrary-insertion attack must add. To
construct a prime target `p` from an old value `h`, it must solve

```text
sum_P(a_v(P)-1)+O_v=p-h,                               (27)
```

not merely exhibit a left-perfect matching saying the left side is
nonnegative. The `7p` lane has the same requirement plus its indivisible
seven-adic carry.

## 8. Finitary many-sorted compactness obstruction

> **Theorem 8.1 (finite-fiber uniform bound).** Let `L` be any finitary
> many-sorted first-order language, `Theta` an `L`-theory, and `E(x,y)` a
> formula with `x` in an object sort and `y` in a witness sort. If every fiber
>
> ```text
> E_M(a)={b:M satisfies E(a,b)}                         (28)
> ```
>
> is externally finite for every `M models Theta` and every `a`, then one
> finite `N` satisfies `|E_M(a)|<=N` for every model and object.

### Proof

Suppose there is no uniform bound. Expand `L` by a constant `c` of the object
sort. For each standard `r>=1`, let `phi_r(c)` say that distinct witnesses
`y_1,...,y_r` satisfy `E(c,y_i)`. Every finite subset of

```text
Theta union {phi_r(c):r>=1}                             (29)
```

has a model: choose a fiber at least as large as the greatest requested `r`
and interpret `c` by its owner. By first-order compactness, all of `(29)` has
a model. Its `c`-fiber has at least `r` elements for every standard `r`, so it
is externally infinite, contradiction. QED.

Repeated order joins of `C_3` have value `3^k`, so finite tournament/path
fibers are unbounded. Theorem 8.1 implies that the class requiring all
intended path fibers to be finite and unbounded is not elementary under
ordinary all-structure semantics.

Adding any ordinary first-order sorts does not remove compactness. An
arithmetic sort and function `H:X->N` may take a nonstandard value whose
internally finite initial segment is externally infinite. The repair must use
external finite-model semantics, a standard finite code/rank convention, an
explicit metatheoretic schema, or stronger categorical semantics.

For each fixed standard `h`, exact fiber size is first-order expressible:

```text
Exact_h(x) := exists distinct y_1,...,y_h
              forall y [E(x,y) iff y is one of y_1,...,y_h].          (30)
```

Thus global completeness is the external schema

```text
{exists x Exact_h(x): h odd, h notin {7,21}}.           (31)
```

Compactness does not refute `(31)`: different sentences may use different
finite objects, and a nonstandard model may contain additional infinite
fibers while satisfying every standard instance.

A sparse matching model shows separately that unboundedness is not
completeness. Let

```text
X={x_k:k>=0},
E(x_k)={(k,i):0<=i<3^k},
x_k ▷ x_l=x_(k+l).                                      (32)
```

Use the canonical product bijection on fibers, and inject `E(x_k)` into
`E(x_(k+1))` by `i->i`. This has finite odd unbounded fibers, exact
multiplicativity, and insertion monotonicity, but support only `{3^k:k>=0}`.
Its elementary theory necessarily has an infinite-fiber model. There is no
finite compactness countermodel: violating external fiber finiteness is
inherently infinite.

## 9. Loss map and mandatory sidecars

The relevant quotient is

```text
(finite tournament code, full old/new path fibers, full deletion incidence)
       |
       v
(old fiber, one selected extension per old path).        (33)
```

It preserves one chosen extension per old path, injectivity, monotonicity,
and equality under neutral source/sink padding. It destroys extra compatible
insertions, orphan target paths, the exact increment, standard finite rank in
arbitrary elementary models, support versus multiplicity, and any target-
address-to-tournament realization rule.

The repairs are typed:

1. **Local scalar response:** retain the total `(10)`.
2. **Composable insertion response:** retain the full deletion fibers, orphan
   paths, and their legal/failed endpoint-cut signatures.
3. **Logical finiteness:** retain an external standard finite code or work in
   explicitly finite semantics.
4. **Global spectrum:** construct both lanes `(21)`--`(22)`, or prove an
   equivalent exact ear/interval realization theorem.
5. **Sequence semantics:** deduplicate the image `(19)` before taking support
   statistics.

This agrees with the typed lessons of the related canon without transferring
their target theorems. THM-2290 keeps the complete selector-indexed matching
kernel and contraction; THM-3456 restores a named branch by an inverse-
boundary word; THM-796 retains path, line, and node sorts plus incidence.
Their common grammar is controlled forgetting, not object identification.

## 10. Exact referee and the legacy count boundary

The companion uses literal permutation enumeration, not a fitted recurrence.
It checks:

1. the two minimal hostiles, including every path;
2. all `5,404` pairs `(T,v)` over every labeled tournament of orders two
   through five, with zero failures of `(10)`;
3. the actual `H=9` support/multiplicity hostile;
4. product closure of all `303` prefix values through target `5,000`, locating
   `613` as the first allowed miss;
5. the arithmetic factor decomposition used in Theorem 6.1's sufficiency
   direction for all `499,998` allowed odd values through one million,
   conditional on the stated lane and carry atoms, with zero failures; and
6. `C_3`-join witnesses for finite subsets of the compactness type.

It separates two counts that the legacy prose of
`THM-012b-insertion-decomposition` conflated before the repair recorded as
MISTAKE-499. At order five:

```text
pairs with at least one orphan:                         3200,
pairs where sum TypeII != O_v:                          3080,
orphan-bearing pairs where the two totals balance:       960.         (34)
```

Thus `3080/5120` is the failure count for the obsolete paper equality, not
the count of pairs carrying orphans. The corrected deficit `(10)` has zero
failures. This bookkeeping correction is independent of the unconditional
proof of Theorem 2.1.

Reproduce from the repository root with

```bash
python3 -B 04-computation/hamiltonian_matching_deficit_two_prime_lane_thm4094.py
python3 -B -O 04-computation/hamiltonian_matching_deficit_two_prime_lane_thm4094.py
python3 -B 04-computation/hamiltonian_matching_deficit_two_prime_lane_thm4094.py \
  | cmp - 05-knowledge/results/hamiltonian_matching_deficit_two_prime_lane_thm4094.out
python3 -B -O 04-computation/hamiltonian_matching_deficit_two_prime_lane_thm4094.py \
  | cmp - 05-knowledge/results/hamiltonian_matching_deficit_two_prime_lane_thm4094.out
```

Both streams byte-match
`05-knowledge/results/hamiltonian_matching_deficit_two_prime_lane_thm4094.out`.
Every executable gate uses `require`; optimization cannot remove a check.
Raw file hashes and the canonical semantic hash are pinned above.

## 11. Scope

THM-4094 does **not** prove or refute `(5)`, and does not assert that `613` or
`623` is absent from the actual spectrum. It proves only that they are not
forced by the stated prefix and multiplication, and that completeness reduces
exactly to `(21)`--`(22)`.

Compactness does not prohibit proving every standard instance of `(31)`, a
constructive realization algorithm, induction in ordinary mathematics,
finite-model reasoning, or stronger categorical encoding. The matching
hostile does not apply when the complete target fiber and incidence remain.
The construction class `C_B` is not a model of all arbitrary insertions. The
finite referee verifies only its stated universes; it is not tail evidence.

No LRC(14), Jacobian, factorial, cellular-automaton, hafnian, or sequence-tail
claim follows. **QED.**
