# Prime-lane asymmetry and the sharp first symmetric `7p` count

**PROVED elementary consequences of the cited Busch strong minimum, plus
FINITE-EXACT tables of 5- and 7-vertex circulants and 192 marked substitutions.**
The global Hamiltonian-path spectrum conjecture remains **OPEN**. No external
novelty or priority claim is made. Source repository pin:
`91cc1210ed2bc94ca52686b8091c0eec2c57a441`.

Let `H(T)` count directed Hamiltonian vertex words, and let `Aut(T)` mean
vertex automorphisms. The new structural statements are:

1. If `H(T)=p` is prime and `p>3`, then `Aut(T)` is trivial.
2. If `H(T)=7p` with prime `p>=11`, then `|Aut(T)|` is either 1 or 7.
3. For every odd prime `p<593`, any tournament with `H(T)=7p` is
   asymmetric. (For `p=3`, no such tournament exists.) The threshold is sharp:
   a strong tournament on ten vertices has `H=4151=7*593` and automorphism
   group of order 7.

These statements concern **every** witness, not merely a specially chosen or
minimum-order witness. They narrow the two inherited existence lanes; they
neither construct all lane values nor prove a missing value exists.

## 1. Inherited frontier and source scope

The source pins in `sources.json` identify the exact blobs used.

* `rotation_owners_spectrum_20261020e/spectrum.md`, Sections 3--5, records
  SCC multiplication, the exact Busch floor, and the already proved
  two-prime-lane reduction. Its public primary source is A. H. Busch,
  *A Note on the Number of Hamiltonian Paths in Strong Tournaments*,
  Electronic Journal of Combinatorics 13 (2006), N3, Corollary 1:
  <https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf>.
* THM-4094 already proves that the spectrum conjecture, using carry witnesses
  49, 63 and 343, is equivalent to strong realization of every odd prime
  `p!=7` and every `7p` for odd prime `p!=3`. This is inherited, not new here.
* THM-4861(c) proves that a **vertex-transitive** tournament with prime H is
  only `C3`. Statement 1 above strengthens vertex-transitive exclusion to
  exclusion of **every nontrivial automorphism** for prime H greater than 3.
* The bounded duplication audit inspected THM-4861 and the named core,
  descent, spectrum-structure, local-closure, constructor and actual-context
  proof notes from Sessions 133--138. None stated the full asymmetry theorem
  or the sharp `7p` threshold. This is not a global literature search.

The only non-elementary all-order input needed below is the cited exact
strong floor. The marked-subtournament and path-product facts are proved
here, avoiding a dependence on an unspecified pancyclicity or OCF theorem.

## 2. Automorphisms, SCCs, and strong order bounds

### Free action on paths

An automorphism fixing a Hamiltonian vertex word fixes every vertex at its
position in that word. Therefore the action on Hamiltonian paths is free and

```
|Aut(T)| divides H(T).                                      (1)
```

### Arbitrarily many singleton SCCs do not create symmetry

The strong components of a tournament are totally ordered by dominance.
An automorphism preserves strong connectivity and the order of that chain,
hence fixes **each component setwise**. In particular it fixes every singleton
SCC vertex individually. On distinct components its restrictions are
independent, giving the direct product of their automorphism groups.
Hamiltonian paths concatenate component paths, so their counts multiply.

If `H(T)=p` is prime, exactly one SCC is nontrivial and has count p. If
`H(T)=7p`, any product with two nontrivial factors would have to include a
factor 7, forbidden by the inherited all-order omission. Thus here too there
is one nontrivial strong core C, of count `7p`. Padding singleton SCCs are
pointwise fixed, and `Aut(T)=Aut(C)`.

### Cited exact floor and elementary consequences

For strong order `n>=3`,

```
f(3k)=3*5^(k-1), f(3k+1)=5^k, f(3k+2)=9*5^(k-1), k>=1.    (2)
```

In particular `f(n+3)=5f(n)` and f is increasing: the four successive values
at orders `3k,3k+1,3k+2,3k+3` are `3,5,9,15` times `5^(k-1)`.
Exact base checks give

```
f(4)=5>4, f(5)=9>5, f(6)=15>6,
f(11)=225>77, f(12)=375>84, f(13)=625>91.
```

The recurrence propagates

```
f(n)>n   for every n>=4,
f(n)>7n  for every n>=11,                                (3)
```

because `5n>n+3` in either range.

For prime `H(C)=p>3`, the core is not a triangle and has order r with
`p>=f(r)>r`. By (1), a nontrivial automorphism group would have order p.
Cauchy's theorem would give an element of order p acting nontrivially on
fewer than p vertices, impossible. This proves Statement 1.

For `H(C)=7p`, `p>=11`, if `r>=p`, monotonicity and (3) would give
`7p>=f(r)>=f(p)>7p`. Thus `r<p`. By (1), the group order divides `7p`, and
it cannot have p as a prime divisor, again by Cauchy's theorem. Its order
is therefore 1 or 7. This proves Statement 2.

A direct orbit variant is useful: an automorphism of prime order q has a
q-vertex orbit whose induced tournament is invariant under a full q-cycle.
For odd q this is a circulant tournament, hence regular and strong. Its
count is at least `f(q)`; induced monotonicity gives `H(T)>=f(q)`. This also
excludes the p factor above without mentioning singleton padding.

## 3. Two elementary tools for the sharp threshold

### Disjoint induced path-product inequality

For disjoint vertex sets A and B,

```
H(T[A union B]) >= H(T[A]) H(T[B]).                       (4)
```

Choose a directed Hamiltonian path in each set. Merge the two path words as
follows: while both are nonempty, output whichever current head points to the
other head and remove it; append the remainder once one word is empty.
Successive outputs form a directed path. If the next output is from the same
old path, their edge was already directed; if from the other, the selection
rule supplied their edge. The output preserves each old word as a subsequence.
Restricting it to A and B recovers both inputs, so this deterministic merge
is injective. This proves (4). Taking one part to be a deleted vertex also
proves the monotonicity used here. Iteration proves the analogue for any
finite disjoint partition.

### Marked strong four-vertex lemma

Every strong tournament Q on at least four vertices, with a specified
vertex v, contains an induced strong four-vertex tournament containing v.

First v lies in a directed triangle C: write its out- and in-neighbor sets
as O and I. Both are nonempty, and some arc goes from O to I, or O would
have no directed exit to `I union {v}`. Such an arc completes a triangle
through v.

If a vertex outside C has both an in- and an out-neighbor in C, adjoining it
to C is strong and finishes. Otherwise every outside vertex is either in A,
dominated by C, or in B, dominating C. Strongness makes A and B nonempty and
forces an edge `a->b` with `a in A, b in B`; otherwise A could not exit to
`C union B`. Choose c in C with `v->c`. Then

```
v -> c -> a -> b -> v
```

is a directed spanning cycle on four vertices. Its induced tournament is
strong, proving the lemma.

## 4. Exact small tables

An orbit of a prime-order cyclic automorphism may be labelled by Z/q so that
rotation is addition by one. A tournament invariant under this action is
specified by choosing one of d and `q-d` for each `1<=d<q/2`. Thus there are
exactly `2^((q-1)/2)` presentations to check, with no isomorphism-census input.

The standalone checker enumerates all these presentations. It cross-checks
counts by literal permutations at orders 5 and 7:

| orbit order | presentations | H values |
|---:|---:|---|
| 5 | 4 | 15 |
| 7 | 8 | 175 (six presentations), 189 (two presentations) |

If exactly one nontrivial orbit M is present, every exterior vertex is fixed
by the automorphism. Its adjacency to M is constant around the orbit.
Therefore M is a homogeneous module, and contracting it gives a quotient
Q. If the whole tournament is strong, Q is strong (unless Q is the single
module vertex, which is harmless). The converse needed in the exact table
also holds for a strong inner module and strong quotient, by lifting paths.

Replacing marked vertex 0 in each labelled strong quotient by each cyclic
orbit gives the following **complete finite** tables. For quotient order 3
there are two strong labelled quotients; for order 4 there are 24.

| orbit order | quotient order | substitutions | H histogram |
|---:|---:|---:|---|
| 5 | 3 | 8 | `135:8` |
| 7 | 3 | 16 | `2163:12, 2205:4` |
| 7 | 4 | 192 | `4151:72, 4221:24, 5789:72, 5859:24` |

Counts after the colon are the numbers of labelled input presentations,
not isomorphism-class counts. These substitutions have at most ten vertices.
The checker evaluates each actual tournament by forward endpoint subset DP
and a separate top-down remaining-vertex recurrence.

## 5. Sharp `7p` asymmetry theorem

First handle p=5 and p=7. If `H=35`, (2) bounds the strong core order by 7.
A nontrivial automorphism has a prime-order element of order 5 or 7.
An order-7 orbit would induce a tournament with at least 175 paths,
contradicting monotonicity. An order-5 element has one 5-orbit and at most
two fixed vertices. The quotient is either the single module (H=15), has
order 2 (not strong), or is a directed triangle (H=135). None gives 35.
Thus every 35 witness is asymmetric. If `H=49`, any nontrivial automorphism
has an element of order 7, whose orbit alone has H at least 175; impossible.
The case p=3 is the known omission 21.

Now let `11<=p<593` be prime and suppose `H=7p`. Statement 2 says that a
nontrivial automorphism group would have order 7. Let g be a generator.
Each nontrivial g-orbit induces one of the eight 7-circulants. Two such
orbits would, by (4), force

```
H >= 175^2 = 30625 > 4151 > 7p.                          (5)
```

Thus there is exactly one nontrivial orbit M; all other core vertices are
fixed. Contract M to a marked vertex in the strong quotient Q.

* If `|Q|=1`, H is 175 or 189, with respective quotients `H/7=25,27`,
  not prime.
* `|Q|=2` is impossible for a strong tournament.
* If `|Q|=3`, the table gives H=2163 or 2205, with respective quotients
  `309=3*103` and `315=3^2*5*7`, not prime.
* If `|Q|>=4`, the marked strong-four lemma supplies an induced strong
  four-vertex subquotient containing the marked module. Lifting it keeps
  **the entire actual M**, along with three fixed exterior vertices, as an
  induced subtournament. The complete 192-entry table bounds its H by
  4151 from below; monotonicity then gives `H>=4151`, contradicting `p<593`.

Every case contradicts nontrivial symmetry. This proves Statement 3 below
its claimed threshold.

### Actual threshold witness

Use bit rows

```
[142,156,184,240,225,195,135,512,255,383].
```

Bit j in row i is 1 exactly when `i->j`. The first seven vertices form the
circulant with forward steps `{1,2,3}`. They replace marked vertex 0 in the
strong four-vertex quotient of arc code 17 (lexicographic pairs,
least-significant bit first, bit 1 means `i->j`).

Direct checking gives strongness and `H=4151=7*593`, with 593 prime. The
permutation `(0 1 2 3 4 5 6)`, fixing vertices 7,8,9, is an automorphism.
Its seven powers are distinct. By (1), the entire group has order dividing
`7*593`; since the graph has only ten vertices, Cauchy excludes factor 593.
Thus its automorphism group has **exactly** order 7. This proves sharpness.

## 6. What this does and does not close

SCC multiplication already reduces any additional least missing allowed odd
value, beyond the known carry atoms, to `p` or `7p`. The argument here adds:

* A prime-lane witness for `p>3` must be completely asymmetric. Searching only
  vertex-transitive, circulant, or other nontrivially symmetric families can
  never settle that lane.
* A `7p` witness at `p>=11` permits at most the precise group orders 1 and 7.
  The symmetric branch really occurs, first at prime quotient 593.
* At a prime-order automorphism with t nontrivial q-orbits, (4) gives
  `H>=f(q)^t`; for q=7 the exact orbit table improves this to `H>=175^t`.

The strong core is of order at most `max{r:f(r)<=h}` in either lane, but that
finite membership cutoff does not assert that a witness exists. Likewise,
an automorphism restriction is not a constructor or a surjectivity theorem.
The still-missing all-order step is an **actual witness-preserving inverse
construction or interval overlap theorem** covering both infinite lanes.
There is no deduction here from modular reachability, a formal state space,
positive Jacobian determinant, local rank, or a finite witness bank to global
integer-tournament realization.

## 7. Replay

```
python3 -B verify.py > verification.json
python3 -B -O verify.py > verification.optimized.json
cmp verification.json verification.optimized.json
```

The checker uses only the Python standard library. It also independently
counts every labelled tournament through order five by literal permutations,
checks the prime asymmetry statement there by enumerating all permutations,
and tests the marked strong-four lemma on every strong labelled tournament
of orders four and five. The all-order parts are the proofs above; the
finite program does not pretend to test infinitely many cases.
