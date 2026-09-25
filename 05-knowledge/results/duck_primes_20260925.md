# A literal divisor-word model for (10^k-7)/3, and the prime-prefix correction at 196

**Status:** PROVED elementary prefix exclusions and the orbit-count construction; INHERITED PROVED divisor classification and prime-bracket classification; FINITE-EXACT checks in the accompanying script. No primality theorem for the orbit counts, preservation of arithmetic graph labels, or Collatz convergence claim follows. No novelty claim is made.

The constructive result is an exact interpretation of the user's expression: for every balanced mixed number `M=p^2 q r`, its ten proper divisors carry a specified order-three action. Delete seven constant words from the length-k word set, then quotient by the action. The number of remaining orbits is exactly

\[
a_k=\frac{10^k-7}{3}=1,31,331,3331,\ldots\quad(k\ge1).
\]

This explains how the divisor counts `10=7+3` can generate that expression. It is a statement about word orbits, with an explicit loss of numerical divisor relations.

The coordinating [duck decoder synthesis](duck_decoder_20260925.md) gives the compatible model by unordered pairs with repetition from `F_2^2`, the natural `K_5` edge realization, and the pair-XOR charge. That charge is an extra specified coordinate of the model; the proof below needs only the group action.

## Inheritance and live board

Closest mechanisms: [THM-2422, operation fibres and twin-center ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md), whose swap quotient retains repeated-pair diagonals; [arithmetic braids, divisors DB1--DB3](arithmetic_braids_20260917_divisors.md), which proves the balance classification; and [seam_prime_balance_20260925](seam_prime_balance_20260925.md), which supplies its exact finite profile controller and support-label obstruction. The prime-removal set is recovered from [decoder_prime_square_20260925, section 2](decoder_prime_square_20260925.md).

Canonical hostile: `4=2*2` is composite despite having no proper factorization with distinct factors. Corrected near miss: `196=14^2` is not the sum of the first twelve primes. Least-used sidecar: the distinction between a divisor value, a count of divisors, and a letter in a quotient alphabet.

The live board is prime cutoffs, filtered prefix length, prime exponent profiles, repeated-pair diagonals, finite group actions, and arithmetic edge labels. The orbit model gives a positive connection between the middle two combinatorial objects. The counterexamples below locate the point at which it stops preserving the arithmetic graph objects.

## 1. Removing primes does not produce a prefix sum of 196

Let `p_i` be the increasing sequence of primes. There are two different operations:

\[
P_k(E)=\sum_{1\le i\le k,\ p_i\notin E}p_i,
\qquad
R_k(E)=\text{sum of the first k primes retained after deleting E}.
\]

Here `k` is an original cutoff in the first expression and a retained count in the second.

| Rule | k=11 | k=12 | k=13 |
|---|---:|---:|---:|
| Unfiltered first k primes | 160 | 197 | 238 |
| Original first k, then delete 11 | 149 | 186 | 227 |
| First k retained, deleting 11 | 186 | 227 | 270 |
| Original first k, then delete 2,3,11 | 144 | 181 | 222 |
| First k retained, deleting 2,3,11 | 265 | 312 | 365 |

**PROVED for every prefix length:** none of the three increasing prime lists—unfiltered, deleting 11, or deleting `{2,3,11}`—has a prefix sum 196. Their consecutive retained-prefix sums straddling 196 are respectively

\[
160<196<197,\qquad186<196<227,\qquad181<196<222.
\]

Strict increase proves the all-prefix statement; it is not an extrapolation from a finite search. An original-cutoff prefix after deletion is also a retained prefix, so the exclusion covers that convention too.

There are nearby genuine square identities:

* The first nine unfiltered primes sum to `100=10^2`.
* Deleting `{2,3,11}`, the first three retained primes are `5,7,13`, summing to `25=5^2`.
* With that same deletion, the first eight retained primes end at 31 and sum to `144=12^2`.

The last identity is equivalently the original first eleven primes, with `2+3+11=16` removed. It gives 144, not 196.

### A finite-filter square theorem

Once the cutoff includes all deleted primes, deletion subtracts the constant `c=sum(E)`. If the unfiltered and filtered sums are both squares, say `u^2` and `v^2`, then

\[
(u-v)(u+v)=c.
\]

Thus a finite list of factor pairs of `c`, with equal parity, gives every possible simultaneous square pair. For `c=11`, the only pair is `36-25`; 36 is not a prime-prefix sum because it lies between 28 and 41. For `c=16`, the only nonnegative square pairs are `25-9` and `16-0`, both below the sum 28 already attained when the cutoff reaches prime 11. Consequently neither deletion can turn a square prime-prefix sum into another square after all its excluded primes are present. This is an exact bridge from a deletion problem to a difference-of-squares factor fibre.

The role of `{2,3,11}` in the inherited bracket theorem is separate. For

`B_m=((2m-1)^2,(2m+1)^2]`,

the integers having a nontrivial multiple in their own bracket are exactly `2,3,4,10,11,12`; restricting to primes gives `{2,3,11}`. For `m>=3`, even twice the smallest bracket integer exceeds the upper boundary. The strict-upper-bound convention removes the pair `3*3=9` but does not change the prime set. This scale property does not privilege the sum of the remaining primes.

## 2. The three meanings of 10 and the repeated-factor coordinate

For an integer `M>=2`, let `D(M)` be its proper nontrivial divisors, `S_set(M)` the squarefree members, and `U_set(M)` the prime members. Write their cardinalities as `F,S,U`. Here `U` is a divisor count, not a Collatz operation. The inherited classification is

\[
F=S+U\quad\Longleftrightarrow\quad M=p,\ p^3,\text{ or }p^2qr,
\]

where `p,q,r` are distinct primes in the last shape. The count triples `(F,S,U)` are `(0,0,0)`, `(2,1,1)`, and `(10,7,3)`.

For completeness, if `M=prod p_i^{e_i}` has `r` distinct prime factors, then

`F=prod(e_i+1)-2`, `S=2^r-1-[M squarefree]`, and `U=r-[M prime]`.

Squarefree composites fail because `F=S` and `U>0`. In the remaining cases balance becomes `prod(e_i+1)=2^r+r+1`. At support sizes 1,2,3 this yields respectively exponent 3, no solution, and profile `(2,1,1)`. At `r>=4`, the nonsquarefree minimum `3*2^(r-1)` already exceeds the right side.

Three assertions involving the integer 10 must now be separated:

1. **A value:** `10=2*5` has `(F(10),S(10),U(10))=(2,2,2)` and is not balanced.
2. **An additive decomposition:** `10=3+7` is its unique unordered decomposition into two distinct primes. If repeated primes are allowed, `5+5` is another.
3. **A cardinality:** any `p^2qr` has ten proper nontrivial divisors, seven squarefree and three prime. For example `M=60`, whose actual prime factors are `2,3,5`, not the displayed counts 3 and 7.

Since every prime divisor is squarefree, `U_set` is a subset of `S_set`; the equation is not a disjoint partition `D=S_set disjoint-union U_set`. The correct disjoint complement is

\[
E=D\setminus S_{set}=\{p^2,p^2q,p^2r\},
\]

and it has a canonical labelled bijection with `U_set={p,q,r}`:

\[
p\mapsto p^2,\qquad q\mapsto p^2q,\qquad r\mapsto p^2r.
\]

This creates the disjoint labelled model used next.

By contrast `196=2^2*7^2` has `(F,S,U)=(7,3,2)` and is not balanced. Its total prime multiplicity is four, as for `p^2qr`, but its exponent profile `(2,2)` is different. Total degree alone loses the balance predicate.

## 3. A genuine order-three action on the divisor alphabet

Fix `M=p^2qr` with three distinct primes, with `p` the uniquely squared one. The squarefree divisors are the seven nonempty subsets of the three-element prime support `P={p,q,r}`. The other three divisors are the labelled copy of `P` under the preceding bijection. Thus

\[
D\simeq\bigl(\mathcal P(P)\setminus\{\varnothing\}\bigr)\sqcup P.
\]

Let the cyclic group of order three rotate the prime labels, on both parts. This group can be specified without a choice of generator: it is the unique order-three subgroup of the permutation group of `P`. Reversing the cyclic order chooses the inverse generator and leaves all orbits unchanged. The displayed chart specifies how the same permutation is transported to the nonsquarefree part.

The action has one fixed letter `pqr` and three cycles of length three: singleton subsets, two-element subsets, and the extra prime-label copy. For `M=60=2^2*3*5`, a generator is

\[
(2\ 3\ 5)(6\ 15\ 10)(4\ 12\ 20)(30).
\]

This is a combinatorial permutation. It is **not** literal prime substitution in the original number `M`: cyclically substituting primes would move the squared exponent to a different prime. Nor is it an automorphism of the divisor lattice. For example `4` divides `12`, but their images `12` and `20` are not divisible in that direction.

## 4. The decimal family counts free word orbits

For `k>=1`, form all words `D^k` and apply the same group element to every letter. Delete the seven constant words

\[
\{(s,s,\ldots,s):s\in S_{set}\}.
\]

This deleted set is invariant. It consists of **one fixed word** `(pqr,...,pqr)` and **six other words in two free orbits**. It is not the set of seven fixed words.

Every nonidentity element of the order-three group fixes exactly the letter `pqr`. Therefore a word is fixed only if all its letters are `pqr`, and that word was deleted. Every remaining orbit has exactly three words. Hence

\[
a_k=\frac{|D|^k-|S_{set}|}{3}=\frac{10^k-7}{3}.
\]

Equivalently, Burnside's formula gives `(10^k+2)/3` orbits on all words; the deletion removes exactly three orbits. This is a **global colour/letter action**, not the usual necklace action rotating the k word positions. At length two, the latter gives 55 necklaces on ten letters, whereas the full colour action gives 34 orbits and the deleted version gives 31.

The first orbit counts are

`a_1=1, a_2=31, a_3=331, ..., a_8=33333331, a_9=333333331`.

The count is meaningful even at `k=1`: the remaining three nonsquarefree letters form one orbit. The formula has no length-zero interpretation with seven distinct deleted constant words.

### An explicit recursion, including its new boundary states

Every surviving length-k orbit has exactly ten surviving length-`k+1` children obtained by appending a letter. Its stabilizer is trivial, so those ten choices remain distinct in the quotient. New orbits also appear whose prefix was a deleted constant word: choose its letter in `S_set` in seven ways, and append a different letter in nine ways. These 63 words form 21 free orbits. Thus

\[
a_{k+1}=10a_k+21\qquad(k\ge1),\qquad a_1=1.
\]

This gives a concrete construction of the family rather than merely checking congruence modulo 3. The word process has ten extensions per existing orbit plus 21 new boundary orbits; the order-three quotient should not be mistaken for three-way branching in word length.

The other nonprime balanced profile also gives a familiar word count: for `M=p^3`, the two proper divisor letters `{p,p^2}`, after removing the single squarefree constant word `(p,...,p)`, leave `2^k-1` words. Thus the two nonprime balance profiles supply binary repunit cardinalities and the decimal `333...31` cardinalities under specified constructions. This does not identify the resulting counts' Collatz dynamics.

## 5. What this quotient preserves and what it loses

Source: the named prime support, the proper-divisor alphabet, a word, and the squarefree-constant exclusion. Target: its orbit under the specified simultaneous prime-label action. Preserved: length, the alphabet orbit channel of every position, the seven-word deleted family, and equivalence under a global cyclic label permutation. The construction is valid for every `p^2qr`, so 60 and 84 produce the same counts despite having different prime supports.

Lost: actual divisibility, square-sum adjacency, and graceful edge differences. For the literal `M=60` action,

* `4|12` becomes `12 not dividing 20`.
* The square-sum pair `3+6=9` becomes `5+15=20`, which is not a square.
* Rank the ten divisors increasingly by labels `1,...,10`. The star centred at label 1 is gracefully labelled, with edge differences `1,...,9`. The transported divisor permutation sends its centre to label 2; the edges to labels 1 and 3 both have difference 1, so the new labelling is not graceful.

These are direct witnesses against treating the word action as an arithmetic graph automorphism. A graph decoder would need additional edge data and would have to prove its transport independently.

Unique prime factorization also does not mean a unique grouping into composite factors. The number 60 has the unique prime factorization `2^2*3*5` but five unordered proper two-factor groupings:

`(2,30), (3,20), (4,15), (5,12), (6,10)`.

The prime-factor multiset, the grouping, and the prohibition on repeated groups are different data. THM-2422 explains their common swap-quotient mechanism and why the repeated-factor diagonal must be retained.

Finally, the orbit-count value need not itself satisfy `F=S+U`. Values `a_2,...,a_8` are prime and therefore balanced with `(0,0,0)`. But

\[
a_9=333333331=17\cdot19607843
\]

has `(F,S,U)=(2,2,2)` and fails balance. The combinatorial construction still works exactly at that index. The [preceding digit audit](ternary_digits_20260925.md) independently proves the phase `17|a_k iff k=9 mod16` and the two-step Collatz descent of this decimal family for `k>=4`; neither property is assumed in the word model.

## 6. Reproduction and stopping boundary

```text
python 04-computation/experiments/duck_primes_20260925.py
python -O 04-computation/experiments/duck_primes_20260925.py
```

The standard-library script writes [duck_primes_20260925.out](duck_primes_20260925.out). Explicit checks survive `-O`. The universe includes the first thirty retained prime prefixes for each specified filter; exact monotonicity witnesses around 196; direct divisor sets at the listed values; literal word partitions through length five for 60 and length three for 84; independent Burnside counts; ten-child/21-new-orbit checks; and the arithmetic and graceful-label hostiles. The formula for arbitrary length and the exclusion of 196 for every prefix are established by the proofs above.

The positive result is the explicit action and its recursive orbit enumeration. Its first unsupported extension would be to infer primality of an orbit count, arithmetic graph incidence, or a terminating number-theoretic rank from that enumeration. Any further decoder must state and retain those extra predicates; no such universal extension is asserted here.
