# First birth by strong factorization, with unique class generation and converse counts

**PROVED all-order finite reduction.** This does not prove that every allowed odd target has a birth. It gives an exact finite algorithm and identifies the precise strong-data obligations.

## Inheritance and conventions

The ordered SCC word and its converse reversal are classical. The repository already explicitly proves their all-order converse-palindrome use in `05-knowledge/results/tournament_converse_context_20261011d.md`, §§2–3. Its ordered/Lyndon distinction is important: do not replace a sequence of SCC types by an unordered multiset. `metagraph_classification_20261021a/README.md` supplies another exact SCC enumeration, but for the distinct statistic H/|Aut|; its numerical classification is not imported for H. The current finite-core bound is `pubsub/publisher/PROOF.md`.

**Additional recovered priority:** `birth-atlas-now/PROOF.md` and `birth_program/README.md` already state the minimum-over-strong-factorizations theorem and a stronger logarithmic finite-core bound using Busch. The recursions below make the multiplicity/converse bookkeeping explicit; `solve.py` implements the exact decision procedure with the Busch cap rather than the weaker elementary cap. No novelty is claimed for the earlier first-birth factorization or decidability statement.

We count ordinary tournament isomorphism classes, not labelled tournaments, Hamilton-word orbits, or converse-merged classes. For h>1 a first creator has no singleton SCC: deleting such a singleton preserves H and lowers its order. Every nontrivial strong component has order at most its own H, by cutting a Hamilton cycle. The whole singleton-free core has order at most h. Hence any attained h>1 first appears by order h, and existence/nonexistence is decidable by finite enumeration through h.

Use an *empty-core* convention at h=1: b(1)=0 and a(1)=1. For actual nonempty tournaments, the first appearance of H=1 is the singleton at order1; the convention is only an algebraic identity for factorizations.

For d>1 let S_d(z) be the polynomial counting strong isomorphism classes with H=d by vertex number. It is finite, of degree at most d. Let s(d) be its least degree (infinity if zero), v(d) the coefficient at that degree, and c(d) the number of those minimum-order strong classes which are self-converse.

## 1. Exact first-birth recursion

Let B_h(z) count all singleton-free tournament classes with H=h by order, with B_1(z)=1. Unique first-SCC decomposition gives

    B_h(z)=sum_(d|h, d>1) S_d(z) B_(h/d)(z).                 (1)

Every factor h/d is smaller than h, so this is a terminating divisor recursion. It is an equality of finite polynomials. Its least-degree and leading-coefficient recursion is

    b(h)=min_(d|h,d>1) [s(d)+b(h/d)],                        (2)
    a(h)=sum_(d|h,d>1, s(d)+b(h/d)=b(h)) v(d)a(h/d).          (3)

If no finite candidate exists, b(h)=infinity and a(h)=0. Otherwise b(h) is the true first tournament order and a(h) its exact isomorphism-class multiplicity. Only strong minimum data are needed, not every coefficient of S_d.

Equivalently, b(h) is the minimum sum of strong birth orders over multiplicative factorizations h=d₁⋯d_k, d_i>1. However **do not use arbitrary binary splits to count classes**: the same SCC word has many parenthesizations. Recurrence(3) is unambiguous because the first SCC is unique.

A component in a first creator must itself be a globally minimum-order realization of its count. Otherwise replacing it with a smaller tournament preserves the product H while lowering total order. In particular every contributing strong component satisfies s(d)=b(d). This restricts the atom catalog without assuming that all first creators are strong.

For a finite h, enumerate strong tournaments only for divisors d of h and only through order d to obtain all required data. When proper-factor constructions give an upper bound U for b(h), the search for a competing strong H=h realization can stop at U: larger strong classes cannot change the birth. This is exact pruning, not a heuristic size cutoff.

## 2. Separate the existential and universal strength questions

Define

    q(h)=min_(d|h,1<d<h) [b(d)+b(h/d)],

with infinity for an empty minimum. Then b(h)=min(s(h),q(h)). For attained h:

- Some first creator is strong iff s(h)=b(h).
- Every first creator is strong iff s(h)<q(h).
- Both strong and decomposable first creators exist iff s(h)=q(h)<infinity.
- Every first creator is decomposable iff q(h)<s(h).

Thus proving s(h)<=q(h) for every attained h would settle the **existence** of a strong first creator, but not the universal class statement. The latter requires strict inequality. These are genuinely different proof obligations, not interchangeable readings of the same conjecture.

For prime h, q(h)=infinity, so every first creator is strong. For composite h the SCC product law alone does not decide the comparison.

## 3. Exact converse multiplicity at first birth

Converse reverses the SCC word and replaces each component by its converse. A fixed word therefore consists of an arbitrary left word, its forced reversed-converse right word, and either no central letter or one self-converse strong central letter.

For d²g=h, the left word has count d and the central letter count g. At first birth the left word must itself be first-born, and the central letter must be minimum-order strong; otherwise the complete palindrome could be shortened. Define s₀(1)=0,c₀(1)=1 for an absent center; for g>1 use s₀(g)=s(g),c₀(g)=c(g). The exact number f(h) of self-converse first-creator classes is

    f(h)=sum_(d²g=h, 2b(d)+s₀(g)=b(h)) a(d)c₀(g).           (4)

There is no division by2 in the left-half count: a chosen SCC prefix determines a unique full palindrome. Distinct prefixes give distinct full SCC words. The converse-pair count is

    (a(h)-f(h))/2.

This is an all-order formula. It does not assert that a fiber consists of one fixed class or one converse pair; arbitrarily many orbit contributions are allowed.

For squarefree h, d=1 is forced in(4). Therefore every self-converse first creator at squarefree h is strong. Any decomposable first creator at squarefree h must lie in a distinct converse pair. For prime h this agrees with the preceding strength statement.

## 4. Unique generation and its boundary

To generate the entire first-birth fiber, recursively choose a divisor d satisfying(2), choose a minimum strong class from S_d, then generate the first-birth suffix at h/d. Join them in ordinal order. Every first creator appears exactly once, because its ordered SCC-type word is unique.

This is not a universal generation theorem for all tournaments from first creators under ordinal joins. A strong tournament that is not first-born remains an additional strong atom for that larger task. Nor is it a proof that the strong atom catalog is irredundant under other operations such as substitution or selective vertex extension.

## 5. Proof-strength consequence for Chan–Pak

The global existence problem reduces exactly to nonemptiness of the strong divisor catalog under recursion(2). Any odd h with a proper factorization into already attained values is automatically attained. Thus unresolved targets can be restricted to those with no such factorization, but a prime target still requires an actual strong witness and cannot be supplied by the SCC product mechanism alone.

This reduction is useful for local-to-global proof design: separate strong prime-like construction obligations from automatic multiplicative closure; retain minimum order and converse class counts as independent data. It does not turn finite strong enumeration into an eventual-threshold proof.

## Replay

The accompanying verifier exhaustively canonicalizes tournaments through order5, recovers their strong birth inputs, and reconstructs every observed first-birth fiber and self-converse multiplicity via(2)–(4). A separate synthetic atom dictionary tests tied strong/decomposable births and the no-parenthesization rule. Synthetic controls illustrate the algebra only, and are not asserted to describe physical tournaments.

## 6. Squarefree symmetry at every later size

The palindrome argument also gives a sharp later-order consequence. For squarefree h>1, a self-converse tournament can have only one nontrivial SCC: every off-center nontrivial SCC would force a squared factor greater than1 in h. All remaining SCCs are singleton pairs placed symmetrically before and after the strong center. If D_h(z) counts all self-converse strong classes of count h by order, then the full self-converse generating function is exactly

    D_h(z)/(1-z²).

Since strong order is at most h, its coefficient is eventually periodic with period2, and bounded independently of total order. In particular, for an attained prime h the ordinary class count is A_h(n+1)-B_h for n>=h, where A_h is the number of all strong H=h classes and B_h is their total order. Its self-converse subclass stays bounded, so the self-converse fraction tends to zero. This is an all-order structural conclusion, not a first-birth census extrapolation. Symmetric padding can create new self-converse sizes while still making symmetry rare within the full growing fiber.
