# Session284: concise rules for entire structural families

**Chan–Pak remains OPEN.** This session answers the structural question with
all-order grammars and exact birth rules, not by checking another isolated
large size. It does not claim a practical explicit list of all classes or
transitions at orders50 or200.

## The main distinction

There already is a complete concise rule generating every tournament layer:
attach one vertex with an arbitrary binary incidence word and identify the
resulting isomorphism classes. It works beyond every tested size. What is
missing is not permission to extrapolate that rule, but a manageable
classification of its enormously many inequivalent outcomes and the
invariants that predict them. [Architecture](architecture/PROOF.md) gives the
exact marked/unmarked orbit identity and illustrates the distinction:
**finite rule syntax is not a bounded state space**. Even one explicit rigid
prime family has exponentially many distinct prime children.

The response to that difficulty is not to give up on structure. Here are
three classes of exact all-order answers, together with their interactions.

## 1. A lossless grammar for every marked growth step

The [signature-run grammar](extension/PROOF.md) handles any old tournament
compressed into its canonical transitive chains and any batch of newborns.
Each old vertex has a binary signature toward the newborn set. Constant
signature runs become variable-length modules; a small skeleton retains
the old core, newborn arcs and run order. Positive run lengths restore the
entire child. The theorem identifies exactly which description symmetries
represent the same marked extension.

For **every fixed template**, one rational generating function counts the
whole infinite family by size and gives its exact first appearance. Child
core rank is determined by the skeleton, independent of positive lengths.
This is not just a way to count: it reconstructs every member and every
marked transition. Iterating with one newborn from the singleton is a
complete generation language. Forgetting construction ancestry requires
the explicit child-automorphism correction or canonical deduplication.

The same language retains Hamiltonian semantics through the inherited
walk-content and path-cover substitution law. Keeping only core rank or a
scalar path count would lose that information.

### A new pattern with a known birthday at every level

A k-vertex boundary whose exterior witnesses realize all2^k orientation
signatures first becomes possible at exactly **k+2^k vertices**. Hence new
external-interface dimensions appear at3,6,11,20,37,70,135,264,521,... .
This is an all-order formula, not an extrapolated list. A Gray-code witness
realizes the pattern using single-coordinate changes along one old chain.
It differs from ordinary neighborhood VC dimension, whose witnesses may
lie inside the boundary and whose earlier repository threshold is2^k.

An exact random-growth rule also emerges. For any parent of size n and
chain-core rank b, its old homogeneous edges survive k independent uniform
newborns as Binomial(n−b,2^−k). Explicit small-error formulas give child rank
and a logarithmic fracture threshold. The law is for labelled independent
arcs, not uniform isomorphism classes. Fracturing the chain compression does
not mean the remaining pair-free core has been classified.

## 2. Complete prime inventories for a whole unbounded family

The [feedback-word grammar](births/PROOF.md) completely describes all
apex-marked tournaments whose apex deletion is transitive. A binary word
determines every module, its unique ordinal/run decomposition, every induced
prime type and its multiplicity.

The number r of adjacent10 transitions is an exact structural birth clock:
when r≥1 the prime quotient is W_(2r+1); r=0 gives a transitive object. New prime arities arise indefinitely, but are
predicted by one rule. In the standard tower W_(2k+1), the number of induced
W_(2j+1) copies is binom(k+j,2j). Its entire prime inventory obeys

    P_k(t)=(2+t)P_(k−1)(t)−P_(k−2)(t),
    sum P_k(t)z^k=(1−z)/(1−(2+t)z+z²).

This is the requested small-object intimacy at all sizes: the same tiny
transfer predicts the exact module and prime-content structure rather than
merely a total number. It also proves that no fixed finite *bounded-arity*
substitution dictionary can be universal. A finite *parametric* grammar is
not ruled out; W itself is the example of one.

This class already has an inherited global omission of path count13. Thus
even complete structural understanding of an infinite class need not solve
Chan–Pak.

## 3. An exact recursive host theory and all rooted symmetry types

The [ternary-host grammar](wildcard/PROOF.md) classifies the whole hereditary
class built from ordinal sums and cyclic triples. It proves that this class
is exactly the induced substructures of balanced C3 towers. Every object has
an unambiguous recursive isomorphism form.

More importantly, **the first tower in which any such object appears is
computed exactly**. Cyclic nodes use1+max of the child costs; ordered strong
components use an alphabetic minimax binary-tree recurrence. An exact greedy
dyadic allocator computes the same birth depth directly from their ordered
child costs, without exploring binary parenthesizations. The actual
order matters: C3→1→1 first embeds in the9-vertex tower, while1→C3→1 first
embeds in the27-vertex tower, despite identical component-cost multisets.
This is a rigorous bridge to ordered coding, not a numerical analogy.

A symmetry-refined cyclic-necklace recurrence classifies every rooted
one-vertex extension of every tower, including its root-fixing automorphism
order. The q-color version handles any named newborn interface. Thus the
family's next layer is not mysterious at large size; the recursive sections
are the sufficient state.

The same tree gives exact pointwise endpoint stabilizers. Arbitrary
attachments may destroy visible automorphisms while still forcing powers of
3 into every Hamiltonian count. The all-level bound extends to cyclic
prime-arity towers using an explicit subgroup. Its sharpness concerns the
pointwise group floor, not an asserted exact path-count gcd in all cases.

## 4. Comparing the grammars rather than merely collecting them

The [intersection theorem](BRIDGE.md) says an apex-word tournament belongs
to the ordinal/C3 class **iff its word has at most one10 transition**.
There are exactly(m+1)+binom(m+1,3) such marked objects with a spine of length
m. A W5 quotient is the first escape, and the same criterion decides every
later size. This is a complete relationship between two infinite families.

Larger deletion interfaces are encoded by words over2^k signatures. The
k=1 inventory is now explicit; the general batch grammar is lossless, but
uniform prime inventories and unmarked quotient simplifications for all k
remain a real research frontier. This is a concrete next parameter, not an
unbounded checklist of unexplained exceptions.

## What is and is not completed

Completed here: exact rules for all marked growth; rational counts of every
fixed signature template; complete module/prime inventories in the one-apex
class; complete ordinal/C3 family isomorphism and host-birth theory; all
rooted ternary-tower symmetry types; and the exact family intersection.

Not completed: a uniformly small classification of every unbounded prime
core; a practical explicitly enumerated extension metagraph at50 or200; or
an inverse selecting a tournament for every allowed odd Hamilton count.
A000568 totals and exact correction certificates do not supply those
objects automatically. No impossibility theorem against future compact
structural theory is claimed.

[Reproduction](REPRODUCE.md) · [Independent review](audit/REVIEW.md).
