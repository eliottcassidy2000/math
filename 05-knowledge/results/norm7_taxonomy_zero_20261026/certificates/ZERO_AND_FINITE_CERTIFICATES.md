# Zero, bottom, and what a finite certificate actually decides

**Status: PROVED elementary reductions and certificate lemmas; CITED classical computability results are identified below.** No claim that Collatz, RH or motivic identification is undecidable follows. Fix an acceptable enumeration φ_e of partial computable functions N⇀N, with decidable bounded execution predicate Halt(e,x,s,y). Write bottom for the empty-domain function and Z(x)=0 for the total-zero function.

## 1. Inheritance and the promised carrier

The closest inherited mechanism is exact finite witness extraction: the TT7 data in `HYP-9391-the-seven-vertex-paley-census-and-four-new-weight-five-motives.md` already contain an exact rank computation, but its four fitted residual motives remain unidentified. MISTAKE-671 correctly warns that the observed statistics do not fix rank4. The present certificate separates finite linear independence from a global representation identification.

The canonical hostile is an evaluator that agrees on every tested input and changes at the next untested input. The repaired near miss is treating “all computed outputs are zero” as “the function is total zero.” The least-used sidecar is the **promise**: totality, candidate membership, represented Galois object, rank, or finite determining theorem must be supplied separately.

Board: graph extension order; positive halting witness; candidate membership; exact trace bank; residual representations; quotient losses. Anchor: a finite TT7 independence certificate. Niche: precise computability levels. Wildcard: a noncubic Frobenius bank that fails when the representation rank changes.

## 2. The exact index-set taxonomy

| Property of index e | Normal form | Complexity |
|---|---|---|
| φ_e is bottom | ∀x∀s, no halting computation on x by s | Π⁰₁-complete |
| φ_e(0) diverges | ∀s, no halting on input0 by s | Π⁰₁-complete, a different set |
| φ_e is total | ∀x∃s,y Halt(e,x,s,y) | Π⁰₂-complete |
| φ_e=Z | ∀x∃s Halt(e,x,s,0) | Π⁰₂-complete |
| whenever φ_e is defined, its output is0 | ∀x,s,y [Halt(e,x,s,y)⇒y=0] | Π⁰₁-complete |

The first Π⁰₁ completeness follows by reducing nonhalting on one fixed input: on every input x, wait for the specified computation, then return0. The resulting function is bottom exactly when that computation never halts. For the last row, return1 after that same trigger; “no nonzero output” then has the same reduction. Blank-input divergence has the standard direct reduction but is not empty domain: a program that diverges at0 and returns0 at every positive input is a witness.

For the Π⁰₂ rows, let an arbitrary Π⁰₂ predicate be ∀x∃s R(e,x,s), with R decidable. The program on x searches for such an s and returns0. It is total, equivalently total zero, exactly when the predicate holds. These reductions prove the stated lower bounds and the displayed normal forms give the upper bounds. No hierarchy classification is inferred merely from Rice's theorem.

In the graph-extension order f⊆g, bottom is the least element. Every total function is maximal; no function is greatest, because constant0 and constant1 are incompatible. Conversely a partial function missing x can be extended by assigning a value at x. This remains a partial computable function when f was partial computable. Thus maximal elements are exactly total functions; Z is one maximal element, not a greatest element.

## 3. Rice's witness need not be Z; promises change the problem

Rice's theorem says every nontrivial extensional property of partial computable functions has an undecidable index set. Its usual reduction chooses a function f whose property value differs from bottom. There is no requirement that f=Z. For the property “the value1 occurs,” Z and bottom are both negative; constant1 is the required positive witness. Alex Simpson's author notes, Theorem8.4, explicitly use an arbitrary member of the property after arranging that bottom is excluded. Theorem8.5 gives the finite positive-information condition for semidecidability (Rice–Shapiro): https://alexksimpson.github.io/Teaching/tinotes.pdf . This source is cited for those classical principles; the reductions and application boundaries here are explicit.

Under the promise φ_e∈{Z,bottom}, running φ_e(0) semidecides the Z side: any halt decides it, by the promise. The bottom side is not semidecidable, and no total decision procedure exists, by the trigger reduction above. Without this promise, one halt with output0 decides neither total zero nor even eventual termination at another input. In the unrestricted index space, the complement of {e:φ_e=Z} is Σ⁰₂-complete, not a merely enumerable set of finite nonzero-output witnesses.

For Collatz, encode the positive starting integer as n+1 for input n∈N. A program that iterates the specified map from n+1 until it reaches1 and then returns0 has total-zero behavior exactly if every positive starting integer eventually reaches1. Its mathematical statement has the form ∀n∃t R(n,t) with bounded simulation decidable. This supplies a Π⁰₂ **form**, not Π⁰₂-completeness for that one fixed conjecture, and not undecidability. A nontrivial periodic cycle would be a finite counterexample certificate; absence of such a cycle does not certify absence of unbounded trajectories.

## 4. Total trace streams have a different quantifier profile

Suppose A(e,n) and B(e,n) are uniformly total, computable, exact values in a decidable set. Their equality is

    ∀n [A(e,n)=B(e,n)],

which is Π⁰₁. A disagreement has a finite witness n. This classification is sharp for unrestricted uniformly total families: let A(e,n)=1 if machine e on input e halts within n steps, and0 otherwise; B(e,n)=0. Equality is equivalent to nonhalting. Both evaluators are total, with explicit finite runtimes.

If the two programs are arbitrary and totality is **not** supplied, equality to a total target includes a termination requirement and can be Π⁰₂-complete. A promise of totality must not be erased when counting quantifiers. An infinite numerical Euler-data identity is Π⁰₁-shaped only after its index domain and exact total evaluators are fixed; that statement does not make an arbitrary trace list a Galois representation.

## 5. Finite-candidate separation theorem

**Theorem C1.** Given a finite list f₁,...,f_m of pairwise distinct uniformly total computable functions N→D, where D has decidable equality, there is a computable terminating procedure, under those hypotheses, that finds a finite set W such that restriction to W distinguishes every pair. It returns pair witnesses n_(ij) with f_i(n_(ij))≠f_j(n_(ij)).

**Proof.** Dovetail over inputs and the finitely many pairs. Distinctness guarantees a disagreement for each pair; totality makes every comparison finish. When every pair has a recorded witness, take their union W. Each pair is separated on W by construction. There are finitely many pairs, so termination follows. QED.

If an unknown total function g is **promised to belong to this list**, evaluating g on W determines its unique member. Without membership, the output is only “one listed function remains consistent.” If D contains two computably specified distinct values, then for every finite W and any f_i, the total computable function equal to f_i on W and changed at one outside input has the same certificate but is not f_i: at that input choose the first specified value unless it equals f_i there, and the second otherwise. The singleton-codomain case has no such hostile. A list of all candidates, together with a proof of target membership or a theorem giving a determining set for the whole carrier, is the missing step.

This is not a solution of arbitrary function equality: pairwise distinctness and target membership are promises, not decidable preprocessing supplied by the theorem. If duplicate functions occur, the pair-search can run forever; a supplied equivalence classification would permit taking a quotient first.

## 6. A compact TT7 application: thirteen exact rows

Let the columns be the exact functions

    1,q,q²,q³,a_q,q a_q,c_q,d_q,|b_q|²,Y02(q),Y05(q),Y09(q),Y18(q)

on primes q≡3 mod4, with the labels and raw Legendre-sum definitions retained from HYP-9391. At the thirteen primes

    7,11,19,23,31,43,47,59,67,71,79,83,103

the evaluation matrix has determinant

    295608050652838265614106593306487371138496200704 ≠ 0.

The full integer matrix is in `tt7_rank13_input.json`. The checker reconstructs each raw Y by direct quadratic-character summation, freshly compares all modular columns with PARI, independently checks a_q and d_q by eta products and c_q by the inherited proved twisted fourth-moment identity, and computes the determinant by fraction-free elimination plus modular elimination. A duplicated-column hostile has determinant0.

**Consequence.** No nonzero rational linear combination of these thirteen exact functions can vanish on every prime of their domain. In particular the four raw Y classes are linearly independent modulo the span of the nine inherited functions. Indeed any hypothetical global linear relation restricts to a null vector of the invertible displayed matrix. This is a global negative conclusion certified by finitely many evaluations. Thirteen rows are cardinality-minimal for a full-rank certificate on thirteen columns; no claim of the smallest possible maximum prime is made.

This compresses the inherited rank observation into a small replayable witness; it is not a new claim of four pure irreducible motives. It does not identify the exact known-part subtraction, motivic ranks, conductors, Hodge numbers or semisimplifications. Four independent functions can be traces of reducible objects, virtual combinations, or objects of different ranks.

For the much weaker candidate-list problem consisting only of the four labelled raw Y functions, q=19 already gives four distinct values (2304,−1280,1536,1280). That separates those four names under membership; it identifies no external candidate motive. A late-perturbation hostile agrees with Y02 at all thirteen certificate primes and differs by1 at107.

## 7. What the certificate adds, and what it leaves open

Source→target: exact evaluators→finite evaluation matrix→independence obstruction. Preserved predicate: absence of a rational linear relation. Lost information: motivic rank, purity, extensions, lattice and geometry. Needed sidecar: the definitions of all thirteen columns and exact evaluation receipts. Cheapest decisive test: the determinant, including a zero-determinant duplicate control.

The representation comparison gates and explicit losses are in `GALOIS_COMPARISON_GATES.md`. The two lanes share a mechanism: a finite test proves a universal claim only when its determining property is proved on the correct carrier. Nothing here reduces RH or Collatz to the four TT7 columns, or makes a one-sided computability obstruction into an arithmetic classification.
