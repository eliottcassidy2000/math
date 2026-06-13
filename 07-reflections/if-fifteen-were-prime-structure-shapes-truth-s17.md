---
source: oracle-2026-05-31-S17
status: reflection (meta-pattern / philosophy of mathematical structure)
tags:
  - counterfactual-mathematics
  - factorization
  - reverse-mathematics
  - lonely-runner
  - paley-tournaments
  - structure-vs-truth
---

# If 15 Were Prime: A Knockout Experiment on How Structure Carries Truth

> Prompt (paraphrased): consider how the structures we study would differ if 15
> were prime instead of 3·5. Think in this absurd way, but use it to understand
> *why and how structure shapes mathematical truth.*

## The premise is incoherent — and that is the method

In any model of arithmetic, `15 = 3·5` necessarily. "If 15 were prime" is not a
possible world; it is a *false antecedent*. But mathematicians use false
antecedents all the time — that is what proof by contradiction is. Here we are
not chasing a contradiction; we are running a **knockout experiment**. In
biology you disable a gene and watch which traits break, and the pattern of
breakage tells you what that gene *was for*. "Disable the factorization of 15"
and watch which theorems break: the breakage pattern tells you which truths were
secretly *about the structure 3·5*, and which only wore the name "15."

A number is not an atom. It is a little program, and its prime factorization is
its source code. `15` is the name; `3·5` is the implementation. Asking "what if
15 were prime" is asking to keep the name and swap the implementation — a
controlled experiment that isolates **which facts compile against the
factorization and which compile against something else.** This is, in effect,
*reverse mathematics localized to a single integer*: not "which axioms does this
theorem need," but "which structural facts about this number does it ride on."

Three layers fall out.

## Layer A — truths that FLIP (they were `3·5` all along)

These are not really facts about "fifteen." They are facts about *a product of
two distinct primes*, and they invert the instant you make 15 prime.

- **No field of order 15.** Finite fields exist only at prime powers. `Z/15` has
  zero divisors — `3·5 ≡ 0` — so it is a ring with holes, not a field. Make 15
  prime and `GF(15)` springs into existence, multiplication becomes a group on
  the nonzero elements, and division is everywhere defined. The *existence of an
  algebraic universe* on 15 symbols is gated entirely by the factorization.

- **No primitive root mod 15.** `(Z/15)* ≅ (Z/3)* × (Z/5)* ≅ Z/2 × Z/4`, order 8
  but exponent 4 — *not cyclic*. There is no single generator of the units. Make
  15 prime and `(Z/15)* ≅ Z/14` becomes cyclic with primitive roots. Whether the
  multiplicative world is "one cycle" or "a product of two cycles" is read
  straight off `3·5` versus prime.

- **No Paley tournament on 15 vertices** (repo-relevant). The Paley construction
  needs a prime power `q ≡ 3 (mod 4)`. Now `15 ≡ 3 (mod 4)` *passes* the
  congruence — the only thing that disqualifies it is that 15 is not a prime
  power. A prime 15 would hand us `QR_15` immediately, a doubly-regular
  tournament built from quadratic residues, with all the circulant spectral
  machinery this project studies. The factorization, not the residue class, is
  the gate.

- **CRT and RSA.** `Z/15 ≅ Z/3 × Z/5` is the whole reason 15 is the canonical
  toy RSA modulus — the smallest number you can "factor" as a lesson. A prime 15
  is unfactorable; `φ` jumps from 8 to 14; the entire pedagogical and
  cryptographic role evaporates. 15 is famous *for being composite*.

The lesson of Layer A: some "facts about 15" are loans taken out against its
factorization. Change the collateral and the loan is called.

## Layer B — truths INVARIANT but RE-ROUTED (the conclusion holds, the proof changes)

These are subtler and, I think, the heart of the matter. The *truth value* is
the same in both worlds, but the *path to it* — the proof technology — is
keyed to the structure.

- **Only one group of order 15, and it is cyclic.** True if 15 is prime
  (trivially: prime order ⇒ cyclic). True if `15 = 3·5` (Sylow: `n_3 = n_5 = 1`,
  both subgroups normal, and `3 ∤ 5−1`, `5 ∤ 3−1`, so no nonabelian
  semidirect product — only `Z/3 × Z/5 ≅ Z/15`). Same conclusion, *utterly
  different reason*. In fact 15 is the **smallest composite number with a unique
  group**, and that distinction exists *because* `3·5` happens to satisfy the
  no-twist condition. The prime world gets cyclicity for free; the composite
  world earns it through a coincidence of `3` and `5`. The truth is robust; its
  *meaning* is not.

- **A doubly-regular tournament on 15 vertices exists either way** — but by
  different builders. The Paley builder fails (Layer A). Yet a homogeneous
  tournament on 15 vertices corresponds to a skew-Hadamard matrix of order 16,
  and order 16 = 2⁴ *does* admit one. So the object survives, rehoused in the
  2-adic world of Hadamard constructions rather than the multiplicative world of
  quadratic residues. **Existence can be anchored in more than one structure**;
  knock out one anchor and the object may still stand on another.

- **The Lonely Runner at n = 15** (direct sequel to the n=14 session). LRC is
  expected true at every `n`; that won't change. But *how you'd prove it* is
  dictated by the arithmetic type of `n`:
  - If 15 were **prime**, the odd-prime polynomial-method theorem (it requires
    `k+1` to be an odd prime) would fire on the near-initial-segment families
    `≡(1,…,14) mod p`, exactly as it did for the genuinely-prime cases n = 11
    and n = 13. The proof would be *algebraic*.
  - Because `15 = 3·5`, that tool is dead, just as it was for the even composite
    n = 14. What you gain instead are **folding symmetries**: the divisors 3 and
    5 give rotations `t ↦ t + 1/3` and `t ↦ t + 1/5` of the time circle (the CRT
    analogues of the `t ↦ t + 1/2` antipodal tool that even-ness handed us at
    n = 14). The proof would be *combinatorial / sieve-theoretic*: cover the
    residues, exploit the folds, descend through the proven sub-cases n = 3 and
    n = 5.

  This is the cleanest possible illustration of the prompt. **The factorization
  of `n` selects the proof technology.** Prime `n` → algebra and polynomials;
  composite `n` → folding, CRT, sieve. Two different "physics" for the same
  conjectured truth, and which physics applies is decided by `n`'s source code.
  Compositeness is double-edged: it *gives* more folding symmetry to exploit and
  *demands* more residue scales to cover. Prime-ness removes both the gift and
  the burden and replaces them with a single sharp algebraic blade.

## Layer C — truths anchored ELSEWHERE (untouched)

Some facts about 15 never mention multiplication at all, so the knockout does
nothing.

- `15 = 1+2+3+4+5 = T₅ = C(6,2)`: a triangular number, the number of edges of
  `K₆`, the number of 2-subsets of 6 things. These are *additive / ordinal /
  combinatorial* facts. They are blind to whether 15 is prime. (The 3×3 magic
  square's magic constant `15 = 45/3` is downstream of the `3`, not of `15`'s
  own factorization.)

- The Conway–Schneeberger **15-theorem** — an integer-matrix positive definite
  form representing `1,2,3,5,6,7,10,14,15` represents every positive integer —
  uses 15 as a *magnitude*, a critical threshold, not as a multiplicative
  object. Its 15 would survive the surgery essentially intact.

Layer C is the control group: facts that depend on 15's *position on the number
line* rather than its *internal structure* don't flinch.

## The synthesis: truth is layered by its anchor

Put the three layers side by side and a picture emerges. A mathematical
statement involving "15" is true *relative to the structures it actually
invokes*:

```
Layer A  multiplicative structure (factorization)   → flips under the surgery
Layer B  the conclusion is over-determined; several  → invariant, but the proof
         structures could carry it                     re-routes to whichever
                                                        structure survives
Layer C  additive / ordinal / combinatorial anchor   → untouched
```

"Structure shapes mathematical truth" is not mysticism. It is the precise
observation that **every theorem silently rides on some structures and not
others, and a number's prime factorization is a *bundle* of structures
(a field-or-not, a cyclic-or-product unit group, a Paley-able-or-not, a
CRT-splitting) that many theorems mount without announcing it.** The
"15 is prime" fiction is a debugger breakpoint: it lets us step into a theorem
and read off which of these structures it dereferences.

This reframes the whole n=14/n=15 lonely-runner frontier from the previous
session. The *difficulty profile of the conjecture across `n`* — prime levels
yielding to algebra, composite levels resisting until folding/sieve methods
catch up — is **literally a portrait of the multiplicative structure of the
integers, drawn by a problem that only ever mentions distances on a circle.**
The runners do not know about primes. The factorization of `n` reaches in
through the back door, through which residue covers are forced and which folding
symmetries exist, and decides how hard the truth is to reach. The same shadow
falls across this project's other objects: the Cayley–Dickson tower
(R → C → H → O, then sedenions lose associativity-of-the-rest at dimension 16)
is the identical phenomenon one level up — *existence and properties gated by
the arithmetic type of a dimension*, not by anything the algebra "wants."

## What this is, methodologically

The "absurd" move is a real tool, and it deserves a name in this workspace:
**counterfactual / knockout reasoning on the structure of a specific object.**
Where reverse mathematics asks *which axioms* a theorem needs, this asks *which
structural facts about the particular numbers in play* it needs — by pretending,
locally and impossibly, that one of those facts is otherwise, and cataloguing
the wreckage. The wreckage is the dependency graph. A theorem that survives the
surgery was never really about 15-as-`3·5`; a theorem that dies was.

The deepest takeaway: we habitually treat an integer as a single thing with a
single name, but it is a confluence of independent structures that usually agree
and occasionally are forced to. Most of the time we never notice which structure
is load-bearing, because the number quietly supplies all of them at once. The
only way to see the seams is to imagine pulling one structure out — and 15,
sitting at the smallest crossing of two distinct odd primes, with no field, no
primitive root, no Paley tournament, a uniquely-cyclic group of its order, and a
lonely-runner case that trades algebra for folding, is an almost perfect
specimen on which to perform the cut.
