# The logical status of LRC(14): it is Π⁰₁, so "impossible to disprove" means "true" — and there is no third way

*kind-pasteur-2026-06-22-S31z. The owner asks: consider that LRC(14) might be both impossible to prove
and impossible to disprove, and how to rigorously prove that indeterminacy. The answer turns on the
arithmetic complexity of the statement — and it collapses the premise in an instructive way.*

## Step 1 — LRC(14) is a Π⁰₁ statement
LRC(14) reduces to integer speeds (the standard density/compactness reduction; the integer case is the
hard core). For a FIXED 13-tuple `S` of integers, `M(S) = max_t min_i ‖v_i t‖` is **computable**: the
max-min is attained at a rational critical point `t = k/(v_i ± v_j)` or `k/(2v_i)` (finitely many), so
`M(S) ≥ 1/14` is a **decidable** predicate of `S`. Moreover the witness, when it exists, has denominator
**polynomially bounded in `max(S)`** (THM-565: `V* = arcCount/c`, a fixed multiple of the speeds). Hence

> **LRC(14) ≡ `∀S ( M(S) ≥ 1/14 )`, a Π⁰₁ sentence** — a universal quantifier over a primitive-recursive
> predicate. (Equivalently `∀S ∃t<poly(S): t is a witness`, the bounded `∃` being decidable.)

## Step 2 — the Π⁰₁ collapse: "impossible to disprove" = "true"
For any **sound** theory `T` (PA, ZFC, …) and a Π⁰₁ sentence `φ`:

> **`T ⊬ ¬φ` ⟺ `φ` is true ⟺ `Con(T + φ)`.**

*Proof.* If `φ` is false, there is a specific counterexample `S₀` with `M(S₀) < 1/14`; that is a finite
computation, so `T ⊢ ¬φ` (every theory proves true Σ⁰₁ facts). Contrapositive: `T ⊬ ¬φ ⟹ φ` true.
Conversely if `φ` is true, sound `T` cannot prove the false `¬φ`. ∎

So for LRC(14): **"impossible to disprove" is not a separate epistemic state — it is logically identical
to "LRC(14) is true."** A disproof would be a finite, checkable counterexample; its non-existence *is*
the conjecture.

## Step 3 — so the premise is "independence," and the truth value is DETERMINED
"Impossible to prove AND impossible to disprove in `T`" therefore means exactly:

> **LRC(14) is independent of `T` = ( true ) ∧ ( `T`-unprovable ).**

There is **no genuine alethic indeterminacy**: in the standard model `ℕ`, LRC(14) is either true or false
(bivalence holds for arithmetic sentences). "Indeterminate" can only mean **`T`-undecidable** — a fact
about the *theory's strength*, not about LRC(14) lacking a truth value. (Contrast set-theoretic
statements like CH, which are genuinely independent of ZFC because they quantify over uncountable sets
and *change truth value across models*; a Π⁰₁ arithmetic statement is **absolute** — true in one ω-model
iff true in all.)

## Step 4 — how one would rigorously prove the independence (and why not in `T`)
To prove "LRC(14) is independent of `T`" one must, **in a metatheory `M ⊋ T`**, establish both:
- **(unrefutability = truth)** `M ⊢ LRC(14)` — i.e. *prove the conjecture* (in `M`). By Step 2 there is
  **no shortcut**: you cannot show `T ⊬ ¬LRC(14)` without proving LRC(14) true.
- **(unprovability)** `M ⊢ "T ⊬ LRC(14)"` — exhibit a model of `T + ¬LRC(14)`: a **nonstandard** speed
  set whose `M < 1/14` *inside that model*, where `T` lacks the infinitary resources to refute it.

And one **cannot** prove the independence *inside `T`*: a `T`-proof of `T ⊬ ¬LRC(14)` would, by Step 2
formalized, yield `T ⊢ LRC(14)`, contradicting `T ⊬ LRC(14)` (this is the second-incompleteness pattern).
So a rigorous indeterminacy proof is intrinsically a **metatheoretic** result: prove LRC(14) above `T`,
and bound `T` below it.

## Step 5 — for LRC(14) specifically, indeterminacy is almost certainly FALSE
The Paris–Harrington/Goodstein independences from PA come from witnesses growing like
PA-non-provably-total functions (Ackermann, ε₀). LRC(14)'s witnesses are **polynomially bounded in the
speeds** (THM-565), and the only known unboundedness (THM-566: no `V`-uniform denominator) is the lcm
family, whose witness denominators grow at most **exponentially** in the construction parameter —
**PA-tame** (PA proves exp-bounded witness existence). The proof input the repo isolated (effective Weyl
equidistribution, three-distance, the comb-teeth count) is **elementary and PA-formalizable**. So:

> **If LRC(14) is true, it is almost certainly PA-provable** (a fortiori ZFC-provable). The premise
> "impossible to prove" is false for any reasonable `T ⊇ PA`; LRC(14) is hard, not independent.

"No finite certificate" (THM-566) is the signature of needing an *unbounded but tame* argument — it would
at most separate LRC(14) from very weak fragments (bounded arithmetic `I∆₀`, which doesn't prove `exp`
total), never from PA/ZFC. That is a technical proof-complexity statement, not the deep undecidability the
premise evokes.

## Step 6 — the payoff: the premise reveals there is NO third way
The most useful consequence is for our own program. By Step 2, **"impossible to disprove" = "true."** So
the apex-7 realizability argument (S31y: a counterexample would realize the forbidden `K_3`, `H=7`) is
**not an indeterminacy escape — it is a proof route**. If we make "the apex-7 over-cover IS the forbidden
`K_3`" rigorous, we have proved `T ⊬ ¬LRC(14)`, which *is* `LRC(14)` true. There is no logical room for a
statement that is simultaneously "shown impossible to disprove" and "left unproven": for Π⁰₁, those are
the same act. The indeterminacy premise, taken rigorously, **collapses onto the conjecture itself** — and
tells us the realizability obstruction, if it works, settles LRC(14) outright.

## Net
LRC(14) is Π⁰₁; "impossible to disprove" is *identical* to "true"; so the premise = `T`-independence =
true-but-`T`-unprovable, with the truth value fully determined (no alethic indeterminacy). Rigorously
proving independence is a metatheoretic task that *contains* proving the conjecture, and cannot be done
inside `T`. For LRC(14) the witnesses are PA-tame, so independence from PA/ZFC is almost surely false — it
is a hard but ordinary true Π⁰₁ statement. The lasting lesson: **proving LRC(14) "impossible to disprove"
would be proving it true** — the realizability/`K_3` route is the proof, not an escape from it.

→ S31y (the `K_3`/`H=7` realizability obstruction), THM-565 (poly-bounded witnesses), THM-566 (no
`V`-uniform certificate), HYP-2899 (the analytic equidistribution input), [[lrc14-thread]].
