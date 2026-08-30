# Higher-order truth versus finitary Heyting certificates -- source audit

**Status:** **CITED PREPRINT v1 + MANUAL PROOF AUDIT, 2026-08-29.**  This
note imports a mechanism and records two repairable typographical errors.  It
does not promote the preprint to proved repo canon, and no open problem below
is solved by the citation alone.

## Primary paper and exact claim

Lingyuan Ye and Yiqi Xu, [*Failure of Higher-Order Truth within
Intuitionistic Propositional Logic*, arXiv:2608.26874v1](https://arxiv.org/abs/2608.26874v1),
submitted 2026-08-27.

Theorem 5.5 claims that the free Heyting algebra `F_2` on two generators is
not isomorphic to `Sub_E(1)` for any elementary topos `E`.  Corollary 5.6
claims the same exclusion for every Heyting algebra admitting a surjection
onto `F_2`.  The direction is important: this is a **quotient obstruction**,
not a subalgebra obstruction.

The sharp free-rank boundary is:

- `F_0 ~= 2` is Boolean and topos-admissible;
- `F_1` is already its profinite completion and hence is complete and
  realizable by a localic topos; see Lingyuan Ye,
  [arXiv:2604.08267v2, Example 4.6](https://arxiv.org/abs/2604.08267v2),
  together with Darniere--Junker below;
- every free `F_n` for `n>=2`, including the infinitely generated cases,
  surjects onto `F_2` and is excluded by the preprint's corollary.

Containment goes the other way: the complete algebra
`widehat(F_2) ~= O^up(K_2)` contains `F_2` and is topos-admissible.  Thus an
unproved identification of a repo certificate lattice with `F_2`, or merely
finding an `F_2` subalgebra, would not activate the theorem.

## Audited proof mechanism

The proof has three load-bearing layers.

1. **Completion escape.**  Bellissima's universal image-finite Kripke model
   embeds `F_2` in the upset algebra of `K_2`, its profinite completion.  A
   finite cyclic recurrence constructs `x_n,y_n,r_n`, and then incomparable
   points `z_n`.  The infinite upset

   ```text
   A = union_n up(z_n)
   ```

   is in the completion but not in `F_2`.  The nonmembership uses finite
   decomposition into join-irreducibles and the fact that the corresponding
   completion elements are downward filtered.

2. **Impredicative orbit closure.**  Principal and coprincipal formulas
   package the recurrence as a global successor `S` on `Omega^6`.  Without a
   natural-numbers object or countable coproduct, higher-order logic defines
   the least `S`-invariant predicate by the Church-style formula

   ```text
   Reach(u) = forall Q subset U,
              ((Q(u_0) and forall v (Q(v) -> Q(Sv))) -> Q(u)).
   ```

   A uniform observable `Z` then makes
   `theta = exists u (Reach(u) and Z(u))` a global truth value.

3. **Finite observers plus separators.**  Restriction to each finite layer
   `K_(2,d)` gives finitely many equivalence classes of the external orbit.
   Crucially, each equivalence is an `S`-congruence.  A finite internal
   disjunction of representatives is therefore inductive, while coprincipal
   formulas separate points outside `A`.  These two directions force
   `theta=A`, contradicting `theta in F_2`.

Darniere--Junker had already observed that the union of cones over an infinite
antichain need not be finitarily definable.  Ye--Xu's new device is to make
such an antichain the orbit of one finite Heyting recurrence, so its union is
higher-order definable.

## Adversarial audit and repairs

No fatal gap was found in the proof steps above.  Two arrows/ranks in Lemma
3.1 need correction.

- The construction gives
  `z_n in K_(2,n+1) \ K_(2,n)`, not
  `z_n in K_(2,n) \ K_(2,n-1)`.
- In the antichain argument the required inequalities are
  `x_n not<= z_m`, and ultimately `x_0 not<= z_m`.  The displayed reverse
  inequalities are not the ones used and the terminal reverse statement is
  generally false.  Reversing those arrows repairs the induction.

In Proposition 3.2, the asserted `z_l<=w` is valid because the chosen common
lower bound also lies in `Q_i subset A`; that membership sentence is omitted.
Lemma 5.2 needs only an externally chosen finite representative list and an
internal finite disjunction.  It assumes neither internal choice nor an NNO.

These repairs should be rechecked against any later arXiv version.

## Positive boundary and related primary sources

Marco Abbadini, Rodrigo Nicolau Almeida, and Igor Arrieta,
[*A topos for etale-finite Heyting algebras*,
arXiv:2606.03861v2](https://arxiv.org/abs/2606.03861v2), construct an
elementary topos for every etale-finite Heyting algebra (Theorem 5.25) and
give an iff characterization for finitely propositional topoi (Theorem 6.5).
This is a positive boundary, not a claim that every non-`F_2` algebra is
admissible.

Other proof inputs and provenance:

- Fabio Bellissima, [*Finitely generated free Heyting
  algebras*](https://doi.org/10.2307/2273952);
- Luck Darniere and Markus Junker,
  [*On Bellissima's construction of the finitely generated free Heyting
  algebras, and beyond*, arXiv:0812.2027](https://arxiv.org/abs/0812.2027);
- Alasdair Urquhart,
  [*Free Heyting algebras*](https://doi.org/10.1007/BF02945107);
- Guram Bezhanishvili and Nick Bezhanishvili,
  [*Profinite Heyting algebras*](https://doi.org/10.1007/s11083-008-9089-1).

## Reusable transfer contract

The useful repo pattern is not the noun "topos".  It is the following typed
closure test.

```text
source:       finitary certificate algebra H and a term-defined successor S;
target:       a completion L containing an orbit-union A outside H;
map:          orbit u_n -> finite observer q_d(u_n) and observable Z(u_n);
preserved:    successor dynamics, fixed-depth observations, separators;
destroyed:    iteration count and cross-depth transition compatibility;
sidecars:     transition maps, owners/phases/addresses, physical admissibility;
decisive test:q_d must be an S-congruence and preserve the target predicate.
```

This is a no-go tool for an allegedly complete **finitary certificate
language**.  It does not directly prove LRC(14), a Rule 30 prize, Mahler
`3/2`, JC(2), or a tournament conjecture.  In each case an actual infinite
state object, definable successor, observable, and separating family must be
constructed first.

## Repo consumers and guardrails

- [THM-4286](../../01-canon/theorems/THM-4286-signature-response-nonfactorization-and-two-deck-surgeries.md)
  uses only the finite-observer audit principle: a response predicate factors
  through a signature quotient iff every signature fibre is pure.  Its exact
  surgery is repo-derived and does not depend mathematically on the preprint.
- [THM-4210](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md)
  supplies a genuine two-successor physical state tree, but its odd channel,
  boundary current, and quadratic admissibility are mandatory sidecars.
- [THM-4072](../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md)
  supplies exact reset/native successor data; liveness is not determined by
  finite terminal flags.
- [THM-4090](../../01-canon/theorems/THM-4090-two-sort-matching-logic-global-completeness-obstruction.md)
  is fixpoint-free.  `Reach` quantifies over predicates and cannot be silently
  encoded in that basic fragment.
- [MISTAKE-397](../../01-canon/MISTAKES.md) is the closest general hostile:
  identical finite objects do not determine an inverse limit without the
  transition arrows.

**Scope firewall:** all paper results in this note are **CITED PREPRINT
CLAIMS**, with the named local repairs.  Only separately audited repo
computations or canon theorems may carry `FINITE-EXACT` or `PROVED` status.
