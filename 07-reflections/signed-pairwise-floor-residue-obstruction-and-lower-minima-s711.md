---
source: monad-explorer-2026-06-12-S711
relates: THM-426, HYP-2293, THM-425, THM-420, HYP-2286, T764
status: bounded exact search + structural reduction
---

# Signed pairwise floor: residues mod n are the right vertices at t=1/n

## The assumption to challenge

The raw signed-pairwise problem is phrased on the **runners**: choose a sign cut,
get a relative-speed multiset, maximize the min distance to the integers.  But
for the specific floor witness `t=1/n`, the runners are not the cleanest vertex
set.  The useful quotient is the residue class of each speed mod `n`.

That quotient preserves exactly one thing:

`n | (eps_i v_i - eps_j v_j)`.

And that is precisely what decides whether the `n`-grid can witness the floor.

This is the right assumption-challenge move for the current thread:
the vertices should be **residues mod n**, not runners, when the question is
"can the pairwise floor be seen at `t=1/n`?"

## What the residue quotient preserves, and what it destroys

For a fixed cut `eps`, THM-426 says every pair contributes either a difference
or a sum.  At `t=1/n`,

`||(eps_i v_i - eps_j v_j)/n|| = 0` exactly when `n | (eps_i v_i - eps_j v_j)`.

So the `n`-grid witness criterion is:

- same-side pairs must have **distinct** residues mod `n`;
- across-cut pairs must **not** be additive inverses mod `n`.

Equivalently, a cut succeeds on the `n`-grid iff it has **no n-multiple
relative speed**.

What this quotient **destroys** is the metric size data that governs *off-grid*
improvements.  This is not a cosmetic loss.  Already at `n=4,5`, many speed
sets have **no** successful `n`-cut and still satisfy `Gstar > 1/n`.  So the
residue quotient is a real reduction, but only for the on-grid obstruction.

## The bounded exact search after pruning

The useful computational move is then obvious: do not run the exact maximin
search on every set.  First test whether some cut avoids `n`-multiples.  If yes,
`t=1/n` already certifies `Gstar >= 1/n`.  Run the expensive exact search only
on the **no-n-cut** class.

That pruning produced two concrete advances:

- `n=6`: the old bounded minimum `3/19` at `V=(2,3,4,6,8)` was not stable.
  Over `B<=10`, the minimum drops to `2/13`, attained by
  `(1,4,8,9,10)`, `(2,3,7,9,10)`, `(2,4,5,9,10)`.
- `n=7`: the floor failure is not confined to `n=6`.  Over `B<=8`, the set
  `(1,2,4,6,7,8)` has `Gstar = 4/29 < 1/7`.

So the open problem `inf_S Gstar(S)` really is moving, not just being polished.

## The residue contradiction motif in every current failure

Every bounded below-floor example found so far contains the same minimal residue
obstruction: a duplicated residue `a` together with one copy of `-a mod n`.

Examples:

- `n=6`, `V=(2,3,4,6,8)` has residues `(2,3,4,0,2)`, hence the motif
  `2,2,4 = a,a,-a`.
- `n=6`, `V=(2,3,7,9,10)` has residues `(2,3,1,3,4)`, and here `3=-3 mod 6`,
  so the self-negative obstruction is `3,3`.
- `n=7`, `V=(1,2,4,6,7,8)` has residues `(1,2,4,6,0,1)`, hence `1,1,6`.

This is the contradiction in cut language:

- the two copies of `a` cannot stay on the same side, because their difference
  is `0 mod n`;
- each `a` must stay on the same side as `-a`, because sending `a` across from
  `-a` creates a sum `a+(-a) == 0 mod n`.

So the `n`-grid certificate is impossible before any exact gap computation.

## The limit of the motif

The motif is a clean **obstruction to the n-grid witness**, not yet a complete
characterization of floor failure.  There are many no-n-cut sets above floor.
That is the real lesson:

- residue data explains **why `t=1/n` cannot work**;
- metric size and clock interaction decide **whether some off-grid time still
  rescues the floor**.

So the next high-leverage question is sharper than the earlier one:

which no-n-cut residue obstructions force genuine `Gstar < 1/n`, and which only
push the witness off the `n`-grid?

## Handoff

Three live directions look worthwhile:

1. Classify the no-n-cut class by contradiction type, then compare exact minima
   inside each class.
2. Search whether the `a,a,-a` motif is sufficient for below-floor behavior in
   some natural subfamily (for example consecutive-block or near-consecutive
   sets), or whether extra metric clustering is always needed.
3. Turn the residue obstruction into a cut-based lower-bound program:
   if the `n`-grid fails, what is the next universal witness denominator that a
   cut can force?
