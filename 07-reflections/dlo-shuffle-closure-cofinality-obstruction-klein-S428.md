# DLO shuffle closure: a cofinality obstruction

**klein-S428, 2026-07-31.** Reflection, not canon. Answers a question the owner
posed outside the repo's programme: for which `n < 1000` is it consistent with
Zermelo set theory that a set carries exactly `n` dense linear orders without
endpoints up to isomorphism?

## The reduction to order theory

Let `(X,R)` be a DLO and `S = {s_1 < ... < s_k}` a finite cut set.  The complement
splits into `k+1` open intervals (each a DLO) and `k` points.  A **shuffle**
permutes the intervals, optionally reverses each, and redistributes the points;
density plus absence of endpoints force strict alternation
`interval, point, interval, ..., interval`.

Every shuffle is **definable from `R` and finitely many parameters**, so in any
model of Z, if `R` is a DLO on `X` then so is every shuffle of `R`.  Hence

    n(X) >= |shuffle closure of R|   for every DLO R on X,               (1)

the closure being the set of isomorphism types reachable by iterated shuffling.
`tau*` is always in the closure (zero cuts, reverse the one piece).

So the question "is any finite `n >= 2` achievable" is bounded below by the pure
order-theoretic question: **which finite sizes can a shuffle closure have?**

## The obstruction

Write `coi` and `cof` for coinitiality and cofinality.  A shuffle's `(coi,cof)`
is `(coi of the first piece as oriented, cof of the last piece as oriented)`, and
first and last are distinct pieces.

**One cut.**  `tau = lambda + 1 + upsilon`; set
`(a,b) = (coi,cof)(lambda)`, `(c,d) = (coi,cof)(upsilon)`, so `tau = (a,d)` and
`tau* = (d,a)`.  The one-cut shuffles realise every pair in
`{a,b} x {d,c}` and `{c,d} x {b,a}`.  Requiring all of them to lie in
`{(a,d),(d,a)}` leaves exactly two equality patterns, and **both** force

    coi(lambda) = cof(lambda)   and   coi(upsilon) = cof(upsilon).        (2)

One cut therefore permits `coi(tau) != cof(tau)`.

**Two cuts, and this is the decisive one.**  `tau = P + s + Q + t + R` with
`(p1,p2), (q1,q2), (r1,r2)` the three pieces' characters.  Exhausting all
`(first,last)` choices and both orientations and demanding every result lie in
`{(p1,r2),(r2,p1)}` leaves exactly **one** equality pattern out of the whole
space: `p1 = p2 = q1 = q2 = r1 = r2`.

**Theorem.**  If a DLO type `tau` has shuffle closure of size at most `2`, then
at every two-cut decomposition all three pieces satisfy
`coi = cof = coi(tau) = cof(tau)`.  In particular

    coi(tau) != cof(tau)   ==>   closure >= 3.                            (3)

**Worked kill.**  `tau = Q . omega_1` (`omega_1` copies of `Q`) is the obvious
non-self-reverse candidate: `coi = omega`, `cof = omega_1`, so `tau !~ tau*`.
Cutting anywhere gives `lambda ~ Q` and `upsilon ~ Q . omega_1`.  The shuffle
`lambda + 1 + upsilon*` has `(coi,cof) = (omega, omega)`, which is neither
`(omega, omega_1)` nor `(omega_1, omega)`.  So its closure is at least `3`, in
agreement with (3).

## What this does and does not settle

It does **not** close the question.  By (3) any closure-`2` witness must be
**cofinality-homogeneous**: every open interval and both ends share one character.
That is a severe rigidity but not an emptiness proof -- a DLO of size `aleph_1`
with `coi = cof = omega` throughout is not excluded by this argument, and the
`(coi,cof)` invariant is by construction blind to it.

Current status of the whole question, stated honestly:

* `n = 0` (finite nonempty `X`) and `n = 1` (countably infinite `X`, Cantor,
  choice-free and Z-legal) are PROVED achievable;
* under AC the spectrum is exactly `{0,1}` (uncountable well-orderable `X` carries
  infinitely many, via `B = W x_lex Q` and `L_n = B + ... + B` separated by
  "`[x,y]` is strictly smaller than `X`");
* the ordered Mostowski model gives `aleph_0`, **not** `2` -- the natural folklore
  reading is REFUTED, with explicit witness `(c,infinity) + {c} + (-infinity,c)`;
* a "trivial filter" permutation model does **not** work: `F = {G}` is not a
  normal filter, no atom is symmetric, so no atom is in the model.  I asserted
  such a construction earlier in the session and **retract it**;
* whether any finite `n >= 2` is consistent remains **OPEN**.  Hence
  `sum >= 1`, with equality conjectural.

The remaining target is exactly: a cofinality-homogeneous DLO type, not
self-reverse, whose every finite cut-reverse-permute shuffle is isomorphic to it
or its reverse.  (3) says nothing else can work.

---

## Addendum (klein-S428, 2026-08-01): one closed door, one completed classification

### A closed door: "make it gapless and the shuffles disappear" -- REFUTED

I tried to weaken the requirement in (3) as follows. Call a cut `(A,B)` of `L` a
**gap** when `A` has no max and `B` has no min. If the shuffle were `B + A`, then
`B + A` is a DLO without endpoints *only* at a gap (otherwise the realized point
becomes a global min or max), and a **Dedekind-complete** `L` has no gaps at all
-- so gapless DLOs would be shuffle-rigid for free, and the open long ray
(dense, no endpoints, Dedekind complete, `coi = omega != omega_1 = cof`, hence
**not** self-reverse) would be an immediate closure-`2` candidate.

**This is wrong, and the reason is exactly the three-block form already used
above.** The shuffle is not `B + A`; it is

```text
   L = A + {c} + B   |-->   B + {c} + A,        A = (-inf,c),  B = (c,inf).
```

`B` has no min and `A` has no max because both are *open*, so `B + {c} + A` is a
DLO without endpoints for **every** point `c`, gap or no gap. Dedekind
completeness buys nothing. Concretely at `L = R`, `c = 0`:
`(0,inf) + {0} + (-inf,0) ~ R + 1 + R`, and `R + 1 + R !~ R` -- take
`S = R + {pt}`, whose upper bounds are the second copy of `R`, which has no
least element, so `R + 1 + R` is not Dedekind complete while `R` is. So even the
continuum carries at least two DLO types in plain ZF, and the long ray is not a
shortcut. Statement (3) and its cofinality-homogeneity requirement stand
unweakened.

Recorded so the "gapless" idea is not re-attempted.

### A completed classification: `Sym(A)` in the ordered Mostowski model

The reflection above refutes the folklore `n = 2` reading with one explicit
witness. The following lemma turns that into a complete classification, because
it identifies *every* isomorphism available in the model.

> **Lemma.** In the ordered Mostowski model (atoms `A` order-isomorphic to `Q`,
> `G = Aut(A,<)`, finite supports), the bijections `A -> A` that belong to the
> model are exactly the **finitary** permutations: those equal to the identity
> outside a finite set.

*Proof.* Let `f` be a bijection with support `F`, and let `a` not in `F`. Every
`g` in `G_F` fixing `a` satisfies `g(f(a)) = f(g(a)) = f(a)`, so `f(a)` is fixed
by the pointwise stabilizer of `F ∪ {a}`. That group is `Aut(A,<)` fixing
`F ∪ {a}` pointwise, and each of its orbits other than the fixed points is an
infinite open interval; hence its fixed-point set is exactly `F ∪ {a}`, giving
`f(a) in F ∪ {a}`. By injectivity at most `|F|` points can map into `F`, so
`f(a) = a` for all but finitely many `a`. `f` is a bijection, so it permutes a
finite set. Conversely every finitary permutation has finite support. `QED`

Two consequences, both needed to make the `aleph_0` count rigorous rather than
suggestive:

1. **The order relations.** A DLO `R` with support `E`, `|E| = k`, is
   `G_E`-invariant, and `G_E` is transitive on ordered pairs inside each of the
   `k+1` open `E`-intervals and on each product of two distinct blocks. So `R`
   is determined by a linear order on the `2k+1` blocks together with one
   orientation bit per interval. Density plus "no endpoints" forces the block
   order to **alternate**, beginning and ending with an interval: with `k`
   singletons and `k+1` intervals there are exactly `k` internal gaps, so each
   receives exactly one singleton. This gives `(k+1)! · k! · 2^(k+1)` orders
   supported by `E`.
2. **The isomorphisms.** By the Lemma an isomorphism is finitary, so it is the
   identity off some finite `F`; enlarging `F` to contain both supports, two
   isomorphic orders must have *literally equal* restrictions to the cofinite set
   `A \ F`. Since each block is infinite, the block arrangement and the
   orientation bits are read off that cofinite restriction. Hence orders with
   different arrangement/orientation data are non-isomorphic, and the count is
   infinite -- the folklore `n = 2` fails by an unbounded margin, not narrowly.

Neither point changes the headline status: `n = 0` and `n = 1` achievable,
`n >= 2` OPEN, `sum >= 1`.
