# P-adic-zeta density: specialization kernel and transverse-jet frontier

> **SUPERSEDED REFLECTION / HISTORICAL PROVENANCE ONLY.** MISTAKE-527 records
> that this session mistook a universal torsor identity for `u->f`. The exact
> graph-kernel algebra survives, but the asserted manuscript gap and downstream
> dependency collapse do not. See THM-4255 and the corrected 2026-08-26 source
> audit before reusing any claim below.

**Session date:** 2026-08-26--27

**Historical status (SUPERSEDED):** THM-4255's algebra remains proved; the
external verdict in this reflection was withdrawn. LRC(14), JC(2), DC(2), and
all Rule-30 prizes remain `OPEN`.

## 1. Outcome first

The user's proposed counterexample identifies a real load-bearing gap. For
the diagonal specialization `ev:k[u][[f]] -> k[[f]]`, `u->f`,

```text
ker(ev)=(u-f),
ev^(-1)((f^n))=(u-f)+(f^n).
```

Thus vanishing after specialization does not imply coefficientwise vanishing
before specialization. Explicit-`f` Cartier channels split `u-f` into the
nonzero pieces `u` and `-1`, so they do not descend through the graph quotient.
This invalidates the proofs of the two narrow window propositions.

The repair is exact, not heuristic. For the genuine digit space `V`, prove

```text
V/(V intersect f^nS) -> k[[f]]/(f^n)
```

injective, and separately prove strictness for the `ell`-adic coefficient
filtration. Alternatively, specialize first and prove a degree-`D+C` bound on
each scalar digit. Neither theorem is present in the audited source.

## 2. Inheritance pass

### Anchor: p-adic-zeta density gate

- **Closest proved mechanism:** THM-4091 separates cumulative LCM transport
  from unsafe literal coefficient depth under coordinate change.
- **Canonical hostile:** `u-f`, zero on the graph but nonzero in the source.
- **Corrected near miss:** the external proofs use valid Taylor/block
  no-backflow lemmas only after an invalid specialization-order promotion.
- **Least-used sidecar:** the graph conormal coordinate and its transverse
  Hasse tower.

### Niche: planar-Jacobian wall entry

- **Closest proved mechanism:** the current `W=0` results recompute residual
  objects directly and do not lift scalar wall identities coefficientwise.
- **Canonical hostile:** a relation can vanish on `W=0` merely because it is
  divisible by `W`.
- **Corrected near miss:** divide by the exact `W`-order before restriction;
  do not infer a universal relation from the wall shadow.
- **Least-used sidecar:** the conormal strict transform
  `W^(-nu_W(R))R|_(W=0)`.

### Wildcard: local packets versus global height

- **Closest proved mechanism:** THM-3488 keeps a Rule-30 transverse Hasse
  channel after its base Cartier value vanishes.
- **Canonical hostile:** LRC local Hasse packets can agree while owner heights
  differ, so local completeness does not imply physical loneliness.
- **Corrected near miss:** a complete local coordinate packet still needs the
  height/owner sidecar.
- **Least-used sidecar:** associated-graded strictness for the consumer's
  actual filtration, rather than for an ambient quotient.

## 3. Live concept board

| Concept | Source | Target | Preserved predicate | Destroyed information | Cheapest decisive test |
|---|---|---|---|---|---|
| Graph specialization | `k[u][[f]]` | `k[[f]]` | scalar germ | multiples of `u-f` | compute the graph-ideal intersection |
| Explicit-`f` Cartier | bidegree channels | residue blocks | source `f` residue | `u`-degree shifts channels | test `u-f` in channels zero and one |
| Transverse Hasse tower | graph neighbourhood | conormal jets | graph multiplicity | nothing with full tower | find first nonzero `D_u^[j]` |
| Arithmetic lattice | restricted digit module | `ell`-graded image | primitivity if strict | hidden `ell` factors | test `u-1 -> ell f^N` |
| Wall specialization | planar-JC family | `W=0` residual | strict transform | wall-normal order if divided too early | compare `gr_W(ker)` and `ker(gr_W)` |
| Fixed-prime compiler | exponent orbit | suffix cylinder | declared suffix predicate | higher bits, tracked explicitly | prove orbit/lift/union before density |

Every proposed transfer was tested against all six concepts. The p-adic
source has the exact graph and channel maps. The planar-JC transfer preserves
only the strict-transform operation. The LRC and Rule-30 links are sidecar
analogies, not theorem transport. The percolation paper preserves no relevant
object or predicate and was rejected as a typed connection.

Incoming THM-4257 supplies the positive density control. Its exponent-to-
suffix map is an exact iff for the declared predicate; it then proves monotone
lifts, identifies the full-height set with the union of finite levels, and
shows that fixed orbits become complete high-bit cylinders. Only after that
loss ledger does it take density one. The p-adic manuscript instead takes a
noninjective scalar shadow before its coefficientwise operation. The shared
word “density” is not the connection; the connection is the proof order of
quotient, lift, and measure.

## 4. Independent hostile and positive controls

Two independent standard-library paths were retained.

1. The primary triangular-linear-algebra path checks degrees one through nine,
   54 graph-plus-tail jet cells, five Cartier moduli, and Hasse multiplicities
   one through six.
2. The clean-room path exhausts all 512 polynomials in a `3x3` monomial box
   over `F_2`. A raw rectangle has short-jet rank `5/9`; graph-normal
   representatives have rank `5/5`. This is the finite model of exactly what
   the missing quotient must accomplish.

The abstract window hostiles are:

```text
Prop. 6.2 syntax:  (D,C,ell,N)=(10,0,23,15), pivot 38;
Prop. 12.3 syntax: (D,C,ell,N)=(12,1,29,17), pivot 46.
```

Neither pivot lies in its asserted `j+ell*r` window. These controls refute the
unqualified algebraic implication, not the modular-geometric conclusions.

## 5. Exact external dependency boundary

The `1<=a<=s_*` pole-grade bound, Lemma 6.1, Lemma 12.2, and the separate
Fitting-length/capacity calculations are not invalidated. What is lost is the
narrow prime window. Therefore the genus-zero and positive-genus density
chains, including the all-prime headline, remain unsupported until one of the
two repair routes is proved.

This is not entered in `MISTAKES.md`: no current repo theorem asserted the
false implication. THM-4255 is a new firewall and the source audit records the
external proof failure. The older 22-value draft is now a targeted audit risk,
not automatically false; its closest decisive case is `(p,s)=(5,5)` because
its terminal numerical margin is narrowest there.

## 6. Generated next tasks

### Anchor tasks

1. Extract the actual finite digit modules in Propositions 6.2 and 12.3 and
   compute the obstruction quotient
   `V intersect ((u-f)S+f^nS) / (V intersect f^nS)` at small weights.
2. Test whether the negative pole-grade symbols descend to scalar base germs
   before specialization. A genuine descent theorem would bypass the graph
   kernel; a single nonzero graph multiple blocks that route.
3. Attempt the scalar-first repair: after the actual section, bound every
   nonzero digit order by `D+C`. Track the exact line bundle and divisor so a
   pre-specialization global-section bound is not silently reused.
4. Audit `ell`-adic associated grades independently of graph normalization.

### Niche tasks

5. For each current planar-JC residual relation `R(W,t)`, record `nu_W(R)` and
   its conormal leading coefficient; test whether entry obligations live in
   `gr_W(ker)` or only in `ker(gr_W)`. Do not subtract any of the current
   1,512 incidence obligations without a target-level exclusion.
6. Revisit the old 22-value p-adic draft at `(5,5)`: identify vertices as
   digit blocks, not scalar pivots, and compute its graph-kernel intersection.

### Wildcard tasks

7. Generalize THM-4255 to a moving section `u->sigma(f)` and then to several
   torsor coordinates. The expected kernel is the full graph ideal; verify the
   row-finiteness hypotheses before canonization.
8. Build a “channel migration matrix” sending bidegrees `(a,i)` to scalar
   Cartier residues `a+i mod ell`. Test whether the actual source lies in a
   subspace on which this matrix is triangular or invertible.
9. Model graph-normal digit spaces as an inverse system and ask whether good
   scalar jets form complete cosets or cylinders under lifting. Use THM-4257
   only as a proof-order template; construct the p-adic objects independently.
10. Apply the same graph-ideal audit to any future density, specialization, or
   boundary proof before using a coefficientwise operation.

## 7. What changed across the main frontiers

- **P-adic density:** a proof gap became an exact obstruction theorem and two
  sharply testable repair programs.
- **Planar Jacobian:** no incidence was removed, but wall entry now has a
  concrete conormal computation instead of a vague “specialization caution.”
- **LRC(14):** no bound changed; the owner-height sidecar gained an independent
  algebraic analogue.
- **Rule 30:** no prize changed; retained Hasse channels became a positive
  control rather than a speculative analogy.

The durable research move is: calculate the kernel of every forgetting map,
quotient by the consumer's own tail, and only then apply coefficientwise or
graded operations.
