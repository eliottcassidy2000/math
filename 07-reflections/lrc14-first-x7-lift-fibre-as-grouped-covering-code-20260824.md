---
status: FINITE-EXACT niche result plus proved arithmetic reduction; LRC(14) OPEN
source: lift_tower_covering / root LRC(14) session, 2026-08-24
tags:
  - lonely-runner
  - lrc14
  - finite-check
  - lift-tower
  - covering-code
  - grouped-sat
  - endpoint-phase
---

# The first `x7` lift fibre is a grouped covering code

## Verdict

There is a faithful future-lift object which the level-one capacity tests do
not see.  If a parent tuple `v=(v_1,...,v_13)` is represented modulo `q`, with
`gcd(q,7)=1`, its first `x7` fibre is

```text
x_i = v_i + q a_i,                 a=(a_1,...,a_13) in F_7^13.       (1)
```

Every genuinely new clock `h/(7q)`, `7 not| h`, supplies one forbidden digit
in each coordinate, except at an exact threshold endpoint, where that
coordinate has no forbidden digit.  Thus the improper children are exactly a
**grouped covering-code / positive-CNF fibre**: choose one symbol in each of
thirteen seven-symbol groups and hit every clock word in at least one matching
coordinate.  The new divisor-`7` gcd alternative is a separate, load-bearing
sidecar.

This changes one parent check from a flat `7^13=96,889,010,407` child scan to
`91` labelled digit literals and at most `6q` clock clauses.  It does not make
the resulting grouped covering problem automatically easy, classify the
parents reaching the lift, or prove LRC(14).

The exact companion is
[`04-computation/lrc14_x7_lift_fibre_clause_code_20260824.py`](../04-computation/lrc14_x7_lift_fibre_clause_code_20260824.py),
with frozen transcript
[`05-knowledge/results/lrc14_x7_lift_fibre_clause_code_20260824.out`](../05-knowledge/results/lrc14_x7_lift_fibre_clause_code_20260824.out).

## Inheritance pass

- **Closest proved mechanism:**
  [THM-573, level-seven lift sieve](../01-canon/theorems/THM-573-lrc14-level7-sieve-sharpening.md)
  chooses among seven physical lifts and uses the fact that one 7-coprime
  speed forbids at most one.  The present object transposes that exact
  one-forbidden-sheet law into the *coefficient fibre* over a fixed paper-sieve
  parent.
- **Canonical hostile:** the AP parent at `q=86=2*43` has an actual improper
  child in its first `x7` ansatz fibre.  Hence neither “one prior `x2` lift is
  enough” nor the clock-code construction itself is an emptiness theorem.
- **Corrected near miss:**
  [HYP-2758](../05-knowledge/hypotheses/HYP-2758-lrc14-is-the-polynomial-methods-composite-c14-lift.md)
  already retracts the idea that the S--T lift factors as a bare CRT product.
  The legal statement here is about the literal fibre `(1)`, not an
  identification of the entire lift algorithm with `Z_2 x Z_7`.
- **Least-used relevant sidecar:** THM-3995's oriented endpoint parity and
  [THM-4002's signed cross-phase](../01-canon/theorems/THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family.md)
  say that boundary locations, not support size alone, are the faithful data.
  Here the exact endpoint flag is what shortens a clock clause.

## Exact derivation

Put `N=7q` and fix `h` with residue `r=h mod 7` nonzero.  At coordinate `i`,

```text
h x_i mod N = h v_i + q r a_i mod 7q.                    (2)
```

Write `h v_i mod 7q = j_i q + u_i`, with `0<=u_i<q`.  As `a_i` ranges over
`F_7`, equation `(2)` visits the seven lifts

```text
u_i, u_i+q, ..., u_i+6q
```

in a permuted order.  At radius `1/14`, danger means

```text
min(z,7q-z) < q/2.                                      (3)
```

Consequently:

- if `2u_i<q`, only layer `0` is dangerous;
- if `2u_i>q`, only layer `6` is dangerous;
- if `2u_i=q`, both endpoints are safe and no layer is dangerous.

Let `d_i(h)` be `0`, `6`, or the endpoint symbol `*` accordingly.  The unique
forbidden digit, when present, is

```text
f_i(h)=r^(-1)(d_i(h)-j_i) mod 7.                         (4)
```

The clock `h/(7q)` is a witness precisely when

```text
a_i != f_i(h) for every non-* coordinate i.              (5)
```

Therefore a child has no new clock witness precisely when

```text
for every h with 7 not|h, there is i with a_i=f_i(h).    (6)
```

This is the promised grouped covering code.  It is equivalently the positive
CNF with variables `X_(i,b)`, exactly one `b` selected for every `i`, and one
clause

```text
OR_(i:f_i(h)!=*) X_(i,f_i(h))                            (7)
```

per distinct clock word.

The old clocks `7|h` reduce to the parent `1/q` grid and must be checked on the
parent; they are not clauses in `(7)`.  This is exactly the future-multiple
discipline required by MISTAKE-194: no branch is pruned by pretending the old
and new clocks are interchangeable.

## The divisor sidecar is not optional

For the first `x7` lift, the new Sungkawichai--Trakulthongchai gcd alternative
fires if at least twelve child coordinates are divisible by `7`.  Since
`gcd(q,7)=1`, coordinate `i` is divisible by `7` at the unique digit

```text
g_i=-v_i q^(-1) mod 7.                                  (8)
```

Thus the exact augmented code adds thirteen negative clauses forbidding any
twelve of the literals `X_(i,g_i)` simultaneously.

At the target control `q=382=2*191`, the clock code alone accepts

```text
g=(5,3,1,6,4,2,0,5,3,1,6,4,2).
```

All thirteen resulting children are divisible by `7`, so this is not an
improper tuple: it is killed by the gcd condition.  Both independent solvers
find the *augmented* code infeasible.  Dropping `(8)` would create an exact but
false survivor.  This is the cleanest lesson of the computation.

## Finite-exact controls

The companion performs the following checks.

1. It exhausts all `7^4=2,401` assignments with four free digits at `q=22`
   and independently compares raw modular witness scans with clause
   satisfaction.  It also compares the new-`7` count with the full gcd
   definition at lift modulus `14`.
2. At `q=86=2*43`, both CaDiCaL and OR-Tools CP-SAT return `SAT`.  The exact
   child

   ```text
   (259,518,175,4,91,350,7,266,525,268,441,356,357)
   ```

   has no old or new ansatz witness and does not satisfy the gcd alternative.
   This is an ansatz-improper child, not an LRC counterexample.
3. At `q=382=2*191`, the smallest prime proposed in the local `k=13` cost
   ledger, both engines return `UNSAT` on `1,267` clauses.  Hence every one of
   the `7^13` children above the literal AP parent is proper in the
   `(13,191,14)` ansatz.  This is one parent at one prime, not the whole
   initial sieve or lift tower.
4. Normal and optimized Python produce the same frozen transcript.  The two
   engines are algorithmically independent, but neither result is promoted as
   a kernel-checked UNSAT proof.

The `q=382` result is a concrete proof-engineering win: the canonical AP parent
after one `x2` lift can be discharged without materializing its `x7` children.
It remains a tiny part of the proposed prime-by-prime computation.

Reproduce the frozen transcript (the environment needs `python-sat` with
CaDiCaL 1.9.5 and OR-Tools):

```bash
python3 -B 04-computation/lrc14_x7_lift_fibre_clause_code_20260824.py
python3 -B -O 04-computation/lrc14_x7_lift_fibre_clause_code_20260824.py
```

## A clean scalar hostile

The obvious capacity statistic is

```text
sum_i max_b #{clock words h : f_i(h)=b}.                 (9)
```

At `q=382`, `(9)` is `1,824`, versus only `968` distinct clock words.  The
union-capacity bound therefore has excess `856` and says nothing, even though
the augmented grouped code is infeasible.  This is the exact future-fibre
analogue of the level-one acceptance-test verdict: independent best-column
capacity cannot see incompatibility among choices in different groups.

Nor does boundary erosion alone suffice.  Both the `SAT` hostile `q=86` and
the `UNSAT` target `q=382` have exactly one distinct length-six boundary
clause and otherwise length-thirteen clauses.  Minimum clause length and the
number of shortest clauses are identical.  What distinguishes them is the
full cross-phase incidence of digit labels across all clocks, together with
the gcd sidecar.  This is consistent with THM-4002 and warns against promoting
the reserved THM-4003 boundary-strip lane using only a scalar strip count.

## Primary-source audit

The exact imported claims from
[Sungkawichai--Trakulthongchai, arXiv:2604.23906v1](https://arxiv.org/abs/2604.23906)
are:

- lifting a surviving set by multiplier `c` checks at most `c^k |S|` children;
- their proved `k=11` diagram uses four `x2` and two `x3` lifts;
- the polynomial shortcut for the canonical parent assumes both `k+1` and
  the auxiliary modulus `p` are odd primes;
- Section 7 explicitly names efficient computation of `I(13,p,1)` as the
  primary bottleneck in extending to `k=13`.

The statement that a future `k=13` implementation must confront a costly
`x7` fibre is a **repo inference** from `14=2*7` and the general lifting
architecture.  The paper does not publish a `k=13` lifting diagram, assert
that one particular unmitigated `x7` step is unavoidable, or price a completed
fourteen-runner calculation.  The local finite-check ledger is useful as a
cost model, but its “two walls” language should not be attributed verbatim to
the source.

Rosenfeld's earlier finite check supplies the covering/backtracking ancestor;
see [arXiv:2509.14111](https://arxiv.org/abs/2509.14111).  The current
`I(k,p,l)` and lift/project formalism used above is the later S--T definition.

## Typed connections

### Covering code

```text
source: one literal first-x7 child fibre over a fixed parent
target: grouped positive CNF / coordinate-hyperplane covering of clock words
map: child -> digit word a; clock -> forbidden word f(h)
preserved: existence/nonexistence of a genuinely new ansatz witness, exactly
destroyed: old clocks, gcd properness, integer speeds beyond the residue child
sidecars: parent-improper flag, equation (8), endpoint-* flags, deck modulus q
hostile: q=86 AP fibre
test: run grouped SAT before enumerating any x7 children of a real survivor
```

### Conductor/normalization theme

There is a useful analogy, but no scheme-theoretic transfer.  Projection from
the `7^13` digit sheets to one parent is the finite covering.  The endpoint set
`2u_i=q` is where the generic one-forbidden-sheet law degenerates to no
forbidden sheet; retaining `*` is the analogue of retaining the conductor
boundary after normalization.  The all-divisible ghost shows that the sheet
code without the divisor sidecar is not the original properness predicate.

```text
preserved: exact sheet incidence at new clocks
destroyed: gcd/divisor ownership and equality provenance
sidecar: special digits g_i plus endpoint-* labels
counterindication: this is a finite-set analogy, not an algebraic
normalization, class-group, or Hopf theorem
```

### THM-4009 square-norm Graver relation

The incoming square-norm-`<=195` relation can in principle add a linear
syndrome to the digit code: for a fixed integer relation `c`, equation `(1)`
gives

```text
c dot x = c dot v + q (c dot a).                         (10)
```

But this is not yet a lawful pruning rule.  THM-4009 concerns an actual global
hypothetical counterexample; an ansatz child is only a residue representative,
and further `14p` deck integers alter `(10)`.  The relation is also existential,
not a fixed `c` attached to every parent.  A transfer would require the global
deck word and a finite coefficient-bank disjunction.  The cheapest honest test
is therefore to compile `(10)` only after a real lift survivor is paired with
its deck constraints; no relation clause is used in the present computation.

## Concept board and next operations

| object | retained predicate | obstruction/loss | next exact test |
|---|---|---|---|
| level-one uncovered-set DFS | improper cover modulo `p` | no future digits | export actual `I(13,191,1)` parents |
| first-`x7` clock code | all new lift clocks | grouped CNF can remain hard | benchmark parentwise SAT versus `7^13` enumeration |
| gcd conductor sidecar | new divisor properness | old divisors live at parent | carry one divisor flag per lift prime |
| endpoint `*` word | threshold equality exactly | scalar length histogram blind | retain labelled word, not strip count |
| Graver syndrome | global short relation | deck and existential coefficient lost | add only after deck realization |

Ranked continuation:

1. Obtain even a small authentic survivor bank after the initial and `x2`
   sieves at `p=191`; compare parentwise grouped-CNF survival rates with the
   ledger's flat `7^13 |S|` cost.
2. Compile clock clauses incrementally during the `x2` tower.  A parent should
   be projected only together with the minimal clause antichain and divisor
   flags; projecting just the residue tuple repeats the scalar-quotient loss.
3. Learn cross-parent clauses only after quotienting by the paper's lawful
   permutation/sign/unit equivalence.  A clause on a chosen representative
   must be transported with its coordinate and digit action.
4. Search for small independently checkable UNSAT cores, not only solver
   verdicts.  A useful mathematical theorem would identify a bounded family
   of physical clocks whose labelled forbidden words plus gcd clauses already
   kill a parent class.
5. Test whether the boundary `*` clauses generated by even `q` systematically
   accelerate the later `x7` solve.  The `q=86/q=382` hostile proves that
   their count alone cannot be the invariant.

## Honest frontier

LRC(14) remains **OPEN**.  The initial `I(13,p,1)` computation remains the only
extension bottleneck explicitly stated by the current primary source.  The
new contribution is a faithful and much smaller object for a possible later
`x7` lift, an exact canonical-parent discharge at `p=191`, and two sharp
failure boundaries: scalar capacity is blind, and clock incidence without the
gcd conductor produces a false survivor.
