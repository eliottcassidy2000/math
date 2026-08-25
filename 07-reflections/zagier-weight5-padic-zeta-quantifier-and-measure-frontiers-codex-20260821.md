# Zagier weight five, p-adic zeta quantifiers, and measure frontiers

**UPDATE -- 2026-08-25.** The singleton `OPEN` ledger below is a correct
literature/source snapshot for this August 20--21 audit, but it is no longer a
complete claim census. Christopher D. Long's August 24 GitHub research draft
now claims 22 individual p-adic-zeta irrationalities. Its terminal interval
certificate is `FINITE-EXACT`; the global theorem remains an **EXTERNAL
PREPRINT CLAIM / UNDER SPECIALIST AUDIT**, not proved repository canon. See
[`padic-holonomy-matching-logic-and-coordinate-depth-frontiers-codex-20260825.md`](padic-holonomy-matching-logic-and-coordinate-depth-frontiers-codex-20260825.md)
and THM-4089/THM-4091. The Zagier weight-five and rejected same-lune measure
conclusions below are unchanged.

**Session date:** 2026-08-20--21  
**Status:** literature statements `CITED`; elementary identities `PROVED`;
rank, recurrence, and coefficient audits `FINITE-EXACT`; the lune exclusion
`FINITE-INTERVAL`; proposed singleton irrationalities and the exceptional-prime
bridge `OPEN` or `RESERVED` as stated below.

Exact replays:

- [`mzv_weight5_rank_transfer_audit_codex_20260820.py`](../04-computation/mzv_weight5_rank_transfer_audit_codex_20260820.py),
  with [stored output](../05-knowledge/results/mzv_weight5_rank_transfer_audit_codex_20260820.out);
- [`zeta13_legacy_modular_approximant_scout_codex_20260820.py`](../04-computation/zeta13_legacy_modular_approximant_scout_codex_20260820.py),
  with [stored output](../05-knowledge/results/zeta13_legacy_modular_approximant_scout_codex_20260820.out);
- [`padic_zeta5_holonomy_energy_scout_codex_20260820.py`](../04-computation/padic_zeta5_holonomy_energy_scout_codex_20260820.py),
  with [stored output](../05-knowledge/results/padic_zeta5_holonomy_energy_scout_codex_20260820.out);
- [`padic_zeta5_lune_jensen_interval_codex_20260821.py`](../04-computation/padic_zeta5_lune_jensen_interval_codex_20260821.py),
  with [stored output](../05-knowledge/results/padic_zeta5_lune_jensen_interval_codex_20260821.out);
- [`zeta2_5_apery_gcd_candidate_audit_codex_20260820.py`](../04-computation/zeta2_5_apery_gcd_candidate_audit_codex_20260820.py),
  with [stored output](../05-knowledge/results/zeta2_5_apery_gcd_candidate_audit_codex_20260820.out).

## 1. Outcome first

Three attractive announcements separate into three different mathematical
types, and none currently proves the strongest advertised conclusion.

1. The Zagier dimension conjecture for actual real multiple zeta values is
   still open.  At weight five it is **exactly equivalent** to

   ```text
   zeta(5) not in Q * zeta(2)zeta(3),
   ```

   but motivic dimension, formalization of the statement, or reduction at
   almost every prime does not supply the missing injectivity of the real
   period map.

2. The number `543606522303979` is prime, and the publicly displayed
   denominator-`p` rational lies strictly above
   `zeta(5)/(zeta(2)zeta(3))`.  No public theorem, proof, code, determinant, or
   certificate explaining why this prime is exceptional was located.  The
   announcement is therefore `RESERVED`, not a theorem available for use.

3. Under the standard diagonal Kubota--Leopoldt convention, all six queried
   individual values

   ```text
   zeta_5(3), zeta_3(5), zeta_7(3),
   zeta_2(7), zeta_3(7), zeta_2(9)
   ```

   remain open as singleton claims in the literature audited here.  The
   published theorems say **at least one member of a set is irrational**.  An
   existential disjunction is not a coordinatewise theorem.  The same issue
   leaves `zeta_13(3)` open.

There is nevertheless a substantive new negative result about the proposed
irrationality-measure candidate `13.806686`.  Holding the published
Calegari--Dimitrov--Tang lune and arithmetic data fixed would require energy

```text
I <= 2.2100946773550234.
```

Two explicit `Gamma_0(2)` companion preimages, combined with Jensen's formula
and Arb interval arithmetic, prove

```text
I_BC > 2.4831774271554007.
```

Since the rearrangement energy dominates this Bost--Charles/Jensen energy,
`13.806686` is `REFUTED` for that unchanged template.  This is not a global
lower bound on every possible method.

The decimal also has a striking finite-arithmetic fingerprint.  In the
Lai--Sprang--Zudilin recurrence, the exact coefficient gcd at `n=320` gives a
naive frozen proxy `13.806057...`, only `0.000629` away.  That shortcut fails
twice: the gcd rate later falls, and its factor `2^63` must be charged against
the same `2`-adic decay.  The corrected finite proxies are `16.42674` (full
gcd charged) and `16.63138` (odd gcd only).  This is strong evidence for a
finite-gcd/asymptotic or local/global conflation, not proof of the decimal's
provenance.

## 2. Inheritance pass and portfolio

### Anchor: actual weight-five MZVs

- **Closest proved mechanism:** Terasoma and Deligne--Goncharov's upper
  dimension bound, followed by Brown's motivic Hoffman `{2,3}` spanning
  mechanism.
- **Canonical hostile:** an integral matrix `diag(1,p)` has full rational
  rank and drops rank only modulo `p`.  An exceptional reduction can measure
  torsion in a chosen lattice without producing a kernel of the real period
  map.
- **Corrected near miss:** motivic/formal independence is not independence of
  actual periods unless the relevant period map is injective.  A Lean theorem
  whose body is `by sorry` is a formalized conjecture statement, not a proof.
- **Least-used sidecar:** saturation and Smith normal form of the integral
  presentation connecting a claimed prime exception to actual periods.

### Niche: isolate `zeta_13(3)`

- **Closest proved mechanism:** Lai--Lupu--Sprang's linear form proving that
  some odd `13`-adic zeta value through weight `19` is irrational.
- **Canonical hostile:** irrational carrier values can cancel in their
  trivial-character sum.
- **Corrected near miss:** Calegari's genus-zero analytic continuation exists
  at `p=13`, but its decay/height exponent is only `0.59894`, below the
  irrationality gate.
- **Least-used sidecar:** a character-DFT selector on the six reflection
  orbits of normalized `13`-adic Hurwitz values.

### Wildcard: the `zeta_2(5)` measure

- **Closest proved mechanism:** the Calegari--Dimitrov--Tang holonomy measure
  and the independent Lai--Sprang--Zudilin Apéry recurrence.
- **Canonical hostile:** mixing the energy of one conformal map with the
  derivative of another fabricates a much better constant while violating the
  proposition's types.
- **Corrected near miss:** a finite common divisor is neither an asymptotic
  liminf nor a free saving when it contains powers of the distinguished prime.
- **Least-used sidecar:** explicit modular companion preimages turn Jensen
  collision mass into a rigorous lower bound on conformal energy.

## 3. The real Zagier conjecture and the exact weight-five gate

Let `Z_w` be the `Q`-span of convergent real MZVs of weight `w`.  Define

```text
d_0=1, d_1=0, d_2=1, d_w=d_(w-2)+d_(w-3).
```

Zagier's actual-number conjecture is

```text
dim_Q Z_w = d_w.
```

The proved general direction is the upper bound `dim_Q Z_w <= d_w`.
Brown's Hoffman theorem gives spanning by words in `{2,3}`; it does not add a
lower bound for the associated real numbers.

At weight five, the two Hoffman words are `(2,3)` and `(3,2)`, and Euler's
reductions give

```text
zeta(2,3) =   9/2 zeta(5) - 2 zeta(2)zeta(3),
zeta(3,2) = -11/2 zeta(5) + 3 zeta(2)zeta(3),
zeta(4,1) =       2 zeta(5) -   zeta(2)zeta(3).
```

The first two rows have determinant `5/2`; after clearing denominators their
Smith factors are `[1,10]`.  Consequently

```text
Z_5 = span_Q{zeta(5), zeta(2)zeta(3)},
d_5=2,
dim_Q Z_5=2  iff  zeta(5)/(zeta(2)zeta(3)) is irrational.
```

This equivalence is unconditional.  Its right-hand side is not presently
known.

### Why the nearby realizations lose the needed column

There is a useful common loss mechanism.  The finite/symmetric MZV model is a
model of the quotient by the ideal generated by `zeta(2)`.  Its predicted
weight-five dimension is

```text
d_5-d_3 = 2-1 = 1.
```

Likewise, the standard `p`-adic **multiple-zeta** realization kills the
motivic `zeta(2)` period.  These are not the one-variable Kubota--Leopoldt
values discussed later.  Both quotient functors erase exactly the second
column `zeta(2)zeta(3)` required to distinguish weight-five real dimension
one from two.

This yields a precise type diagram:

```text
motivic/formal weight-five space (two columns)
        | real period map                 | finite or p-adic-MZV realization
        v                                 v
actual real periods (injectivity open)    zeta(2)-quotient (one column)
```

An all-primes result in the right branch cannot by itself prove injectivity in
the left branch.

## 4. The exceptional-prime announcement

The public source located is Christopher Long's 19 August 2026 X post.  It
reports that an unspecified theorem holds at every prime except

```text
p0 = 543606522303979.
```

No exact theorem statement or audit artifact accompanies the report.

The following pieces are reproducible.

1. `p0` is prime (`FINITE-EXACT`; SymPy and an independent Pocklington chain).
2. The displayed rational is

   ```text
   285075330345178 / p0.
   ```

3. Arb at 80 decimal digits certifies

   ```text
   285075330345178/p0 - zeta(5)/(zeta(2)zeta(3))
   = 4.7388330086937421767... * 10^-17 > 0.
   ```

This excludes one rational with denominator exactly `p0`.  It does not
exclude other denominators, powers `p0^e`, or rationality itself.

Nor can `p0` come from the canonical weight-five basis change above: that
cleared determinant is `10`, with primes only `2` and `5`.

### Necessary audit packet

Before an exceptional-prime result can cross to the real weight-five claim,
one needs all of the following.

1. The exact integral module and matrix whose reduction is being tested.
2. Its characteristic-zero rank and its saturated lattice.
3. Smith factors showing whether `p0` survives saturation rather than a poor
   choice of generators.
4. A source-to-target map from that module to the two actual weight-five
   periods.
5. The direction in which rank information travels.
6. An injectivity statement for the real period map, or a denominator theorem
   reducing rationality to the single displayed rational.

Until that packet exists, the exceptional-prime bridge is `RESERVED`.

## 5. The `p`-adic singleton ledger

Fix the diagonal Kubota--Leopoldt convention

```text
zeta_p^Delta(s)=L_p(s,omega^(1-s))
               = lim zeta(k),
```

where negative integers `k` tend `p`-adically to `s` in the congruence class
`s mod p-1`.  This convention matters: some sources write `zeta_p(s)` for the
principal branch `L_p(s,1)`, and Diamond--Hurwitz values have an additional
variable.

| Individual cell | Audited status | Strongest nearby theorem |
|---|---|---|
| `zeta_5(3)` | `OPEN` | some odd value `3,...,13` is irrational |
| `zeta_3(5)` | `OPEN` / no singleton theorem found | direct modular exponent below one |
| `zeta_7(3)` | `OPEN` | some odd value `3,...,13` is irrational |
| `zeta_2(7)` | `OPEN` | one of `7,9,11,13` is irrational |
| `zeta_3(7)` | `OPEN` / no singleton theorem found | no applicable small singleton result |
| `zeta_2(9)` | `OPEN` | one of `7,9,11,13` is irrational |
| `zeta_13(3)` | `OPEN` | some odd value `3,...,19` is irrational |

The exact failure is logical, not numerical:

```text
exists i in I: Irr(zeta_p(i))
```

does not imply `Irr(zeta_p(i0))` for a named `i0`.  A Boolean countermodel may
set the named cell false and use another member of `I` as the witness.

The positive singleton controls are Calegari's `zeta_2(3)` and `zeta_3(3)`,
and Lai--Sprang--Zudilin's `zeta_2(5)`.

### Hurwitz carriers and the equality boundary

Beukers proves irrationality of normalized Hurwitz carriers

```text
u_a = omega(a)^(-2) H_p(3,a,p).
```

Reflection gives `u_(p-a)=u_a`, while the diagonal zeta value is twice the sum
over reflection orbits.  At `p=3` there is one orbit, so carrier irrationality
descends to the sum.  At `p=5,7,13` there are respectively `2,3,6` orbits.
Irrational summands can cancel; componentwise irrationality does not imply an
irrational sum.

## 6. A quantified attempt at `zeta_13(3)`

### 6.1 Legacy genus-zero modular route

For

```text
f=(Delta(13 tau)/Delta(tau))^(1/12),
```

the direct exponent at odd weight three is

```text
theta_p = 4 log(p)/(2 log(p)+p-1),
theta_13 = 0.5989409278 < 1.
```

In natural-log units, the analytic decay is `A=log(13)=2.564949...`,
the closest-branch growth is `G=(1/2)log(13)`, and the Eichler-integral
denominator exponent is `d=3`.  Thus

```text
B-A = G+d-A = 1.717525...
```

must be recovered.  Perfectly cancelling the closest branch is still
insufficient because `3>log(13)`; even then the denominator exponent must fall
by more than `0.435051`.

The exact `f`-coefficient scout reproduces Calegari's printed `p=2` rows.
At `p=13`, the tested denominators through `n=80` show no saving remotely
large enough.  This is a finite hostile, not an asymptotic exclusion.

### 6.2 Refined holonomy gate

For the packet `1,H,H',H''`, the unrefined denominator column is
`(0,3,3,3)` and

```text
tau = 45/16 = 2.8125 > log(13).
```

Without changing the analytic radius, the effective denominator slope must
satisfy

```text
b_eff < (16/15)log(13) = 2.735946... .
```

Without denominator saving, the usable `13`-adic radius must instead grow
from `13` past

```text
exp(45/16) = 16.651... .
```

These are concrete construction targets.

### 6.3 Lai--Lupu--Sprang support gate

At `p=13`, the smallest working template parameter is `r=8`.  The resulting
linear form involves all nine odd values from `3` through `19`; the theorem
proves their disjunction.  Numerically, `r=7` misses its scalar decay gate and
`r=8` clears it, but improving that scalar gate does not isolate the first
coordinate.  A selector problem remains.

### 6.4 Generated directions

Each direction is stated with its lost information and cheapest decisive
test.

1. **Hurwitz DFT selector.**  Source: the six reflection-orbit carriers.
   Target: their trivial-character sum `zeta_13(3)`.  Map: character Fourier
   transform.  Preserved: exact carrier relations.  Lost: irrationality under
   arbitrary summation.  Sidecar: a simultaneous Padé/Casoratian determinant.
   Cheapest test: exact minors through order `30`, including the height cost
   of the cyclotomic norm.

2. **Multi-template annihilator.**  Source: shifted, differentiated, or
   moment-weighted Lai--Lupu--Sprang forms.  Target: a row supported only at
   weight `3`.  Map: exact row combinations.  Preserved: coefficient
   rationality.  Lost: decay and nonvanishing under a large combination.
   Sidecar: coefficient height.  Cheapest test: rational rank and Smith form
   of the unwanted columns `5,...,19` before asymptotics.

3. **Atkin--Lehner/`U_13` continuation.**  Source: the genus-zero modular
   function.  Target: an overconvergent disk of radius greater than `16.651`.
   Preserved: modular functional equations.  Lost: control at intervening
   poles.  Sidecar: divisor/Newton polygon.  Cheapest test: exact pole divisor
   and `U_13` slopes.

4. **Denominator thinning.**  Source: the exact `f`-coefficient recurrence.
   Target: effective exponent below `2.735946`.  Preserved: the analytic
   packet.  Lost: none if the common divisor is prime to `13`; a `13`-primary
   divisor also weakens local decay.  Sidecar: prime-by-prime valuation
   classifier.  Cheapest test: `300--1000` exact rows with `p=3` positive and
   `p=5` hostile controls.

5. **Character sidecars on `Gamma_1(13)`.**  Source: character-projected
   Eisenstein components.  Target: a larger differential module with some
   zero-denominator columns.  Preserved: `13`-adic analytic continuation.
   Lost: independence if a sidecar is rational over the existing module.
   Sidecar: exact differential-module rank.  Cheapest test: enumerate
   characters and recompute the holonomy denominator functional.

The selector directions are conceptually prior to further scalar
optimization: a better nine-value disjunction is still not a singleton.

## 7. The irrationality measure of `zeta_2(5)`

### 7.1 Published packets

Calegari--Dimitrov--Tang prove an effective bound reported as

```text
mu(zeta_2(5)) < 19.7439,
```

whose displayed constants replay to `19.743897509211...`.  The relevant data
are

```text
m=6, gamma=1/6, L=12 log(2), tau=175/36,
phi'(0)=14/29, I<=3.92881.
```

The independent Apéry-like result of Lai--Sprang--Zudilin gives

```text
mu(zeta_2(5)) <= 16 log(2)/(8 log(2)-5)
               = 20.342651...
```

and supplies a different order-five recurrence useful for arithmetic
experiments.

### 7.2 Why `13.806686` is impossible in the unchanged lune

Solving the CDT equality for the energy while holding every other datum fixed
shows that `13.806686` requires

```text
I <= 2.2100946773550234.                         (7.1)
```

The published lune is

```text
psi(z)=(29/63)(1+z-sqrt(1-(82/841)z+z^2)),
psi^(-1)(q)=q(63q-58)/(2(29q-14)).
```

Put `q=psi(exp(it))=exp(2 pi i tau)`.  The two matrices

```text
[[1,0],[+2,1]], [[1,0],[-2,1]] in Gamma_0(2)
```

give companion points

```text
q_+ = exp(2 pi i tau/(1+2 tau)),
q_- = exp(2 pi i tau/(1-2 tau)).
```

The Hauptmodul `x(q)=q product_n(1+q^n)^24` is invariant under these
transforms.  Pulling `q_+` and `q_-` back through `psi^(-1)` yields explicit
additional preimages of the same boundary value.  Jensen's formula therefore
counts

```text
max(0,-log|z_+|)+max(0,-log|z_-|)
```

inside the energy.  Arb interval subdivision, discarding every uncertain
cell, proves the normalized contributions

```text
A_- > 1.6057464673146107,
A_+ > 1.6056694602120054.
```

Including `log|psi'(0)|=log(14/29)` gives

```text
I_BC > 2.4831774271554007,
```

which exceeds (7.1) by `0.2730827498`.  The inverse-sheet condition is
certified on every retained cell using the exact quadratic

```text
63q^2-(58+58z)q+28z=0,
```

whose two roots multiply to `4z/9`.  Uncertain cells contribute zero, so the
enclosure is one-sided in the safe direction.

This closes the **same-lune, same-arithmetic** route.  It does not close other
domains or local systems.

### 7.3 The finite-gcd fingerprint

In the Lai--Sprang--Zudilin packet, obtaining `13.806686` by a pure
Archimedean common-divisor improvement would require the asymptotic saving

```text
sigma = 0.2580822869982988.
```

At the single index `n=320`, the exact recurrence gives

```text
log(gcd(A_320,B_320))/320 = 0.258118873235879,
```

whose naive frozen measure proxy is `13.80605717265`.  This near match is too
close to ignore as provenance evidence, but it is not mathematics of the
required asymptotic kind.

The rate subsequently falls to `0.1442417` at `n=1000`.  Moreover

```text
gcd(A_320,B_320)
=2^63 3^6 5^3 7^3 11^2 13^2 17^2 19 23.
```

The factor `2^63` reduces the `2`-adic decay along with the coefficient
height.  Charging that local cost changes the finite proxy to `16.42673966`;
using only the odd part gives `16.63138364`.  Hence the attractive decimal can
arise only by making both an asymptotic extrapolation and a local/Archimedean
accounting error in this packet.

### 7.4 Live measure directions

1. Apply CDT holonomy to the distinct order-five LSZ local system and recompute
   `gamma`, `tau`, radius, and energy from scratch.
2. Search multi-parameter Rhin--Viola families for a proved prime-saving step
   function of total Chebyshev weight at least `0.258082`.
3. Search conformal domains only as complete packets: derivative, energy,
   branch geometry, and arithmetic must travel together.
4. Automate modular-companion enumeration.  Each known deck transformation
   supplies a Jensen collision lower bound and can cheaply rule out an
   overoptimistic energy target before heavy optimization.
5. Separate the prime-to-`2` gcd from the `2`-primary gcd in every recurrence
   search.  The first is a possible free height saving; the second carries a
   local decay price.

## 8. Cross-lane synthesis

The live concept board was:

```text
period injectivity | lattice saturation | quotient/character loss
existential support | local-vs-global height | collision energy
selector determinants
```

The strongest merges are these.

1. **Saturation and prime-primary accounting are the same warning at two
   scales.**  In the MZV claim, a prime can be an artifact of an unsaturated
   integral lattice.  In the measure claim, a `2`-primary gcd is not a free
   global saving because it changes local decay.  Both require prime-by-prime
   transport through the actual map, not a scalar determinant or gcd.

2. **Quotient loss and existential loss both erase coordinates.**  Finite or
   `p`-adic-MZV realization removes the `zeta(2)zeta(3)` column.  A multi-value
   irrationality theorem retains a row but forgets which column witnesses it.
   Both problems demand a selector sidecar before a statement can be pulled
   back to a named real or `p`-adic value.

3. **Jensen preimages are analytic analogues of hidden carrier orbits.**  In
   the Hurwitz problem, six orbit values may cancel in their sum.  In the
   holonomy problem, modular deck transforms create extra preimages that add
   unavoidable energy.  Character projection is needed to control the former;
   collision counting exposes the latter.

4. **A decimal is not a theorem object.**  The exceptional prime and the
   measure candidate both come with striking numerical sidecars.  Each becomes
   meaningful only after identifying the exact module, recurrence, map, and
   asymptotic or injectivity obligation that would let the number cross types.

## 9. Procedural next-session generator

The following loop generates directions without merely resampling constants.

1. Choose one object from each column:

   ```text
   object:       motivic word | Hurwitz carrier | linear form | conformal map
   representation: lattice | character DFT | recurrence | modular cover
   invariant:    SNF | selector rank | valuation polygon | Jensen collision
   operation:    saturate | project | shift/differentiate | deck transform
   quotient:     mod p | reflection orbit | support disjunction | branch sheet
   scale:        weight | template order | recurrence index | analytic radius
   ```

2. State the source, target, map, preserved predicate, destroyed information,
   needed sidecar, and cheapest hostile test.
3. Run the hostile before optimizing constants.
4. If it survives, demand either an exact determinant/divisor identity or a
   validated interval certificate.
5. Revisit every other lane: ask whether the new sidecar repairs a coordinate
   loss, a prime-primary accounting loss, or an asymptotic gap elsewhere.

The highest-value next computations are therefore:

1. an exact small-order DFT/Casoratian census for the six `p=13` Hurwitz
   reflection orbits;
2. an exact row-span/Smith census for several Lai--Lupu--Sprang templates,
   targeting annihilation of weights `5,...,19`;
3. a prime-classifier audit for the LSZ recurrence rather than another raw gcd
   plot;
4. automatic enumeration of `Gamma_0(2)` companion preimages for alternative
   domains, so energy impossibility is detected before expensive searches.

## 10. Reference surface

- F. Brown, *Mixed Tate motives over Z*:
  <https://arxiv.org/abs/1102.1312>.
- T. Terasoma, *Mixed Tate Motives and Multiple Zeta Values*:
  <https://arxiv.org/abs/math/0104231>.
- P. Deligne and A. Goncharov, *Groupes fondamentaux motiviques de Tate
  mixte*: <https://arxiv.org/abs/math/0302267>.
- F. Calegari, *Irrationality of certain p-adic periods for small p*:
  <https://arxiv.org/abs/math/0408214>.
- F. Beukers, *Irrationality of some p-adic L-values*:
  <https://arxiv.org/abs/math/0603277>.
- L. Lai, *On the irrationality of certain 2-adic zeta values*:
  <https://arxiv.org/abs/2304.00816>.
- L. Lai, C. Lupu, and J. Sprang, *On the irrationality of certain p-adic
  zeta values*: <https://arxiv.org/abs/2505.23088>.
- L. Lai, J. Sprang, and W. Zudilin, *A note on the irrationality of
  zeta_2(5)*: <https://arxiv.org/abs/2505.05005>.
- F. Calegari, V. Dimitrov, and Y. Tang, *Arithmetic holonomy bounds and
  effective Diophantine approximation*: <https://arxiv.org/abs/2510.04156>.
- The exceptional-prime announcement (no theorem artifact located):
  <https://x.com/octonion/status/2090089114154504208>.

No claim in this reflection promotes the actual Zagier conjecture, an
individual `zeta_13(3)` irrationality theorem, any of the six queried singleton
values, or the decimal `13.806686` to proved status.
