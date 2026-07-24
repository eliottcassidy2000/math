---
source: opus-2026-07-24-puzzle-atlas
status: SYNTHESIS. The proved outputs are THM-2159 and THM-2162. Every other
  bridge below is explicitly typed as proved input, live proposal, or no-go.
related:
  - THM-2022
  - THM-2047
  - THM-2051
  - THM-2054
  - THM-2145
  - THM-2159
  - THM-2161
  - THM-2162
  - THM-2163
  - MISTAKE-245
---

# Puzzle atlas -- expose first, choose scale second

## 1. What the supplied sources actually were

The source list was intentionally heterogeneous. The first validity task was
to identify the objects rather than force a common vocabulary.

- `arXiv:2506.24088` is Brittenham--Hermiller, *Unknotting number is not
  additive under connected sum*, not the Shunia paper. Its reusable move is
  a common intermediate knot which beats the sum of separately optimized
  local costs.
- Shunia's quotient-ring paper is `arXiv:2404.00332`. Its Section 6 states
  the finite modular root formula as Conjecture 6.1 after proving only a limit
  representation. THM-2159 now proves the finite formula.
- The Annals paper is Ben Green, *Roth's theorem in the primes*. The relevant
  technique is a restriction estimate plus a structured/uniform split, not
  “primes make every sparse set random.”
- OEIS A005670 is Mrs Perkins's quilt (minimum orders of primitive squared
  square dissections). OEIS A360665 is the array
  `T(n,k)=k((2n-1)k+1)/2`, with column increment `k^2`. Neither sequence
  arrives with a proved map to LRC.
- the Zbl item `1260.11001` is Bugeaud's book on distribution modulo one and
  Diophantine approximation.

The remaining named topics -- Bonse, Siegel--Walfisz,
Elliott--Halberstam, Baker, Ramanujan sums, partial cubes, Robbins,
Frobenius--Koenig, Steinhaus, Worpitzky, Asano contraction, graceful trees,
and N-body regularization -- are technique prompts. None is a theorem about
LRC merely because an analogous noun occurs in the repo.

## 2. Inheritance board

The closest proved mechanisms at the start were:

| lane | proved mechanism | hostile control | least-used sidecar |
|---|---|---|---|
| LRC anchor | THM-2051/2054 relation-free exit and relative Fejer line; THM-2047 phase-height carrier | AP/GW isolated unit witnesses; THM-2161 fixed-bank mimic | signed endpoint primitive with Euler/BV split |
| Shunia niche | quotient-ring recurrence and roots-of-unity filter | `(a,n)=(3,2)` and perfect-power exclusion | explicit algebraic distance to adjacent integer |
| adaptive-prime wildcard | THM-2022 NC2 Frobenius lowest-face prime | ramified/bad primes and `A_0=0` | prime chosen after the face/carry is exposed |

The live concepts were:

```text
common intermediate;
positive smoothing;
signed endpoint current;
arithmetic quantum;
adaptive good prime;
owner-labelled topology;
relation carry.
```

Every proposed connection below was compared with all seven.

## 3. Two proved outputs

### 3.1 Shunia: contraction plus an arithmetic quantum

THM-2159 rescales the coefficients of

```text
(1+y)^t mod (y^n-a)
```

by `rho^r`, `rho=a^(1/n)`. The recurrence becomes a lazy directed random
walk on the `n`-cycle. Its nonconstant Fourier modes contract by an explicit
factor, and `K=2a^n` makes the last-coefficient ratio lie within half the
algebraic distance

```text
min(rho-floor(rho),ceil(rho)-rho)>=1/(n a^(n-1))
```

of `rho`. Kronecker evaluation at `X=a^K` is lossless and the lower digits
perturb the ratio by less than the other half. Natural division therefore
extracts the exact integer root.

The mechanism is:

```text
dominant spectral mode + explicit contraction + arithmetic quantum
  => finite exact extraction.
```

Baker's theorem is not needed: the relevant separation is between an
algebraic root and adjacent integers, and the factorization of `d^n-a`
already gives the required effective gap. This is a useful assumption
challenge: invoking a stronger transcendence tool would hide the elementary
quantum which makes the formula finite.

### 3.2 LRC: integrate the comb before paying its boundary

THM-2162 uses the continuous periodic primitive `H` of the centered danger
indicator. For a good set with positive intervals `[l_i,r_i]`,

```text
epsilon_W
 =1/W sum_i [H(Wr_i)-H(Wl_i)].
```

At danger density `1/7`, `osc(H)=6/49`. Isolated weak-safe points contribute
one left and one right term at the same phase and therefore cancel exactly.
This distinction halves the binding-core tail: its fourteen Euler entries
are eight positive BV intervals plus six isolated unit witnesses.

The result proves the full infinite one-swap neighborhood of the sharp
drop-6 `7/858` core. The unique proper-neighbor minimum is

```text
3859/420420=7/858+1/980
```

at replacement `10->20`. It also supplies MISTAKE-245: the incoming claim
that `426/35035` was the global second value came from a bank ending at 18
and missed this structured speed-20 row.

This is the unknotting paper's legitimate transfer. Separate interval costs
overcharge. Passing through the common intermediate -- the whole core union,
with its signed endpoint current intact -- reveals cancellation.

## 4. Typed connection ledger

| source | target | map | preserved predicate | destroyed information | status / cheapest test |
|---|---|---|---|---|---|
| nonadditive unknotting | LRC one-far discrepancy | common intermediate knot -> whole core safe union | total post-operation cost/measure | attribution to separately optimized summands | **REALIZED**, THM-2162 |
| Shunia quotient ring | exact integer root | coefficient scaling -> cyclic Markov walk | last-coordinate ratio | raw modular appearance | **REALIZED**, THM-2159 |
| Green restriction | LRC relation dichotomy | uniform Fourier mass -> relation-free branch | structured/uniform split | owner labels if scalarized | **ALREADY REALIZED**, THM-2051/2054 |
| Ramanujan sums | endpoint cocycle | group owner endpoints by exact period after integrating | signed harmonic current | individual boundary owner unless retained | **LIVE**, exact finite form is THM-2162; arbitrary-body form THM-884 |
| partial cube | phase-height cells | wall cells -> even-cycle/Theta carrier | adjacency and component topology | interval length and signed current | **NO-GO alone**; THM-2162 shows Euler/BV split |
| Robbins | safe-cell graph | bridgeless graph -> strong orientation | cycle connectivity | existence/length of a selected safe component | **NO-GO alone**; selected safe paths have bridges |
| Frobenius--Koenig | endpoint cancellation | zero comparison matrix -> Hall-deficient phase block | matching obstruction | signed magnitudes unless matrix is weighted | **LIVE SCOUT**, must retain `H(Wx)` weights |
| Asano contraction | trigonometric safe product | multiaffine contraction | zero-free regions under valid hypotheses | signed/complex channel sums | **NO-GO raw**; positivity is orthogonal to the target |
| Steinhaus | positive safe body | difference set contains a neighborhood | qualitative positive measure | quantitative `7/858` floor | **TOO WEAK** for OPEN-Q-108 |
| graceful trees | bounded relation templates | edge labels -> coefficient differences | tree-supported relation address | existence of a graceful labelling is itself open | **PROMPT ONLY** |
| N-body regularization | wall collisions | time/height blow-up | collision ordering | LRC owner and threshold selection without sidecar | **PROMPT ONLY**; THM-2047 is the faithful existing blow-up |
| Worpitzky triangle | signed descent channels | binomial-basis transform | ordered descent statistics | LRC phase height | **NO DIRECT MAP** |
| quilt dissections | component packing | squared tiles -> safe-component ledger | area partition | comb periodicity/owner phase | **NO DIRECT MAP** |

## 5. Adaptive primes: why Bonse points in one direction only

The NC2 proof THM-2022 and the LRC no-go THM-2161 form a clean pair.

THM-2022 first exposes the lowest balanced face, its nonzero constant term,
and the factorial valuation `A_0`. Only then does it choose an unramified
prime larger than the finite face data. Kummer's carry theorem isolates the
whole minimum-valuation layer.

THM-2161 constructs, for every fixed `Q`, a primitive defect-seven row which
copies the AP modulo every `q<=Q` but escapes at a rational phase whose
reduced denominator exceeds `Q`. A fixed prime/modulus prefix is therefore
provably blind.

Bonse's primorial inequality can help with an **availability lemma**: a
bounded nonzero determinant or carry cannot be divisible by every prime in
a sufficiently long initial bank. Siegel--Walfisz or Elliott--Halberstam can
help only after a residue class and a quantitative modulus range have been
legally produced. They do not turn a fixed bank into an LRC certificate.

The correct architecture is:

```text
expose a face/relation/carry
  -> list its finitely many bad primes
  -> choose a good prime or scale adaptively
  -> descend or isolate.
```

This is the natural interface between THM-2145's `6+7` crossing carry and
the reserved THM-2163 radix descent. The prime must be selected from the
actual carry state, not prescribed before the row is seen.

## 6. Why the raw analogies fail together

Several prompts fail for the same controlled-forgetting reason.

- A partial-cube or Robbins carrier sees whether cells connect, but not the
  signed endpoint current deciding a measure change.
- A Frobenius--Koenig zero submatrix sees failure of a matching, but an
  unweighted matching can pair a large positive `H` value with a large
  negative one and lose the needed sum.
- Asano contraction preserves a zero-free domain under separately affine
  hypotheses, but LRC's exact integral is an alternating, owner-colliding
  channel sum.
- Steinhaus turns positive measure into a qualitative neighborhood of zero;
  it supplies no uniform measure quantum.
- Baker supplies spectacular lower bounds for linear forms in logarithms,
  but neither the Shunia integer gap nor rational LRC endpoint separation is
  naturally such a form.

The recurring error would be to keep the theorem's conclusion while dropping
the coordinate on which its hypothesis acts.

## 7. Next decisive tests

1. **OPEN-Q-108 swap distance two.** THM-2162 is now the exact boundary
   condition. Test the two-swap structured collar, but first derive a
   two-new-comb signed primitive which retains their joint overlap; do not
   sum two absolute one-comb errors.
2. **Relation carry.** Combine THM-2145 with THM-2163. At each radix step,
   choose a prime avoiding the current coefficient/carry minors. The cheapest
   hostile control is THM-2161's AP-mimic family, which defeats every
   preselected prime bank.
3. **Endpoint Hall carrier.** Build the weighted Ferrers comparison matrix
   between left and right endpoint `H`-values. Frobenius--Koenig is legal only
   if the weight threshold implies the target signed sum. Test first on the
   exact `10->20` equality row and on the AP/GW isolated controls.
4. **Shunia exponent sharpening.** The proof uses deliberately coarse
   contraction. Determine the least uniform `K(a,n)` which still beats both
   the algebraic quantum and the Kronecker digit error. This is separate from
   the proved sufficiency of `2a^n`.

## 8. Meta-pattern candidate

Evidence now comes from four distinct threads:

- NC2: expose the lowest balanced face, then choose the prime;
- Shunia: rescale to the spectral walk, then compare with the integer quantum;
- LRC one-swap: form the whole core union, then integrate its signed boundary;
- relation carry: expose the crossing frequency, then choose the radix.

The candidate rule is:

```text
EXPOSE FIRST, CHOOSE SCALE SECOND.
```

Trigger: local costs, fixed moduli, or unsigned norms are losing large
cancellations. Counterindication: the scale is part of the problem statement
and cannot adapt, or the exposure map forgets the target predicate. Evidence:
THM-2022, THM-2159, THM-2161/2162, and THM-2145/2163.

This should be promoted to META-PATTERNS only after THM-2163 is audited:
at present its fourth instance is a live provisional proof, not canon.
