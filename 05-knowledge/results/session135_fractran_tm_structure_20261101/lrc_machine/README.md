# Session 135: an equality-aware LRC certificate machine

**Status:** PROVED elementary machine/certificate construction and pointwise
corollary of inherited canon; FINITE-EXACT tests. **LRC(14) remains OPEN.**
This is a representation simplification, not a new solution or a claim of
mathematical novelty for the harmonic/Farey geometry.

## Useful result first

The harmonic compression for the interval core is exact **as a measure**, but
needs a small equality sidecar to become an exact membership machine. At the
LRC(14) deep well, `{1,...,12}` has twelve positive safe intervals and also
four isolated safe points:

```
3/14, 5/14, 9/14, 11/14.
```

The points contribute no measure or Fourier integral. A machine accepting
only positive intervals loses genuine witnesses. For `{1,...,13}`, all safe
points at threshold `1/14` are the six primitive residues modulo 14: the
safe set has measure zero but is nonempty.

We preserve both strata and supply actual executable machines: a five-fraction
FRACTRAN two-comparison primitive, a 53-fraction verifier for thirteen-speed
rational witnesses, and a twelve-transition unary Turing comparator. Modular
arithmetic preprocessing remains an explicit sidecar, not hidden in the
claim that these are entire universal LRC solvers.

## 1. Inheritance and concept board

Read-only inherited source baseline: `2a18212bfe0e` (session 134). Bounded
startup used AGENTS, START-HERE, CURRENT-FRONTIER, ACTIVE-GUARDRAILS,
RESEARCH-PROTOCOL and targeted canonical files. The coordination peek was
blocked by the source checkout lacking `.machine-id`; no identity invented.

- Closest proved mechanism: THM-819,
  `01-canon/theorems/THM-819-primitive-harmonic-law-for-interval-cores.md`,
  plus THM-826,
  `01-canon/theorems/THM-826-farey-profile-theorem-interval-core-measure.md`.
- Canonical hostile: THM-762,
  `01-canon/theorems/THM-762-small-denominator-signed-pair-deck.md`, and
  THM-764,
  `01-canon/theorems/THM-764-covering-small-period-signed-pair-deck-and-q25-refutation.md`:
  `26*{1,...,12} union {339}` has first good period **27**, not at most 25.
- Corrected near miss: THM-819 explicitly repairs a pointwise Dirichlet
  confinement claim to confinement **up to a finite null boundary**.
- Underused sidecar: all equality points and endpoint owner signs; compare
  THM-4331,
  `01-canon/theorems/THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape.md`,
  and THM-4335,
  `01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md`.
- Finite/unbounded guard: THM-763,
  `01-canon/theorems/THM-763-strict-finite-height-for-tight-lrc-instances.md`.

The five live concepts are: closed safe sets; harmonic/Farey compression;
labelled owner/equality events; prime-exponent comparison; finite reduction
versus an actually discharged finite search. They share exact inequalities,
but measure discards singleton events and unary prime encodings inflate size.

## 2. Exact pointwise completion of the primitive harmonic law

Let `k>=1`, `q=k+1`, `lambda=1/(q+1)`, and `epsilon=1/(q(q+1))`.
For each `a in {1,...,q-1}` coprime to `q`, let

```
u_a = a^(-1) mod q, with 1<=u_a<=q-1,
J_a = [a/q-epsilon/u_a, a/q+epsilon/(q-u_a)].
```

Then the following is an equality of **closed sets**, not merely almost
an equality:

```
G_k = {t in [0,1]: ||jt||>=lambda for every 1<=j<=k}
    = disjoint union_a J_a
      disjoint union {b/(q+1): 1<b<q and gcd(b,q+1)=1}.       (1)
```

There are `phi(q)` nondegenerate components and `phi(q+1)-2` isolated
components. The total measure is the inherited primitive harmonic formula

```
mu(G_k) = 2/(q(q+1)) * sum_{1<=u<q, gcd(u,q)=1} 1/u.
```

### Proof, including the exceptional set

Put the `q+1` circle points `0,t,...,qt` in cyclic order. If a cyclic gap
is smaller than `1/(q+1)`, its two endpoints differ by `jt` with
`1<=j<=q` and `||jt||<lambda`. Goodness for `j<q` forces `j=q`.
Consequently `t=a/q+s` with `|s|<epsilon`. Nonprimitive `a` is impossible:
if `g=gcd(a,q)>1`, speed `q/g<q` has distance at most `|qs|/g<lambda`.

At a primitive `a/q`, the unique residues +1 and -1 have speeds `u_a` and
`q-u_a`. Their exact safe inequalities are
`-epsilon/u_a <= s <= epsilon/(q-u_a)`. Every other residue has distance
at least `2/q`; its change on `|s|<=epsilon` is at most `(q-1)epsilon`,
which leaves it strictly above `lambda`. This yields precisely `J_a`.
These intervals are disjoint because each lies within `epsilon` of a
separate `q`-grid point and `2epsilon<1/q`.

If no cyclic gap is smaller, all `q+1` gaps are equal. Thus
`t=b/(q+1)` with `gcd(b,q+1)=1`, and every such point is good. The points
`b=1,q` are already the outer endpoints of `J_1,J_(q-1)`. For
`1<b<q`, the distance to the nearest `a/q` is

```
min(b,q+1-b)/(q(q+1)) >= 2epsilon,
```

so these points lie outside every `J_a` and are isolated. This proves (1),
including `k=1`, when there are no exceptional points.

**Not an independent canon upgrade:** THM-826's exact Farey-gap remnants
also yield (1) if zero-length remnants are retained. The gain here is the
explicit equality-aware normal form and its executable machine interface.

## 3. A total exact finite algorithm for every supplied row

For positive integer speeds `V=(v_1,...,v_n)` and reduced rational
`delta=c/d` with `0<delta<=1/2`, each single-speed safe set is

```
G_v = union_{0<=k<v} [(dk+c)/(dv), (d(k+1)-c)/(dv)].        (2)
```

Intersect these sorted closed interval lists. Retain an intersection iff
`max(left endpoints) <= min(right endpoints)`; **equality is retained**.
A two-pointer merge exhausts each finite list, so the algorithm halts and
returns exactly `G_V`. Applying the integer residue predicate independently
to every arrangement wall and every open-cell midpoint gives a second
implementation of the same set.

A nonempty result has a least point `a/b` equal to an entering wall of
some speed `v`; hence `b<=d max(V)`. Since `dk+c` is coprime to `d`, its
reduced denominator satisfies **`d|b`**. Every endpoint is accompanied by
all labelled tight owners and signs, so a tie never becomes an arbitrary
winner. The endpoint denominator bound is **row-dependent**; it does not
revive the false universal `q<=25` claim, and it concerns the first safe
wall, not the smallest denominator of any rational witness.

The certificate predicate, requiring no real arithmetic, is

```
r_i = (a v_i) mod b, 0<=r_i<b,
d r_i >= c b  AND  d(b-r_i) >= c b    for every i.          (3)
```

A positive witness has polynomial binary encoding length in the input
magnitudes, and (3) is an ordinary polynomial-time integer verifier. The
full interval construction has a magnitude-dependent, pseudo-polynomial
cost; no polynomial-time algorithm for its binary input is asserted here.

## 4. Actual FRACTRAN, not only an appeal to universality

The five ordered fractions are

```
[1/6, 1/35, 13/33, 13/77, 17/11].                       (4)
```

On the promised input

```
N = 11 * 2^A * 3^B * 5^C * 7^D,   A,B,C,D nonnegative,
```

the machine halts with marker prime 17 exactly when
`A>=B AND C>=D`, and otherwise with marker prime 13.

**Invariant and proof.** The first fraction removes one matched `(A,B)`
pair, preserving `A-B`, until one exponent is zero. The second does the
same for `(C,D)`, preserving `C-D`. Thus the first two rules run exactly
`min(A,B)+min(C,D)` times and then become forever ineligible. If residual
`B>0`, rule 3 replaces the sole control factor 11 by 13 and consumes one
3. Otherwise residual `D>0` triggers rule 4. If neither exists, rule 5
replaces 11 by 17. There is then no eligible rule. The exact time is

```
min(A,B)+min(C,D)+1.
```

For one speed in (3), the explicit arithmetic preprocessing is

```
(A,B,C,D) = (d*r_i, c*b, d*(b-r_i), c*b).
```

This preprocessing computes the residue, products and prime powers. It is
**not** performed by (4). Confusing verified exponent comparisons with a
self-contained modular-arithmetic compiler would overstate this result.

### A single whole-row FRACTRAN program

For `m` inequalities `A_j>=B_j`, choose distinct register primes `p_j,q_j`
and distinct control/reject/accept primes `s,r,h`. Use the **ordered** list

```
1/(p_1 q_1), ..., 1/(p_m q_m),
r/(s q_1), ..., r/(s q_m), h/s.                         (5)
```

Input is `s product_j p_j^A_j q_j^B_j`. The same invariant proves totality,
exact time `sum_j min(A_j,B_j)+1`, and acceptance iff every comparison
holds. There are `2m+1` fractions; for an `n`-speed rational point,
`m=2n` and there are **`4n+1` fractions**. `machines.json` contains the
literal 53-fraction program, all primes, exponents, input bit lengths and
big-endian integer hashes for four thirteen-speed runs (two input encodings). `check.py` builds
and executes the actual integer at every FRACTRAN step.

### Exact rounding is a genuine representation simplification

Define `K=ceil(cb/d)=(cb+d-1)//d`. Since `r` is integral, (3) is equivalent
to `r>=K AND b-r>=K`. Thus the same machine accepts the smaller registers

```
(A,B,C,D) = (r,K,b-r,K).                               (6)
```

This erases the denominator multiplier from every register without changing
the predicate. On an accepted n-speed point, its exact FRACTRAN runtime is
`2nK+1`; on any point it is at most that. For `V26=26*{1,...,12} union {339}`:

| Time | Encoding | Result | Steps | Input bits |
|---|---|---|---:|---:|
| `2/27` | cross-multiplied | accept | 703 | 33,888 |
| `2/27` | rounded | accept | 53 | 2,440 |
| `1/26` | cross-multiplied | reject | 353 | 33,607 |
| `1/26` | rounded | reject | 28 | 2,431 |

The program is unchanged; only the mathematically equivalent input encoding
shrinks. This is exact rounding with a proof, not an approximate threshold.
At `r=0`, the first inequality rejects; using `b-r=b` in the second causes
no false acceptance.

A small program is not a compressed input: encoding arithmetic quantities as
prime exponents can be exponentially larger than their binary encoding.
The prime-exponent computation is useful as an exact operational model,
not as a claimed speedup or evidence of undecidability.

**Ordering hostile.** Moving `17/11` to the first position falsely accepts
`(A,B,C,D)=(0,1,0,0)`. First-eligible order is part of the mathematical
program, not cosmetic formatting.

## 5. Actual unary single-tape Turing machine

Input promise: `1^A#1^B`, head on the first symbol (on `#` when `A=0`),
blank symbol `_`, work symbol `X`; all other cells blank. Start state `H`;
`ACCEPT` and `REJECT` halt. This literal twelve-transition table decides
`A>=B`:

| State | Read | Write | Move | Next |
|---|---|---|---|---|
| H | 1 | 1 | R | H |
| H | X | X | R | H |
| H | # | # | R | B |
| B | X | X | R | B |
| B | 1 | X | L | C |
| B | _ | _ | S | ACCEPT |
| C | 1 | 1 | L | C |
| C | X | X | L | C |
| C | # | # | L | D |
| D | X | X | L | D |
| D | 1 | X | R | H |
| D | _ | _ | S | REJECT |

The source uses `A,R` for the halt-state names. One completed cycle marks
one 1 on each side; the difference of unmarked counts is invariant.
If the right side is empty, accept. If it supplies another 1 and the left
is already empty, reject. Each cycle strictly reduces the right count, so
the promised input always halts. The implementation checks the loose bound
`4(A+B+2)(min(A,B)+2)`. No transition-table implementation of the full
interval search is asserted: that algorithm is separately explicit Python;
this TM is the exact inequality primitive used by the certificate.

## 6. Finite decision versus an open theorem

For a fixed supplied integer row, Section 3 is a total exact decision
procedure. For a finite declared row box, it gives a total exhaustive search.
Unbounded enumeration of all thirteen-speed rows can discover a failure but
will not certify success merely by continuing to find witnesses.

It would nevertheless be wrong to call fixed-13 LRC inherently undecidable,
or to say no finite global bound is known: inherited THM-763, conditional on
its lower-dimensional LRC citation, places every primitive thirteen-speed
failure or tight row inside

```
sum_i v_i <= 91^12.
```

Under that lower-dimensional assumption, a finite scan of the bound,
combined with the cited reduction, would be a terminating global decision. **We have not run that scan.** Structural
compression/discharge of the remaining bounded locus is still needed for a
practical proof. No large-machine, small-program, or universality statement
supplies the missing proof. We inherit the repository's conditional/cited
status; this session performs no external literature audit.

## 7. Connection contract and verification

```
source: rational loneliness predicates and interval-core geometry
target: owner-labelled closed cells, prime-exponent comparisons, unary TM
map: intervals (2); residues (3); exponent input and ordered fractions (5)
preserved: closed safety, equality, tied owners, every comparison's verdict
lost: measure loses isolated points; bare verdict loses owners and location
sidecar: exceptional primitive grid, owner signs, t, delta, residues, registers
cheapest test: tight AP; half-threshold singleton; q25 hostile; reversed rule order
```

`python3 -B check.py` and `python3 -O -B check.py` produce byte-identical
JSON. Every check uses explicit exceptions, never an optimization-elided
`assert`.

- Five-fraction primitive: all `0<=A,B,C,D<=8` (6,561 cases).
- Literal unary TM: all `0<=A,B<=64` (4,225 cases).
- Rational-to-FRACTRAN interface: 4,860 cases for each of the two exact
  encodings (9,720 machine runs).
- Exact intersection versus independent wall/midpoint oracle: 1,020 cases.
- Harmonic closed-set formula versus two independent realizations: `1<=k<=80`.
- Complete bounded-search smoke census: all 560 increasing thirteen-speed
  subsets of `{1,...,16}` (primitive filter explicit; no other pruning),
  each returning an independently integer-checked witness. This revalidates
  the instrument on a known small universe; it is not a new frontier bound.
- Canonical hostiles: coherent V26, gcd-incoherent S*, zero-free saturated S0;
  first good periods respectively 27, 26, 27.
- Negative/equality controls: `{1,2}` at `1/2` empty; `{1}` at `1/2`
  singleton; owner tie at `t=1/4`; zero-measure tight AP; reordered fractions;
  six malformed inputs rejected.

Artifacts: `check.py`, `machines.json`, `results.json`,
`results.optimized.json`, `MANIFEST.sha256`. The `inherited/` directory holds
read-only source extracts for local review, not new mathematical content.

### Independent audit

The separate audit lane reports **PASS** on the serialized program artifacts,
not just this file's helper functions: 2,401 literal five-fraction cases,
1,681 runs of the serialized twelve-transition TM, all four whole-row
FRACTRAN runs and integer hashes, and an independent Farey construction for
`1<=k<=120` with all component endpoints/midpoints checked against every
speed. The theorem/quantifier review accepts the harmonic equality sidecar
and retains the lower-dimensional assumption on the global finite bound.
Reproducible independent artifacts are packaged as siblings at
`../independent/audit_lrc.py` and `../independent/lrc.normal.json`; these
counts are separate from the primary test universes above.
