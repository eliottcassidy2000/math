# From Berggren edge addresses to Collatz block budgets

**Status:** PROVED elementary statements with independent audit: exact
repetition counters, residue representatives, explicit reset families, and
the finite linear valuation-rank obstruction. FINITE-EXACT controls are
separately quantified below. Proposed global descent mechanisms are OPEN.
**Collatz remains OPEN.** No literature-priority or Lean claim.

The [pasted source](../reference/COLLATZ-KURATOWSKI-2026-09-25-SOURCE.md)
is preserved with a warning. Its edge dictionary is useful; its claims about
complete orbits, completeness of cycle obstructions, and the scope of a
finite-bit no-go require the repairs in section 1.

## Inheritance and research board

Closest proved mechanism: [THM-3756, odd-square ordinal Berggren affine
descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),
the source's explicit edge-address formula, and the exact guarded affine
ports in [creative_descent](creative_descent_20260925.md). Canonical hostile:
the legal edge `7->11` has a Berggren parent that is not a plus edge.
Corrected near miss: choosing only the smallest sibling loses a usable
certificate at `241`, while `483` needs an intermediate ladder height;
see [creative_sibling](creative_sibling_20260925.md). Least-used sidecar in
this session: the rational fixed point of an entire orbit block, together
with the integer's dyadic distance from it and its ordinary height.

The anchor is orbit-block composition. The niche is a precise repair of
finite-obstruction language. The wildcard is the transfer of dyadic
proximity between different rational centers.

| Live object | Retained predicate | Cheapest hostile test / outcome |
|---|---|---|
| Berggren edge address | exact endpoints, sign, valuation | ancestry only groups same-target fibres |
| Affine orbit block | exact input cylinder and endpoint | arbitrary finite exponent words occur |
| Periodic center | exact number of consecutive block copies | finite escape does not imply descent |
| Center transition | carry determinant and dyadic cancellation | `27->41->31` replenishes the counter |
| Integer boundary | canonical representatives of nested cylinders | minus cycles stabilize at 5 and 17 |
| Candidate rank | ordinary height plus dyadic counters | every finite linear center bank fails |

Relevant prior work was recovered before deriving these statements:
[THM-4474, strategy cube](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md)
already rules out its precisely defined bounded-lookahead certificates;
[the valves](collatz_guards_20260921_valves.md) and
[the run transducer](creative_transducer_20260925.md) already contain the
rise/fall macro and its descent threshold; the
[incoming synthesis](collatz_procgen_20260925_incoming_synthesis.md) already
demotes Berggren antichains as divergence certificates. The negative
preperiodic point `-13/9` occurs in
[the choice ladder](collatz_procgen_20260922_choice_ladder.md).
The contribution here is an explicit repetition/reset mechanism and its
rank obstruction, not the discovery of those underlying objects.

The meta-patterns used are **Separate unbounded local support from a
height-bounded modular cover** and **Correct the object before sharpening
the technique**, from [META-PATTERNS](../../00-navigation/META-PATTERNS.md).
No new universal method card is promoted from this one thread.

## 1. Repairs needed before using the pasted argument

### 1.1 Complete-orbit incomparability is false as stated

Under `U_-(n)=(3n-1)/2^v2(3n-1)`,

    23 -> 17 -> 25 -> 37 -> 55 -> 41 -> 61 -> 91 -> 17.

The distinct edges `23->17` and `91->17` have Berggren addresses
`ABCC` and `ABCCAA`, computed directly with the parent map. Both targets
are outside `{1,5}`, yet one address is a prefix of the other.
The failed implication in source 3.3 is "an orbit visits each odd value
once." That assumes away recurrence at cycle entry.

The survivor is source 3.2: for distinct same-sign edges with nonexceptional
targets, comparability is exactly equality of targets. Apply incomparability
only to orbit segments whose targets are distinct. This counterexample
does not refute the edge-address formula or the separate consecutive-edge
statement outside its specified exceptions.

Also, a lossless edge address can recover the forward dynamics. It is
**bare ancestry**, not the full encoding, that supplies no forward rank.
Infinite antichains in the tree are compatible with an escaping orbit.

### 1.2 Cycles are not all possible obstructions

An orbit on a fixed reduced positive rational sheet either eventually
cycles or, if it never repeats, eventually leaves every finite height
interval. This follows because the numerator is an integer and bounded
sets are finite. Thus source 7.2/7.4 omit aperiodic escape. Excluding
nontrivial positive integer cycles addresses only the cycle half of
Collatz; it is not equivalent to coverage without a no-escape theorem.

For a reduced sheet `q` coprime to 6, retain `gcd(a,q)=1` in `a/q`.
The raw `3a+q` census can include lower-denominator cycles. For example,
`a=q=5` represents the ordinary integer fixed point 1, not reduced sheet 5.

### 1.3 The finite-bit conclusion needs its exact quantifier

For every `L`, `n=-1 mod 2^(L+1)` has `L` initial strict rises. This
excludes a uniform bounded first-descent horizon. It does not exclude a
finite-description proof using variable-length blocks, unbounded counters,
or induction on additional arithmetic data.

More strongly, section 2 shows that **every finite positive exponent word
occurs on infinitely many positive integers**, for each sign. Consequently
there are no universally forbidden finite valuation blocks to exploit.
This does not prohibit a finite rule set with unbounded parameters.

Kuratowski-style reasoning needs an operation preserving the relevant
property and an appropriate ordering theorem. Neither is supplied by
matching girths or by an embedding of single edges into a tree. Compressing
an orbit block while keeping its exact affine boundary data is a valid
operation; forgetting those data is the failed transfer.

## 2. Exact arithmetic block compression

Let `q>=3` be odd and `sigma=+1 or -1`. On positive odd integers use

    U(n) = (q n + sigma) / 2^v2(q n + sigma).

For a nonempty exponent word `w=(k_1,...,k_r)`, all `k_i>=1`, put

    S = sum k_i,   A=q^r,   B=2^S,
    C = sigma * sum_(j=0)^(r-1) q^(r-1-j) 2^(k_1+...+k_j).

The empty inner sum is zero. The exact affine block is

    F_w(n)=(A n+C)/B.                                      (1)

Its source is one odd residue class modulo `2B`:

    r_w = A^(-1)(B-C) mod 2B,  1<=r_w<2B.                 (2)

**PROVED.** The first `r` exact valuations equal `w` iff `n=r_w mod 2B`.
Necessity follows by multiplying the step equations and retaining the
oddness of the terminal value. Sufficiency follows backwards through
the congruence: the first numerator is divisible by exactly `2^k_1`,
and its odd quotient satisfies the remaining word's congruence. Equivalently,
the inverse branch `(2^k y-sigma)/q` maps odd 2-adics bijectively to the
appropriate cylinder. This also proves the all-finite-words assertion.

Concatenation is exact: applying `w` then `v` gives coefficients

    (A_v A_w, B_v B_w, A_v C_w+B_w C_v).                    (3)

Preserved: source legality, endpoints, sign, exact valuations, and their
composition. Lost if only `(A,B,C)` is retained: the displayed internal
vertices; keep the word when those matter. The source cylinder is mandatory.
The branches are actual forced steps, never a choice of convenient exponents.

## 3. A complete countdown for repeated blocks

Define the odd integer `D=A-B`, rational center `alpha=-C/D`, and defect

    E_w(n)=D n+C.

**PROVED repetition lemma.** The maximum number of complete consecutive
copies of `w` at the beginning of the actual orbit of positive odd `n` is

    R_w(n) = floor((v2(E_w(n))-1)/S),       if E_w(n)!=0;   (4)
    R_w(n) = infinity,                    if E_w(n)=0.

The zero case is a genuine periodic orbit with this block, not a discarded
singularity. The value in the nonzero case is nonnegative: `D` and `C`
are odd, so `E_w(n)` is even.

*Proof.* Compose the inverse branches of `w`. They contract the odd
2-adics by `2^-S` and have a unique fixed point, necessarily `alpha`.
Its forward itinerary repeats `w`. By (2), the cylinder for `t` copies is
therefore `n=alpha mod 2^(tS+1)`. Since `D` is odd, this is equivalent to
`v2(E_w(n))>=tS+1`. The maximum integer `t` gives (4). Alternatively,
direct algebra gives the exact consumption law

    E_w(F_w(n))=(A/B) E_w(n).                              (5)

Every copy consumes exactly `S` units of valuation. After the last complete
copy, the residual valuation lies in `{1,...,S}`. This proof applies to
both signs and every odd multiplier `q>=3`.

For plus Collatz and an expanding block (`A>B`), `D n+C>0`, so infinite
repetition is impossible at positive integers. Its duration obeys

    R_w(n) S + 1 <= log_2(D n+C).

This explains, with an exact counter, why positive integers can follow
a negative periodic center for a long time but must eventually leave it.
It does not bound how many different centers an orbit can follow.

### 3.1 A long sequence of short rise/fall runs

For `q=3,sigma=+1,w=(1,2)`,

    F_w(n)=(9n+5)/8,   alpha=-5,
    R_w(n)=floor((v2(n+5)-1)/3).

For every `N>=1`, the start `x_N=2*8^N-5` has exactly `N` copies, with
boundary values `2*9^j*8^(N-j)-5`, `0<=j<=N`, all strictly increasing.
In maximal rise/fall notation each is the same passport `(a,b)=(2,1)`.
Thus bounded lookahead in maximal runs fails even when each rise count
is only two. On the minus sheet replace `-5` by `+5`; the limiting
center is then the genuine positive cycle point 5.

At the plus family's first escape the next exponent is 3, giving

    e_N=(3*9^N-7)/4,        v2(e_N+1)=1+v2(N).             (6)

For the valuation, `e_N+1=3(9^N-1)/4` and
`v2(9^N-1)=3+v2(N)`. Nevertheless `e_N>x_N` exactly for `N>=9`.
Indeed four times the difference is `f_N=3*9^N-8*8^N+13`;
check `1<=N<=9`, then `f_(N+1)=9 f_N+8*8^N-104>0` for `N>=9`.
Finite escape from the negative shadow has not necessarily paid back its
ordinary-height growth.

## 4. Where the countdown is replenished

For two blocks set `Delta=D_w C_v-D_v C_w`. At any seam `y`,

    D_w E_v(y)=D_v E_w(y)+Delta.                           (7)

Both D's are odd. If `v2(E_w(y)) != v2(Delta)`, the new valuation is
their minimum. It can exceed that minimum only when the two valuations
agree and cancellation occurs. If `Delta=0`, the centers agree and no
valuation reset occurs. This is an exact arithmetic description of
the transition, not an assumed loss of information.

**PROVED unbounded reset.** For `H=6j+5`, `j>=0`, put

    n_H=(2^(H+3)-13)/9.

It is a positive odd integer, and exactly one `(1,2)` block sends it to
`y_H=2^H-1`. The next `(1)` block repeats `H-1` times. In detail,

    v2(n_H+5)=5,   v2(y_H+5)=2,   v2(y_H+1)=H.

The first examples are

    27 -> 41 -> 31,           followed by 4 rises;
    1819 -> 2729 -> 2047,     followed by 10 rises;
    116507 -> 174761 -> 131071, followed by 16 rises.

Both block types expand in ordinary size. For this switch `Delta=-4`;
the residual valuation 2 meets precisely the cancellation condition in (7).
A rank that counts only the current remaining repetitions therefore fails.

The missing proximity was already present at the source, around the
preperiodic rational center

    alpha_*=-13/9 -> -5/3 -> -1,   exponents (1,2),
    9 n_H+13=2^(H+3).

Every periodic rational center for `3x+1` has reduced denominator coprime
to 3, since it divides `3^r-2^S`. Hence `alpha_*` is not periodic.
For any fixed finite bank of periodic centers `beta`, the source counters
`v2(n_H-beta)` eventually stabilize: `n_H` converges 2-adically to
`alpha_* != beta`. Meanwhile `v2(y_H+1)=H` diverges.
This proves a precise finite-bank limitation; it does not exclude a
finite transition formula with unbounded arithmetic registers.

The general pullback of a center `beta` is

    alpha=(B_w beta-C_w)/A_w,
    v2(F_w(n)-beta)=v2(n-alpha)-S_w.                       (8)

So pulling centers backwards accounts for a reset exactly, but may
introduce growing powers of 3 in their denominators. A useful next target
must price that extra arithmetic complexity as well as orbit height.

## 5. A sharp obstruction to a natural candidate rank

**PROVED.** Fix finitely many negative rational numbers `beta_i`, real
coefficients `c_i`, and `a>0`. No function

    V(n)=a log n + sum_i c_i v2(n-beta_i)                   (9)

strictly decreases at every plus odd step `n->U_+(n)`, `n>1`.
Here rational valuations are numerator valuation minus denominator
valuation, and repeated centers are combined. Negative centers keep
the function finite on every positive input.

*Proof.* Set `n=2^H-1`. Its successor is `3*2^(H-1)-1`. For every
center except `-1`, the two valuations eventually agree and are constant.
At `-1` the valuation drops from `H` to `H-1`. Descent therefore forces
the combined coefficient `c_(-1)>a log(3/2)>0`; equality also fails because
the actual growth ratio is strictly greater than `3/2`. If the center is
absent this is already impossible.

There are infinitely many distinct negative centers

    alpha_d=1-2(4/3)^d,    d>=1.

Choose one absent from the finite bank. For
`H=1+2*3^(d-1)t`, `t>=1`, the number

    n_H=1+(4/3)^d(2^H-2)

is a positive odd integer. Divisibility follows from
`3^d | 2^H-2`; its first `d` valuations are exactly 2, ending at
`y_H=2^H-1`. The intermediate points are
`1+(4/3)^(d-i)(2^H-2)`, so integrality and exact valuations are explicit.
All these sources exceed 1.

As `t` grows, `n_H` tends 2-adically to the absent `alpha_d`, so every
source counter eventually stays constant. At the target only the `-1`
counter is unbounded and equals `H`. Meanwhile
`log(y_H/n_H)` tends to `d log(3/4)`. Thus

    V(y_H)-V(n_H)=c_(-1) H+O(1) -> +infinity,

contradicting the sum of the `d` strict descent inequalities. QED.

Scope: (9), finitely many fixed negative rational centers, linear
combination of their valuations, positive logarithmic coefficient, and
descent at every odd step. This does not exclude nonlinear couplings,
variable center families, common-future dependencies, or descent at
an adaptively chosen later time. It identifies an obstruction to a
particularly natural repair of ordinary log drift.

## 6. A second reframe: which symbolic rays contain an integer?

For an infinite exponent word `(k_1,k_2,...)`, let `K_m=sum_(i<=m) k_i`
and `C_m=sum_(j=0)^(m-1) 3^(m-1-j) 2^K_j`, with `K_0=0`. Its canonical
positive odd source representative is

    r_m^sigma = [3^(-m)(2^K_m-sigma C_m)] mod 2^(K_m+1).   (10)

Nested cylinders imply

    r_(m+1)^sigma = r_m^sigma + b_m 2^(K_m+1),
    0<=b_m<2^k_(m+1).                                    (11)

**PROVED.** The infinite word is realized by a positive integer iff
these nondecreasing representatives eventually stabilize. If realized
by `n`, the representative equals `n` once the modulus exceeds `n`.
Conversely, an eventual constant belongs to every cylinder and realizes
every finite prefix. Therefore boundedness and stabilization coincide.

The plus orbit reaches 1 iff its exponent word is eventually all 2:
the reverse implication follows from the repetition lemma for `w=(2)`,
whose only indefinitely repeating integer is 1. Hence

    Collatz iff every word not eventually all 2
               has r_m^+ tending to infinity.             (12)

This reformulation includes both cycles and aperiodic escape. It is not
a proof or a claim of easier complexity. It states the exact missing
ordinary-integer boundary condition inside the 2-adic coding.

For the same word, the two signs satisfy

    r_m^+ + r_m^- = 2^(K_m+1).                            (13)

Thus the known minus cycles are hostile controls, not a structural
contradiction: their minus representatives stabilize at 5 or 17, while
the plus representatives of those same words grow by (13). Every finite
word is possible on both signs; the infinite bounded representative is
where the sign-sensitive integer question lives.

## 7. Concrete proof-attempt reframes worth retaining

1. **Arithmetic block reductions.** Use the graph-theoretic inspiration
   as a demand for legal reductions. A reduction must carry `(w,A,B,C)`,
   source congruence, and endpoint/common-future evidence. Search for
   parameterized reduction families; the supplied graph triple does not
   itself furnish one. The first test is a sound rule covering an
   additional infinite residual family of the existing sibling grammar.

2. **Coupled height and center complexity.** Equations (5), (7), and (8)
   turn repeated loops and their changes into exact arithmetic. Seek a
   certificate that charges cancellation resets to a decreasing resource
   retained across the seam. A fixed linear bank is ruled out by section 5;
   the residual candidates must involve nonlinear coupling, an unbounded
   symbolic family of centers, or a different proof-dependency order.
   First decisive tests: `x_N`, `n_H`, the pullback family `alpha_d`,
   minus cycles, and the `5x+1` cycles. Merely proving finite residence
   near each negative center does not meet the obligation.

3. **Integer exclusion on bad symbolic rays.** Try to force infinitely
   many positive extension digits `b_m` in (11), for a rigorously defined
   family of nonterminal words. Periodic and highly repetitive words
   are the easiest controls; arbitrary changing words are the hard part.
   This connects to the repository's existing
   [Diophantine hard-class work](collatz_procgen_20260922_hard_class.md),
   rather than bypassing it. A successful theorem for a new explicit
   infinite word family is a useful partial result, not universal coverage.

No tournament is imposed: the intrinsic objects here are directed orbit
blocks and their boundary identities, not pairs with a natural tournament
orientation. The `(3,4,5)` geometry contributes a precise dictionary;
the new proof obligation concerns composition and well-foundedness.

## 8. Computation, independent audit, and scope

Reproduce from the repository root:

    python3 04-computation/experiments/kuratowski_reframe_20260925.py
    python3 04-computation/experiments/kuratowski_reframe_20260925_audit.py
    python3 04-computation/experiments/kuratowski_reframe_20260925_centers.py

Outputs: [main](kuratowski_reframe_20260925.out) and
[independent audit](kuratowski_reframe_20260925_audit.out), plus the
[finite-center family](kuratowski_reframe_20260925_centers.out).
The main controls compare direct odd iteration with (4) for multipliers
3 and 5, both signs, all 340 exponent words of lengths 1--4 over `{1,2,3,4}`,
and all 1024 positive odd starts through 2047: **1,392,640 cases**, zero
mismatches, including 23 zero-defect controls. There are also 4080 exact
cylinder-lift checks, the explicit source correction, both shadow families
through `N=64`, and the reset family through `H=197`.

The independent implementation was written by another agent before seeing
the main code: **1,361,360** repetition checks (odd starts through 2001),
23 zero-defect controls, **13,644** switch identities, and 21 reset-family
checks. A separate logical audit checked the congruence modulus, zero case,
switch cancellation, finite-bank scope, and section 5 proof.
The center-family probe checks 56 exact instances, `d=1..7, t=1..8`,
including integrality, every prescribed step, the endpoint, the reduced
denominator, and the dyadic source proximity. File hashes are recorded in
[the manifest](kuratowski_reframe_20260925_manifest.json).
The infinite assertions rest on the proofs, not on these finite ranges.
No general `5x+1` divergence claim is used: its exact cycles at 1,13,17
are sufficient hostile controls and are retained by both verifiers.

Primary literature context only: Lagarias,
[*The set of rational cycles for the 3x+1 problem*](https://eudml.org/doc/206298),
Acta Arithmetica 56 (1990), 33--53, studies the rational-cycle viewpoint.
The elementary derivations here do not require a theorem from it.
No novelty claim follows from the targeted prior-work search.
