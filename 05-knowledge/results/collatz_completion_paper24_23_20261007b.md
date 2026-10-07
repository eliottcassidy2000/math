# Completion premises and an exact deletion request

**PROVED:** the source-specific request interval, sharp uniform numerical
budget, and authenticated transport below. **FINITE-EXACT:** the declared
implementation controls. The supplied paper headlines are accepted premises
for mechanism comparison; none is a dependency of the elementary Collatz
proofs. **OPEN:** supplying the requested actual deletion for every source.

The productive form of “assume completion, then recover the missing rule”
is to specify the certificate that would make the completion usable. This
package computes that obligation at the unchanged source. It also detects
sources where no positive odd bit-deletion target can pay, even though a
stronger supplied ROOT proof still transfers successfully.

## 1. Inheritance, papers, and the comparison map

The closest proved algebra is
[guarded child normal forms](collatz_child_normal_forms_20261007.md): the
three native inverse maps below retain the ordered carry and an actual
forward word back to their input. The six-bit limit and its carry-dependent
height cut were established in
[mixed child rules, section 5](collatz_mixed_child_rules_20261007.md).
[Grounded child closure](collatz_grounded_child_closure_20261007.md) already
closes a genuinely supplied ROOT family under compatible inverse words;
its nonzero-phase reseeding hostile shows why the original parameter must
remain attached. The present task keeps a single supplied source fixed.

The canonical hostile is `35 -> 23 -> 15` under two inverse-1 operations:
its only positive odd bit-deletion target is 17, which is larger than 15.
The corrected near miss is to treat a sufficiently large numerical depth
as an available, proved Collatz deletion. The least-used sidecar is the
finite interval of legal depths, alongside the actual common-future words.
The live board comprises fixed source, native guard, affine carry,
numerical payment, authenticated deletion, and completed parent proof.

The local PDFs were read as supplied mathematical premises. Page numbers
below are PDF pages; complete pages, including the indicated diagram, were
inspected.

* [Paper 24](C:/Users/Eliott/Downloads/paper-24.pdf), **Semialgebraic universal
  covers of normal projective varieties**, OpenAI, September 24, 2026,
  65 pages. Its accepted headline characterizes the stated semialgebraic
  universal covers by a bounded-symmetric-domain, Euclidean, and compact
  projective factorization. Pages 4–5 and Figure 1 separate local chart
  control, quotient maps, and the additional inverse-chart argument needed
  for a global identification. Remark 2.3 on page 7 keeps the base parameters
  fixed during correction: constants and thresholds need not be uniform
  merely because a power estimate is. Page 39 adds boundary and nonvanishing
  guards to a limit argument. Page 47 uses properness to upgrade a countable
  image to a locally finite one. The transferable discipline is to preserve
  the parameter being repaired and explicitly prove a uniform bound.
* [Paper 23](C:/Users/Eliott/Downloads/paper-23.pdf), **Ramanujan–Arthur
  Decompositions of Cuspidal Functions at Full Finite Level**, OpenAI,
  September 24, 2026, 64 pages. Its accepted headline is a full-finite-level
  cuspidal decomposition over global function fields with the stated
  Ramanujan–Arthur behavior. Pages 5 and 42–44 retain ordered Bernstein
  weights before quotienting and prove detector/ideal statements with
  finite algebraic expressions. In particular, page 43 does not replace a
  strong grading by an assumed group-of-units splitting. Page 47 keeps
  dimension and operator norm in a uniform bound on powers of possibly
  nonnormal operators. The transferable discipline is that a detector or
  support statement and its executable finite witness are separate data.

| Source mechanism | Target and map | Preserved predicate | Lost information and needed sidecar | Decisive test |
|---|---|---|---|---|
| Fixed-parameter correction in paper 24 | Hold parent `m` fixed while asking for a deletion | Same supplied arithmetic problem | Replacing `m` by another native/CRT member changes the problem | `m=35`, source 15 has an empty request interval |
| Uniform bounds with retained scale/norm | Maximize `(m+1)/(h+1)` on one native progression | Numerical payment for every member | Binary deletion guard and actual receipt are not supplied by the maximum | Seven bits suffice numerically for `G5^36`, but must still be legal |
| Finite detector expressions in paper 23 | Consume an explicit common-future receipt | Exact equality of actual odd endpoints | Counts, slope, and an assumed completed rule set omit the word witness | The authentic `35 -> 17` receipt is rejected after transport to 15 |

These are mechanism transfers, not a functor from the papers' geometric or
automorphic objects to Collatz orbits. No compactness, Hecke, spectral-gap,
or positivity theorem from either paper is asserted for Collatz here.

## 2. The fixed-source theorem

Let `U(n)=(3n+1)/2^v2(3n+1)` on positive odd integers. The input guards are

\[
G_1(n)=\frac{2n-1}{3},\quad n\equiv2\pmod3;
\qquad G_5(n)=\frac{8n-5}{9},\quad n\equiv4\pmod9;
\qquad G_{17}(n)=\frac{2048n-2363}{2187},\quad n\equiv2170\pmod{2187}.
\]

For a nonempty chronological legal word at parent `m`, write its child as

\[
h=\frac{qm-C}{p},\qquad q=2^A,\quad p=3^B.
\tag{1}
\]

The inherited exact normal form proves `3 <= h < m`,
`p-q <= C <= 17(p-q)`, and supplies an actual positive-odd word from `h`
back to `m`. Final integrality on the native cell enforces every prefix
guard. These assertions concern the whole ordered word, not just `(A,B)`.

Consider only bit-deletion targets

\[
z_D=\frac{m+1}{2^D}-1,
\qquad D\ge1.
\tag{2}
\]

**Exact interval.** Define

\[
D_{\min}=\operatorname{bitlength}\!\left(
 \left\lfloor\frac{m+1}{h+1}\right\rfloor\right),
\qquad D_{\max}=v_2(m+1)-1.
\tag{3}
\]

Then `z_D` is a positive odd integer smaller than `h` if and only if

\[
D_{\min}\le D\le D_{\max}.
\tag{4}
\]

Indeed, positivity and oddness are exactly `D <= v2(m+1)-1`.
The payment inequality is exactly `2^D(h+1)>m+1`. For a rational `R>1`,
the least integer `D` with `2^D>R` is `bitlength(floor R)`, including the
strict boundary when `R` is an exact power of two. The smallest available
target is `2*oddpart(m+1)-1`, so (4) is nonempty exactly when

\[
2\operatorname{oddpart}(m+1)-1<h.
\tag{5}
\]

This is an iff about numerical targets. It asserts no Collatz connection
between `m` and `z_D`. For a legal Mersenne parent `m=2^E-1`, the interval
is always nonempty because its maximal deletion has target 1 and `h>=3`.
Demanding an authenticated deletion to 1 can therefore demand the very
ROOT theorem still missing; arithmetic feasibility is not its proof.

## 3. A sharp uniform numerical budget

Fix the whole word. Its positive odd native pairs are exactly

\[
m=m_0+2pt,\qquad h=h_0+2qt,\qquad t\ge0,
\]

with `(m0,h0)` the least native pair. The ratio

\[
R(t)=\frac{m_0+1+2pt}{h_0+1+2qt}
\]

is nonincreasing: the numerator determining the sign of its derivative is
`2(p(h0+1)-q(m0+1))=2(p-q-C)<=0`. Consequently

\[
D_* = \operatorname{bitlength}(\lfloor R(0)\rfloor)
\tag{6}
\]

is the least uniform integer for the strict numerical comparison
`2^D(h+1)>m+1` on the entire progression. It does **not** imply that depth
`D_*` satisfies the separate binary guard at every member.

Two finite upper bounds make the obligation predictable without a search
through supplied sources. Every legal letter obeys
`G(n)+1 >= (2/3)(n+1)`. For `G1` equality holds; for `G5` the difference
is `2(n-1)/9`; for `G17` it is `(590n-1634)/2187`, positive on its native
positive odd domain, whose least input is 4357. Thus for length `ell`,
`R(t)<=(3/2)^ell`. Alternatively, using `h>=3` and the carry envelope,

\[
R(t)=\frac pq+\frac{q+C-p}{q(h+1)}
\le \frac{5p}{q}-4 < \frac{5p}{q}.
\]

Hence either `2^D>(3/2)^ell` or `2^D*q>=5p` is a sufficient uniform
numerical budget. The exact calculation gives `D_*=7` for both `G5^36`
and `G17^64`; six bits fail at the least native pair. This repairs the
fixed-six-bit *comparison*, while leaving the required seven-bit receipt
as an explicit unproved obligation if one has not been supplied.

## 4. Consume a completion premise rather than assume its contents

`request(m,word)` authenticates (1) and returns (3). `audit_request`
recomputes the fields and rejects malformed exact types. `consume_deletion`
requires an independently supplied, replayed receipt

\[
m\xrightarrow{u} e \xleftarrow{v} z_D,
\]

with actual positive odd valuation words. It checks the original parent,
the exact power-of-two deletion identity, and (4). Prepending the normal
form's actual return word gives
`h -> m -> e <- z_D`, now a paid dependency because `z_D<h`.
No path is discovered by this API. The inherited receipt validator retains
the first-ROOT boundary; a supplied proof of the smaller target may then
ground the dependency. Finding a different source where the guard works
is deliberately not an operation of the compiler.

At parent 11, `G1` gives source 7, and the exact interval is `{1}`.
The supplied receipt `11 --(1,2,3)--> 5` transports to
`7 --(1,1,2,3)--> 5`. At parent 35, `G1^2` gives 15, but `(Dmin,Dmax)=(2,1)`.
There really is an actual common future
`35 --(1,5)--> 5 <--(2,3)-- 17`; transporting it cannot pay 15 because
17 is larger. This is a failed comparison, not a failed Collatz orbit.

A stronger supplied premise still works at exactly that source:
`35 --(1,5,4)--> 1` grounds
`15 --(1,1,1,5,4)--> 1` through `consume_parent_root`. The actual inverse
prefix reaches `m>1` and cannot pass through ROOT first, so concatenation
with a first-hit ROOT word remains first-hit. Thus a genuine completed
parent proof closes all its legal inverse children. Assuming a modified
system is completed supplies such a proof only after each modified edge
used in it has been replaced by an authenticated original-system receipt.

## 5. Exact finite controls and what remains open

Reproduction, from the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_completion_paper24_23_20261007b.py
python -B -O -X utf8 04-computation/experiments/collatz_completion_paper24_23_20261007b.py
```

The independent paths compare the affine carrier, literal return words,
all allowed depths, and the closed-form interval. The universe is all
363 nonempty generator words of length at most five, with 17 native pairs
each: 6,171 pairs, of which 1,882 have a nonempty arithmetic request and
4,289 do not. There are 31 additional literal Mersenne phase controls with
exponents at most 4,096, two long-word seven-bit controls, exact positive
and hostile receipts, and malformed type/source/depth/ROOT controls.
Both modes perform **51,610 checks** and produce the saved `.out` exactly.
Normalized-LF SHA256:
`4bc376b92ba0fcbd0f2477ae00fd305a3ff589cd5cd267eb7c105c3cb9a8f6ad`.

The surviving positive theorem is a computable request and a sound
certificate consumer. No finite controls establish arbitrary-source ROOT
coverage. The open obligation is concrete: supply a legal receipt at one
of the depths in (4), or supply another paid rule or a stronger completed
parent proof. An empty interval is a reason to change the rule family,
not to silently increase the depth or change the supplied source.
