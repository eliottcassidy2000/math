# Child closure after the ternary reduction leaves Mersenne form

**PROVED:** the guarded mixed-word paid exits, exact normal forms, and
closure of the explicitly grounded seed families described below.
**FINITE-EXACT:** the selected first-hit ROOT certificates and declared
computation universes. **OPEN:** closure of the entire arbitrary-source
obligation graph, and therefore a Collatz proof for every integer.

**Forward route:** [Completion-first research and paid exponent phases](collatz_completion_synthesis_20261007b.md)
adds explicit exits in every long reset shell, exact fixed-source deletion
requests, and a fully solved digit-factorial comparison model. Arbitrary-source
Collatz coverage remains open.

The new progress is a way to continue after a child changes shape. Its
return to a larger Mersenne parent can be composed with a stronger existing
deletion, giving an obligation smaller than the changed child itself.
A complementary construction keeps complete ROOT certificates through
arbitrary finite child-word refinements of generated families.

## 1. The stopped child now has a paid exit and a grounded instance

Start with the previous parent `m=2^1457-1`, writing `X=2^1457`. The
native child chain is

\[
m\xrightarrow{G_5}\frac{8X-13}{9}
 \xrightarrow{G_1^5}\frac{256X-2315}{2187}
 \xrightarrow{G_5}z=\frac{2048X-29455}{19683}.
\]

The simple inverse rules stop at `z`. Its known return word
`(1,2),1^5,(1,2)` goes back to the larger `m`; that alone is not payment
against `z`. Composing with the preceding six-bit common-future deletion
at `m` gives the new obligation

\[
q=2^{1451}-1<z.
\]

The exact comparison is

\[
z-q=\frac{111389X-625408}{64\cdot19683}>0\qquad(X\ge6).
\]

Both sides retain actual valuation words and their common future. The
already supplied first-hit certificate at `q` completes the new receipt:

| Changed child | Odd steps to first 1 | Total valuation |
|---|---:|---:|
| `(8*2^1457-13)/9` | 7349 | 13105 |
| `(2048*2^1457-29455)/19683` | 7356 | 13113 |

No new orbit discovery was used to ground these points. The same stopped
chain and paid exit hold on an explicit infinite mixed exponent phase,
where grounding the arbitrary smaller Mersenne child remains an obligation.
See [the paid-exit proof and code](collatz_mixed_child_rules_20261007.md).

## 2. The full mixed word has an exact payment budget

The inherited generators are

\[
G_1(n)=\frac{2n-1}{3},\quad
G_5(n)=\frac{8n-5}{9},\quad
G_{17}(n)=\frac{2048n-2363}{2187}.
\]

Each has its own actual ternary guard. Their guards are disjoint, so a
supplied integer does not get to choose whichever map has a convenient
cost. A chronological mixed word has the exact normal form

\[
h=\frac{2^A m-C}{3^B}
  =\frac{2^A X-D_0}{3^B},\qquad D_0=2^A+C,
\]

where the carry `C` retains the order. Given an authenticated deletion at
`m` with smaller obligation `q=X/2^D-1`, the composite pays against `h`
exactly when

\[
\boxed{(2^{D+A}-3^B)X>2^D(D_0-3^B).}             \tag{1}
\]

This is a rule for the actual source and its guarded receipt, not a
comparison against a reset intermediate height. For this generator
alphabet `D_0>=3^B`, so payment requires `2^(D+A)>3^B`. The carry fixes
the remaining exact height threshold. Equal slopes alone lose that data.

For the existing `D=6` deletion:

- Pure `G5` depth is paid exactly through 35 blocks when `X>=8192`;
  36 or more fail at every positive height.
- Pure `G17` depth can reach 63. Depth 64 fails the slope test.
- Every funded mixed word has at most 63 generators. The large native
  deep-deletion phase satisfies the height cutoff for all of them.

Every finite word beginning with `G5` or `G17` has one nonempty Mersenne
exponent phase; words beginning with `G1` have none. The compatible odd
phase can be intersected with the existing deep binary deletion phase.
Thus **every mixed word passing the six-bit budget yields an infinite
guarded family of paid exits back to a smaller Mersenne obligation**.
This quantifier is over constructed word families; it does not say that
every Mersenne exponent, or every arbitrarily long child word, is covered.

[The normal-form package](collatz_child_normal_forms_20261007.md) proves
the phase classification, freeness of full affine carriers, and exact
word decoder. The router terminates within `11*bitlength(n)` steps,
but its terminal is a proof obligation: `31 -> 27` is the least Mersenne
hostile. No generator can reach 1 from a nonroot positive source.

## 3. Genuine grounded closure is possible with a retained seed coordinate

Take a supplied seed `a>1` with a first-hit ROOT word of length `L` and
cost `A`. Put `P=3^L`, `Q=2^A`. The inherited complete family is

\[
N_s=a+\frac{4Q(4^{Ps}-1)}{3P},\qquad s\ge0.
\]

Every member already has a finite explicit certificate. The new exact
identity is

\[
v_3(N_s-N_t)=v_3(s-t)\qquad(s\ne t).              \tag{2}
\]

Therefore every finite contracting inverse valuation word of length
`ell` selects **exactly one** class `s mod3^ell`. Its child is positive,
odd, legal, smaller, and grounded by prepending that actual return word
to the known seed-family certificate. Further inverse words refine that
class, keeping the entire certificate.

This gives a closed family language with source residues and terminal
words computable without expanding its enormous integers. Stabilization
of `(4^(3^L*s)-1)/3^(L+1)` modulo `3^j` removes dependence on huge `L`
once `L>=j-1`.

The source coordinate is essential. `refine_family` constructs a
subfamily; `extend_point` retains the exact supplied parameter and fails
if its new guard is false. Seed 5 at parameter 0 can extend to 3, but
another `G1` extension of that same point is rejected. A legal different
parameter is not a certificate for the original integer.

The two selected children in section 1 are parameter-0 points of this
grounded language, with phases `s=0 mod9` and `s=0 mod19683`. Their
positive-parameter companions are also fully grounded. This is stronger
than merely producing smaller unproved obligations, while still leaving
arbitrary-source membership open.

For these **zero-phase** children, there is exact recursive reseeding:
the whole child family becomes the same canonical completed-lift
construction based at its newly certified child, with parameter `t=s/3^ell`.
Every source and ROOT word is preserved. This permits literal recursive
reuse of the small construction after its surface form has changed.

The boundary is sharp. At a **nonzero phase**, rebuilding the canonical
completed lift using only the child's ROOT word shares just one integer
with the old child family: their common base point. A modulo-4 argument
proves that all other members differ. The original phase and exponential
coefficient must remain as sidecar data. This is an explicit reason to
retain history, beyond merely saying that it might be useful.

The parameter density of a fixed word is `3^(-ell)`; it is positive at
every finite stage but tends to zero along unbounded refinement. This
does **not** yet give the sought source-specific positive lower bound
that survives all refinement. See [the complete grounded-family
proof, implementation, and hostile controls](collatz_grounded_child_closure_20261007.md).

## 4. What the seven supplied papers changed in the construction

The paper headlines were accepted as supplied. The useful transfers were
specific proof interfaces, each given its own elementary Collatz proof:

| Mechanism extracted | Resulting arithmetic requirement or construction |
|---|---|
| Campana–Peternell overlap witnesses and Deligne–Drinfeld projection injectivity | Keep the ordered affine carry and its native cylinder; prove the word decoder rather than infer it from a slope. |
| Cubic-torus endpoint propagation with a retained inequality | Prove a common height bound for every legal child step; distinguish termination at a frontier from ROOT grounding. |
| Hilbert–Smith fixed integral lattice | Keep the source fixed while using valuation depth; changing the source changes the bound. |
| Homotopy exact boundary gluing | A reversible prefix transports a supplied terminal witness and can also transport an empty witness set. |
| Lipschitz-height conditioning and six-vertex fixed-word limits | Intersect guards with the retained boundary; separate a fixed receiver bank from adaptively growing receiver depth. |

Exact filenames and inspected pages are in the
[three-paper normal-form note](collatz_child_normal_forms_20261007.md#6-three-supplied-papers-exact-mechanisms-and-transfer-boundaries)
and [four-paper receiver note](paper_child_closure_transfers_20261007.md#1-inheritance-and-the-four-paper-interfaces).
No geometric, stochastic, or representation-theoretic hypothesis has
silently become a theorem about integer Collatz orbits.

The receiver note also gives a lossless total decoder of every maximal
initial `(1,2)` run:

\[
r=\left\lfloor\frac{v_2(n+5)-1}{3}\right\rfloor,
\qquad m=\frac{9^r(n+5)}{8^r}-5.
\]

It retains the original source through receiver composition. The actual
path `3067 --(1,2)^3--> 4369 --(2)--> 3277` shows why a decrease from
the changed core is insufficient. A fixed finite bank of actual forward
receivers cannot pay arbitrary forced-prefix depths. Equation (1) supplies
the broader common-future interface instead of discarding that obstruction.

## 5. Precise continuation targets

The concept board after the new computations is: immutable source;
ordered carrier; ternary phase; deletion budget; grounded seed coordinate;
terminal obligation. Each new construction now has to preserve the
coordinates needed by the others.

1. **Adaptive paid receiver:** for a forced child word with carrier
   `(A,B,C)`, construct an authenticated deletion of depth `D` whose
   native phase includes that same supplied parent, and prove (1).
   Merely choosing a numerically adequate `D` does not construct a receipt.
2. **Grounded membership:** recognize or route an uncovered supplied
   source into one of the completed families without changing its parameter
   or assuming its ROOT certificate. Ternary isometry solves guard matching
   inside a family, not entry of every source into one.
3. **Well-founded union:** join the paid Mersenne exits and grounded mixed
   families with a rank that decreases across changes of representation.
   The full graph needs rules at its remaining terminals. Finite-word
   closure and residue frequencies alone do not establish that coverage.

These are the missing implications toward a universal proof. The newly
grounded points, exact family compilers, and finite checks do not hide them.

## 6. Reproducibility checkpoint

The four packages reproduce **574,517 explicit exact checks** in normal
and optimized Python modes, with matching saved outputs:

| Package | Checks | Principal independent control |
|---|---:|---|
| [Mixed paid exits](collatz_mixed_child_rules_20261007.out) | 23,593 | Actual two-sided receipt replay, funded/unfunded native cells, exact phase arithmetic |
| [Child normal forms](collatz_child_normal_forms_20261007.out) | 87,614 | Generator-by-generator replay, full-carrier decoding, complete finite residue universes |
| [Grounded closure and reseeding](collatz_grounded_child_closure_20261007.out) | 322,863 | Literal ROOT replay, independent old recognizer, direct modular exponentiation, reseeding identity/hostile |
| [Prefix/receiver interfaces](paper_child_closure_transfers_20261007.out) | 140,447 | Literal decoder on all positive odds below20,000 and gcd-aware guard comparisons |

Each linked note specifies its universe, positive and hostile controls,
normal/optimized reproduction commands, and saved output hash. Independent
peer review covered the proofs and code. These finite controls support the
implementations; the quantified family statements use the displayed proofs.
The repository documentation gate and whitespace check pass.
