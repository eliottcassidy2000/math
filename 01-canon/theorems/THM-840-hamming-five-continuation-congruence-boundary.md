---
id: THM-840
title: Exact endpoint geometry is insertion-Markov while the Hamming-five handoff quotient is not
status: PROVED (operation-kernel congruence, minimal linear continuation-sidecar rank, and exact insertion/deletion boundary) + FINITE-EXACT (handoff, residual, maximin, and tournament replay)
source: codex-2026-07-15-S10 continuation
depends_on: [THM-822, THM-828, THM-832]
related: [THM-837, THM-2230, THM-2237, THM-2240, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_continuation_congruence_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_continuation_congruence_codex_S10.out
---

# THM-840 — the Hamming-five continuation-congruence boundary

This theorem separates static reconstruction, present evaluation, and
deterministic continuation.  The distinction is operation-relative: an
observation may be an exact Markov state for adding a comb and fail completely
for deleting that same comb.

## 1. The operation-kernel criterion

For an observation `O:X->Y` and a family of operations `T_a:X->X`, an update
map `F_a` satisfying

```text
O(T_a x)=F_a(O(x))
```

exists if and only if

```text
ker(O) subset ker(O after T_a).                           (1)
```

Here the update need only be defined on the reachable image `O(X)`.

### Proof

If `F_a` exists and `O(x)=O(y)`, then

```text
O(T_a x)=F_a(O(x))=F_a(O(y))=O(T_a y),
```

which proves the kernel inclusion.  Conversely, define

```text
F_a(O(x)):=O(T_a x).
```

Condition (1) says exactly that this definition is independent of the chosen
representative `x`.  Thus it is a well-defined update on `O(X)`.  This proves
the criterion. ∎

Static injectivity or gluing is a different question.  It asks whether some
family of observations reconstructs `x`; (1) asks whether one observation's
equivalence classes are congruences for the named operations.

## 1a. Exact rank of a linear continuation sidecar

There is a quantitative linear form of the same criterion. Let

```text
R:V -> Y,                    N:V -> Z                 (1a)
```

be linear maps of vector spaces over one field. Here `R` is the present
response and `N` is one named next observable. A linear sidecar

```text
S:V -> W                                               (1b)
```

makes `N` recoverable from `(R,S)` if and only if

```text
ker(R) intersection ker(S) subset ker(N).             (1c)
```

When (1c) holds, the induced recovery map on `im(R,S)` is automatically
linear. Define the continuation-defect rank

```text
d_R(N)=dim N(ker R).                                   (1d)
```

Then the minimum possible rank of a repairing linear sidecar is exactly

```text
min_S rank(S)=d_R(N).                                  (1e)
```

This statement includes infinite dimensions: if `d_R(N)` is infinite, no
finite-rank linear sidecar can repair the response for this named
continuation.

### Proof

Condition (1c) is the operation-kernel criterion applied to the joint
observation `(R,S)`. If it holds, restriction to `ker R` makes `N` factor
through `S|_(ker R)`. Therefore

```text
dim N(ker R)<=dim S(ker R)<=rank(S),                  (1f)
```

which proves the lower bound.

Conversely, restrict `N` to `ker R` and identify

```text
ker R/(ker R intersection ker N) ~= N(ker R).        (1g)
```

Use this quotient map as `S` on `ker R`, and extend it linearly from
`ker R` to all of `V`, keeping codomain `N(ker R)`. Its rank is
`d_R(N)`, and its kernel on `ker R` is exactly
`ker R intersection ker N`. Thus (1c) holds and proves (1e).

The ambient size of `V/ker R` is irrelevant. Only the image of the present
kernel under the named future observable measures the required memory.
Accordingly:

```text
d_R(N)=0       means the present kernel is a genuine gauge for N;
d_R(N)>0       means it contains control/memory directions.             (1h)
```

Two current applications delimit the scope.

- In THM-2237, truncated Boolean moments leave one affine direction, and
  the even/odd parity laws prove that the top atom is nonconstant on it.
  Hence the missing top-Walsh scalar has continuation-defect rank exactly
  one.
- In THM-2240, take `R=d_6` on pure grade-six Ore corrections and let `N`
  be the change in the grade-seven residual. The arbitrary-`q` syzygy axis
  maps injectively under `N`. Thus `d_R(N)` is infinite. On the span of
  parameters

  ```text
  C(q,u)=sum_(r=0)^R sum_(s=0)^n c_(r,s)q^r u^s,
  ```

  the independent `q` layers and the injectivity of each layer already
  force sidecar rank at least `(R+1)(n+1)`.

This does not prove Ore or Weyl nontermination: THM-2240 permits
state-dependent next-rung corrections. It says only that no fixed
finite-rank linear statistic of the unrestricted current fiber makes that
specific one-step response exact. Conversely, THM-2230's target-shear
kernel is a true gauge for every downstream predicate proved invariant
under target shears; for nonlinear predicates this is the set-theoretic
kernel criterion, not a use of the linear rank formula.

## 2. Exact endpoint geometry is Markov for monotone addition

For pure LRC comb insertion, the ordered exact open-endpoint word is such a
state because

```text
E(S union {u})=E(S) intersect {t:||ut||>1/13}.            (2)
```

Indeed, the endpoint word determines the literal open set `E(S)`, and the
right side of (2) is obtained by an exact merge scan with the known teeth of
the incoming speed `u`.  Therefore equality of endpoint words is preserved
by every common insertion.  The update decides the new endpoint word and
whether the strict residual is empty, at every height and without owner data.

This positive statement is deliberately narrow.  The endpoint word does not
record the labelled combs which formed it, and threshold geometry need not
determine the exact maximin value.

## 3. `H0=H1` is not insertion-Markov

THM-822's coarse `H0=H1` observation is not a congruence even for one named
insertion.  Take missing labels `(1,2,3,4,5)`.  The rows with
replacement speeds `(14,15,16,17,18)` and `(27,28,29,30,31)` share `H0=H1`.
Their common complete decorated live relation is

```text
{(2->1,k=2),(3->1,k=3),(4->2,k=2)}.                      (3)
```

Now insert the same labelled replacement `(q,u)=(6,19)`.  In the low row the
new incident handoffs are

```text
{(6->2,k=3),(6->3,k=2)},                                 (4)
```

whereas in the high row they are

```text
{(5->6,k=3),(6->3,k=2)}.                                 (5)
```

These are exact consequences of the half-open rule

```text
-u_owner < z u_owner-2u_provider-13m u_owner <=u_owner.  (6)
```

Thus the two inputs lie in one `H0` and `H1` fibre, but their images under the
same insertion do not.  Criterion (1) proves that neither quotient has a
deterministic insertion update.  Adding the bounded-bank integer centre does
not help because it is already a function of the live labelled edge there.

## 4. Endpoint geometry is not deletion-Markov and does not determine `M`

The limitation of the positive endpoint result is exact even at fixed
cardinality.  Put

```text
S={1,2,...,13},
T={1,2,3,4,5,6,7,8,9,11,12,13,20}.                       (7)
```

Both strict `1/13` residuals are empty, so their exact endpoint words agree.
Nevertheless,

```text
M(S)=1/14,                 M(T)=2/27.                     (8)
```

Delete the common speed `13`.  Then

```text
E(S minus {13})=empty,        M(S minus {13})=1/13,       (9)
```

while

```text
E(T minus {13})
 = (79/260,4/13) union (9/13,181/260),
M(T minus {13})=2/23.                                   (10)
```

Hence the kernel of the endpoint observation is not a congruence for common
deletion, and the threshold endpoint word does not determine exact `M` even
before deletion.  A deletion-, replacement-, or transport-ready state must
retain the labelled active tooth bank or equivalent proof-obligation
ancestry.  This is precisely the extra coordinate retained by THM-837's
recursive carrier.

## 5. Static Cech data versus acted state

THM-822's three exact `H2` faces are injective on the height-at-most-two bank,
and THM-828/832 exhibit a constant rank-four Cech core with trivial positive
cohomology in the separate n=9 codec.  Neither fact implies (1).  Static
gluing can reconstruct a present input while a proposed quotient still fails
to update; conversely, literal endpoint geometry can update under addition
without reconstructing the labelled history that produced it.

The operation-indexed hierarchy is therefore

```text
H0/H1:                    not insertion-Markov;
exact residual endpoints: insertion-Markov, not deletion-Markov;
endpoints + tooth bank:    recompute additions and deletions;
owners/transport sidecars: needed when the teeth themselves move in a deck.
```

This is a theorem about sufficiency for named operations, not a claim that
the last line is a globally minimal LRC state.

## 6. Tournament Analysis

Put replacement labels at vertices, use the directed half-open handoff as the
pairwise observable, switch from the low-lift to the high-lift branch after
the common insertion, and complete silent pairs along the declared tie
Hamiltonian path `(1,2,3,4,5,6)`.  The child fingerprints are

| gauge | score histogram | triangles | SCCs | Hamiltonian paths |
|---|---|---:|---|---:|
| low | `{1:1,2:1,3:4}` | 7 | `(6)` | 33 |
| high | `{1:2,3:3,4:1}` | 5 | `(6)` | 23 |

The live supports differ on `(6,2)` versus `(5,6)`, although the completed
tournaments differ on only one edge.  Tournament completion preserves one
chosen pairwise orientation but destroys live-versus-tie provenance, integer
centres, residual geometry, and threshold truth.  The exact continuation
carrier is not a tournament.

## 7. Verification

The verifier uses only integer inequalities and `fractions.Fraction`.  It
reconstructs (3)--(6), all four endpoint/maximin claims (8)--(10), the
operation-kernel truth table, and the tournament fingerprints.  Its payload
and certificate hashes are

```text
b5c04163a21209f82231228b6fcb8c4eb5beb96fe5d254157253d21f21d14f27,
b9580c2c4793555ae76ae2c5bea5a2cc0b22a1f6bc2900feff6acf1735ce63d2.
```

Run

```bash
python3 04-computation/lrc13_hamming_five_continuation_congruence_codex_S10.py
```

and compare byte-for-byte with the stored output in the frontmatter.
