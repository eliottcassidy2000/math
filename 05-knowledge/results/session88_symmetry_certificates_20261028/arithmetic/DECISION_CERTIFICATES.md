# Keep the state, shrink the proof: source-faithful dual certificates

**Session 88, 2026-10-28. PROVED elementary theorem and applications;
FINITE-EXACT controls.** No external-priority claim. This is not a general
short-NO-certificate theorem for integer programming, primality-constrained
families, graceful-tree existence, or arbitrary relational joins.

## Result first

A state sufficient for **every named future operation** and a certificate
for **one decision about one supplied input** solve different problems.
Here is a precise useful replacement for “keep every guard”:

> For a system of integer linear equalities of rational rank `r`, failure
> has a dual proof using at most `r+1` source equations. Rational failure
> has a row-dependence witness; purely integral failure has a cyclic
> congruence witness using at most `r` source equations. The latter may
> require a prime power, not a prime. Either witness composes by adding
> boundary covectors, without listing every feasible boundary state.

The proof is classical linear algebra/Smith theory, organized here as a
source-audited certificate interface. Its uses below are concrete:

* fixed aliquot core words emit modular or positivity NO receipts without
  searching for the unknown primes;
* signed incidence distinguishes cycle inconsistency from integral torsion;
* one proposed graceful signed permutation has a short path NO witness,
  although rejecting one proposal does not reject its existential family;
* square-`sigma` YES and NO receipts can use authenticated source
  factorization and adjacent squares, not just a binary parity bank.

The certificate is checked **against the original input**. Its size does
not include that input; it does not replace it for future operations.

## 1. Inheritance, corrections, and the live board

Sources were read at repository pin
`7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51` (current materialized navigation
was used only for routing). Canon IDs below include their exact slugs.

| Role | Correction-aware source |
|---|---|
| Closest module mechanism | `05-knowledge/results/relation_gluing_7_15_20261028/THEOREMS.md`, T1, T4, T6: exact image-intersection gluing, torsion under coefficient change, weighted assembly |
| Least-used old sidecar | `05-knowledge/results/aliquot_guarded_cycles_20261014e.md`, §2: unreduced positive affine carry for a fixed core word |
| Original-choice/LCA discipline | `05-knowledge/results/session56_rmq_lca_connections_20261025/correlated-gluing/THEOREM.md`, §§1–6: compatible original assignments, once-only source weights, lossless-join hypotheses |
| Signed integration | `01-canon/theorems/THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion.md`, §§1–2: graceful iff every nonempty path has signed sum with absolute value in `1..m` |
| Signed holonomy | `05-knowledge/results/free_bits_signed_incidence_20261028/THEORY.md`, P4: chords impose root compatibility and integral parity |
| Exact future-state contrast | `05-knowledge/results/session84_guarded_unification_20261028/GUARDED_OBSERVERS.md`, Proposition 3 and results §4: aliquot observer with named coprime multiplication distinguishes exact `N`; the defect observer alone instead retains `(rad(N),tau(N))`, and named exponent promotions refine it to exact `N>=2` |
| Representation correction | `01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md`, front correction; `05-knowledge/results/free_bits_seven_fiftysix_20261028/INCOMING_REPAIR.md` |

In particular, we do **not** infer a signed-graph model from regularity, a
TU signing of a prescribed matrix from regularity of its unaugmented
column matroid, or equality of rational and binary dependencies. We do not
inherit the obsolete OPEN dyadic-three-cycle status from the body of the
old frontier: its front update routes to
`05-knowledge/results/dyadic_aliquot_three_cycle_exclusion_20261019/DYADIC_THREE_CYCLE_EXCLUSION.md`.
Our examples use the older compiler, not that obsolete status.

Live board: **source equations / boundary covectors / torsion / unique
candidate / native predicate / source weight**. The useful change of
vertices is from objects being assembled to the equations certifying an
obstruction. Canonical hostiles are the false `46 -> 26 -> 16 -> 46`
aliquot cycle, the signed path `(+1,-3,+2)`, and the non-lossless three-bit
parity projection. The new hostile `4x=2` separates prime-field tests from
integral decision.

## 2. Sparse dual alternative (PROVED)

Let `A` be an integer `m x n` matrix, `b in Z^m`, and
`r=rank_Q(A)`. Rows are the actual source equations, with their identities.

**Theorem.** Exactly one of the following two sides is possible:

1. **YES:** `x in Z^n` with `Ax=b`.
2. **NO:** at least one of these explicitly verifiable certificates:
   * **rational obstruction:** an integer row vector `lambda` with
     `lambda A=0`, `lambda b != 0`, and at most `r+1` nonzero entries;
   * **torsion obstruction:** an integer `q>=2` and row vector
     `y mod q` with `yA=0 mod q`, `yb!=0 mod q`, and at most `r`
     nonzero entries.

If `Ax=b` is rationally consistent but has no integer solution, the second
kind always exists. Choose any `r` independent rows `S`. The modulus can
be taken to divide the order of

    Z^r / im_Z(A_S),

and hence to divide the absolute determinant of every nonzero `r x r`
column minor of `A_S`. It may be chosen a prime power. When `r=0`,
rational consistency already gives an integer solution, so this clause is
vacuous. The two NO kinds are not asserted mutually exclusive.

### Proof

Soundness follows by applying the displayed linear combination to `Ax=b`.
For the rational obstruction, select a row basis of `A`. If its equations
are inconsistent with the full system, one additional row disagrees with
the corresponding linear combination of the basis rows. This uses at most
`r+1` rows; clear denominators. On a support-minimal choice, it is a row
circuit whose right-hand sides do not satisfy that circuit relation.

Otherwise the system is rationally consistent. Every omitted equation,
including its right-hand side, is a rational combination of the equations
on `S`. Thus these `r` equations have exactly the same rational solution
set, and consequently the same integer solution set, as the full system.
The map `A_S:Z^n -> Z^r` has full row rank, so its cokernel is finite.
A nonzero class of `b_S` in that cokernel is detected by some cyclic
coordinate of a Smith decomposition. Equivalently, a homomorphism to
`Z/q` kills every column of `A_S` but not `b_S`. Any homomorphism from
`Z^r` to `Z/q` is dot product with a row `y`; extend it by zero outside
`S`. This gives the second certificate. A nonzero element of a cyclic
group remains nonzero in at least one prime-power component, so a
prime-power `q` suffices.

For any invertible column minor `C`,
`im_Z(C) subseteq im_Z(A_S)`. Hence the order of the latter cokernel
divides `[Z^r:im_Z(C)]=|det C|`; this gives the stated modulus bound. QED.

Over an arbitrary field `K`, the same basis argument gives the familiar
alternative: `Ax=b` or `lambda A=0, lambda b!=0`, with at most
`rank_K(A)+1` source rows. Coefficients and rank must be taken in **that
field**; a binary certificate is not a rational certificate.

### Length, verification, and construction are different invoices

If `|A_ij|<=H` and `r>=1`, Hadamard gives
`q <= r^(r/2) H^r` for the torsion case. A certificate can store source
row indices and residues in `[0,q-1]`, with bit size

    O(r log(m+1) + r^2 log(rH+1)).

Rational row circuits likewise admit minor-based integer coefficients
of polynomial bit length. This is a bound for **unrestricted integer
linear equalities**, not for bounded variables or nonlinear predicates.
Sparse means few named equations, not bounded coefficient magnitude or
geometric locality.

Verification requires the indicated weighted sums against **every column
of the original selected rows**, plus the right-hand side. For explicitly
stored rows it uses `O(nnz(A_S)+|S|)` ring operations, with operand bit
sizes charged separately. Standard exact elimination/Smith methods
construct witnesses; this note proves no runtime bound for a particular
implementation. The test program deliberately uses exhaustive small
character search as an independent path, not as a claimed scalable Smith
algorithm. Compiling the native problem into a complete matrix may itself
be expensive, and the native problem may not be linear.

Two indispensable hostiles:

* `4x=2` is rationally solvable and solvable modulo **every prime**:
  modulo 2 it is `0=0`, and modulo odd primes 4 is invertible. But modulo
  4, `y=1` gives `0=2`. Checking all prime fields is not integer feasibility.
* For `A=[2 3]`, `b=1`, the chosen column `2` has determinant 2, but its
  putative mod-2 annihilator does **not** kill column 3. The full system has
  integer solution `(-1,1)`. A selected minor bounds a legitimate modulus;
  its adjugate is not automatically a full-lattice character.

Also `x1+x2=3` has an integer solution but no solution with both variables
in `{0,1}`. The theorem provides no NO witness for that additional domain.

## 3. A dual gluing law, not another full state table (PROVED)

Put source rows at the leaves of any assembly tree. For a certificate
`lambda` (integer, modular, or field-valued), let

    C_v = sum_(i below v) lambda_i A_i,
    d_v = sum_(i below v) lambda_i b_i.

At a join, `(C_v,d_v)=(C_left+C_right,d_left+d_right)`. At the root a NO
certificate has `C_root=0` and `d_root!=0` in the chosen coefficient
ring. This is a locally checkable proof by cancellation.

For a valid full certificate, `C_v` is supported only on variables
occurring in rows both below and outside `v`: a variable occurring only
below `v` has its entire global coefficient there, and it must be zero.
Thus the message may be just a **boundary covector and a scalar**, with
vanishing of every internal coordinate checked locally. The last
contributions to a variable cancel at the LCA of its source-row support.
No running-intersection hypothesis is needed if this actual crossing
boundary is retained. A smaller purported separator needs its own proof.

A table of all `p^k` possible boundary assignments is not required to
verify one supplied covector on `k` boundary coordinates. This does **not**
assert that every faithful representation needs such a table: affine
relations themselves often have compact matrix descriptions. Nor does it
make the whole certificate constant-sized; source multipliers, local
checks, coefficients, and the original input remain.

### Signed incidence exposes both kinds of obstruction

Use equations `x_u+x_v=w_e`. The triangle with edge data `(4,4,1)` has
rational solution `(1/2,7/2,1/2)`. Adding all three rows modulo 2 gives
`0=1`, a torsion certificate. This is the inherited odd-cycle factor 2.

Now take two triangles sharing vertex 0, with rows in order

    01,12,20,03,34,40,
    w=(2,2,2,2,2,4).

Each triangle alone has an integral solution: the left pins `x0=1`, the
right pins `x0=2`. All individual odd-cycle parity tests pass. Nevertheless

    lambda=(1,-1,1,-1,1,-1)
    lambda A=0,    lambda w=-2.

The row rank is 5; all six equations are needed for this rational
inconsistency circuit, so the `r+1` bound is sharp here. Every deletion of
one row is rationally feasible. Splitting at the two triangles produces
boundary covectors `+2 e0` and `-2 e0`, with scalar charges 2 and -4.
Their LCA cancels the shared variable and leaves `0=-2`.

In contrast, the valuation matrix for source primes `(7,19,89)` is

    [3 2 1]
    [0 0 2]    rows v2,v3,v5; determinant -6.
    [0 1 1]

It has no rational kernel but `(1,1,1)` is a binary kernel vector. This is
exactly the torsion distinction inherited from relation-gluing T4, not a
new Fano-code assertion or a claim that the same dependencies lift to Q.

## 4. Aliquot core words become proof-producing templates (PROVED)

For an ordered cyclic word `c_i>=2`, put
`S_i=sigma(c_i)`, `a_i=S_i-c_i>0`, `d_i=c_(i+1)`. The marked prime model
has the source equations

    d_i p_(i+1) - a_i p_i = S_i.

For the actual aliquot interpretation one additionally needs `p_i` prime,
`gcd(p_i,c_i)=1`, and distinct positive source nodes `c_i p_i`.

Let

    D = product c_i - product a_i,
    w_i = (product_(j>i) a_j) (product_(j<i) d_j),
    B = sum_i w_i S_i > 0.

Here empty products are 1. In the original row order,

    w^T A = D e_0^T,       w^T S=B.                     (1)

**Proof:** every intermediate `p_i` coefficient cancels between its
incoming and outgoing equations; the remaining coefficient of `p0` is
`product d_i-product a_i=D`. This is exactly the old compiler's positive
carry, now read as a proof multiplier rather than only a candidate formula.

Rotating the source word gives the same identity for each coordinate:
`D p_j=B_j`. Therefore:

* `D<=0` gives a positivity NO receipt, since `D p0=B>0` contradicts
  `p0>0`. At `D=0` it is also a rational row certificate; at `D<0` it is
  **not** rational inconsistency and the positivity premise is essential.
* If `D>0` and some `D` does not divide `B_j`, reducing its source
  multiplier modulo `D` gives the theorem's integer-lattice NO receipt.
* Otherwise the entire rational candidate is uniquely forced. A failed
  prime, coprimality, or distinct-node test rejects that unique candidate,
  hence the **fixed template**. This does not certify all core words.

For `(c0,c1)=(2,4)` the equations are

    -p+4q=3,   2p-3q=7.

The row combination `(3,4)` gives `5p=37`. Modulo 5 this is the two-row
certificate `0=2`; no search cutoff on primes is involved.

For `(2,2,16)`, the unique candidate is `(23,13,1)`. All affine/integer
checks pass, but the last marked prime is 1. The proposed source loop
`(46,26,16)` fails natively because `s(16)=15`. Uniqueness is what permits
one failed native factor to reject the whole fixed template. If the
linear solution space had several points, one failed candidate would not
suffice. A positive control is `(20,4)` with marked primes `(11,71)`, giving
`220 <-> 284` after actual divisor-sum verification.

The parameter word is still an exact source object. This bridge does not
prove general aliquot termination or exclude all three-cycles; the old
special dyadic exclusion is a separate theorem with its own hypotheses.

## 5. Graceful signed candidates: a path is a NO certificate

This application inherits, rather than reproves as new, **THM-2787** at
the path in §1. On a rooted tree with `m` edges, choose one signed
permutation `d_e` of `1,...,m`, and integrate it to root-zero potentials
`x_v`. Tree incidence is unimodular: every such `d` integrates to integer
potentials, so the linear system itself has no obstruction.

The chosen candidate is graceful after translation iff the potentials
are distinct and their span is at most `m`. Failure therefore supplies
vertices `u!=v` and their unique tree path with

    sum_(e on path u->v) signed d_e = 0,
    or absolute value of that sum > m.

This is a telescoping dual evaluation: if `h` is the oriented path row,
then `h partial_T=e_v-e_u`. It uses at most `m` edges and is verified by
one signed sum and the displayed failed bound. Integration plus duplicate
and extrema detection constructs a witness, without checking all paths.

Examples are `P3` with `(+1,+2)` (span 3 exceeds `m=2`) and `P4` with
`(+1,-3,+2)` (the endpoint path sums to zero). At an LCA, the shared root
prefix cancels; branch signs and order are not dropped.

This is a short NO certificate for a **supplied signed permutation**.
Gracefulness of an unlabeled tree asks whether **some** signed permutation
succeeds. One path only rejects one candidate; rejecting all candidates
may require far more work. Neither the Graceful Tree Conjecture nor the
frozen-fiber/global-connectivity distinction is resolved here.

## 6. Square-sigma receipts beyond a parity bank

The predicate here is **`sigma(N) is a square`**, not squarehood of
` s(N)=sigma(N)-N`. For any supplied positive integer `N`, provide

1. distinct source primes and exponents `N=product p_i^a_i`;
2. a primality certificate for each `p_i`;
3. either `r^2=S` or `t^2<S<(t+1)^2`, where

       S=product_i (1+p_i+...+p_i^a_i)=sigma(N).

The first equality is a YES receipt and the strict adjacent-square
interval is a NO receipt. For `N=1`, use the empty factorization and
`S=1`. No factorization of `S` is needed. Examples: `N=70` has `S=144`
and `r=12`; `N=30` has `S=72`, between `8^2` and `9^2`.

One fully elementary prime-certificate verifier is the recursive
Lucas/Pratt-type order proof. For `p>2`, supply a complete prime
factorization of `p-1`, recursively authenticate its prime factors, and
a value `g` such that

    g^(p-1)=1 mod p,
    gcd(g^((p-1)/q)-1,p)=1 for every distinct prime q | p-1.

For any prime `ell|p`, these conditions force the order of `g mod ell`
to be exactly `p-1`, hence `p-1 | ell-1`; therefore `ell>=p` and `p` is
prime. Base case 2 is direct. Such certificates exist for primes because
the multiplicative group of a finite prime field is cyclic. This is
classical prime certification, not a new primality theorem.

Verification is polynomial in the **written input-and-certificate bit
size**, using repeated-squaring arithmetic. A verifier first bounds each
source exponent by `bit_length(N)` and every recursively claimed exponent
by `bit_length(p)` and checks strictly smaller child primes; absurdly huge
binary-encoded exponents must not trigger enormous unguarded expansions.
Constructing the source factorization and certificates is a separate
cost. The included tiny generator uses trial division and makes **no**
polynomial-time factorization or generation claim. Nor do we claim these
receipts are size-minimal or that a decision procedure must use them.

For a free squarefree subset bank, the older binary alternative is still
useful. If prime 2 is forced and only primes 3 and 7 are optional, every
possible source has odd `v3(sigma(N))`; that one parity functional
annihilates both optional columns and rejects all four completions.
Allowing primes 5 and 7 instead gives `N=70`, a positive completion.
Exact one-hot exponent choices or repeated-prime fusion are not free
binary variables; their original constraints/corrections remain needed.

Session84's exact-state theorem concerns the **aliquot observer** with
named coprime multiplications, not this square-sigma predicate. Its proof
recovers `N` from `s(N)` and `s(Np)` at one common fresh prime. Consequently
there is no conflict between full contextual separation in that language
and these decision-specific receipts. We make no claim that the
square-sigma predicate itself has that same contextual quotient.

## 7. Source weights and the boundary of the transfer

Proof multipliers `lambda_i` are not probability weights. On original
assignment relations, retain the inherited compatible-union rule and
assign each original source weight once. A shared choice with weights
`(2,3)` has total 5; joining two occurrences and multiplying their scalar
weights gives the incorrect total `2^2+3^2=13`.

Boolean assembly decides existence. With strictly positive real weight
on every feasible global assignment, total weight is positive iff the
join is nonempty. With zero weights, total zero excludes only
positive-weight assignments; with signed weights, `(1,-1)` can sum to
zero although two feasible assignments exist. A zero weighted answer is
not automatically a NO certificate.

The dual theorem applies when the complete retained constraints are
linear equalities over the stated field or over unrestricted integers.
Affine local blocks may be stacked with literal shared variables and
then certified; arbitrary local tables, hidden global factors, native
physical realization, positivity, all-different conditions, and
primality are not thereby converted into linear algebra. Session56's
cyclic pairwise-consistent but globally empty parity triangle is exactly
linear over F2 and does admit a three-equation dual contradiction;
arbitrary cyclic constraints need not.

## 8. What the exact checker establishes

Run from this directory:

    python3 -B check.py
    python3 -B -O check.py

`normal.json` and `optimized.json` are byte-identical. All checks use
explicit failures, not removable assertions. Current replay: **360,819
active checks**.

* All 15,625 integer `2 x 2` systems with entries and RHS in `[-2,2]`:
  7,313 integer YES; 2,648 rational NO; 5,664 torsion NO. Independent
  determinant/gcd membership is compared to exact rational elimination
  and exhaustive modular witness search, checking every original column.
* All 256 binary `2 x 3` systems: exhaustive primal and dual alternatives.
* The signed triangle, the minimal six-row/rank-five bowtie circuit, and
  the integer-versus-binary valuation lift; the selected-minor and
  prime-power hostiles.
* All 3,615 ordered core words in `{2,...,16}`, lengths 1–3 (no rotation
  quotient): 3,415 integrality NO, 82 positivity NO, 98 native-guard NO,
  16 repeated-node NO, and four marked YES templates. Rotated source
  multipliers are independently compared with rational elimination;
  accepted native edges use direct divisor-pair sums.
* All 100,724 signed-permutation candidates on rooted paths and stars
  with 1–6 edges: 1,820 YES; 24,108 zero-path NO; 74,796 span-path NO.
  These are two tree families, not an all-tree census.
* Every `1<=N<=4096`: authenticated source factorization, recursive prime
  certificates, product sigma versus an independent divisor-pair sum,
  and square/adjacent-square receipts (138 YES, 3,958 NO).
* Once-only source weighting, signed cancellation, and zero-weight
  nonemptiness controls.

Finite tests do not prove the general sparse-dual theorem or any universal
native exclusion; the proofs above establish exactly their stated scopes.

### Independent review

A separate Session88 auditor accepted the sparse-dual proof, the
prime-power/minor bound, boundary-covector cancellation, aliquot multiplier
identity, and authenticated square-sigma receipt argument. Its independent
HNF/determinant implementation imports no producer code and rechecks all
15,625 `2 x 2` systems above, plus all 19,683 `3 x 2` systems with entries
and RHS in `[-1,1]` (5,065 integer YES, 14,330 rational NO, 288 torsion NO).
It constructs a modular certificate for **every independent selected row
set** in each integral-NO case, within the smallest nonzero maximal-minor
bound. The independent checks used a separate implementation.
The auditor inspected the native-application code but did not independently
exhaustively replay those three native censuses; their stated independent
calculation paths are the producer's internal controls, not a separate-agent
reproduction.

**Outcome:** full future state, solution set, current decision, proof
object, verifier, and proof construction are six separate interfaces.
The new useful operation is to pass from a potentially large primal
compatibility problem to a source-authenticated dual contradiction—not
to pretend that one short contradiction stores every future state.
