# Session75 incoming reconciliation (bounded, read-only)

2026-10-27. Reviewed incoming commit `235333940cab42a6b4ef7126a742c7cf040e91b0`
and close-out messages `c8f2940bdbce85ac2367f5c490ac1892f6aabf9f`, against the
frozen Session75 producer packets and their separate independent audits.
Local publisher HEAD observed at `5c570e1a0ff60f5672480cdbaee6cd185ea18c39`.
No shared edits, rebase, remote action, or new research experiment was made.
This is a targeted reconciliation, not an audit of the incoming enumeration.

## 1. HYP-9422(3): resolved in its fixed-bank substitution scope

Exact source: HYP-9422, lines 42–46; incoming README §3, lines 135–155;
THM-4814(C). Full paths, git blobs and SHA-256 hashes are in `SOURCE_PINS.json`.
The hypothesis expressly asks about **individual disagreeing named pairs** of
Paley versus bounded-arity substitution parses, including join/C3 and fixed
banks. It is not a claim about block-exception occurrences or native moves.

For effective maximum arity `b>=2`, Session75 `session75-spectral/RESULT.md`
§§1–4 proves, uniformly over labels, irregular children and all parse shapes,

    e >= max(0,ceil((q(q−1)−C_b q^(3/2) log₂q)/4)),
    C_b=2sqrt(α(1−α))/h₂(α), α=1/(b+1).

The exact input is `A Aᵀ=qI−J`, `||A||op=sqrt(q)` and
`<A,B>=q(q−1)−4e`; operator/nuclear duality gives the bound. Transitive
baselines give `e_min<=floor(q(q−1)/4)` by averaging over total orders.
Thus for **fixed b**, `e_min=q²/4−o(q²)` for join/C3 or a fixed bank augmented
by ordered joins. No claim that it equals the finite transitive optimum is
needed. A bank lacking baselines at some orders needs its admissible orders
specified; the transitive-upper proof assumes transitive baselines are allowed.

Scope checks:

- Children are actual disjoint modules with all internal tournaments retained.
  Any irregular/unbalanced modules with lawful recursive arity-b parses are
  covered. Arbitrary undecomposed large modules are not free leaves.
- Reused templates expand to disjoint occurrences. Repeated port *references*
  are harmless; physically identifying vertices across modules is not licensed.
- `b=b(q)` keeps the finite formula but not the fixed-b conclusion. A sufficient
  asymptotic condition is `C_b log q=o(sqrt(q))`; no uniform arbitrary-growing-bank
  assertion is made. A bank containing Paley as one quotient defeats such an assertion.
- Paley means the finite-field construction for every prime power `q=3 mod4`,
  not modular arithmetic modulo a composite integer. This includes q=27.
- Edge flips, rank-two additive corrections, block reversals, and fixed-size
  path complements have distinct tariffs; only the first is HYP-9422(3).

**Suggested C2 status:** `(3) PROVED in the fixed-b disjoint-substitution,
individual-pair scope by Session75 spectral §§1–4, independently audited;
(2) and (4) unchanged.` Preserve the original filing as dated history and add
this resolution rather than retroactively making the incoming proof stronger.

### Attribution and constants

The incoming README explicitly cites an earlier *live*, uncommitted
`dyadic_module_spectrum_20261026` all-tree cap `(2b−3)n ceil(log₂ n)` by a
heavy-child argument. No file of that name is in the inspected incoming tree.
Credit that cited all-tree result; Session75 supplies a self-contained
balanced-module/entropy proof with a different constant and the sharp-leading
flat-target dual lower bound. Do not claim invention of the all-tree O(n log n)
phenomenon. The incoming triangle/Frobenius inequality itself gives leading
`q²/8`; it is consistent but weaker, not refuted.

Its displayed universal subleading coefficient `1/π` does not follow from
its cited b-dependent all-tree cap. If restating that old estimate for all
trees, use `q²/8−O_b(q^(3/2)log q)` or its explicit cited cap. This is a
proof-support limitation, not a demonstrated counterexample to that coefficient.

## 2. HYP-9422(1): correct all-order closure, distinct named-start guard

Exact source: incoming README §4 lines 159–182; THM-4814(D); HYP-9422 lines
28–35. The source proves only the lower bound and endpoint/parity law and
records finite witnesses at all regular classes n=7,9 plus P11 and 12 sampled
regular n11 tournaments. Its reflection expressly says it has no all-order
construction. Commit/message titles abbreviating this are not stronger proofs.

Session75 `session75-black-core/BLACK_CORE_GEODESICS.md` §3 supplies the
all-order construction for n=2m+1>=7, with exactly the same literal native
move and target `1→strong-core` or its converse. Its reversal-gauge proof is
converted explicitly into the labelled native complement word. The independent
`session75-black-audit/ADVERSARIAL_REPORT.md` passes both proof and certificates.
Thus HYP-9422(1)'s unprescribed-singleton theorem is resolved once that audited
packet is published. Its additional prescribed-singleton question is **not**
resolved by this construction: the first black path chooses its terminal.

The incoming phrase “two distinct path starts” means the singleton alternates
as native start/end. This is not the load-bearing guard in our proof. We keep
one fixed reversal-gauge endpoint v and choose the second *other* start to be
a different named vertex via two successors on one actual deleted-core cycle.
That guard both excludes later blue paths and rules out a terminal
`1→strong-middle→1`. The supplied n7 same-start route is a hostile to the
unguarded argument, not a counterexample to the source lower bound or witnesses.
The even balanced extension is additional work, not part of HYP-9422(1).

## 3. ROOT finite-norm sentence: still unlocated, no conflicting theorem found

Neither new THM-4814, HYP-9422, incoming README nor close-out supplies a Hilbert
weight/operator for the quoted zero-forcing finite-norm ROOT mode. The new
README §5 cites the paid-ROOT identity and live dyadic rows; it does not prove
a spatial Hilbert-domain statement. Targeted searches found no such later
committed theorem. The live dyadic-module source is not in this local Git tree.
Status remains **UNLOCATED TRANSCRIPT EVIDENCE**, not false.

Session75 ROOT `REPORT.md` retains the native primitive source law, integer
clock and first-hit ROOT boundary (ROOT self-edge removed). Its statement is:
`h(s)=epsilon(s)4·2^T/(3^d(3s+1)) z` on the ROOT basin, zero elsewhere;
finite norm for a supplied positive weight w iff `sum_ROOT w(s)B_s^(-2)<∞`.
Counting ell2 fails; finite-cutoff `mu_N=Q(A_N)` gives finite seminorm but native
pullback fails to descend to the quotient for every N>=2. This is **not** a
refutation of a different strictly positive Hilbert weight or operator.

The independent ROOT audit passes these scoped claims with two required
exposition guards: signed `F=0` in `e'=epsilon*b*e+F` differs from positive-b
`f=0`; a single common ROOT boundary value differs from source-labelled
transports with different desired terminal charges. Keep both guards.
The older weighted-renewal zero-forcing sentence concerns the fully odd map,
`b=3n/(3n+1)`, and physical clock `sum(1+k)`; it is not the missing native
mode definition. The older Collatz-amplitude zero forcing concerns Terras-map
Dirichlet characters, also a different operator and source space.

## 4. 433/79 and pairing: numerical recovery, not historical identification

Incoming README §5 lines 197–205 proves/cites only the four-child mean
`sum_p M_p=E+2S`, its averaged-path positive-part identity, and the integrality
obstruction to attaining the exact mean. It expressly says the quoted numbers
433/79 are not in a committed packet and are not reproduced there. This is
consistent with, not superseded as historical provenance by, our recovery.

Session75 dyadic revision 2 independently recovers at N=512, complete first
arrival (maximum depth42): 768 original labels; full source norm1244;
terminal energy1313; positive labelled-pair sum433; positive parent-bin sum293;
**whole-level positive-part debt79**. The independent audit reconstructs the
source law and all routes without importing producer code. The actual original
historical computation remains unlocated; the matched numbers are FINITE-EXACT
recovery, not proof of source identity.

No pairing contradiction: an averaged middle energy lies in the endpoint
interval, but some actual middle energy need not. Our original-native examples
show that failure while retaining original charges, terminal tags, clocks,
ports, empty bins and low-bit ownership. Mean equality/integrality does not
prove zero debt for a realized pairing. The alternative compatible-tree
selection theorem changes the partition tree; it does not close native R3.

## 5. One demonstrated incidental correction: “maximal order”

Incoming README §5 line190 calls `Z[(1+sqrt(−q))/2]` the maximal order for
arbitrary doubly regular q. This adjective is false in that scope: q=27 has a
finite-field Paley tournament, and with `η=(1+sqrt(−3))/2`, `η²−η+1=0`,

    Z[(1+sqrt(−27))/2] = Z[3η−1] = Z + 3Zη ⊊ Z[η].

The element η is integral but absent from the displayed order; the index is3.
Smallest correction: replace “maximal order” by “quadratic order”, and use that
explicit order instead of ambiguous `O_{−q}` when q is not squarefree. The
polynomial identity, integral action, rank statement and q=3 Eisenstein case
remain intact. THM-4814(E) already avoids the adjective and needs no theorem
retraction. `MINIMAL_ORDER_CORRECTION.patch` is proposed, not applied.

## 6. Evidence and audit labels

Incoming packet explicitly says **not independently audited**. Its normal/-O
receipt equality is not an independent audit. The four separate Session75 audit
reports are genuinely independent of their producer agents; the dyadic producer
file `AUDIT.md` and black producer `independent_audit.py` alone are not proof of
agent independence. This reconciliation relies on the separate audit reports.

Spectral auditor binding is to candidate hash `ac9ce92c...`; the present
`RESULT.md` hash `09872984...` differs only in its approved status line, with
exact status-only diff and promotion manifest retained. Black proof bytes match
the independent frozen binding directly. No new proof status is inferred from
a headline or file name. `SOURCE_PINS.json` records current exact hashes.

Stopping condition met: the incoming claims have been compared and narrow
integration/correction proposals supplied. The only missing evidence is the
uncommitted dyadic-module original and its finite-norm mode definition; no
publication or access action was attempted to obtain it.
