# Reconciliation with Session315: three different dyadic questions

**Status:** incoming proof mechanisms are **INHERITED PROVED** in their
stated scope; the low-target inventory is **INHERITED FINITE-EXACT** with a
**CITED** all-order cutoff; this note adds a typed reconciliation, not a new
proof of those results. Chan–Pak remains **OPEN**.

The incoming surface is pinned at
`ae6c9620d8718e3ce5dd710ef6ecb1e00b93436a` (Session315 synthesis).
[Raw-source fingerprints](incoming-sources.json) give SHA256 and byte lengths
for every incoming document used below. Source proofs were read; their
large or full paired-site computations were not rerun for this note. Our
independent low-H census and diagram checks retain their separately stated
universes. No chronology-based priority claim is made.

## The incoming recurrence is exact and retains its cycle variables

[Session315, pointed-cycle recurrence](../session315_dyadic_pairing_20261113/SPECTRUM_AND_SIEVE.md)
uses odd permutation-cycle types and the full inventory of cycle variables:

```text
n P_n^(a,b)(x) = a sum_(odd d<=n) x_d b^((d-1)/2)
                    P_(n-d)^(a,b)((b^gcd(d,e) x_e)_e).
```

At unit cycle variables, `(a,b)=(1,2)` gives A000568; `(2^epsilon,4)` at
site order m gives the self-converse count at vertex order `2m+epsilon`.
The source already distinguishes this exact recursion from a bounded-state
scalar recurrence: each future insertion changes the gcd substitutions.
The construction and evaluations are inherited here, not claimed as an
additional contribution of our packet.

Its SCC Dirichlet-weighted identities are the same ordered-SCC and
square-central mechanisms used in [pairings P1–P2](pairings/PROOF.md),
and already credited there to Session24 and the November2 spectrum packet.
They isolate the unknown strong inventories. They do not enumerate those
inventories merely by writing the generating functions. Nor does the
ordinary Burnside scalar recurrence determine individual H fibres or
physical complement-line partners.

Our P1 fixes a **finite H cutoff B and precision k first**, then proves
periodicity modulo2^k of every ordinary and self-converse coefficient with
H<=B as the vertex order n varies, from n>=B. For B>=3 the period bound is

```text
2^(k-1+ceil(log2(floor(log3 B)+1))).
```

P2 computes each fixed-H least eventual parity period from the reduced
pole at x=1 over F2. These results use multiplicities of actual strong
cores, including cancellations. They are not periodicity claims for all
of F_n, for A000568 itself, or for integer multiplicities before reduction.
Evaluating F_n at1 after dropping H>B terms need not give A000568(n).

## All odd residue classes coexist with exact self-converse holes

[Session315, dyadic spectrum §4](../session315_dyadic_pairing_20261113/DYADIC_SPECTRUM.md)
constructs a strong self-converse tournament in every odd residue class
modulo every power of two. Its quantifiers are

```text
for every k and odd c modulo2^k, there exists a strong self-converse T
with H(T)=c modulo2^k.
```

The witness, its order, and its exact H may depend on both k and c.
The four rooted charts plus mirrored wrapper are the incoming source's
construction, not a new construction in this packet. Its fixed wrapper
bank still has only `O((1+log X)^2)` distinct exact counts up to X. The
positive wrap correction vanishes modulo the chosen modulus without
vanishing as an integer.

[Session315, exact self-converse gaps](../session315_dyadic_pairing_20261113/SPECTRUM_AND_SIEVE.md)
gives the missing positive odd values **strictly below75**:

```text
7,19,21,23,35,39,43,47,53,55,59,63.
```

This is a globally quantified nonexistence claim over all tournament
orders, within a finite target range. Its finite paired-site inventory
through order8 becomes complete by Busch's strong minimum of75 from
order9, followed by the square-central SCC decomposition. The result
concerns the self-converse sector; it does not exclude ordinary19,23,63
and makes no claim about the boundary value75. This note inherits that
larger inventory rather than promoting our through-order6 census into it.

There is therefore no conflict at19. For each k some strong self-converse
witness has count `19+2^k q` for an integer q, yet no self-converse witness
has count exactly19. The witnesses cannot be replaced by one finite
fixed tournament satisfying all precisions. Our P4 instead concerns the
**ordinary exact-H19 fibre**: for every n>=6 it has `2(n-5)` classes, all
non-self-converse, and compatible free binary rank `1+v2(n-5)`. That rank
is unbounded and not eventually periodic, while every fixed cap is.

Likewise P3's `Delta_21(n)` is the minimum number of fixed classes of an
H-preserving, converse-commuting involution on the order-n carrier H<=21.
Its eventual pattern `2,4,3,5` is a count of residual classes at fixed n,
not a list of realized H values or a modular realization statement.
The exact self-converse gap19 contributes no residual class because its
entire self-converse fibre is empty, not because its ordinary fibre is empty.

## Full half-cuts repair H; actual partner wiring is a further sidecar

[Session315, half-cut fold](../session315_dyadic_pairing_20261113/HALF_CUTS.md)
retains endpoint counts `E_S(i)`, actual junction arcs, and both transversal
and nontransversal cuts. Mirror-nonfixed cuts contribute twice their cut
value. Fixed cuts contribute quadratic expressions in even order and
squares in odd order; only the diagonal in the actual path-choice product
counts fixed paths. Thus the source supplies an exact per-chart Hamilton
count together with its fixed-path count. It explicitly does not turn the
full black incidence into a two-state transfer.

This agrees with our [exact-wiring obstruction](PAIRING_WIRING.md).
At order5 an abstract optimal H-preserving involution commuting with
converse fixes4 classes. Requiring preservation of the literal merged
complement-line dart matrix raises the minimum to8: the weighted merged
graph is vertex-rigid. The H9 exchange loses a specific transitive-node
incidence, and the H15 exchange loses degrees5 versus3. A fibre cardinality,
its parity, or its total H is therefore insufficient for a diagram action.
This finite obstruction does not assert that the complete half-cut table
is a minimal carrier or that it can never support a richer reconstruction.
It identifies what a proposed lift must actually reproduce: the physical
partner map and its dart multiplicities.

The incoming [paired-site construction](../session315_dyadic_pairing_20261113/PAIRED_SITES.md)
already retains physical blue partners and distinguishes their marked
centralizer from the full automorphism group. That is stronger than a
blue count alone, and still narrower than full black wiring. We preserve
all three levels of information in this reconciliation.

## Connection contract and inherited method

The method card used is [Controlled forgetting and unlabeled quotients
require a sidecar](../../../00-navigation/META-PATTERNS.md#controlled-forgetting-and-unlabeled-quotients-require-a-sidecar).
Its instruction to inspect fibres and retain off-diagonal seams applies
literally: nontransversal cuts repair the H observer, while endpoint/dart
incidence controls whether a class pairing respects the actual diagram.
Its quantifier-order boundary separates a varying modular witness from
one exact target witness. No new meta-pattern is proposed.

**Source:** Session315's pointed cycle inventory, marked half-cuts, and
rooted dyadic constructions. **Target:** our fixed-H coefficient clocks,
optimal abstract pairing defects, and actual diagram-compatible pairings.
**Maps:** evaluate the cycle inventory for totals; project marked path
counts to `(n,H)`; reduce multiplicities modulo2^k; or demand the physical
partner map for the last target. **Preserved:** precisely the named
coefficient or incidence under its own map. **Lost:** cycle substitutions,
integer carries, witness order, nonfixed cut contributions, or dart partners
when their corresponding sidecars are removed. **Decisive hostiles:** the
exact19 gap despite residue coverage; the H9 multiplicity cancellation;
and the order5 defect increase from4 to8. The incoming work improves exact
interfaces without supplying Chan–Pak completeness or a compressed full
black-wiring recurrence.
