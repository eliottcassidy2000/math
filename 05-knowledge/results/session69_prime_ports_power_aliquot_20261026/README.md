# Session69 — prime actions in field coordinates and exact literal-target composition

**Status:** independently audited elementary proofs relative to explicitly
**CITED** canonical decomposition/realization inputs, with separately scoped
**FINITE-EXACT** controls. No formal-assistant verification, unrestricted prime
census, or small-state/efficient general metagraph algorithm is claimed.

Two concrete advances go beyond matching the small numbers. The full Paley7
action now has a complete field-coordinate classifier and marked move law.
More generally, ordered multi-rooted complemented-path-cover responses form
an **exact recursively closed ordinary, uncoloured target-incidence carrier**
under canonical tournament substitution. The full response state is large;
the theorem is closure, not an economical classification of all prime inputs.

## The prime-action atlas obtained

[Native theorem/proof report](session69-tournament/REPORT.md), F1–F5, uses
`F=F2[alpha]/(alpha^3+alpha+1)` and the actual arc predicate on `F*`:

    x -> y  iff  Tr(y/x)=0,       x != y.

This is Paley7 in exponent coordinates. The full group is
`GammaL1(8)=C7 ⋊ C3`, not merely the cyclic subgroup. Every binary decoration
is uniquely

    w(x)=a+Tr(bx+c x^3),       a in F2, b,c in F8.

For `b!=0`, its normalized parameter is `r=c/b^3` modulo Frobenius. Together
with the two `b=0` cases this gives six types per parity, hence **12 ordinary
classes**, with root stabilizer strata `2,6,4` of orders `21,3,1`. The 12-class
census is recovered prior work; the full field classifier and marked law are
the continuation. The cyclic-only necklace quotient has 20 classes and loses
the extra action.

For `T_w=Paley7[A_(w_i)]`, `A0=1->C3`, `A1=A0^op`, the complete ordinary
cyclic-triple-realization fibre has those 12 types. The full labelled fibre
has `2^15` members; standard palette words and labelled realizations differ.
The group is `C3^7 ⋊ Stab(w)` with its actual slot action retained. These
fibre assertions use the explicitly cited reconstruction theorem.

Marking a quotient slot z gives `(a,B,C)=(a,bz,cz^3)` modulo simultaneous
Frobenius. A modular block flip is exactly

    (a,B,C) -> (a+1,B+1,C+1).

There are **48 block-pointed classes**, 24 unoriented flip-port orbits, and
96 actual-vertex-pointed types. This modular flip is not a literal
Hamilton-path-complement line. The function-module dimensions are **1+3+3**;
the exponent orbit `{1,2,4}` is not a representation-theoretic 1+2+4 split.
The trace rule also fails as a universal binary-field recipe: the report's
F32 inverse-pair hostile would orient a pair both ways.

## An exact, recursively closed literal-target interface

For an ordered spanning c-path cover C of a child S, retain

    rho(C) = [ S^op xor E(C) ;
               (first(P1),last(P1)),...,(first(Pc),last(Pc)) ].

The brackets mean the full multi-rooted ordinary target type. Endpoint
coincidences from singleton paths and the number of source-cover orbits
giving each response are retained. They cannot be replaced by scalar counts.

**P2:** a Hamilton word in `Q[S1,...,Sq]` splits uniquely into maximal block
runs. Its quotient walk orders the child-cover occurrences. Start with all
cross arcs reversed, then restore precisely the inter-run endpoint hinges.
Multiplying child response weights and dividing the invariant aggregate by
the actual child-type stabilizer `|K|` computes each ordinary target-class
incidence, not only H or triangle flux. The retained walk is the free carrier;
a bare response marginal cannot simply be assumed to carry a free `|K|`
action.

**P3:** use an ordered d-tuple of quotient walks and keep the ordered outer
endpoint pair of each global path. Hinges occur only inside each walk.
There is no `d!`. The same construction computes every parent cover response
from child responses, with singleton bases and the analogous ordered linear
node rule. Thus it is recursively closed on the canonical strong-module tree
given the prime/linear quotient relations and their actions. It does not
assume that source modules survive in the literal target.

Two exact state-loss controls make this closure substantive:

- In `C3[S4,1,1]`, two covers with the same local target and same two outer
  endpoints produce distinct literal target codes **81 and 83**. Two outer endpoints alone fail;
  retaining all ordered component marks repairs this witness.
- In `C3[S5,1,1]`, replacing positive response weights by one changes the
  correct **93** ordinary source ports to **88**. Equal responses genuinely
  aggregate multiple covers; their multiplicity is not disposable.

## What even the full flux observer loses

All 128 Paley28 decorations have
`H=5307704039306577787515` and 903 cyclic triples. The complete raw literal
partner flux histogram depends only on `k(7-k)`, where k is the number of
reversed blocks. In particular, two different weight-three Fano classes have
the same H, full automorphism order **6561**, and entire raw flux histogram,
yet are not isomorphic or converses. This is a same-group-order hostile, not
merely a collision of scalar H.

For the prime source Paley7 itself, **nine ordinary source ports reach seven
ordinary target classes**. Two prime target classes are reached twice, while
some other targets are nonprime. This does not refute a conjecture restricted
to sources whose every partner is prime, nor the distinct carousel conjecture.
The complete target codes and paths are in the native report and receipt.

[Independent tournament audit](session69-tournament-audit/AUDIT.md) accepts
F1–F5 and P2–P3: **4,464 active checks**, 31 literal-target substitution
controls, and seven all-d cases totaling **7,523 ordered covers / 3,044
response categories**. It imports no producer code. Its Paley28 endpoint
calculation uses determinant/cofactor rational-walk expansion, independently
of the producer's bounded-content dynamic program.

## What the small numbers actually say

The two appearances of six are real, but concern different predicates.

- In the **odd-root doubling forest**, six is first with its square a positive
  doubling of a later odd root: `36=4·9`, `9>6`. For `n=2^a u`, u odd,
  the exact condition is `a>=1` and `u>2^a`. If arbitrary roots are allowed,
  four is first (`16=2·8`). Nine is already in the power row of three, so
  “later odd root” does not mean “later new power component.”
- For the aliquot map `s(n)=sigma(n)-n`, six is the first **finite** target
  fibre with more than one source, and the first positive fixed point:
  `F(6)={6,25}`. The fibre `F(1)` is already infinite (all primes).
  The transcription correction is `s(9)=4`, not 3.

| target m | complete inverse fibre F(m) |
|---:|:---|
| 0 | {1} |
| 1 | all primes |
| 2 | empty |
| 3 | {4} |
| 4 | {9} |
| 5 | empty |
| 6 | {6,25} |
| 7 | {8} |
| 8 | {10,49} |
| 9 | {15} |
| 10 | {14} |

The collision family is exact: for odd twin primes `p,p+2`,
`s(2p)=s((p+2)^2)=p+3`. Only `p=3` gives a self-preimage in this family.
This supplies witnesses, not every source of every such target, and assumes
no infinitude of twin primes.

**A stronger global uniqueness theorem now explains the rare conjunction:**
for positive N and prime p,

    s(N)=N and s(p^2)=N  iff  (N,p)=(6,5).

Thus six is the **unique perfect target with a prime-square preimage**,
not just the first such target or a special case within the twin-prime family.
[Elementary proof](session69-root/RARE_SIX.md) uses the even-perfect-number
structure and a modulo-three obstruction; it assumes nothing about odd
perfect numbers or the infinitude of Mersenne primes.
[Separate independent addendum](session69-arithmetic-corollary-audit/ROOT_ADDENDUM.md)
accepts the proof and repaired [connection firewalls](session69-root/CONNECTIONS.md),
with a different sigma sieve through one million and 196 active checks.
The separate square-later predicate still has density one half.

## Exact arithmetic classifiers

[Power report](session69-power/REPORT.md) distinguishes all-power rows
`S_n={2^k n^r:k>=0,r>=1}` from repeated-square rows
`Q_n={2^k n^(2^j):k,j>=0}`. Their operation law, earliest seeds, exact row
intersections and indexed occurrence counts are proved. Their address counts
first disagree at 8; earliest seeds first disagree at 27. A new seed need not
be prime: 15 is already an odd composite new seed in both models.

The odd-root square-later set has density **1/2**, not the rarity associated
with perfect numbers. Retaining its positive difference gives a bijection to
primitive triples with distinguished binary term:

    (u-2^a)+2^a=u,     n=2^a u,     a>=1, u>2^a odd.

Six gives `1+2=3`. The first source with `C>rad(ABC)` is 72, giving
`1+8=9` with radical 6. This is a finite minimality check, not an ABC result.

[Aliquot report](session69-aliquot/REPORT.md) proves the complete-fibre source
bound `n<=(m-1)^2` for `m>1`. After removing prime-square possibilities the
stronger `4n<=(m-1)^2` holds. The fixed-core branch `s(ap)=s(a)p+sigma(a)`
requires `p` prime and `p` not dividing `a`. With
`d=gcd(s(a),sigma(a))`, `U=s(a)/d`, `V=sigma(a)/d`, its normalized primitive
triple obeys `C<=(U+V)rad(ABC)`. The core-dependent constant is unbounded;
this is not the uniform ABC conjecture. Normalizing `(n,s(n),sigma(n))`
also sends every perfect number to `1+1=2`; the removed gcd is indispensable
for reconstruction.

## The genuine seven bridge, and the false shortcut

The exponents `{1,2,4}` sum to seven and form a doubling orbit modulo seven.
Their residue difference predicate defines Paley7; exponent coordinates in
`F8*` carry residue doubling to Frobenius. This preserves actual tournament
arcs, unlike the passage from positive integers to residue labels, which loses
order, factorization and birth addresses.

Keep the three image questions separate:

- **PROVED examples:** 2 and 5 are aliquot-untouchable. Their completeness
  among all prime targets is **not proved here**; a distinct-prime Goldbach
  hypothesis is sufficient. The packet independently checks prime targets
  through one million and carefully distinguishes larger cited computations.
- **Inherited raw-H theorem:** 7 and 21 are omitted as directed Hamiltonian
  vertex-word counts. This is not the same as the orbit-normalized target
  `tc=H/|Aut|`, and does not imply an omitted `7·3^k` tower (`H=63` occurs).
- **Decisive hostiles:** `s(8)=7`, `s(18)=21`, while raw tournament `H=5`
  occurs. Paley7 itself has `H=189`, `|Aut|=21`, and `tc=9`.

No target-preserving map identifying these missing-fibre predicates is given.
The useful common method is certified image analysis with its sidecars kept.

## Arithmetic audit and reproduction

[Independent proof and scope audit](session69-arithmetic-audit/REPORT.md)
accepts these statements without a requested mathematical correction.
Its independent checker imports no producer code and uses a different
aliquot sieve. It performs **3,173,300 active checks**: power origins/address
counts through 20,000; the odd-root predicate through one million; complete
aliquot fibres for every target 2..1000; 27,827 guarded core/prime branches;
and independent order-four/Paley7 controls. Producer checks include 57,332
power checks and 8,785 guarded aliquot branches. Finite checks supplement,
rather than establish, the all-order elementary proofs.

From this directory:

    python3 -I -B verify_packet.py --replay
    python3 -I -B -O verify_packet.py --replay
    python3 -I -B verify_mutations.py

The envelope checks the exact file universe, bytes and hashes, then replays
isolated copies in normal and optimized modes against immutable receipts.
Negative controls exercise both mathematical gates and corrupted/missing/extra
packet bytes. The manifest is an integrity ledger, not an external trust anchor.
Native source reports, scripts, snapshots and receipts remain byte-for-byte
unchanged. Their research pin is `61d7cca5d8df80169f939268269239c2f67ea724`;
this is not the publication-parent pin.

## Provenance, concurrent work, and the remaining boundary

The research pin is `61d7cca5d8df80169f939268269239c2f67ea724`, not the
publication parent. Native producer and independent-audit files are frozen
byte-for-byte; [source snapshots](source-snapshots/SOURCE_SNAPSHOTS.json)
identify the exact inherited inputs.
No source snapshot is a blanket endorsement of its historical claims.

Concurrent reservation `fa62e88c5f36159c4c5c2c5c9d9f78600afe3ffc` contains
empty **RESERVED / UNPROVED** THM-4805/HYP-9415 stubs, not proved dependencies.
Concurrent `563022f89374a78270832c30fa8b165d340911da` supplies the separately
published [joint-port/reflection-fold and run-forest packet](../rooted_lines_frobenius_sheets_20261026/README.md).
The [incoming scope addendum](session69-tournament-audit/INCOMING_SCOPE.md)
accepts its R1 joint-port ordered-sum rule, R2 canonical-SCC blue folding,
and R4 full explicit run-data reconstruction of arbitrary-substitution colour
and legal re-marking. Incoming code and finite witnesses, including H=289/253,
were **not replayed here**. This is scoped proof review and explicit credit,
not a blanket computational endorsement.

Concurrent `52fde1af0df4543960465978c7788726304d4208` publishes the
[aliquot/power/1–2–4 packet](../aliquot_power_124_20261026/README.md), including
an independent proof of the **same unique perfect prime-square target** and
the compatible distinctions between first occurrences at 3, 4 and 6. We retain
our own frozen proof/audit and explicitly credit the overlapping result; no
priority claim is made. Its [scope repairs](../aliquot_power_124_20261026/CORRECTIONS.md)
restore S57's distinct-Goldbach antecedent and THM-4748's minimum-over-a-fibre
quantifier. Those repairs are preserved; this packet never imports the
stronger discarded assertions. Its computations were not replayed here.

**Closed here:** sufficient full multi-rooted ordinary, uncoloured response
composition, recursively from canonical children with their quotient data.
**Still unestablished here:** an unrestricted structural prime atlas,
minimal/economical state compression, efficient full metagraph generation,
an economical universal colour-aware, complement/re-marking-compatible
aggregate algebra,
prime-untouchable completeness, uniform ABC, or any LRC(14) consequence.
The separate coloured ordered-join rule and lossless full-run reconstruction
must not be erased by an overly broad claim that all colour handling is open.
P3 itself proves ordinary uncoloured aggregation; no target-only colour-
sufficiency probe or small universal closure theorem is claimed.
