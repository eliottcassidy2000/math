# Session69 — distinct first sixes, exact occurrences, and guarded bridges

**Status:** independently audited elementary proofs and separately scoped
**FINITE-EXACT** checks. This checkpoint recovers the arithmetic interfaces
needed by the continuing prime-action and literal-metagraph investigation;
it does not claim an unrestricted prime atlas or full metagraph closure.

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

## Independent audit and reproduction

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

Concurrent mainline reservation `fa62e88c5f36159c4c5c2c5c9d9f78600afe3ffc`
contains empty **RESERVED / UNPROVED** THM-4805/HYP-9415 stubs for overlapping
questions. It is preserved and is not a proved dependency of this packet.
The inheritance/audit reports credit prior work by full slug, with source
snapshots and pins; no global novelty or priority claim is made.

**Still unestablished:** unrestricted prime-tournament classification, an economical closed
general metagraph grammar, prime-untouchable completeness, uniform ABC, or any
LRC(14) consequence. The continuing tournament lane must retain the full
prime permutation action and literal path/port targets rather than infer them
from arithmetic seeds, scalar counts, or shared small labels.
