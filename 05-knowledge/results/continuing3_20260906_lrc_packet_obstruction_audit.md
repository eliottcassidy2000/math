# Independent audit of both actual packet rows and their bounded scale class

**Status: independent full proof/source audit PASS.** No mathematical repair
is requested in the
[producer report](continuing3_20260906_lrc_packet_obstruction.md).
This review covers both scales 412164 and 412183, the entire declared bounded
coprime scale class, actual equality entry, all inherited profiles, and the
first-denominator and compatibility claims. It supersedes the earlier
bounded arithmetic review only in audit coverage, leaving its chronology intact.

## Actual entry and the full scale quantifier

The cofactor construction is correct. The six denominators are pairwise
coprime and their product is lcm(1,...,28)=80313433200. Their cofactors have
collective gcd one. The normalized shapes are exactly the displayed V and
U, and every gcd(V_i,U_j) is one. For

    412164<=t<=25880486141010, gcd(t,305613336829)=1,

it follows that gcd(tV_i,U_j)=1 at every cross pair. The physical row has
13 distinct positive coordinates and is primitive because U contains one.
The upper endpoint is exactly floor((Q^2-sum U)/sum V), retaining the actual
physical box; this is a bounded class rather than an asymptotic construction.
The gcd restriction is part of the class and is not silently dropped at
either endpoint.

The complete inert-atlas graphs on V and U are the displayed connected
graphs. A separate multiplicative construction of all eligible sums verifies
every edge and nonedge. The physical cross ratios already exceed 355 at
the class's lower endpoint and increase with t, so there is no cross edge
anywhere in the class. Internal edge coefficients have height at most 355,
and a connected weighted edge graph spans its full component relation
kernel over Q. For the two displayed rows, independent rational elimination
of the literal 14-by-13 weighted edge matrices confirms rank 11 directly.

The stronger equality assertion does not follow from that rank alone.
The report correctly excludes every bounded crossing relation:

* For two V labels and one U label, any nonzero coefficient on the U label
  is divisible by t gcd(v,w), by cross-coprimality. This is strictly larger
  than Q throughout the class. There are 15*7=105 such supports.
* For two U labels u,w and one V label tv, put D=gcd(u,w). The coefficient
  on tv must be divisible by D. After division, a nonzero target has magnitude
  at least tv, whereas the remaining signed sum has magnitude at most
  Q(u+w)/D. The lower endpoint already gives a strict inequality in the
  opposite direction. There are 21*6=126 such supports.

Both statements allow a redundant coefficient to vanish, so support-two
crossings are included. The exact extrema are dV=1738800 and sU=236;
at t=412164 the margins are

    t*dV=716670763200>Q,
    t*min V/(Q*sU)=4786336882800/4786326552917>1.

Both inequalities increase with t. The stated lower endpoint is the first
integer satisfying the stronger of these two uniform sufficient gates; no
claim is made that smaller scales cannot have equality by another mechanism.

All height-at-most-Q, support-at-most-three relations are therefore internal.
Conversely the actual bounded atlas edges span both internal kernels.
These two inclusions prove W=V_dec under the stated THM-3818 actual decoder
definition, not an abstract rank surrogate. The proof covers the entire
bounded class analytically; two numerical representatives are not being
used to infer that quantifier.

## Hereditary profiles and denominator loading

Every body of at least seven labels either is exactly U, whose gcd is one,
or meets both components and contains a coprime cross pair. Thus its gcd is
one for every scale in the class. Each complementary gcd is also one.
The six all-one signatures are explicitly present in the pinned full
joint-shadow JSON, not merely allowed by the scalar gcd ceilings. The referee
independently reconstructs all 4095 body/complement profiles in each of the
two physical rows and tests their actual membership in that JSON.

For every denominator 2<=d<=28, d divides P and at least one of the six
disjoint prime groups q_i is coprime to d. Such a d cannot meet all six
groups, since that would require at least six distinct prime factors.
Consequently d divides P/q_i for some i. This creates a zero-distance
coordinate at every numerator k/d, regardless of whether k is a unit. It
blocks the primitive V row and every scaled physical row in the class.

The complete primitive V packet at denominator 29 is the stated 14-element
set, so its first safe rational denominator is exactly 29. This does not
give a small-packet obstruction for U: its phase 1/2 remains safe. The
producer correctly limits the failed forcing claim to each-component and
full-row small-grid requirements, keeping adaptive component selection open.

## The two located scales and exact physical first denominators

The denominator-29 U packet and all 29 multiplier overlap counts are correct.
Physical compatibility is precisely k in the U packet and tk mod 29 in
the V packet. Its value depends on the located multiplier, not only the
two packet cardinalities or the component/profile data.

At t=412164, the multiplier is 16 and the full physical packet is
{5,12,17,24}. All denominators 2 through 28 are blocked, so the first
physical denominator is 29. Phase 5/29 has exact clearance 3/29. The
additional phase 166010/332021 has clearance 7/29, independently checked.

At t=412183, the multiplier is 6 and the denominator-29 packet is empty.
The two primitive packets themselves remain nonempty. The complete overlap
profile has only 6 and 23 as nonzero multipliers with zero overlap.
Denominator 30 has a divisible physical coordinate, and the denominator-31
packet is exactly {3,28}. Thus the physical first denominator is 31, and
phase 3/31 has exact clearance 3/31.

Both rows are safe and already satisfy the inherited larger-unit closure
criterion min V<=floor(3Q/28). The result is an actual-entry stopping object
for packet and compatibility forcing, not a new LRC closure or a new CRT
join theorem. A zero overlap rejects only that selected clock.

## Independent arithmetic and frozen artifacts

The [referee source](../../04-computation/continuing3_20260906_lrc_packet_obstruction_audit.py)
imports no producer module or external algebra package. It generates the
inert atlas multiplicatively from primes, uses union-find for components,
performs exact rational elimination on weighted edge rows, and constructs
safe packets as complements of unions of unsafe numerator masks. These
provide distinct arithmetic paths from the producer's factor predicate,
component propagation, and direct packet conjunction. It separately checks
all 231 mixed supports per displayed row, both complete body universes,
the exact scale thresholds, and every relevant denominator numerator.

```text
python -B 04-computation/continuing3_20260906_lrc_packet_obstruction_audit.py --repo .
python -B -O 04-computation/continuing3_20260906_lrc_packet_obstruction_audit.py --repo .
```

Both runs pass **8842 always-active gates**, with byte-identical LF
[output](continuing3_20260906_lrc_packet_obstruction_audit.out).
Canonical-LF normalization applies only to the explicitly inherited JSON
input pin, accounting for checkout line-ending conversion. The primary
source and output are separately pinned by their literal bytes.

```text
audit source  2e5b3ca4327e97180e9295e0825ee934ec2ba39d8315228ea439057e0b63b534
audit output  2f0003f0e3259fb7db17db35942a8a655029ebea7c56a2f49fdee762a5dbaa91
primary source 66eacc9647b765c392fb080ed41a031fd50a307b6d58b49629d2cdf342832080
primary output b3e340d2c5731d9e9b0bf48b78ff8f8b84db887cf1b5c375e129a364a63ee0f4
inherited JSON 935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
```

No primary artifact, repository file, or Git state was changed.
