# A common-future proof keeps its clock barrier; a charged loop can certify a period

**Status:** PROVED elementary certificate rules, with FINITE-EXACT controls.
The clock algebra is the classical bicyclic monoid, not a new algebra.
The useful addition here is an executable, source-checked interface from
asynchronous join banks to periodic-tail and ROOT certificates. Existence
of sufficient receipts for every positive integer remains OPEN.

## 1. Inheritance and the shift in target

The [previous finite-bank construction](collatz_finite_bank_certificate_20261007.md)
gave a legal merge cylinder for every bounded integer translation. Its
canonical hostile is the positive pair(1,2), which never merges at equal
Terras time while both sources reach ROOT. The
[Vitali/selector certificate note](vitali_selector_certificate_20261005.md),
sections4-5, already requires finite edge receipts for a component label.
The corrected near miss is treating a component minimum or a common-future
class name as a ROOT proof. The underused coordinate is the pair of actual
clocks, including the time before a join becomes valid.

Board: **source / time pair / overlap seam / cycle charge / eventual period /
positive domain**. The anchor is a sound rule for reusing asynchronous
Collatz joins; the niche is exact monoid composition; the wildcard is using
a nonzero cycle charge as positive evidence instead of discarding every
dependency cycle as circular.

No novelty claim is made for the underlying semigroup or groupoid facts.
[Hines, arXiv2206.07412, section3.1](https://arxiv.org/abs/2206.07412)
gives the bicyclic normal form as partial translations of natural numbers.
[Armstrong--Brix--Carlsen--Eilers, section2.2](https://arxiv.org/abs/2105.00479)
uses common-future triples and their integer cocycle in the
Deaconu--Renault groupoid. The finite rules below need only a total map;
no topological groupoid theorem is required for their proofs.

## 2. Composition keeps the point where an equality becomes valid

Let F be any total map on a set. A checked receipt

    R=(x,y;a,b) means F^a(x)=F^b(y), with a,b>=0.

Receipts must include independently replayable endpoint evidence. If
R=(x,y;a,b) and S=(y,z;c,d), their literal composition is

    R*S=(x,z; a+max(c-b,0), d+max(b-c,0)).             (1)

Proof: advance both equalities until their y clocks equal max(b,c).
The source seam y is essential. Reversing a receipt swaps both its sources
and its clocks. Its integer charge is a-b; charges add under(1).

Associate to(a,b) the partial translation t -> t+a-b on t>=b. Composition
of these partial maps is exactly(1), proving associativity. This is the
bicyclic monoid. Its idempotents are(a,a), the identity only on the tail
t>=a. In particular,

    (a,b)*(b,a)=(a,a),

which is generally not(0,0). The integer-difference quotient is useful for
cycle arithmetic but discards the domain barrier. The concrete Collatz
receipt(5,4;3,3) has zero charge and is valid, yet5!=4; it cannot be replaced
by the empty-clock assertion(5,4;0,0).

This normalization is the smallest alignment forced by the supplied two
receipts. It is not a claim that the underlying trajectories have no earlier
join: that would require additional orbit evidence.

## 3. Closed walks with nonzero charge force a periodic tail

Compose a closed walk of checked receipts at x to obtain(x,x;a,b).
If a!=b, put m=min(a,b) and q=abs(a-b). Then

    z=F^m(x),  F^q(z)=z.                              (2)

Thus a nonzero time discrepancy is a finite periodicity certificate. A
zero-charge loop provides no such consequence: for the nonperiodic map
F(n)=n+1 on the natural numbers, (n,n;4,4) is always true.

For finitely many nonzero loops at the same source, let their entries be
m_i and positive charges q_i. At M=max(m_i), z=F^M(x) is fixed by each
F^q_i. Euclidean subtraction preserves this property: if q>=r and
F^q(z)=F^r(z)=z, then F^(q-r)(z)=z. Therefore

    F^g(z)=z,  g=gcd(q_1,...,q_r).                    (3)

Both the entry M and the loop witnesses remain part of the certificate.
A gcd calculated from unverified proposed joins is not evidence.

For a finite connected join graph, choose a spanning tree and transport
each non-tree edge to a loop at the root. The gcd of these fundamental
cycle charges is the gcd of all closed-walk charges, because the integer
charge is additive and every graph cycle decomposes in that basis.
This provides a finite certificate target. It does not assert that every
component contains a nonzero-charge cycle.

## 4. Positive Terras periods turn some circular dependencies into ROOT proofs

Use T(n)=n/2 when n is even and(3n+1)/2 when n is odd. The four parity
words of length2 give the only rational candidates for T^2(n)=n:

| Parity word | Affine second iterate | Fixed candidate |
|---|---|---:|
| 00 | n/4 | 0 |
| 01 | (3n+2)/4 | 2 |
| 10 | (3n+1)/4 | 1 |
| 11 | (9n+5)/4 | -1 |

All candidates obey their indicated parity words. Consequently the only
positive solutions are1,2. If a positive-source join bank supplies loops
whose nonzero charges have gcd2, (3) gives a route to ROOT by time M+1.
The code extracts and independently replays its first-hit ROOT path.

This is a sufficient inference rule, not a universal construction of its
premises. The demonstration rows intentionally construct known finite
loops of charges4 and6; they test the inference, not a new discovery of
convergence for those inputs. Gcd1 is impossible for a positive Terras
source because T has no positive fixed point. A zero-charge loop must
remain unresolved. Signed or zero sources cannot use the positive-period
classification: -1 and0 are exact hostile fixed points.

More generally, a finite certified classification of the positive points
fixed by T^g can replace the g=2 classification. The complete word census
through g=12 here is a small control, not an all-period statement. For a
parity word of length g with r odd steps, T^g(n)=(3^r n+C)/2^g, so a fixed
point must be C/(2^g-3^r). Checking all2^g words and replaying their parity
guards is finite and exhaustive.

## 5. The 223/233/322 triangle has a real join but zero cycle charge

Direct integer replay gives

    T^7(223)=T^15(233)=T^25(322)=425.

The receipt(233,322;0,10), composed with(223,233;7,15), gives exactly
(223,322;7,25). Closing these three particular receipts gives charge zero.
It certifies a shared future, not ROOT by itself. The
[prime/partition join-fibre note](collatz_prime_partition_223_20261007.md)
retains their full ordered words and arithmetic progressions.

The exact finite source rows in this checker reach1 after Terras times
0,1,5,70,46,54,64 for1,2,3,27,223,233,322 respectively. They are seven
literal certificates, not a density estimate or a universal bound.

## 6. How the new paper mechanisms connect

| Input idea | Actual map into this interface | Preserved target | Needed information |
|---|---|---|---|
| Shortest common superstring | Store receipt words on a shared tape; keep occurrence intervals | Exact word retrieval | Source, interval, both clocks, seam compatibility |
| Kervaire equations for groups | Compare adjoining a formal relation with checking a concrete join | Existing relations only when the representation is faithful | Actual integer realization; clock partial maps are not a group action |
| Finite algebraic colouring obstruction | Extract and verify a finite consequence object | Sound finite rejection or proof | Search termination premise; no uniform size bound follows |
| Finite-bank Haar positivity | A guarded cylinder supplies a valid edge when its source matches | Conditional positive mass and literal equality | Prescribed source bits; positivity alone does not supply a charged loop |

The practical next search is for a second independently verified path whose
charge differs from a first one by2, or for a small finite set of loops with
gcd2. A successful search provides a ROOT certificate immediately. A failure
at a chosen depth is only failure of that finite search. The root-cycle
phase obstruction is removed by the clock-pair representation; universal
receipt coverage is still the principal open obligation.

## 7. Reproduction and exact controls

    python -B 04-computation/experiments/collatz_clock_holonomy_20261007.py
    python -B -O 04-computation/experiments/collatz_clock_holonomy_20261007.py

The declared universe is all15,625 triples of clock pairs in[0,4]^2;
625 pair compositions independently realized as partial translations at
17 arguments; all288 labelled maps on sets of sizes1..4, with209,982 pairs
of actual closed receipts using clocks0..6; all8190 Terras parity words
of lengths1..12; the seven positive source rows; and ten malformed-type,
source-seam, zero-charge, signed-domain, or unsupported-period hostiles.
Normal and optimized output must agree. The mathematical all-map and
all-clock quantifiers follow from(1)-(3), not the finite census.
