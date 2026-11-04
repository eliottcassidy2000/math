# Closing the two nontrivial-child circulant-nine branches

**Status: exact finite certificates with arbitrary child sizes; independent review pending.** This separate extension leaves the previously frozen nine-family packet unchanged. It covers only its two displayed prime circulant isomorphism classes, not every prime tournament of order nine.

## 1. Statement, physical labels and inherited family reduction

The quotient representatives, on physical child labels0 through8, are

    C0=(480,449,387,263,15,30,60,120,240),
    C1=(226,452,393,275,39,78,156,312,113).

For arbitrary nonempty source children A_i and target children B_i on the same fixed physical sets, assume q_h(A_i)=q_h(B_i) for each i and h. Profiles may differ freely between distinct children. Compare Ck[A_i] and Ck^op[B_i]. The ports satisfy x->y; a_i,c_i mean all physical vertices in child i point to x,y, respectively. The same physical masks are used on both sides.

The separate nine-family packet proves that every singleton-containing assignment is separated in every fixed ambient, and gives explicit common physical isomorphisms transporting these representatives to all twelve prime labeled cyclic step choices on Z/9Z. We now prove that every all-nonsingleton assignment is also separated. Consequently each of those twelve displayed prime circulant orientations is separated for arbitrary matching children in every fixed ambient. Together with the separately proved U9 and W9 branches, this closes the bounded family studied there. It is not a prime-nine census or an advance of the global highest-node frontier beyond nine.

## 2. The complete nontrivial-child domain modulo eight

We use the separately reviewed S188 complete normalized profile domain and Session189 prime-seven profile proof. A nontrivial child has one of the normalized profiles

    P=(1,2,0), E=(1,6,4), O=(1,6,6), C=(1,2,2) modulo8,

with later coordinates zero. Put b_i=1 for E,O and c_i=1 for O,C. The c_i bit is the parity of the actual child order. These four pairs(b_i,c_i) exhaust F2^2. Let L be the common odd product of the child Hamilton counts.

For each actual quotient-level port query z, write K0 for the run-word count with every quotient child once, K_i for child i twice, L_i for child i three times, and K_ij for children i and j twice. The ports themselves occur once. Occurrences of a repeated child are indistinguishable. All arcs and port orders are included. The complete converse contrast is

    Delta H/L = Delta K0 +2 sum_i Delta K_i
                  +4 sum_(i<j) Delta K_ij
                  +4 sum_i b_i Delta K_i
                  +2 sum_i c_i Delta L_i modulo8.        (1)

The fixed-point-free return-segment interchange makes each L_i even. Every selected baseline contrast in(1) is divisible by four. Accordingly the common odd L fixes the result modulo eight, and dividing(1) by four produces an exact affine equation in the eighteen bits(b,c), not a truncation to low port degree.

The packet retains all56 exact integers per side at each query:55 run channels and the rooted Hamilton-cycle count. It uses206 explicit actual queries. They are a finite separating test set; no completeness of this query selection is assumed. A finite subset of their equations is enough for each asserted conclusion.

## 3. Cycle-safe actual-query certificates

For odd total child order, a local scalar gap four modulo eight requires the complete higher-cover check. At each query, count Hamilton cycles of the bare quotient extension, rooted at x. Call the query safe if the source and target counts have the same parity. The inherited child-uniform cycle-collapse theorem makes this test valid after arbitrary child substitution.

When a safe query has Delta H=4 modulo8, the complete augmented-profile contrast is(4,0,...) modulo8. For even total order this follows from universal cover residues. For odd total order the inherited formulas give

    Delta q2=Delta q3=2 Delta H-4 Delta Z modulo8,

where Delta Z is the cycle-count parity difference; safety makes it zero. All h>=4 vanish by factorial divisibility. Thus every nonzero safe query survives every fixed ambient, multiplying by an odd Hamilton count. No unsafe scalar witness is used for the odd-total conclusion.

There are126 safe retained queries for C0 and119 for C1. Their affine systems have these exact certificates:

- C0: eighteen independent safe equations, with unique simultaneous zero(b,c)=(all1,all0), the profile type E^9.
- C1: an explicit collection of six safe affine equations whose left sides sum to zero and whose constants sum to one. No profile assignment can annihilate all of them.

For C1 the six actual queries are(0,7),(3,1),(3,384),(17,320),(32,129),(12,64). The source and target cycle parities agree at each one. Their six normalized contrast bits sum to one for every assignment(b,c).

`linear_certificates.json` lists the actual physical quotient masks and the precise equation indices, not just matrix ranks. The verifier checks the binary sums directly. Therefore every nontrivial C1 assignment is covered. Every nontrivial C0 assignment except E^9 is covered, including every odd-total assignment. The only remaining type E^9 has all children of even order at least four.

## 4. Two exact modulo-sixteen probes close E^9

For an actual even child of order at least four, the complete normalized profile modulo sixteen is

    (1,14+8e_i,4+8f_i,8,8g_i,0,...).

The independently reviewed all-order even-transfer theorem proves that its quotient-converse contrast is obtained from the geometric profile g_k=(-2)^(k-1) by only the q2 parity corrections in its Hamilton coordinate:

    Delta H/L = delta_g +8 sum_i e_i Delta K_i modulo16.

The q3 and q5 scalar corrections are even and vanish, and the complete augmented-profile contrast always satisfies

    Delta q2=2 Delta H,  Delta q_h=0 for h>=3 modulo16.  (2)

This uses childwise matching profiles and allows different even profiles and orders at different positions.

For the displayed C0, take the actual queries

    z0=(a,c)=(0,7),      z1=(a,c)=(8,264).

Exact finite word counts give geometric contrasts0 and8 modulo sixteen. In both queries, the q2 derivative contrast is odd precisely at child positions4 and7. Consequently

    Delta H(z0)=8(e_4+e_7),
    Delta H(z1)=8(1+e_4+e_7) modulo16.                 (3)

The common odd L fixes residues0 and8. Exactly one of(3) is eight modulo sixteen. By(2), its full contrast is(8,0,...), so it survives every fixed ambient. This is a finite two-query identity for C0, not an unproved all-order carousel extrapolation.

The geometric calculation is finite modulo sixteen: a run shape of total excess r contributes(-2)^r, so every shape with r>=4 vanishes. All220 weak-composition shapes of excess at most three are retained on each of the four query/orientation sides. Counts are exact integers, not only parity residues. An independent physical transitive-child interpolation recovers all880 run counts from actual named graphs of at most fourteen vertices by exact multivariate Stirling inversion. This checks both the geometric constant and the marked q2 derivatives without relying on the same multiset recurrence.

## 5. Ambient and original-module scope

Both ports inherit the original core's incidences to every exterior physical name. The augmented core is a common module. The modulo-eight safe gaps from Section3 and modulo-sixteen pure gap from Section4 therefore survive every fixed ambient.

For the all-nonsingleton carousel branch, its eighteen safe basis queries together with z0,z1 give a fixed twenty-query baseline battery independent of the child profiles. For C1, the listed six-query inconsistent safe subsystem is a fixed baseline battery. Singleton-containing sources use the separate nine-family battery. These finite batteries are source-relative separators, not all-pairs reconstruction of arbitrary tournaments.

For original-module XOR families supported inside one such prime core and selecting the whole core oddly, pair each baseline query with its core converse-complement J(a,c)=(511-c,511-a), preserving exterior masks. Complete augmented-profile duality supplies a negative baseline member; hereditary profile monotonicity and nonnegative run coefficients prevent target child increases from canceling it. This gives the same restricted root-family consequence as the reviewed prior packets. Exterior operations, arbitrary newly created modules and a complete order-nine quotient classification are not included.

## 6. Evidence, complete query data and reproduction

`analysis.json` and `run_values.txt` preserve both orientations' exact full run-channel and cycle counts for every selected query. Queries with more than one marked bit are actual queries, retaining all contributions, including cubic and quartic terms. Nothing is replaced by an endpoint product or a quadratic response truncation. The separately frozen family packet retains complete named Boolean banks for both orientations of C0 and C1, with every cubic and quartic coefficient and every zero.

`even_bank.json` preserves every exact shape count and both actual quotient extensions for the two modulo-sixteen probes. `physical_even_records.json` retains each interpolation graph's original quotient rows, actual augmented rows, child sets, physical port masks and exact Hamilton count. The interpolation neither relabels source and target independently nor treats a coarse query as if it had only two physical marked vertices. In addition, `physical_mod8_records.json` retains2,640 actual physical graphs of order at most thirteen. Exact transitive-child interpolation recovers every one of the55 run channels on both sides of all eighteen carousel basis queries and all six noncarousel contradiction queries. Thus every equation actually needed by the certificates has an independent literal physical-count reconstruction.

Run `python3 verify.py` and `python3 -O verify.py`. Only Python standard-library tools and the installed C++17 compiler are required. Every check uses explicit exceptions, not assertions removed by optimization. Retained exact physical counts have order at most fourteen and fit unsigned64-bit integers.

The domain statements are unbounded in child size because of the reviewed complete profile reductions. The quotient scope remains exactly two displayed nine-name circulants and their explicit common isomorphisms. Independent review of this extension is required before treating the family closure as established.
