# Synchronized graft windows, exact gap filters, and the primorial boundary

**Status: PROVED elementary consequences of inherited physical graft identities; FINITE-EXACT controls are separate. No prime-infinitude, Polignac, Brocard, or Chan–Pak conjecture is proved.**

## 0. Inheritance and scope

The live board is: actual source graphs; ordered path-cover products; fixed versus growing prime support; endpoint/count distinction; fixed-gap versus prescribed-target quantifiers. The anchor is a concrete prime-gap candidate family. The niche is an exact unit-intercept obstruction. The wildcard is the endpoint difference compiler in [GAP_UNIVERSALITY.md](GAP_UNIVERSALITY.md).

Recovered dependencies and hostile controls:

- Session229 `05-knowledge/results/session229_combined_constructions_20261109/inverse/PROOF.md`, §§1–3, supplies the proper-graft scalar row, parity, and physical predecessor restriction. These are inherited, not new.
- Session227 `05-knowledge/results/session227_chan_pak_constructions_20261108/multiroot/PROOF.md`, §3, already supplies scalar diagonalization in path-cover coordinates and their ordinal-product law. The two-coordinate diagonal form below is its singleton-root specialization, with a direct proof retained for reproducibility.
- Session237 `symmetry_grafts/PROOF.md` in `session237_chanpak_synthesis_20261109`, and THM-5092-four-test-tail-content-and-invisible-reroot-escape.md, warn that even four coordinates lose reroot information.
- Session241 `session241_prime_relations_20261110/constructions/PROOF.md` supplies actual Fermat-product graphs, but explicitly does not preserve addition or primality.
- Session242 `session242_congruent_descent_20261110/exception/PROOF.md` proves the content-seven restriction and exhibits the nonphysical normalized pair (7,116). Its descent packet supplies actual small prime roots, not a general square-root operation.
- Incoming `kps_s20_dyadic_cones_pillai_20261109/PREREGISTRATION.md` studies dyadic-port chaining and bounded odd-response closures. Here no port-chaining persistence is claimed: we attach a transitive gadget once, or flatten fixed-root words into one ordered gadget.
- The full-count exclusions 7 and 21 are inherited from THM-1370-h-spectrum-omits-7-21-all-n.md. THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md gives the stronger one-vertex +2 obstruction.

The contribution is the explicit four-source synchronization, its exact constant-gap/twin criterion, a simplified complete small-prime wheel, and the unit-intercept/primorial application. No literature-priority claim is made.

## 1. Literal construction and inherited diagonal coordinates

A row's bit j is 1 exactly when its vertex beats j. Let H count directed Hamiltonian paths, with H(empty)=1. For a nonempty gadget B set

    H = H(B),
    Q = sum_{X subset V(B)} H(B[X]) H(B[V(B) minus X]).

For rooted source (A,r), let D=H(A), E count paths ending at r, d=H(A-r), and U count ordered two-path covers whose distinguished nonempty path starts at r; the other path may be empty. Put

    a = D + 2d - E - U = d-f,       b = U-d,
    f = U-D-d+E.

For a nonsink root, the inherited physical constraints give a odd, b>=1, and f>=0. A proper graft retains both tournaments and directs

    A-r  --->  B  --->  r,

while retaining all old arcs between r and A-r. Its count is

    H(G(A,r;B)) = Ha+Qb.                                      (1)

Keeping the old root, its full state is exactly

    (D',E',U',d') = (Ha+Qb, Hd, Qb+Hd, Hd),
    (a',b') = (Ha,Qb).                                        (2)

Equation (1) is the inherited proper-graft row. Deleting the root leaves the ordinal sum (A-r)=>B, proving d'=Hd. A path ending at the root must visit all of A-r and then all of B, proving E'=Hd. For U', independently choose an old distinguished-root cover and an ordered cover of B, concatenating the respective paths. The only invalid combinations have old distinguished path (r) alone and its B-part nonempty, because r loses to B. There are d(Q-H) such combinations. Therefore U'=UQ-d(Q-H)=Qb+Hd. This proves (2).

For ordered join B1=>B2, intersect each of the two paths with each block and concatenate in block order. This bijection proves

    (H,Q)(B1=>B2) = (H1 H2,Q1 Q2).

Moreover two successive grafts at the same root are literally the one graft whose gadget is B1=>B2, not merely scalar-equivalent. Hence for a fixed-root word containing at least one nonempty gadget

    H(final) = a product_i H_i + b product_i Q_i.              (3)

This is explicitly an inherited diagonal-product mechanism, not a new operation. Its useful new application here is retaining the two coefficients while comparing physical source graphs.

## 2. Four strong sources with a common slope

Use the two actual tournaments

    A4 = [8,1,3,6],
    A5 = [20,1,2,7,14].

Their four marked states are:

| Source | Root | (D,E,U,d) | (a,b) |
|---|---:|---|---|
| A4 | 3 | (5,1,7,1) | (-1,6) |
| A4 | 2 | (5,1,9,3) | (1,6) |
| A5 | 1 | (9,3,9,3) | (3,6) |
| A5 | 2 | (9,3,11,5) | (5,6) |

These small finite facts are certified both by literal permutation enumeration of all induced path sets and by a separate subset endpoint dynamic program. Strongness is visible from Hamiltonian cycles: A4 has 0->3->2->1->0; A5 has 0->4->3->2->1->0.

**Four-window theorem.** For every nonempty actual gadget B, all four values

    6Q-H, 6Q+H, 6Q+3H, 6Q+5H                              (4)

are counts of actual strong tournaments, on respectively |B|+4, |B|+4, |B|+5, |B|+5 vertices.

Proof: apply (1) to the four literal sources. To check strongness without assuming B strong, its vertices all reach r; r reaches all old vertices because the source is strong; each old nonroot vertex reaches every gadget vertex directly. The root reaches a nonroot vertex and thus the gadget. All ordered pairs communicate. The source-order distinction is essential: no same-order strong claim is made for all four outputs.

Taking B=TT_m, m>=1, gives H=1,Q=2^m and four consecutive odd counts

    6*2^m-1, 6*2^m+1, 6*2^m+3, 6*2^m+5.                    (5)

The middle third value is divisible by 3 and exceeds 3. This is never a prime quadruplet. The first two are twin-prime candidates; if both are prime they are necessarily consecutive. The second and fourth, if prime, are consecutive primes of gap 4, because the intervening odd integer is the third value. First and fourth can be prime without being consecutive: m=1 gives 11,17 with 13 between. Thus a gap-6 Polignac instance needs the additional compositeness check of the second value.

For arbitrary fixed B, replace it by B=>TT_t, t>=0 (TT_0 empty). This preserves H and sends Q to 2^t Q. Since Q>0, (4) supplies infinitely many distinct strong count pairs at gaps 2H,4H,6H. These are differences of attainable counts, not necessarily adjacent attainable counts or prime gaps. The global all-even difference theorem is already known; see the inheritance in GAP_UNIVERSALITY.md.

## 3. Exact constant-gap criterion in this architecture

Fix two physical nonsink sources with coordinates (a1,b1),(a2,b2), and one actual gadget (H,Q). Pad the gadget identically by TT_t. The ordered difference is

    Delta(t)=H(a2-a1)+2^t Q(b2-b1).                           (6)

It equals a fixed integer g for infinitely many nonnegative t iff

    b1=b2 and H(a2-a1)=g.                                    (7)

Indeed if b2-b1 is nonzero the strictly monotone exponential term takes each value at most once. If it is zero, (6) is constant. Since the a_i are odd, a positive prescribed gap 2k requires H dividing k and a2-a1=2k/H. In particular a twin gap requires

    H=1, b1=b2, a2-a1=2.                                    (8)

Every tournament with H=1 is transitive: a nontransitive tournament contains a directed triangle, and adjoining vertices one at a time cannot lower H. To see the latter, insert a new vertex into each old path at a fixed first admissible position; deleting the inserted vertex recovers the old path. A triangle already has three paths. Thus (8) forces a transitive common gadget. Conversely a transitive gadget and a physical coefficient pair satisfying (8) do give a twin-gap family. This is an iff in the declared common-gadget/synchronized-padding architecture, not a classification of all tournament twin pairs.

For the symmetric pair 6Q*2^t +/- H, let g=gcd(H,6Q). Every output is divisible by g. If g>1, the lower count is at least 6Q-H>=11H>g (Q>=2H), so neither output can be prime. When g=1, each output is coprime to the entire product 6QH: modulo a divisor of H it equals the center, and modulo a divisor of 6Q it equals +/-H. Consequently every prime factor is outside the prime support of 6QH. This is a genuine Euclid-style prime-factor escape statement, not primality of the output. The gadget E6 has (H,Q)=(17,224), and yields 1327,1361,1395,1429; the directed triangle has (3,12) and yields 69,75,81,87, all composite.

## 4. Exact local sieve and a closed wheel through 31

For F(t)=A+B*2^t and odd prime q, its zero-residue set is determined exactly:

1. If q divides B, every t is forbidden when q divides A, and none otherwise.
2. If q does not divide B but divides A, no t is forbidden.
3. Otherwise put d=ord_q(2). There is a zero iff -A/B lies in the subgroup generated by 2 modulo q, and then exactly one residue t modulo d is forbidden.

The proof is the equation 2^t=-A/B in the multiplicative group. The powers 2^0,...,2^(d-1) are distinct. For finitely many primes, take the union of their forbidden classes in the least common multiple of their periods. This is an exact periodic divisibility test. Divisibility proves compositeness only when |F(t)|>q, so equality exceptions must be retained. A nonempty surviving residue set is not a prime-infinitude theorem.

For 6*2^m +/-1, the combined forbidden residues are:

| q | ord_q(2) | m giving at least one zero |
|---:|---:|---|
|3|2|none|
|5|4|0,2|
|7|3|0|
|11|10|1,6|
|13|12|1,7|
|17|8|none|
|19|18|4,13|
|23|11|2|
|29|28|8,22|
|31|5|none|

Each row is a finite modular-power identity; the verifier reproduces both signs separately. Combining q=5,7,13 says exactly m=5 modulo 6. The q=11 condition becomes m!=1 modulo 5; q=23 gives m!=2 modulo 11. The q=19 forbidden class is 1 modulo 3, incompatible with m=5 modulo 6; q=29's forbidden classes are even. Thus the complete no-divisor-through-31 condition is exactly

    m=5 mod6, m!=1 mod5, m!=2 mod11.                           (9)

It has period 330 and exactly 40 classes, hence density 4/33 in the exponents. We do not assert this period is minimal. For m>=2, simultaneous primality forces (9): m=2 fails at 25; for m>=3 both values exceed31. At m=1, 11 and13 are prime equality exceptions to their own divisibility tests.

Exact trial division for m=1,...,30 finds twin-prime pairs precisely at m=1,5,17: (11,13), (191,193), (786431,786433). These are FINITE-EXACT instances. The surviving m=23 is a mandatory hostile:

    50331647 = 6563*7669,
    50331649 = 61*825109.

Thus the wheel does not certify primality. The same finite run lists gap-4 prime pairs at m=1,4,5,7, all consecutive. Gap-6 prime pairs occur at m=1,2,3,5,6, but only m=2,3,6 are consecutive. No unbounded claim follows from these lists.

## 5. Why this does not manufacture prime Euclid numbers

Here P_n is the product of the first n ordinary primes, n>=1, and a prime Euclid number means a prime value of P_n+1. Do not conflate this sequence with iterated Euclid–Mullin choices or with arbitrary smooth-product-plus-one values.

**Exact disjointness.** None of the four values (5), for any m>=1, equals any P_n+1. For n=1 the latter is3, below the window. For n>=2 it is1 modulo3. Among offsets -1,1,3,5, only 1 has that residue. Equality would force P_n=6*2^m, but the right side has 2-adic valuation m+1>=2 and every primorial has valuation1. This proves disjointness at all orders, not merely disjointness from prime members.

**Unit-intercept fixed-bank obstruction.** Fix finitely many physical sources (a,b), and use arbitrary nonempty fixed-root gadget words as in (3). A symbolic unit intercept, a * product_i H_i = 1, forces a=1 and every H_i=1. Therefore all gadgets are transitive and every Q_i is a power of2. Such outputs have exactly the form

    1 + b*2^M,

with b from the fixed source bank and M>=1 for a nonempty word. In fact equality with P_n+1 forces v2(b)+M=1, hence b odd and M=1: exactly one singleton gadget and P_n=2b. This is a sharper squarefree gate. Their center's odd prime support belongs to a fixed finite set. If this is to equal P_n+1, then all odd primes in P_n already divide that fixed b, bounding n. The statement is about this exact unit-intercept realization, not a prohibition on incidental numerical coincidences from nonunit-intercept formulas or arbitrary tournament constructions. Allowing the source bank to grow removes the fixed-support argument, but requires supplying the new physical sources and still satisfies the squarefree gate. A concrete nonunit-intercept hostile prevents overextension: A5 rooted at3 has (a,b)=(-5,18), and a singleton gadget produces the strong count -5+2*18=31=30+1, a genuine prime Euclid value. Thus the unit-intercept and four-window restrictions are not unrestricted Euclid-output exclusions.

As additional elementary type barriers, every H is odd, while P_n is even; and any gadget on at least two vertices has 4 dividing Q. The latter follows by pairing complementary subsets: Q/2 is a sum of 2^(n-1) odd terms, an even number for n>=2. Thus a primorial with exactly one factor2 cannot itself be a nontrivial gadget Q. Formal multiplication of integers is not a physical +1 mechanism.

Contrast with a growing-primorial symmetric pair P_n +/-d for fixed odd d>1: a prime factor of d eventually divides both center and offset, forcing both terms composite once they exceed it. The centers 6Q*2^t here have fixed odd support, so that particular growing-wheel obstruction need not apply when gcd(H,6Q)=1. Stabilized support removes one obstruction; it does not prove that prime outputs occur infinitely often.

## 6. Brown roots, squares, and the invalid reverse map

The verifier independently supplies physical root counts 5,11,71, by the following adjacency rows:

    5:  [6,12,8,1],
    11: [24,17,19,6,8],
    71: [206,252,248,240,225,65,128,32].

The last comes from the edge-root word OOZJOZ, with O a sink, Z a fixed-root bridge and J a bridge that becomes the new root. All three counts are checked both by subset DP and by literal full permutation enumeration. Joining each graph with a second copy supplies counts 25,121,5041, checked again on the actual joined graphs (orders8,10,16). The arithmetic identities 4!+1=25,5!+1=121,7!+1=5041 are the known Brown instances; this is a physical supply of their roots and squares, not a theorem that all factorial-plus-one values are counts or squares.

The forward square map is H(A=>A)=H(A)^2. Its inverse is unavailable:49 is physically realized but7 is forbidden. The roots above are supplied independently and never extracted from square counts. The associated joined graphs need not be strong.

## 7. Reproduction and loss ledger

Run from this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > optimized.json
    cmp optimized.json CHECK.json

The script is standard-library-only and imports no repository mathematics. Explicit failure gates survive -O; --fail must exit nonzero. The complete census is all1099 labelled tournaments of orders1 through5 and all5044 nonsink rooted instances. It finds43 full source states and39 (a,b) profiles. All172 combinations of the43 physical states and four gadget pairs represented through order3 are checked by literal graph DP against the full-state graft law. It also checks48 displayed whole graphs,36 iterated/flattened grafts,2964 constant-gap cases,57600 local-sieve cases, all27720 residues of the general wheel, and the Brown graphs above. These finite facts validate the implementation; the unbounded statements rest on the proofs.

Map: actual source and gadget rows -> literal graft tournament -> integer count. Preserved: specified graft orientation, strongness under the proved guard, count formula, coefficient gap. Lost under taking the scalar: endpoint data, prime status, consecutiveness of primes, arbitrary-target placement. Required sidecars: source/root bitrows, gadget pair, and separate arithmetic primality/consecutiveness certificates. Cheapest decisive hostiles: the third count divisible3; m23 surviving the wheel but composite; primorial2-adic valuation;49 versus7.

There is no universal physical +2 operation on all tournaments: a graph with count5 would require the forbidden7. Even excluding that scalar hostile, the inherited one-vertex result rules out any +2 extension from H>=5. Our four-window compares different source roots/graphs and does not contradict that obstruction. Chan–Pak, Polignac and prime Euclid infinitude remain OPEN.
