# Endpoint arithmetic and aliquot ancestry of the fourteen-step relay

This uses the earlier [lonely-runner endpoint criterion](../../../01-canon/theorems/THM-357-lonely-runner-endpoint-protection-trichotomy.md), its [necessary endpoint cycle](../../../01-canon/theorems/THM-365-lrc-endpoint-protection-cycle-necessity.md), the [divisor phase clock](../../../01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md), and the [retained aliquot guards](../aliquot_guarded_cycles_20261014e.md). No implication between the conjectures is assumed. The resulting advances are a faster exact phase clock, a bounded exact inverse for one relay, and a literal changed-interface reversal of both aliquot edges in a twin-prime fork. The universal selector is still open.

## 1. Lonely-runner lesson: keep the endpoint and its arithmetic label

In the old endpoint work, full interval measure did not exclude a boundary witness. The relevant clock was finite but marked: to protect a candidate endpoint, an actual integer inequality had to hold. We use the same distinction between a real band and its integer cells for the [single-arc relay](../fourteen_axis_20261117/SINGLE_ARC.md). Its lower read is

    H(z)=C+A(z+1)+gF(z),   (C,A,g)=(-5884,13190,37), z≥0,

where z=0 means no further continuation, and z>0 encodes a word ending in a near-sink. The other read is precisely H(z)+14.

For prescribed N, take L=2^k≥2, the first power with A L²≥N−C. Split the code z=Lq+r. Write X_r for the k-bit reverse complement of r and δ_r=1 if the top low bit is zero. The exact low-bit carry is

    F(Lq+r)=F(q)+X_r q+δ_r(q mod 2)+F(r)−1.

Put B_r=AL+gX_r and Y_r=N−C−A(r+1)−g(F(r)−1). Then an actual integer in the r-cell must obey

    ceil((Y_r−2g)/(B_r+g)) ≤ q ≤ floor((Y_r−g)/B_r)
    B_r q ≡ Y_r (mod g),   q≥0.

The interval follows from 1≤F(q)≤q+1. For any live r its width is at most g/A+2g/(AL)<1, by Y_r≤N−C and A L²≥N−C. Thus a cell has at most one candidate integer, and evaluating F there is an exact accept/reject test. The modular label is also reversible: gcd(A,g)=1 and B_r≡AL mod g. Tables for F and reversed low addresses take O(L)=O(√N) arithmetic operations, since q≤(N−C)/(AL)≤L. The distinction between a live real interval and actual integer landing is crucial even when every residue could be steered separately.

This is a complete inverse of this relay only. Both 30,011 and 30,013, and both billion twins, fail as lower quotients in it. Those values already have other witnesses. The relay succeeds at lower prime quotients 3,889,007 and 5,417,681 with addresses 2063 and 2870; their +2 quotients are prime too. Also, to find two addresses inside the **same lower** relay whose outputs differ by 14, the exact modulus g forces Δz≡9 mod 37. Neighboring counts do not correspond to neighboring addresses; the direct one-arc edit moves between two parallel reads rather than being a tiny address change.

## 2. A logarithmic phase compiler

The earlier slope/phase witness used O(m) elementary moves for the odd part m of a modulus. Retain the exact clock instead of choosing a long walk one residue at a time. For T_k on the low address u, put v=u+1 and h=2F(u)−v. Modulo odd m the effects with e=2^(k+1) are

    (v,h) ↦ (e v−1, h+v+1).

We use e=1 or e=2; each stands for a legal k at least the current dyadic cutoff, found within one period ord_m(2). Before this odd step, program F mod 2^a by the old bit toggles. Higher k preserve it.

Let d=ceil(log₂ m). On w=v−1, applying e=2 and optionally e=1 gives w↦2w−bit. Thus d bits steer v from anywhere to any prescribed v, by retaining the required residue in [0,m). First steer to 1; e=2 fixes 1 and adds 2 to h. From 1 steer to 4. At digit j the common v will be 3+2^j. The words (1,2) and (2,1,1) have the same map v↦2v−3; the second adds 3v−1 more to h. Choosing a d-bit integer z therefore adds 3z+8 popcount(z), compared to choosing all first words. The final v is independent of the digits, so steer it to target V by a fixed word.

To hit h_target, calculate the baseline h₀ with all first words and no loops at v=1. Choose r∈{0,1,2} so that h_target−h₀−8d−2r is divisible by gcd(3,m), then solve

    3z ≡ h_target−h₀−8d−2r (mod m),       0≤z<m.

Back at v=1 insert r+4(d−popcount(z)) loops before the steering to 4. This compensates the popcount exactly. The result programs both v and F using O(log m) elementary moves, not O(m).

Thus the fixed low suffix for joint magnitude/phase in the [fourteen-axis inverse](../fourteen_axis_20261117/README.md) uses O(a²+(a+ord_m(2)) log m) bits. For the odd clock m=2^j−1 the order is precisely j, so the cost is O(j²) with dyadic cutoff a=O(j). The existing magnitude placement then works with an odd modulus as large as exp(c_ε√log x), giving error O(x^(1/2+ε)) at quotient scale x. For a fixed odd sieve m₀ take j a multiple of ord_(m₀)(2), so that m₀ divides this clock. This is a quantitative upgrade of the congruence selector, not a modulus larger than the error or an exact landing.

## 3. Aliquot ancestry really divides a fourteen-gap by seven, with a guard

Write s(n)=σ(n)−n. For 7∤q,

    s(7q)=q+8s(q).

Consequently s(7q)=q+8 for q>1 iff q is prime. An odd composite misses this affine read by at least 24, or at least 40 when coprime to 3. The exact form of the prime guard is retained, not a drift heuristic. Also on a fixed 3-free core m the retained σ(m) yields s(3^(k+1)m)=3s(3^k m)+σ(m). With m=7 this gives s(7·3^k)=5·3^k−4 and the first three images 1,11,41. Arithmetic descent is exact; it is not automatically graph descent.

The twin branch has a strict renewal obstruction in the old aliquot compiler. A twin pair q,q+2>3 with 7-free quotients sends the two lifts to q+8 and q+10, which are coprime, so there is no common nonunit next cofactor. The upper successor q+10 is always divisible by 3, so using a shared unit and prime factor at both ends also fails. Conversely a fixed aliquot template c·p read forward at twin primes p,p+2 coprime to c changes by 2s(c). To get 14 requires s(c)=7, whose only solution is c=8: even c≥14 already have proper sum ≥1+c/2; checking smaller evens gives 8; an odd prime or prime square fails, while any other odd composite has two distinct nonunit proper divisors and sum at least 9. The resulting read is 7p+15, even and 1 mod 7. So neither direction gives a simultaneous fixed-core compiler of tournament witnesses.

### A literal changed-interface inverse of the twin fork

For the seed 7·1049 and 7·1051, the aliquot successors are 1057=7·151 and 1059=3·353. Both admit literal one-vertex inverse returns, but with different attachments:

- From rows (2,36,313,83,355,457,135,287,75) of count 1057, append a vertex with outgoing mask 226 into those nine vertices; all its other incident arcs are incoming. The new count is 7343.
- From rows (160,13,73,97,15,86,147,318,127) of count 1059, append a vertex beating only vertex 8; the marked E=3149 gives H'=1059+2·3149=7357.

These are existential physical lifts at both targets, not the same underlying graphs as the earlier single-arc pair. A universal root of this form would have to produce E=3H−28 at a smaller H=p+8, since 7p=7H−56. The upper seed meets it; it is not guaranteed by H alone.

### Why a fixed sequence of uniform attachments cannot enforce that inverse along a dyadic family

Start with a normalized dyadic rooted family, H and 2E having equal leading coefficients, O−E and g fixed at that scale. Set r=leading(2E)/leading(H). For a fixed appended block of count u and two-slot count v, put k=v/u. The old one-/two-cover rules for a uniformly oriented ordinary block or one beating only the root give the slope maps

    ordinary: r ↦ k r;
    near-sink: r ↦ k r / (1+(k−1)r).

For the second formula, with Z=E+O−g the exact H' is u(H+2E)+(v−2u)Z: the no-repeat skeleton contributes the first term, and two nonempty block runs contribute the second. The exact E'=vZ/2 and O'=E'+ug, because the block two-slot split gives equal even/odd contributions. In the ordinary case H'=uH and E'=O'=(v/2)(E+O). The maps therefore preserve leading equality of E and O and keep g bounded with respect to dyadic height.

Every nonempty physical block has k≥2, and the only k in [3,6] is 4. For orders through four the pairs are obtained directly (the strong four-vertex case has 5,32). At order at least five, splitting a Hamiltonian path gives k≥n+1≥6, and equality 6 is impossible since v is divisible by 4 at order at least two (each term of the odd-cycle two-slot expansion has that factor) while u is odd. Thus the inverse of r=6 after the last near-sink, which has r≤2, requires aggregate ordinary factor 4 and prior r=3/2. The same holds if there was no near-sink, taking prior r=1. Now the only way to attain 3/2 is a near-sink singleton (k=2) from r=3; k≥4 has ceiling at most 4/3. The only way to attain 3 is an ordinary singleton from r=3/2. Neither can occur first from r=1. This is a labelled backward cycle in the sense of the endpoint-protection analogy, not a proof about the circle of runners. It excludes this entire fixed-block-word mechanism for E=3H−28 at unbounded dyadic heights, without excluding genuinely changing cores.

Replay with `python inverse.py` and `python phase.py`; `-O` gives the same checks. No finite replay is used to justify the all-order compiler or the slope-cycle obstruction.
