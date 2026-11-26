# Local power obstructions and a two-hole Newton defect

2026-11-26. This follows [independent responses and the two unfilled branches](../independent_responses_and_two_holes_20261126/README.md). It transfers two exact tests from the linked papers, not their main conclusions. Chan–Pak and the degree-57 Moore existence question remain open.

## 1. Ordowski: a local obstruction can defeat every common spoke depth

[Kohlmeyer and Kruer, *A negative answer to Ordowski’s second question*](https://jenwin.io/papers/ordowski-conjecture-2.pdf) use a prime attaining the maximal valuation and an explicit local witness. Their paper treats a **universal** base congruence, not the distinct fixed-base-2 question. Here the analogue is a fixed-base orbit: a right shift $\sigma$ is possible at bit $i$ when $R2^i\equiv\sigma\pmod H$ and the suffix formula is used. Orders do not determine the common clock without also keeping each phase.

For an odd prime $p$, write $d=\operatorname{ord}_p(2)$ and $t=v_p(2^d-1)$. If $d\nmid k$ then $v_p(2^k-1)=0$; otherwise the binomial lifting argument gives $v_p(2^k-1)=t+v_p(k/d)$. Thus the lifted order is $d p^{\max(0,a-t)}$ modulo $p^a$. Given a residue in the cyclic subgroup, the exponent is one congruence class modulo that lifted order. Across prime powers you must choose the **same sign** for a given core and check the exponent phases by generalized CRT, not independently solve each modulus and discard its phase. A concrete fixed-base failure is $p=1093$: $d=364$ and $t=2$, so $1095=2+1093$ has the right residue mod $p$ but lies outside even the signed power orbit mod $p^2$. This is an arithmetic example, not a claimed tournament profile.

For the actual seven cores and a unit response, the exact phase table is:

| $(H,R)$ | phase for $+2$ | phase for $-2$ | period |
|:--|--:|--:|--:|
| $(23,24)$ | 1 | none | 11 |
| $(25,28)$ | 14 | 4 | 20 |
| $(29,31)$ | 0 | 14 | 28 |
| $(59,87)$ | 39 | 10 | 58 |
| $(71,72)$ | 1 | none | 35 |
| $(103,129)$ | 3 | none | 51 |
| $(143,173)$ | none | 53 | 60 |

The $(29,31)$ and $(71,72)$ cores are each unit-capable, but they require $i\equiv0$ and $i\equiv1\pmod7$ respectively. There is no common bit index. A second obstruction is the odd phase of $(143,173)$ versus the even phases of $(25,28)$ or $(29,31)$. Exact CRT over all subsets of the seven profiles gives a maximum of **five**, for the unique subset $H=23,25,29,59,103$, at phases $1{,}091{,}454$ and $2{,}230{,}284$ modulo $2{,}277{,}660$. Two vertex-disjoint edits on one tournament remain possible; this excludes a stronger demand across these separate recipes.

Changing the common amplitude instead of the shared depth can restore synchronization. For a fixed index $i$, the seven path-count moduli are pairwise coprime, so imposing either sign independently gives exactly $2^7$ residue classes for the half-amplitude modulo their product. At **bit 1**, take the all-negative signs. The replay obtains half-amplitude $736{,}403{,}033{,}622$. For each core $\beta=(2R+2m)/H+2$ is even and positive (zero suffix), and at $c=4\beta+1$ the right spoke has shift $-2m=-1{,}472{,}806{,}067{,}244$. The independent left edge can give $+2$. Thus the shared-depth obstruction is genuinely about the prescribed unit, not about the existence of a common nonzero response.

Another exact lift: for $H=11^a13^b$ with $a,b\ge1$, $\pm\langle2\rangle$ is the whole unit group. The two lifted orders are $10\cdot11^{a-1}$ and $12\cdot13^{b-1}$, with gcd 2. Adding a common negative sign flips the parity of the first discrete logarithm and preserves the second, so the signed CRT always succeeds. This is an arithmetic extension of the $H=143$ lane, not a claim that the other lifted profiles exist.

## 2. CFH: exact low-overlap counting, not a real-rootedness shortcut

The [Davis counterexample and Adler rank-threshold supplement](https://github.com/1337isnot1337/solved-conjectures/tree/main/results/caiolo-ferroni-hoster) for [Caiolo–Ferroni–Hoster](https://arxiv.org/html/2605.28474v1) use controlled intersections of replacement blocks, exact flag arithmetic and a discriminant or Newton inequality. The transfer here is the verification method: isolate what the overlap contributes, then test the coefficients directly instead of importing real-rootedness. The following polynomial is a collision enumerator on the marked Moore board, **not** a Chow polynomial.

Keep the same 50 reference contacts and collect the two-hole score by $z^{\text{squares}-3}$. Complete enumeration of the $32\times33$ physical pairs gives

    P(z) = 3 + 63z + 117z² + 811z³ + 62z⁴.

It fails even ordinary log-concavity: $117^2-63\cdot811=-37{,}404$. Exact Sturm counting gives two real roots; its discriminant is $-386{,}162{,}726{,}057{,}424$. Removing all 71 hole pairs that share a center slot still gives

    P_no-mutual(z) = 3 + 60z + 114z² + 808z³,

which has one real root and negative discriminant. Thus the failure is **not** caused solely by the hole–hole interaction. It survives in the sector where the two hole scores add independently against the references. A real-rooted surrogate for this repair distribution would miss the actual coefficient behavior.

For each hole pair store its old-owner hit counts $h_j$ and mutual collision. Releasing a subset $D$ of old contacts **without supplying replacements** subtracts exactly $\sum_{j\in D}h_j$ from its score, reproducing the induced enumerator without reclassifying pairs. This is the board analogue of counting which terms a low-overlap replacement can affect. A released board is a relaxation; it cannot certify that a removed contact can be rerouted.

| Released references | Sets tested | Fail log-concavity | Relaxed minimum 3 / 2 / 1 |
|--:|--:|--:|:--|
| 0 | 1 | 1 | 1 / 0 / 0 |
| 1 | 50 | 39 | 42 / 8 / 0 |
| 2 | 1,225 | 434 | 796 / 414 / 15 |

The failure of matroid exchange is literal on the marker-filtered contact set (before fixing the 50 references). Contact $(33,19)$ with center slots $(19,42,45)$ clashes separately with $(12,19)$ at center 56 and $(38,24)$ at center 0. The latter contacts have triples $(19,7,55)$ and $(24,42,31)$ and are mutually compatible. The singleton cannot augment toward the compatible pair. This is an intersection of partition constraints, not a matroid, so matroid exchange and matroid root theorems cannot be used directly. The quartic's derivative discriminant is also strictly negative, $-3{,}164{,}395{,}667{,}256$; together with $P(-1)=-692$, $P(0)>0$ and its positive leading term, this gives a second exact check on the count of real roots.

For a bounded solver these relaxed minimums give a safe first filter: a specified edit support with relaxed minimum at least three cannot yield a two-square board. The existing exact marked search is stronger at its completed edit limits. These polynomial tests neither fill branches 13 and 33 nor resolve the pending split caps at depth 12.

## 3. Three disjoint physical edits, with two opposite units

There is a concrete prime-square application of Section 1. Take the existing outer $(H,E,R)=(23,8,24)$ core and take inner masks $(14,28,24,16,33,15)$ with root 5. Reverse the inner arc $1\to4$. The exact inner rooted profiles $(h,e_-,e_+,r,v,Q)$ are $(23,1,8,24,75,296)$ and $(23,1,10,24,87,296)$. So the inner edit leaves $A=529$ and $B=7656$ unchanged, and increases the offset by $23\cdot2+24\cdot12=334$. The outer edit increases the offset by two since the terminal multiplier remains one; the inner edit does not alter that multiplier. Both hold for every odd binary address $c$. With edit indicators $a,b$ for the outer and inner arcs, the count is

    529 F(c−1)+7656c−1992−2a−334b.

For the third arc choose $i=191+253k$, $\ell=2^i$ and $\beta=(8185\ell+2)/529$; set $c=2\ell\beta+1$. The congruence $7656\cdot2^{191}\equiv-2\pmod{529}$ gives an even $\beta$. The right spoke at bit $i$ then lowers the count by two from address $c$ to $c+\ell$. Start at $c+\ell$ and reverse it for an increase. Across all eight orientations the count is **$N-2a-334b+2d$**. This proves independence including the third mixed difference. The three arcs lie in separate vertex blocks. All orientations are strong by the marked bridges, the inner unique terminal path, and the highest spine bit. The smallest unit instance of this explicit index progression has 400 vertices; the replay independently counts all eight orientations of a smaller 16-vertex control.

Moreover $v_{23}(2^{11}-1)=1$ and $\operatorname{ord}_{23^2}(2)=253$, an odd index-two subgroup. Consequently $\pm\langle2\rangle$ covers all 506 units modulo 529. Since $7656$ is a unit, the same three physical channels realize magnitudes **$(2,334,2m)$ for every $23\nmid m$**, at unbounded heights. One may set the signs independently by choosing the starting corner.

An ordered incoming singleton on the inner core kills both root-start counts but preserves the two starting-second cover counts 75 and 87. It doubles the other two cover counts and leaves the ordinary and terminal counts unchanged. This gives $B=15312$ and middle amplitude **288** instead of 334, with forward unit spoke bit $i\equiv190\pmod{253}$. The same signed orbit realizes $(2,288,2m)$ for every $23\nmid m$.

The separate complete inner-core census through order nine checks this specific pure-offset architecture. It preserves $h,e_-=1,r,Q$ in both orientations and divides the changes into $(\Delta e_+,\Delta v)$. The numbers of nonzero candidates at orders five through nine are 2, 6, 8, 32 and 32. None has terms of opposite sign. With the seven outer cores considered, the smallest nonzero middle amplitude is 144 without imposing a unit spoke; the smallest that also admits the third unit is 288. The order-nine calculation does not exclude three units in another architecture or at higher order.

## Replay

- `node three_arcs.js` verifies both the literal and large-index nested cases; compile `pure_offset_census.cpp` with `g++ -O3 -std=c++17` and pass it order `9`.
- `node ordowski_clocks.js` computes the phase table, obstruction and common-amplitude examples, using exact integers.
- `python collision_polynomial.py board.txt` enumerates the two-hole polynomial, release counts, and Sturm certificates. Generate the same board with `python nested_core_check.py --board > board.txt` from the [pinned Moore chart](../moore57_quadratic_tensors_20261125/).
