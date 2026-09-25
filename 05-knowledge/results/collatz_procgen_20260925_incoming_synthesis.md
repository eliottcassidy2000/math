# Incoming synthesis: codex session `collatz-bugs-20260925` (decoders, seams, sewing, triples), audited against the procgen session

**Status: SYNTHESIS + INDEPENDENT AUDIT** of the 11 codex commits `3ef7145d74..8e475b4f3d` (2026-09-25): 40 changed result notes, the MISTAKES diff, and the edits to THM-060 and THM-4473. Twenty independent spot checks (own code, `04-computation/experiments/procgen_incoming_20260925_audit_*.py`, none of codex's scripts imported) found **no numerical or logical disagreement**. The flags are interpretive:
* the C9 "plus nonlinear / minus linear" result is a digit-frame artifact, not a sheet invariant;
* the "8a+1 descent family" is the classical stopping-time-2 class;
* no codex lane runs a DRIFT (5n+1) control.

Four observations are new here:
* **Q_N planarity (PROVED by two finite certificates).** `Q_N` is planar iff `N <= 24`, and its first K_{3,3} passes through the ear `11-25-24` that every Hamiltonian path of `Q_25` must use.
* **Pfaffian and THM-4472 (FINITE-EXACT).** `|Pf| = H = 2+sgn(b s_o)` on the THM-4472 quadruples.
* **Berggren sign content (FINITE-EXACT).** The sign content of codex's Berggren theorem is exactly the minus 2-cycle `{5,7}`.
* **Exchange law (PROVED, one line).** The adjacent-exchange law does not depend on q, so THM-4469's adjacent carries are {2,3}-unit equations.

Collatz, E-SCC, the Periodicity Conjecture and HYP-9120..9136 all remain **OPEN**; nothing here changes their status. No HYP or THM file was created.

Session: collatz-procgen-20260922 (incoming lane, 2026-09-25). Scope: the 40 result notes changed in `3ef7145d74..8e475b4f3d`, read in full. Our side: synthesis §§2f–2g and THM-4469..4473.

## 0. Bottom line

1. **What codex built.** Codex produced a large, carefully labelled body of *representations*: tournaments, marked Pythagorean triples, elliptic points, Fano/E8 lattices, ten-letter alphabets and affine word matrices. For each one it gives a proved "decoder" and a list of the sidecar data it must retain. It also proved several limitation theorems:
   * forward closure;
   * the finite-bank obstruction;
   * no fixed-residue polynomial potential (F5);
   * no bounded-low-digit rank (D9);
   * no finite-state itinerary encoder (C10).

   Every open obligation is stated honestly as OPEN. In the terms of our foundry, the decoder/creative program is a detailed build-out of **one of our two unblocked mechanism types, "sound certificate searches"**. Its no-go theorems confirm that such searches need unbounded parametric schemas. Beyond the lookahead-type obstructions already known (the all-ones prefixes `2^a-1`, i.e. the 2-adic point `-1`), they do not locate the missing mechanism.
2. **Genuine repairs to the repository.** Three are legitimate and verified:
   * the square-sum degree-2 forcing error, with a new endpoint-aware `Q_24` proof that I re-derived by hand;
   * THM-060 Type A (a transitive triple contributes 0, not 2);
   * THM-4473's `k >= 2` boundary (single-digit primes).
3. **Relevance to our objects.** Most codex constructions answer owner prompts numerologically: exact identities with no map to Collatz dynamics, such as the fruit triple, the several "36"s, `(10^k-7)/3`, Mills, the Zenodo symbols and the E8/Clifford frame. The pieces that do touch our program:
   * the Mahler exchange law and root-basin sign transport (THM-4469);
   * C10, which closes the automata route to PC;
   * the Berggren consecutive-edge theorem, a sign-aware local law whose minus exceptions are exactly the `{5,7}` cycle entries;
   * the Pfaffian switching invariant (THM-4472);
   * the signed inverse-fibre half-step `F(n,sigma)=(2n+sigma,-sigma)`, which interleaves the two sheets' inverse fibres.
4. **Seeds.** Section 5 extracts 7 ranked proof-angle seeds and pre-screens three of them. The Berggren antichain seed turns out to be DRIFT-blind: divergent 7n+1 orbits satisfy it. The 6-vertex Pfaffian is local. The automata route to PC fails by counting. Four further directions are recorded as dead on arrival.

## 1. Per-family digest

### 1.0 Glossary (codex terms in plain mathematics)

| term | meaning |
|---|---|
| **bug** (collatz_bugs) | An even "identity" `I=2^m q` with the reversed-arrow motif `I/2 -> I -> 2I -> 4I`, plus the odd leaves `a=(I-1)/3` (child of `I/2`) and `4a+1` (child of `2I`). It is admissible iff `I = 4 mod 6`. |
| **f(X,Y)** | The number of edges traversed *against* the reversed arrows on the path from X to Y in the inverse tree of 4 (the root basin C with `{1,2}` removed). It equals `h(X)-h(LCA)`, where h is the height to 4. |
| **decoder** / **sidecar** | A procedure that recovers a target operation (halving, a Collatz step, a route to root 4) from a compressed representation plus explicitly listed extra data (the sidecar). It is always sound. It terminates only where that is proved. |
| **E(n)** (decoder_halving) | `E(2^k q)=H_q[TT_(2^k)]`: a regular odd tournament core `H_q` with transitive fibres. **Pair module**: a 2-vertex module. |
| **four-core** | The marked quotient Q of the substitution `F_Q(H)=Q[H,H,H,{r}]`, which has order `3|H|+1`. |
| **seam** | The "doubles seam": the layer `v2(n)=1` between the odd chain (`v2=0`) and the "sea" (`v2>=2`). Adding 2 keeps odd numbers odd and alternates between the seam and the sea. In seam_mills it is also the prime-shell gap digit. |
| **port** | A parity word w with sign `sigma`, its exact cylinder `n=r+qk` and affine exit `s+pk`, plus the threshold K above which the exit descends. |
| **guarded join** | `(x,y,a,b)` with `0<y<x` odd and `U^a(x)=U^b(y)`. It transports a root certificate from y to x via `h(x)=a+max(h(y)-b,0)`. |
| **run passport** | `(a,u,b,m)` with `n+sigma=2^a u` and `3^a u-sigma=2^b m`, i.e. the maximal block `1^a0^b`. |
| **sibling ladder** | `S_sigma(n)=4n+sigma`, with `U_sigma(S_sigma n)=U_sigma(n)`. **Signed half-step** `F(n,sigma)=(2n+sigma,-sigma)`, with `F^2=S_sigma`. |
| **creation certificate** | `P=2P'+(bits)`, where the bits are read intrinsically from P, together with a terminal bank and a strictly decreasing integer rank. Reading it backwards "creates" P. |
| **Fano carry decoder** | The depth-3 parity-address map `q_(+,3)`, which is the transposition `(1 5)` on `Z/8`. Its polar defect is `4(x0y1+x1y0)`. The transported law `x*y = x XOR y XOR 4(x0y1+x1y0)` makes `q_(+,3)` a group isomorphism. On E8 = Construction A of RM(1,3), the swap `(1 5)` moves E8 to a neighbour E8' sharing a D8. |
| **thirty-six phase** | The `36=2*18` relative-phase square roots of the minimal square-root completion (cycle type `2^2 5 18^2`, 45 states) of the ordinary-clock minus cycles. It is also the 36-cycle of the global half-step R on the 18-cycle through 17. |
| **sewing** | Joining a forward prefix value `(A_u n+B_u)/D_u` in `Z[1/2]` to a backward suffix value `(D_v r-B_v)/A_v` in `Z[1/3]`. They agree only at an integer, and agreement holds iff the word `uv` is a legal route. The residual `K_(n,r)(u,v)` has rank 2. |
| **marked triple** | `Phi(s/t)=(st,(s^2-t^2)/2,(s^2+t^2)/2)`, with the even leg signed: a negative leg marks the reciprocal state. |
| **fruit** | The positive primitive solution of `a/(b+c)+b/(a+c)+c/(a+b)=4`. This is `9G` on `y^2=x^3+109x^2+224x` (Bremner–Macleod), with digit lengths 81/80/79. |
| **color / charge** (duck) | `Q(n)=(n-2b(n),b(n))` with `b(n)=floor((n+1)/phi^2)`. The color is `c(n)=Q(n) mod 2` in `F2^2={0,R,G,B}`. |
| **ten letters** | `D=Sym^2_set(F2^2)`, the 10 unordered pairs with repetition. It is in bijection with `E(K5)` and with the proper divisors of `p^2qr`. |
| **duck**, **zenodo** | Lane labels. "duck" has no mathematical definition in the notes (my reading: an owner prompt name; uncertain). "zenodo" is the owner-supplied DOI 10.5281/zenodo.22800071 (He–Jing–Li, six-gluon MHV symbols through nine loops). |

Status labels below are codex's own, as written. My reading follows the arrow `->`. "Numerology" means an exact identity with no map to Collatz dynamics.

### 1.1 collatz_bugs (1 note)

**Prompt answered (inferred).** The owner's "bug" construction:
* a directional count of backwards traversals between bugs;
* links to square-sum Hamiltonian paths ("negative Collatz"), prime encodings and graceful labelings ("positive Collatz");
* a follow-up hypothesis: square-sum and graceful constructions work in both directions, parallel to 3n+1 and 3n-1.

**Claims.**
* **PROVED:**
  * bugs exist iff `I=4 mod 6`;
  * `f(X,Y)=h(X)-h(LCA)`, and `d=f(X,Y)+f(Y,X)`;
  * universal finiteness of `f(I,4)` is equivalent to Collatz;
  * the ladder `a_(j+1)=4a_j+1` has coprime neighbours, and `gcd(a_j,a_(j+r))` divides `(4^r-1)/3`;
  * **square-sum/graceful:**
    * the sign gauge `z=+-x` turns square sums into square differences (THM-2761);
    * `Q_n` never contains a spanning star, since its maximum degree is at most `floor(sqrt(2n-1))-1`;
    * paths are graceful while `Q_3` has no spanning path;
    * `T_-(-n)=-T_+(n)`;
    * no bijective conjugacy exists between the positive plus and minus systems (`T_-` fixes 1, `T_+` has no positive fixed point).
* **OPEN / proposed:** the square-sum ↔ graceful and 3n+1 ↔ 3n−1 "isomorphisms". No map is established.

-> PROVED-elementary. The equivalence is a restatement. The star obstruction is a clean negative that strengthens our "graceful: ANALOGY" verdict. The sign identity duplicates our transport theorem.

### 1.2 collatz_tournament_core (1 note)

**Prompt answered (inferred).** Encode integers by tournaments made of three copies of an odd-order H plus one vertex (order `3A+1`). A four-vertex core is meant to certify routes to root 4.

**Claims (PROVED).**
* The edge count is `C(3A+1,2)=3C(A,2)+3(A-1)^2+9(A-1)+6`; the owner's `9A` term overcounts by 9.
* Substitution keeps the quotient Q.
* Hamiltonian paths exist but give no rooted arithmetic certificate: the root can be a sink.
* `Q[C3,C3,C3,1]` (source over C3) has **no two-vertex module**.

-> PROVED-elementary; no Collatz content. It explicitly refuses to identify its core with THM-4472's quadruple, which is correct.

### 1.3 decoder_* (5 notes)

**Prompt answered (inferred).** Decode tournament structure into Collatz moves; odd chain / doubles seam / sea; intrinsic halving; links to Mahler/Catalan, graph minors (K5, K_{3,3}, Robertson–Seymour), square sums and `{2,3,11}`.

| note | claims (label) | reading |
|---|---|---|
| decoder_halving | **PROVED**: <br>• nested regular odd cores grow by 2; <br>• the pair-module graph of `E(2^k q)` is `q P_(2^k)`, so every even order has a unique intrinsic halving matching with quotient `E(n/2)`; <br>• every raw `Q[H,H,H,1]` (H regular, odd, `>=3`) has no pair module. <br>**FINITE-EXACT**: guarded D/O root words for all `n<=2000` (max length 179). <br>**OPEN**: generating root words from structure. | The halving is intrinsic *to an encoding built to contain* `v2(n)`, so the dyadic address is circular by design (acknowledged). No Collatz content. |
| decoder_mahler_catalan | **PROVED**: <br>• exact root verifier `3^a n+R(w)=4*2^L`; <br>• minimal clock-loss witness: `01` vs `10` give 14/3 vs 5; <br>• exchange laws `R(u01v)-R(u10v)=2^j 3^b` and `M(u01v)-M(u10v)=-2^j 3^h`; <br>• carry extremes `3^a-2^a <= R <= 2^(L-a)(3^a-2^a)`; <br>• shared-word image of root 4 is `-4/15`, every root-certified n maps to a negative rational, and the suffix `10101` violates Mahler's tail (`266/243>1`). <br>**CITED**: Mihailescu; Andrieu–Eliahou–Vivion. | **Substantive for THM-4469** (section 3). |
| decoder_minors | **PROVED**: the graph G_N (`{x,x+2}`, `{x,2x}`) is planar iff `N<=15` (rotation at 15, K_{3,3} subdivision at 16); row-scaled `+2` gives the planar grid; strong tournament minors (Kim–Seymour) differ from pair quotients; guarded words are not subword-closed (`DDO` legal, `DO` not). **CITED**: Kuratowski, Wagner, Robertson–Seymour. | Exact; relevant only to the owner's Kuratowski theme (section 4). |
| decoder_pair_repair | **PROVED** cost formula and lower bound 3. **FINITE-EXACT**: unrestricted minima `{3:24, 6:24, 7:12, 8:4}` over the 64 cores; canonical `H5[TT2]` minima `(15,13,12,7,10)` by class. | Pure tournament combinatorics. |
| decoder_prime_square | **PROVED**: two-sheet halving `Q_(2N)^(c)[even]/2 = Q_N^(3-c)`; dyadic layer rule; `{2,3,11}` bracket escape (inherited from us); endpoint-aware nonexistence proof for `Q_24`; every `Q_25` Hamiltonian path uses `11-25-24`. **FINITE-EXACT**: path census for N = 14..25. | Correct. It closes the slot left OPEN by codex's own retraction. |

### 1.4 seam_* (4 notes)

**Prompt answered (inferred).** The relations between `4A, 4B, A^4, B^4` and `x^y=y^x`; "what can be ignored above cubic"; divisor balance with values 0, 2, 10; Mills' constant; the thresholds 15/16; "triangular ten" versus a regular five-tournament.

| note | claims (label) | reading |
|---|---|---|
| seam_mills | **PROVED**: <br>• prime-shell chains ↔ unique nested-interval constants; <br>• the fractional carry `theta` is essential (`23.66` vs `24.389`); <br>• `v2((p+2)^c-p^c)=1` for odd c, and `v2(c)+v2(p+1)+1` for even c (LTE); <br>• `p^3+3` is prime only for `p=2`, giving `11=2^3+3`. <br>**FINITE-EXACT**: greedy prefixes, and the cubic bracket `1.30637788386308 < A < 1.30637788386948`. <br>**CITED**: Mills, Caldwell–Cheng, Saito, Matomäki. | Numerology relative to Collatz; correct. |
| seam_power_coordinates | **PROVED**: <br>• `R(A,B)=A^B/B^A` identities; <br>• the only integer tie is `(2,4)`; <br>• power comparison stabilizes from squares; <br>• the gcd-normalized cycle `16->81->256->16`; <br>• odd log-gap series with coefficients `(4j+1)/(j(2j+1))`; <br>• `H(H5)=15` (the quintic layer contributes 4); <br>• Sophie Germain factorization. | Correct; no Collatz content. |
| seam_prime_balance | **PROVED** (re-derived): <br>• `F=S+U` iff `p, p^3, p^2qr`; <br>• any prime-fourth-power divisor safely rejects balance; <br>• an 8-profile controller with an absorbing reject state; <br>• the filter is not Collatz-closed (`16->8`, `81->61`). | Correct; no Collatz content. |
| seam_threshold | **PROVED**: <br>• `Q14` is a 3-leaf tree, and the ear `1-15-10` repairs it; <br>• in G15, the ear `8-16-14` completes K_{3,3}; <br>• `45=30+15` = edges of `L(K5)` plus edges of Petersen; <br>• `144=24+120` perfect matchings of `L(K5)`, with the 24 regular 5-tournaments as K5-minor decoders; <br>• K5 cannot be split into three triangles plus one edge. | Exact; the owner's K5/Petersen/K_{3,3} theme appears here (section 4). |

### 1.5 creative_* (4 notes)

**Prompt answered (inferred).** Replace the tournament decoder by something that really produces root certificates.

| note | claims (label) | reading |
|---|---|---|
| creative_decoder | **PROVED**: <br>• join transport (J1–J2); <br>• forward and inverse port formulas (P1); <br>• **forward-closure obstruction**: no inverse-orbit join enlarges the full-sibling family; <br>• finite-bank obstruction via `2^a-1` (a even, `> L+2`); <br>• the unbounded schema `N=2^a(2^b+1)/3^a-1`. <br>**FINITE-EXACT**: <br>• local coverage 19/24 (15 residual classes mod 144); <br>• 953 inverse words with `L<=8` (I checked the count); <br>• the adaptive compiler certifies all odd `n<=10^4` with 406 learned seeds. | Correct. The compiler's seeds come from orbit search (acknowledged). |
| creative_descent | **PROVED**: <br>• ports and thresholds; <br>• adaptive cover with an explicit residual; <br>• the balanced counter grammar `1^a0^a`; <br>• exact root families (P7) for both signs; <br>• the phase of b modulo `2*3^(a-1)`. <br>**FINITE-EXACT**: depth 14 gives 142 ports and 734 residual classes. | Terras-level plus a classical landing family; correct. |
| creative_sibling | **PROVED**: <br>• signed sibling normal form; <br>• common-future transport; <br>• the infinite family `(32*64^l-5)/9`; <br>• `Q+(n)<n` iff `n=1 mod 4` or `3 mod 16`; <br>• forward closure (10). <br>**FINITE-EXACT**: closure counts 252/274/541/640 up to `10^5`; first uncertified source 7. | Correct. The grammar certifies about 1.3% of odd `n<=10^5`: an infrastructure result, not a mechanism. |
| creative_transducer | **PROVED**: <br>• the 3-state carry machine (D1); <br>• signed passports (D4); <br>• the cylinder compiler; <br>• the bounded-prefix hostile `n=2^(k+1)-sigma` (D9). | Correct; inherited mechanisms. |

### 1.6 ternary_* (4 notes)

**Prompt answered (inferred).** The owner's pasted `(10^k-7)/3` primes (31, 331, ...), `F5`, the cubic `2n^3+4n^2+n`, "12 primes sum to 196", and Pythagorean/Berggren "ternary" recursion versus Collatz and repunit clocks.

| note | claims (label) | reading |
|---|---|---|
| ternary_berggren | **PROVED**: <br>• the Calkin–Wilf ordinal chart (inherited); <br>• the fixed-target ray with first-return height `2^(k-1)`; <br>• the ternary odometer `v3(H(j)-H(i))=v3(j-i)`; <br>• **consecutive nondegenerate plus edge-triangles are never Berggren-comparable at any distance**, while the minus orbit `27->5->7` has `(27,5)=B1^2(7,5)`; <br>• the projective invariant `J=tr^2/det` is an integer on Berggren words and never on Collatz words. | The comparability theorem is the one new sign-aware law in the batch. Section 2 identifies its sign content. |
| ternary_bridge | **PROVED**: <br>• the gcd splitter `g=gcd(3,t)` is 1 or 3; <br>• the strata `P_r` enter the integers in exactly r steps; <br>• the hypotenuse contracts, `C'<C/3`; <br>• the intrinsic triple update (7); <br>• fixed points `1/(2^k-3)`; <br>• a repunit residue-tree isomorphism; <br>• the decimal digit audit. | Correct. The contraction stops exactly at the integer boundary, which is Collatz (acknowledged). |
| ternary_digits | **PROVED**: <br>• `N_2..N_8` are prime and `N_9=17*19607843`; <br>• `17 | N_k` iff `k = 9 mod 16`; <br>• `U^2(N_k)<N_k` for `k>=4`; <br>• the rotation valuation carry; <br>• the circular-prime `k>=2` boundary; <br>• the first 12 primes sum to 197. <br>**CITED**: repunit certificates. | Correct; audit-quality digit work. It corrects THM-4473 legitimately. |
| ternary_triples | **PROVED**: <br>• the marked-triple bijection; <br>• the splitter; <br>• the `P_3` filtration; <br>• contraction `C'<C/3` (plus) and `C'<C/4` (minus); <br>• inverse guards (T10); <br>• the low-bit obstruction; <br>• the counter family `s_r=(2*8^(r-1)+3^r)/5`. | Correct. Exact reduction: plus convergence on `P_3` is equivalent to Collatz. |

### 1.7 duck_* (4 notes)

**Prompt answered (inferred).** Prime removals (`{11}`, `{2,3,11}`) and "196", `(10^k-7)/3` as a count, four colors / ten letters, a "three-color Zeckendorf decomposition", and tournament decoders without pair modules.

| note | claims (label) | reading |
|---|---|---|
| duck_decoder | **PROVED**: <br>• `(10^k-7)/3` counts the free C3-orbits on length-k words over `Sym^2_set(F2^2)` minus 7 constant words; <br>• the charge-orbit splits; <br>• the Fibonacci carry cocycle; <br>• the color odometer; <br>• **the rank-2 Hamiltonian-matching kernel is `K_(a,b)` with a,b odd**. <br>**REFUTED**: the hope "W nonsingular". **FINITE-EXACT**: ranks 2/4/6 at n=6 with counts 1680/17520/13568; K_{1,5}: 960; K_{3,3}: 720. | Numerology relative to Collatz, except the K_{3,3} kernel, which matters for the owner's theme. |
| duck_primes | **PROVED**: <br>• 196 is not a prefix sum of any of the three prime lists; <br>• the finite-filter square theorem; <br>• the divisor C3 action; <br>• `a_(k+1)=10a_k+21` via 21 new orbits. | Numerology; correct. |
| duck_tournament | **PROVED**: <br>• D1: the three blocks of `Q[H,H,H,1]` are exactly the size-q modules; <br>• D2: no cut switch creates a pair module; <br>• **D3: `|Pf|` is a switching invariant, equal to 1 on {TT, strong} and 3 on the two diamonds**; <br>• D4: triangle contraction with ports preserves reachability. | D3 links to THM-4472 (section 3). |
| duck_zeckendorf | **PROVED**: <br>• the charge `Q(n)` is invariant over distinct representations; <br>• the carry `delta` is in `{-1,0,1}` and satisfies the cocycle identity; <br>• the three colors R, G, B carry an order-3 action `M`; <br>• the color odometer. | Numerology; correct. |

### 1.8 thirtysix_* (4 notes)

**Prompt answered (inferred).** "36" (from `36=18+18` at level 11), `8*9/2=36`, the `Q_15` endpoints 8 and 9, half-steps, the ratios 4.2/8.4, the 17-phase, and multipartite graphs with a K4 quotient.

| note | claims (label) | reading |
|---|---|---|
| thirtysix_bridge | **PROVED**: <br>• seven distinct carriers of 36; <br>• half-step clocks; <br>• `F(n,sigma)` with `F^3=(8n+3sigma,-sigma)`; <br>• the six-edge trade between `K_(1,3,3,3)` and `K_(4,2,2,2)`; <br>• a mod-17 decimal/signed-Pell conjugacy (two 16-cycles on the 32 points of norm `+-1`); <br>• if a Pell trace `x_l` is prime then l is a power of 2. | Numerology except F (section 3). |
| thirtysix_digits | **PROVED**: <br>• **the signed germ map `(n,sigma) -> (y,k)` is a bijection, and F shifts k by 1**; <br>• `a_18` is prime (certificate); <br>• forced growth after the two-step descent; <br>• the cross-sheet half-height is non-integral for `y>1`. | F is the useful object. |
| thirtysix_multipartite | **PROVED**: <br>• the rank formula `rank A(K_(n1..nr))` = r (r even), r−1 (r odd); <br>• all twelve 36-edge complete multipartite types, exactly two of them planar; <br>• `C3[C3,C3,C3]`; <br>• exact doubling. | Graph theory; owner's theme only. |
| thirtysix_signed | **PROVED**: <br>• the square-root criterion (even cycles must pair); <br>• 36 roots of `2^2 5 18^2` and 216 of `2^2 5^2 18^2`; <br>• the global suspension `R^2=G`; <br>• minus cycle lengths by clock: odd-only (1,2,7), shortcut (1,3,11), ordinary (2,5,18). | Correct; the clock table matters for our cycle triple (section 4). |

### 1.9 creation_* (4 notes)

**Prompt answered (inferred).** The owner's giant "fruit" integers (with a typo), the ratios about 4.2 and 8.4, a "creation certificate", the Fano plane / E8 / Bott periodicity, and three carry colors.

| note | claims (label) | reading |
|---|---|---|
| creation_numbers | **PROVED**: <br>• literal vs repaired transcription (`d=10b+9`); <br>• the dominance bound `2+sqrt3 < a/(b+c) <= (7+sqrt65)/4`; <br>• the Nesbitt bound `>=3/2`. <br>**FINITE-EXACT**: trajectories; the 9G triple hits the three minus basins `(17,5,1)`. <br>**REFUTED**: "one coordinate per basin", since `13G+T` gives `(5,5,17)`. | Numerology (self-refuted pattern). |
| creation_elliptic | **PROVED** on `L=<G,T>`: <br>• the Kummer square class `alpha=(-1)^m 14^k`; <br>• rational halving by quadratic tests; <br>• a terminating creation decoder, `9G->4G->2G->G->O`; <br>• no positive fruit point is a real double; <br>• `|L/bL|=b*gcd(b,6)`; <br>• the actual Collatz lift `mG -> T(m)G`. | A binary-expansion (2-descent) decoder; the Collatz lift is a re-encoding (section 5, D1). |
| creation_fano | **PROVED**: <br>• F1–F2: octonionic `L_i` and the `Cl_(0,8)=M16(R)` generators preserve `E8+E8`; <br>• F3: lattice admissibility of `Psi(n)` is equivalent to endpoint integrality; <br>• F4: the swap `(1 5)` exchanges the E8 glue classes over D8; <br>• F5: no residue-indexed coercive polynomial potential decreases every L steps. <br>**CITED**: Bott, ABS. | F3/F4 are decoration of the affine action; F5 is a valid no-go. |
| creation_decoder | **PROVED**: <br>• C3–C7 (the Fano carry decoder); <br>• **C9: the top output bit of `q_(+,d)` has ANF degree `d-1`, and `q_-` has degree at most `d-2`**; <br>• **C10: no finite-state synchronous transducer computes `q_+`**; <br>• the dominance bound. | C9 is a frame artifact (section 2). C10 is a genuine no-go for the automata route to PC. |

### 1.10 zenodo_* (4 notes)

**Prompt answered (inferred).** The Zenodo DOI (He–Jing–Li symbol "sewing"), the chain 3, 4, 5 -> 9, 14, 7, 11, 25, square-sum paths, and the D3 representation.

| note | claims (label) | reading |
|---|---|---|
| zenodo_bridge | **PROVED**: <br>• the `9->14->7` block on triangles; <br>• the cylinder `8a+1 -> 6a+1`; <br>• the root family `(7*4^r-1)/3`; <br>• sewing S3a–S3b (rank 2, the seam lies in `Z[1/2] ∩ Z[1/3]=Z`); <br>• the free word-monoid decoder S4; <br>• the square-path map (S8); <br>• the four-reflection macro `+2` (S9). <br>**CITED**: the dataset. | Correct. The 8a+1 "family" is the classical stopping-time-2 class (section 2). |
| zenodo_sewing | **PROVED**: <br>• Z1: an integral lattice for the source's D3, whose reduction mod 2 is `GL2(F2)`, the Fibonacci color action; <br>• Z2: parity is lost; <br>• Z3: the invariant pairing is `diag(3,1)`. | Representation theory; no Collatz content. |
| zenodo_square | **PROVED**: <br>• a primitive triple gives a 6-vertex square-sum path with sums `c^2,(c+1)^2,c^2,(c-1)^2,c^2`; <br>• an extra chord exactly on negative-Pell legs; <br>• maximal repetition bounds (repeat versus escape). | Exact; no Collatz map (section 5, D2). |
| zenodo_triples | **PROVED**: <br>• the marked chart `Psi(n)` for all n; <br>• the spine jump; <br>• the `n=1 mod 4` depth contraction `j'<=3j/4`; <br>• the middle Berggren child `(3p-1,p)`; <br>• the elliptic macro `(3P+G)/4=6Q+G`; <br>• the Lorentz similarity `M^T J M=(9/4^k)J`; <br>• no fixed linear height. | Re-encodings of `n'<n`; correct. |

### 1.11 Repaired inherited notes (square-sum w6, THM-060, THM-4473, synthesis 20260917)

* **Square-sum Hamiltonicity (w6) and summand-closure filter.** RETRACTED: the unguarded rule "a degree-2 vertex has both edges in any Hamiltonian path", together with the certificates it produced at 20–22 and 24.
  * Minimal hostile: a triangle. In-lane hostile: the valid `Q_23` path ends at 22 and omits `22-14`.
  * The cut certificates at 18–22 and the leaf/nonendpoint forcing at 18–19 survive.
  * Nonexistence at 24 is re-proved endpoint-aware.
  * -> **Correct and important.** I checked every adjacency list in the `Q_24` proof and the forced 11-cycle `1-8-17-19-6-10-15-21-4-12-24-1` (all sums square; the six degree-2 vertices force it).
* **THM-060 Type A.** A zero-backbone triple contributes 0 (transitive) or 2 (cyclic) under full reversal, not always 2; the parity argument survives. -> Correct.
* **THM-4473.** "Repunit primes are exactly the prime fixed points" now requires `k>=2`, since 2, 3, 5, 7 are fixed. -> Correct.
  * Cosmetic defects introduced: the missing spaces in "base10", "through1100" and "on2026-09-25".
  * The body bullet "Every other prime with k digits has a free orbit of exact size k" still lacks `k>=2`. It is harmless (orbit size 1 = k).
* **Repunit note §4.6.** "Circular primality is rotation-invariant; ordinary primality of one member is not (19 prime, 91=7*13)." -> Correct wording repair.

## 2. Audit

### 2.1 Independent spot checks (own code; all in `04-computation/experiments/procgen_incoming_20260925_audit_*.py`)

| # | claim (source) | my method | result |
|---|---|---|---|
| A1 | `q_(+,3)=[0,5,2,3,4,1,6,7]`, `q_(-,3)=[0,7,6,1,4,3,2,5]`; `q_-` linear; polar defect `4(x0y1+x1y0)`; C9 degrees; C9 signs (creation_decoder, creation_fano) | own parity-vector map, Möbius-transform ANF, `d<=13` (`audit_01`) | **agree**. Top-bit degree: plus `d-1`, minus `d-2` for `3<=d<=13`; `sign(q+)=-1`, `sign(q-)=+1`. **Also:** `q_-(n)=q_+(-n mod 2^d)` for all `d<=13`, and negation's own top bit has degree `d-1` (see F1). |
| A2 | consecutive plus edge-triangles are never Berggren-comparable; `(27,5)=B1^2(7,5)` (ternary_berggren §5) | own ancestor test with accelerated B1-runs, validated against an explicit depth-7 tree (3280 nodes) (`audit_02`) | **agree**. Odd `x<=100001`: 49,992 pairs; kinds rise-rise 12,500, rise-fall 12,500, fall-rise 12,496, fall-fall 12,496 (codex's exact counts); 0 comparable. Extended to `x<=2,000,001`: 999,990 pairs, 0 comparable. Minus: 20 comparable to `2*10^6`, **all of them entry edges into the 2-cycle `{5,7}`**. |
| A3 | closure counts 252/274/541/640; 483 certified via 181; 241 loses its canonical certificate; first uncertified 7; inverse ports add nothing; residual classes mod 48/144 (creative_sibling/decoder) | own inductive closure (`audit_03`) | **agree** on every number, including 99 new sources (first 483), forward closure in range, and 19/24 = 57/72. |
| A4 | Hamiltonian-matching kernel census at n=6 (duck_decoder §4) | own permutation-completion construction, F2 rank (`audit_04`) | **agree**: `{2:1680, 4:17520, 6:13568}`; rank-2 kernels K_{1,5} (960) and K_{3,3} (720); all kernel degrees odd; orders 2 and 4 full rank; hostile has `H=9`, rank 4, null vectors `e0+e4`, `e2+e5`. |
| A5 | G_N planar iff `N<=15`; K_{3,3} paths at 16; complement of `L(K5)` is Petersen; 144 = 24 + 120 matchings; 651 rank-2 graphs on 6 vertices (6 + 10 with odd degrees) (decoder_minors, seam_threshold, thirtysix_multipartite) | networkx 3.4.2 LR test (codex used 3.5) on my own graph construction; direct path checks; own matching recursion (`audit_05`) | **agree** on all. New observation F2. |
| A6 | `Q_N` Hamiltonian-path counts `N=14..25` = `0,1,1,1,0,0,0,0,0,3,0,10`; every `Q_25` path uses `11-25-24` (decoder_prime_square) | own DFS with reachability prune, `N=2..27` (`audit_06`) | **agree** (and `N=26,27`: 12, 35 = A090460). |
| A7 | `U_(-sigma)(2n+sigma)=U_sigma(n)`, exponent +1, `F^2=S_sigma`; germ bijection; Jacobsthal root fibre; forward-conjugacy hostile (thirtysix_digits) | direct, odd `n<=2*10^5`, both signs (`audit_07` A) | **agree**. Also `T_-(2n+1)=2T_+(n)` and `T_+(2n-1)=2T_-(n)` for odd n. |
| A8 | `a_k=(10^k-7)/3`: primes for k=2..8, the listed divisors at 9..17, `a_18` prime with `n-1=2*3*5*2071723*5363222357`; `17 | a_k` iff `k = 9 mod 16` (ternary_digits, thirtysix_digits) | deterministic Miller–Rabin (13 bases) plus my own Lucas/Pocklington check (`audit_07` B) | **agree**. |
| A9 | four-core Pfaffian classes (duck_tournament D3) | all 64 labelled tournaments, all 16 switchings (`audit_07` C) | **agree**. Also `|Pf| = H (mod 4)` on all 64 (F3). |
| A10 | `C'<C/3` (plus), `C'<C/4` (minus), fixed points `1/(2^k-3)`, exact r-step entry (ternary_triples/bridge) | exact fractions over 16,320 plus and 14,958 minus cases (`audit_07` D) | **agree**: max `C'/C` = 0.2946 (plus), 0.2494 (minus). |
| A11 | exchange laws; cylinder identity; `-3/5`, `-4/15`, `266/243` (decoder_mahler_catalan) | all 4,097 adjacent swaps with `L<=10`; cylinders `L<=8` (`audit_07` E) | **agree**. Also q-independence (F4). |
| A12 | Mills cubic prefix `2, 11, 1361, 2521008887` and bracket (seam_mills) | own next-prime search and exact rational 81st powers (`audit_07` F) | **agree**. |
| A13 | `F=S+U` iff `p, p^3, p^2qr` (seam_prime_balance) | `N<=20000`; closed formulas against direct divisor lists `N<3000` (`audit_07` G) | **agree**. |
| A14 | residual C3-orbit counts `(10^k-7)/3` (duck_decoder/primes) | brute force over words `k<=4` (`audit_07` H) | **agree** (1, 31, 331, 3331). |
| A15 | 36 roots (`2^2 5 18^2`), 216 roots (`2^2 5^2 18^2`), 2 (odd-only), 1 (shortcut) (thirtysix_signed) | square-root count by cycle type, formula validated on all permutations with `n<=7` (`audit_07` I) | **agree**. |
| A16 | depth 14: 142 accepted ports, 734 residual classes (creative_descent) | own first-acceptance DP (`audit_07` J) | **agree**. |
| A17 | 9 -> 4 sewing certificate `M_u=(3,1,4)`, `M_v=(243,347,512)`, seam 7 = 7; residual rank 2 (zenodo_bridge) | own matrices; rank over Q on 24 cases (`audit_07` K) | **agree**. |
| A18 | pair repair `{3:24, 6:24, 7:12, 8:4}` and canonical minima `(15,13,12,7,10)` (decoder_pair_repair) | own enumeration of 945 matchings × 24 regular quotients × 64 cores (`audit_08`, `audit_11`) | **agree**. |
| A19 | `Q_24` endpoint-aware proof (decoder_prime_square §4) | by hand from the displayed adjacency lists; script check of the forced cycle | **agree**: all four endpoint cases saturate 5 and force the 11-cycle. |
| A20 | THM-4469 census `(10,7):4, (16,11):5, (19,13):1, (20,14):8` (our own, as a by-product) | own enumeration (`audit_12`) | **agree**. Extended to `L<=24` in section 5 (seed 1). |

### 2.2 New observations made during the audit

* **F1 (PROVED, one line; label clarification for C9).** The parity-address maps satisfy `q_- = q_+ o neg` on every `Z/2^d`. This is the residue form of our transport theorem. Negation `x -> -x mod 2^d` is itself a permutation whose top output bit has ANF degree `d-1` (A1).
  * Algebraic degree is not invariant under composition with a non-affine permutation. So the plus/minus degree gap (`d-1` versus `d-2`) and the depth-3 "nonlinear plus / linear minus" dichotomy are properties of the pair *(map, binary digit frame)*, not of the side of 0.
  * In the complemented frame the roles swap.
  * The statistic is side-blind: it holds equally for `3n+1` on negative integers, where the cycles −5 and −17 live. So it cannot serve as a sign-aware input.
  * Codex states `q_-(n)=q_+(-n)` inside the C9 proof and says the difference "does not separate convergent from nonconvergent dynamics". But MISTAKES ("Three-bit plus parity is nonlinear") and the concept board ("plus parity is nonlinear, minus is linear at depth three") read as if it were a sheet invariant. **Label drift; flag, do not promote.**
* **F2 (PROVED by two finite certificates, like codex's G_N theorem).** The square-sum graph `Q_N` is **planar iff `N<=24`**.
  * `Q_24` (24 vertices, 30 edges) has a planar embedding with `F=8`.
  * `Q_25` contains a K_{3,3} subdivision with shores `{3,5,12}` and `{4,11,13}`. Its nine paths are `3-13`, `5-4`, `5-11`, `12-4`, `12-13`, `3-6-10-15-21-4`, `3-22-14-11`, `5-20-16-9-7-2-23-13` and `12-24-25-11`; I checked that the interiors are disjoint and every edge sum is a square.
  * Since `Q_24` is planar, **every** Kuratowski subgraph of `Q_25` uses vertex 25, hence the ear `11-25-24`. That ear is exactly the one every Hamiltonian path of `Q_25` must use (A6).
  * So in the square-sum family the first size after the exceptional range (24) gains traceability and loses planarity through **one and the same degree-2 ear**. This is sharper than codex's cross-family analogy (seam_threshold §1: an ear repairs Q in one family and breaks G in another).
  * Monotonicity (`Q_N` is an induced subgraph of `Q_(N+1)`) gives the all-N statement. For the Kuratowski/Tutte lanes this is raw data, not a theorem about Hamiltonicity.
* **F3 (FINITE-EXACT, exhaustive).** On the 64 four-vertex tournaments, `|Pf| = H (mod 4)`, and `|Pf|=3` iff `c3` is odd (TT: `|Pf|=1, H=1`; strong: `|Pf|=1, H=5`; diamonds: `|Pf|=3, H=3`).
  * On all 596 AM-fair quadruples of THM-4472 (both sheets, `5<=|s_o|<=301`), `|Pf| = H = 2+sgn(b s_o)`.
  * So codex's switching invariant is exactly THM-4472's sign indicator: diamonds cannot be switched to TT.
* **F4 (PROVED, one line).** The exchange law does not depend on the multiplier. For `R_q(w)=sum_i q^(a-i) 2^(j_i)` (ones at positions `j_1<...<j_a`) one has `R_q(u01v)-R_q(u10v)=2^j q^b`. Hence an adjacent-carry pair of THM-4469 (`R_B'=R_B+1`) is exactly a {2,3}-unit equation `sum_i 3^(a-i)(2^(j'_i)-2^(j_i)) = 1`.
  * For the smallest instance: `-3^6 - 2*3^5 + 2^6*3^2 + 2^7*3 + 2^8 = 1`.
  * Seed 1 builds on this.

### 2.3 Flags: overclaims, label drift, circularities, conflicts

1. **C9 label drift (F1).** Not a conflict with the transport theorem as long as it is read correctly. It would become one if promoted as "the first sheet-specific carry".
2. **Triviality under a PROVED label (zenodo_bridge §2, zenodo_triples §3/§5).** The "guarded descent family" `8a+1 -> 12a+2 -> 6a+1` is the special case of the classical fact `T^2(n)=(3n+1)/4<n` for `n = 1 mod 4` (stopping time 2, Terras/Everett).
   * The accompanying Berggren-depth, hypotenuse and canonical-height decreases are re-encodings of `n'<n`, since each height is quadratic in n.
   * "closes a real infinite portion of that task" and "repairs the local difficulty at 9G" overstate its significance.
   * Continuing "residue family by residue family" is DEFECT-blind by our foundry (THM-4470 §4, exceptional dimension 0.95). Codex elsewhere acknowledges that "local residue coverage is not root coverage".
3. **Decoration with no dynamical content.** Two cases, both correct and both acknowledged by codex as not giving descent:
   * **The E8/Clifford carrier.** `Psi(n)=(n rho, rho)` with `A(a,b,d)` is the affine action on `(n,1)`; "faithful lattice model" is literally true but adds no information.
   * **The elliptic Collatz lift.** On pure multiples `mG -> T(m)G` is a relabelling, and canonical height is `m^2 h(G)`.
4. **Circular by design (acknowledged).** Tournament halving encodes `v2(n)` into `E(n)`. The adaptive compiler learns its seeds by orbit search. Counts that "match" owner numbers (`(10^k-7)/3`, the 36s) are exact constructions made to hit them. No *hidden* circularity was found.
5. **No DRIFT controls anywhere.** Every codex lane uses the minus cycles 5 and 17 (SHEET) and the all-ones prefixes / −1 (lookahead) as hostiles. None runs 5n+1 or qn+1. By foundry rules no codex mechanism is DRIFT-tested. The one DRIFT test I ran (seed 5) changed the verdict.
6. **Minor wording (C6).** "q exchanges the other four [Fano lines]" means that the four lines missing 4 are sent onto the four **non-lines** `{127,136,235,567}` (the lines of the transported `*`-plane). Minimal witness: `q({1,2,3})={2,3,5}` and `2 XOR 3 = 1 != 5`. Fano incidence is preserved only on the pencil through 4.
7. **Conflicts with THM-4469..4473 or the no-go results.** None found. Codex:
   * scopes its Mahler transfer as not contradicting THM-4469;
   * refuses to identify its four-core with THM-4472;
   * treats density and coverage fractions as non-proofs (consistent with THM-4470's DEFECT result);
   * makes no contraction claim in `|x-y|` (consistent with THM-4471 (D)).

   Its F5/D9/finite-bank obstructions are further instances of our lookahead no-go (THM-4470 §5, the −1 obstruction).

## 3. Connections to our objects

| our object | codex families | relation | details |
|---|---|---|---|
| **Sign law and transport theorem** (mod-192 note Thm 6 / Cor 7: residue data are side-blind; a proof needs the sign law, then gate integrality, then null-set avoidance) | collatz_bugs §6; thirtysix_signed; creation C3–C9; ternary_berggren §5; ternary_triples; decoder_mahler_catalan §4 | **duplicates** (`T_-(-n)=-T_+(n)`); **refines** (clock-explicit gluings; `q_-=q_+ o neg` is the residue-level transport); **consistent** (C9 is a frame artifact, F1) | ternary_berggren's comparability law is the only new sign-aware law. Its minus failures are exactly the `{5,7}` entries (A2), so its sign content is gate integrality of the word `(1,2)` (`2^3-3^2=-1`). The Mahler shared-word map sends the whole plus root basin to negative rationals: a sign transport between conjugate maps. |
| **Pairing ladder THM-4470** | collatz_bugs §4/§6; thirtysix_digits §3; creative_* no-gos | **strengthens** the "graceful: ANALOGY" verdict (spanning-star obstruction, sign gauge); **new, consistent** (F interleaves the sheets' inverse fibres; `T_-(2n+1)=2T_+(n)`); **duplicates in spirit** (F5, D9, the finite-bank obstruction and forward closure are all the −1 obstruction of THM-4470 §5) | THM-4470 says the sheets are the two consecutive pairings. F says their inverse fibres over each target are alternate levels of one ray (the Jacobsthal trunk for `y=1`). |
| **Mahler bridge THM-4469** | decoder_mahler_catalan; seam_mills | **refines** our §5 remark "Mahler's own 3/2 problem does not transfer" with an explicit obstruction (root 4 goes to `-4/15`; every root word has a suffix violating the tail bound) | The exchange law is q-independent (F4). Adjacent carries are {2,3}-unit equations; this is seed 1. |
| **E-SCC (Q1, Q2)** | ternary_triples T9–T10; thirtysix_digits §3; ternary_berggren §3–4 | **duplicates** | Rational inverse edges and the level rise duplicate the Q2 guard structure (our mod-27 lemma and kappa formula are finer). The joint signed fibre is the additive sign-choice relaxation (deck #4: empty exceptional set). The ternary odometer of fibre heights is the inherited triadic odometer. |
| **HARD class / Periodicity Conjecture** | creation_decoder C9–C10; creative_descent (residual language) | **new no-go** (C10) | The parity-vector map (THM-4473's Q, the inverse Bernstein–Lagarias conjugacy) is not finite-state. Finite-state maps preserve ultimate periodicity, so an automaton proof of PC is impossible; section 5 extends this to counter machines by counting. `W_14` (734 classes) is the finite-level exceptional set (duplicate). There is no contact with Theorems S/D/Y or HYP-9127. |
| **Gates `2^p-3^a`** | decoder_mahler_catalan §5; zenodo sewing; ternary_triples T12; creative_descent §5 | **duplicates / reformulates** | The −139 gate with carry `2363=17*139`, and `1/(2^k-3)` as the Banach points of the words `1 0^(k-1)`, are inherited. The sewing seam `Z[1/2] ∩ Z[1/3]=Z` contains cycle-gate integrality as the case `r=n`. |
| **Banach fixed-point chain** (THM-4471 §4) | ternary_triples; creation_elliptic; thirtysix_digits; thirtysix_signed | **analogies** | Terminating creation ranks (`m -> floor(m/2)`) and `F^(-1)` (k decreasing) are Banach-type well-founded descents on encodings. Permutation square roots (even cycles must pair) sit beside THM-4472's free involution on periodic points. |
| **Square-sum / graceful verdicts** (square-sum:brackets REAL; square-sum:Collatz NUMEROLOGY; graceful ANALOGY) | square-sum repairs; decoder_prime_square; collatz_bugs §6; zenodo_square | **corrects inherited material, supports all three verdicts** | The degree-2 forcing retraction concerns the 2026-09-22 mod6 note, not our session. `{2,3,11}` is re-proved (duplicate). The two-sheet halving law and the S8/S9 square paths are exact with no Collatz map. F2 (planarity) is new. |
| **Four-vertex tournaments** THM-4472 | duck_tournament D3; collatz_tournament_core; decoder_pair_repair | **refines** | `|Pf|` separates {TT, strong} from {diamonds}. On AM-fair quadruples, `|Pf| = H = 2+sgn(b s_o)` (F3). A 6-vertex extension is local (section 5). The four-core is correctly kept distinct. |
| **Repunits and digit chains** THM-4473 | ternary_digits; thirtysix_digits; creative_transducer; duck_primes | **corrects** (the `k>=2` boundary); **refines** (the decimal family shadows `-7/3 -> -3 -> -1` like repunits shadow their 2-adic limits); **duplicates** (the run length `v2(n+1)`; `T^k(2^k-1)=3^k-1`) | `(10^k-7)/3`: seven primes, then 17 at the phase `k = 9 mod 16`, then the prime `a_18`. This is numerology relative to Collatz. The two-step descent for `k>=4` is a genuine but local carry law. |

## 4. Triples inventory (raw, for the Kuratowski/Tutte lanes)

**4a. Our session's triples (as given).**
* sheets `b` in `{-1,0,+1}`;
* the controls SHEET / DRIFT / DEFECT;
* the necessary proof inputs: sign law, gate integrality, null-set avoidance;
* the places `{inf, 2, 3}`;
* the swap words Sturmian / square / cube;
* the negative cycles `{-1,-5,-17}` with lengths 1, 3, 11 (shortcut clock);
* GM, AM, QM with `QM^2+GM^2=2AM^2`;
* Banach / Brouwer / Sharkovskii;
* the 4-tournaments TT / diamonds / strong.

**4b. Codex-side triples, with exact definitions.** In the "relation" column, *map* means an explicit structure-preserving map exists, *analogy* means a shared shape with no map, and *numerology* means equal numbers with no map.

| # | codex triple | exact definition | source | nearest session triple | relation |
|---|---|---|---|---|---|
| T1 | **K5, K_{3,3}, Petersen** | K5 arises as the contraction of `L(K5)` along the pairing `P_v` of the two out-arcs of v in a regular 5-tournament (24 decoders among 144 matchings), and as `Sym^2_set(F2^2)` ↔ `E(K5)` (C3 fixes 0 and ∞). K_{3,3} arises as the subdivision in G_16 (shores `{5,8,12}`/`{6,10,14}`), as the rank-2 odd-degree Hamiltonian-matching kernel at n=6 (720 tournaments), and as the witness in ten of the twelve 36-edge multipartite types. Petersen is the complement of `L(K5)` (`45=30+15`; THM-261). | seam_threshold, duck_decoder, decoder_minors, thirtysix_multipartite | 4-tournaments / controls | raw data. Plus F2 (mine): `Q_N` is planar iff `N<=24`, with a K_{3,3} (shores `{3,5,12}`/`{4,11,13}`) through the forced ear `11-25-24`. |
| T2 | three nonzero colors R, G, B (plus neutral 0) | `R=(1,0), G=(0,1), B=(1,1)` in `F2^2`, `R+G=B`, cycled by `M(a,b)=(b,a+b)` (order 3 mod 2). D3 = `GL2(F2)` acts faithfully. | duck_zeckendorf, zenodo_sewing | sheets `{-1,0,+1}` | analogy |
| T3 | Fibonacci carry `delta` in `{-1,0,+1}` | `delta(a,b)=b(a)+b(b)-b(a+b)`, `b(n)=floor((n+1)/phi^2)`, with the cocycle identity. Realized by `(1,1),(1,2),(2,2)`. | duck_zeckendorf | sheets `{-1,0,+1}` | analogy |
| T4 | Berggren branches B1, B2, B3 | `(s+2t,t)`, `(2s+t,s)`, `(2s-t,s)`, with determinants `1,-1,1`. On ratios `r=s/t` they act as `r+2`, `2+1/r`, `2-1/r`. | ternary_berggren | swap words (three growth modes) | analogy |
| T5 | three parent regimes | `s>3t`, `2t<s<3t`, `t<s<2t` (the odd Gauss map) | ternary_berggren | — | map (tree inverse) |
| T6 | three dyadic regions | odd chain `v2=0`, seam `v2=1`, sea `v2>=2`. `+2` fixes the odd chain and alternates seam and sea. The pair-module graph of `E(n)` is isolated points, edges, or even paths accordingly. | decoder_halving/minors, seam_mills | places `{inf,2,3}` (the 2-adic place split three ways) | analogy |
| T7 | minus cycles in three clocks | cycle lengths: odd-only `(1,2,7)`, shortcut `(1,3,11)`, ordinary `(2,5,18)`; locus sizes 10/15/25; minimal square-root completions 12/15/45 with 2/1/36 roots | thirtysix_signed | negative cycles `{-1,-5,-17}`, lengths 1, 3, 11 | **map** (the same cycles; codex shows the length triple is clock-dependent, and ours is the shortcut clock) |
| T8 | carry states `{0,1,2}` | the multiplier `M_c` computes `3x+c`; plus uses `c=2`, minus `c=1`; flushes `[0]#->#`, `[1]#->1#`, `[2]#->01#` | creative_transducer | sheets | partial map: `c=(3+b)/2` for `b=+-1`; `b=0` would need `c=3/2` |
| T9 | three blocks and a singleton `Q[H,H,H,1]` | three copies of a regular odd H plus r (order `3A+1`) | tournament_core, decoder_halving, duck_tournament | 4-tournaments | numerology (explicitly not THM-4472) |
| T10 | Fano pencil through 4 | `q_(+,3)` fixes the lines `{1,4,5}, {2,4,6}, {3,4,7}` and sends the other four lines to non-lines | creation_decoder | — | map (on `Z/8`) |
| T11 | E8, D8, E8' | E8 = Construction A of RM(1,3); E8' is its image under the swap `(1 5)`; `E8 ∩ E8' = D8`, with 112 common roots and 128 others each | creation_fano | — | map (lattice) |
| T12 | three "eights" | Fano/cube `F2^3` labels; Bott period 8 (`Cl_(0,k+8)=Cl_(0,k)⊗M16(R)`); `F^3(n,s)=(8n+3s,-s)` | creation_fano, thirtysix_digits | — | numerology |
| T13 | three half-steps | `F(n,s)=(2n+s,-s)` (inverse fibres); the suspension `R` with `R^2=G` (clock); permutation square roots (cycles) | thirtysix_* | Banach/Brouwer/Sharkovskii (period structure) | analogy |
| T14 | primitive triple `(A,B,C)` with signed B | `Phi(s/t)`. B<0, B=0 (terminal `(1,0,1)`) and B>0 mark x<1, x=1, x>1. Exactly one leg is divisible by 3. | ternary_triples | sheets `{-1,0,+1}` | analogy (orientation sign) |
| T15 | integer-parent residue cases | for an integer target u: `u = 1 mod 3` needs k even; `u = 2 mod 3` needs k odd; `3 | u` has no integer parent (plus; swapped for minus) | ternary_triples T10 | gates / E-SCC Q2 | map (duplicate of the Q2 guard) |
| T16 | primes `{2,3,11}` | the primes with a nontrivial multiple in their own odd-square bracket (inherited); also `11=2^3+3` (parent, correction, child) | decoder_prime_square, seam_mills | our brackets result | duplicate / numerology |
| T17 | balance values `(F,S,U)` | `(0,0,0)`, `(2,1,1)`, `(10,7,3)` for `p`, `p^3`, `p^2qr` | seam_prime_balance | — | numerology |
| T18 | mean-type triple | `A+B=2c`, `AB=c^2-h^2`; Nesbitt `sum a/(b+c) >= 3/2`; the fruit dominance window `(2+sqrt3, (7+sqrt65)/4]` | seam_power_coordinates, creation_numbers | GM, AM, QM (`QM^2+GM^2=2AM^2`) | analogy |
| T19 | Mills exponents 2 / 3 / >=4 | `xi_2` (existence of selected branches, Matomäki), `xi_3` (irrational; transcendental or Pisot), `xi_c` for `c>=4` (transcendental; Saito) | seam_mills | swap words square / cube | analogy: what is provable changes between square and cube levels in both |
| T20 | two-register sewing places | forward value in `Z[1/2]`, backward value in `Z[1/3]`, seam in Z | zenodo_bridge | places `{inf,2,3}` | analogy (product-formula shape) |
| T21 | three content levels | `gcd(3s+t,t)` in `{1,3}`; `P_r` strata; `C'<C/3` | ternary_bridge | places (3-adic) | analogy |
| T22 | 3-4-5 / 7-24-25 / 8-15-17 | the root triangle; its Gaussian square; the Pell triangle of `T_8=36` (`17^2-8*6^2=1`) and the marked `(15,-8,17)` of the edge `3->5` | zenodo_*, thirtysix_bridge | — | numerology |
| T23 | switching classes | `|Pf|=1` (TT, strong) versus `|Pf|=3` (the two diamonds) | duck_tournament D3 | TT / diamonds / strong | **map** (F3: `|Pf|=H mod 4`; on quadruples `|Pf|=2+sgn(b s_o)`) |
| T24 | `36=` three multipartite graphs on 10 letters | `K_(1,3,3,3)`, `K_(4,2,2,2)`, and `K9=K_(3,3,3)+3` triangles; `C3[C3,C3,C3]` | thirtysix_multipartite | — | numerology |

## 5. Proof-angle seeds (ranked)

The ranking weighs three things: contact with the divergence half or a registered hypothesis, survival of the controls, and cheapness of the falsification test. No seed attacks the foundry's single missing mechanism: a non-integrality or irrationality statement reaching positive-entropy words.

1. **S-unit anatomy of the exact Mahler instances (HYP-9134 lever).** *Source:* the decoder_mahler_catalan exchange law made q-independent (F4), together with THM-4469.
   * *Statement.* Every adjacent-carry pair with `3^a>2^(L+1)` is a {2,3}-unit equation `sum_i 3^(a-i)(2^(j'_i)-2^(j_i))=1` with m nonzero terms.
   * By the Evertse–Schlickewei–Schmidt finiteness theorem for non-degenerate solutions (CITED; not re-derived), each m admits only finitely many pairs.
   * The HARD slice therefore carries a complexity ladder in m, and the Mahler exclusion `Z(p/q,I)` becomes a statement about an explicit S-unit identity.
   * *Pre-screen* (`audit_12`):
     * the census `L<=24` reproduces THM-4469's rows;
     * it adds `(21,14):75`, `(21,15):6`, `(22,15):54`, `(23,16):126`, `(23,17):10`, `(24,16):198` and `(24,17):20`;
     * term counts are at least 5 (`m=5` only at `(10,7)`), the minimum has fallen to 7 at `L=23,24`, and no vanishing subsums occur.
   * *Falsification test.*
     * Extend to `L<=30`. A recurring parametric family with bounded m would contradict the non-degenerate reading.
     * Prove or refute "the minimal-m instance at each L has m bounded by f(L)".
   * *Controls.*
     * SHEET: minus carries are −R, so the minus instances need `xi<0`; the argument must use `xi>0` (the sign law).
     * DRIFT: rerun with `q=5`. A Mahler-side exclusion that also covers 5n+1 slices must be checked against actual long `{B,B'}` segments of 5n+1 orbits.
     * DEFECT: explicit words, so not blind.
2. **Sheet-switch density ladder via the signed half-step (diagnostic, E-SCC/choice ladder).** *Source:* the thirtysix_digits bijection `(n,sigma) <-> (y,k)` and `F=(2n+sigma,-sigma)`.
   * *Statement.* Profile the exceptional set of the relaxation in which the orbit may switch sheet on at most a fraction c of odd steps, or only on a residue set S.
   * c unconstrained gives empty (deck #4); c=0 gives dimension 0.95.
   * *Falsification test.* Profile at `2^20..2^26` as in the choice ladder. If every c>0 already empties the set, there is no gradation and the seed dies.
   * *Controls.* DRIFT: the same profile for 5n±1 must differ, since choice ladders are known to be drift-blind (synthesis §4 correction). SHEET is symmetric by design.
3. **Bruhat monotonicity of carries at cycle gates (cycle half).** *Source:* the exchange law.
   * *Statement.* R is strictly increasing along `01->10` moves, with extremes `3^a-2^a` and `2^(L-a)(3^a-2^a)`. So for fixed `(L,a)` the cycle candidates `n=R(w)/(2^L-3^a)` are ordered along the Bruhat order of fixed-weight words. Integrality at the gate plus monotonicity may localize cycles to Bruhat intervals.
   * *Falsification test.* It must recover the minus cycles (`110`, gate −1; `11110111000`, gate −139) and the 5n+1 cycles as the only integral points at the tested clocks, and it must beat brute force in the size of the search. Otherwise it is classical (Böhm–Sontacchi, Eliahou) and dead.
   * *Controls.* SHEET: the sign of the gate must be used. DRIFT: the 5n+1 cycles (13, 17 families) must be found.
4. **Rational 3-content height across the integer boundary (places `{inf,2,3}`).** *Source:* ternary_triples/bridge (`C'<C/3` while `3 | t`; `P_r` enters Z in r steps).
   * *Statement.* Look for a height `h = log C + lambda*log|den|_3 + mu*(2-adic term)` on positive odd-denominator rationals that decreases along plus steps off the boundary *and* extends across it.
   * *Falsification test.* Any such h restricted to the integers is a function of n alone, so F5 and D9 kill every version without a 2-adic sidecar. The test is whether an extra 2-adic term (for example `v2(n+1)`, the run length) can be added without violating F5.
   * *Controls.* SHEET: the minus contraction `C'<C/4` also holds, so h must see the sign. DRIFT: `5x+1` on `P_r` (denominator `5^r`?) must fail.
5. **Berggren ancestor law (pre-screened; demoted).** *Source:* ternary_berggren §5.
   * *Statement.* Along any plus orbit, edge triangles off the terminal spine (edges into 1) are pairwise Berggren-incomparable.
   * *Pre-screen* (`audit_09*`):
     * plus: 0 comparable pairs among all later edges, for every odd start up to 200,001;
     * minus: 70 comparable pairs (to 200,001), **all on cycle edges**;
     * 5n+1 (to 100,001): 18 comparable pairs, of which 15 involve cycle edges and 3 form one sporadic family, `35->11->7->9` with its siblings 563 and 9011;
     * **7n+1 (divergent orbits): 0 comparable pairs**.
   * *Verdict.* DRIFT-blind for divergence, since infinite antichains are cheap and divergent 7n+1 orbits satisfy the law. For cycles it is equivalent to "no repeated edge". What remains is a sign-aware local law whose sign content is small gates. Useful only as a generator of exact local laws: extend codex's case analysis to distance 2 and check whether it closes.
6. **Seam-height of two-register certificates ("sound certificate search").** *Source:* zenodo sewing (rank-2 residual, seam in Z).
   * *Statement.* Measure `H(n)` = the least seam height `max(|A_u n+B_u|, |D_v 4-B_v|)` over successful splits. A bound `H(n) <= n^C` for all n would be a Collatz-strength quantitative target.
   * *Falsification test.* Compute H for `n<=10^6` and fit the growth.
   * *Controls.* DEFECT: sound verifiers pass automatically, but any bound proof must be pointwise. DRIFT: H is infinite for divergent 5n+1 orbits, so the bound must fail there. SHEET: H is infinite on the cycles 5 and 17 with target 4.
7. **Pfaffian invariants of AM-fair orbit tournaments (pre-screened negative).** *Source:* duck_tournament D3 with THM-4472.
   * *Pre-screen* (`audit_10`, 7,988 cases).
     * At 4 vertices, `|Pf|=H=2+sgn(b s_o)`.
     * At 6 vertices (two steps), `(|Pf|, H mod 4)` is a function of `(b, sgn(b s_o), parity of T(s_o))`: local.
     * It is side-aware (it distinguishes minus-positive from plus-negative through the order arcs) but carries no global content.
   * *Falsification test.* 8- and 10-vertex versions; if still local, it is dead.
   * *Controls.* DRIFT is built in: AM-fair pairs exist only for `q=3` (THM-4470).

**Dead on arrival (recorded so that nobody pursues them).**
* **D1. Elliptic creation-rank transfer** (creation_decoder proposal 3). On `{mG}`, `h(mG)=m^2 h(G)`, so every group-embedding height is a function of `|m|` and carries no information beyond integer size.
* **D2. Square-reflection dynamics as a Collatz carrier** (zenodo_square). Compositions of `x -> s^2-x` are `x -> +-x + const`: there are no dilations, so no word encodes `x -> 3x+1` or `x -> x/2`.
* **D3. Fano/E8 frame changes** (creation proposal 2). The plus/minus difference is a digit-frame artifact (F1), and the carrier is the affine action on `(n,1)`.
* **D4. Automata route to PC beyond finite state** (the C10 extension; FINITE-EXACT counts, asymptotic sketched and not proved). A prefix p of length k gives the residual map `m -> q_+(3^(a(p)) m + T^k(p))`.
  * A deterministic synchronous transducer needs one configuration per distinct residual map at time k. `audit_13` counts the distinct residual maps on the next k bits, for `k=1..8`:
    * plus: `1,3,5,11,21,37,69,127`;
    * minus: identical from `k=5` on;
    * 5n+1: 229 at `k=8`.
  * The counts grow like `2^(k-1)`. A one-counter machine has at most `|Q|(k+1)` configurations at time k, so this growth excludes one-counter and similar polynomial-configuration machines.
  * Moreover `neg` is finite-state, so `q_-` lies in the same class: the route is SHEET-blind.

## 6. Reproduction

All audit scripts are standard library, except networkx 3.4.2 in `audit_05*`. Each runs in seconds with peak RSS at most 35 MB:

```text
cd 04-computation/experiments   # scripts prefixed procgen_incoming_20260925_
python3 procgen_incoming_20260925_audit_01_parity_degree.py            # A1, F1
python3 procgen_incoming_20260925_audit_02_berggren.py 2000001         # A2
python3 procgen_incoming_20260925_audit_03_sibling_closure.py 100000   # A3
python3 procgen_incoming_20260925_audit_04_matching_kernel.py          # A4
python3 procgen_incoming_20260925_audit_05_planarity_petersen.py       # A5
python3 procgen_incoming_20260925_audit_05b_q25_kuratowski.py; python3 procgen_incoming_20260925_audit_05c_paths.py   # F2
python3 procgen_incoming_20260925_audit_06_squaresum_paths.py 27       # A6
python3 procgen_incoming_20260925_audit_07_small_checks.py             # A7-A17, F3 (first part)
python3 procgen_incoming_20260925_audit_08_pf_quadruples_and_repair.py # F3 (quadruples), A18
python3 procgen_incoming_20260925_audit_11_canonical_repair.py         # A18 (canonical)
python3 procgen_incoming_20260925_audit_12_adjacent_carry_sunit.py 24  # A20, seed 1
python3 procgen_incoming_20260925_audit_09_seed_berggren_antichain.py 20001 400; python3 procgen_incoming_20260925_audit_09b_exclude_root.py
python3 procgen_incoming_20260925_audit_09c_classify.py 200001; python3 procgen_incoming_20260925_audit_09d_controls_large.py   # seed 5
python3 procgen_incoming_20260925_audit_10_seed_pfaffian6.py           # seed 7
python3 procgen_incoming_20260925_audit_13_residual_complexity.py      # D4
```

The scripts were moved from the lane scratch directory into `04-computation/experiments/` (prefix `procgen_incoming_20260925_`; imports rewritten) by the orchestrator. The infinite claims above (F1, F2, F4) rest on the stated one-line proofs and finite certificates, not on the ranges.
