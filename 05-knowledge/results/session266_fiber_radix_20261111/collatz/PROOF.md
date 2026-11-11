# Punctured Collatz graphs, path-completion fibers, and the missing height coordinate

**Session266. PROVED elementary statements below; FINITE-EXACT controls in `probe.py`. No global Collatz convergence or cycle-exhaustion claim.** H always counts directed Hamiltonian **paths**, not cycles. This lane supplies new synthesis and small constructive controls, not a retrospective verification campaign.

## Inheritance and concept board

Session255's signed central-table atlas and inverse ladder are accepted inputs, pinned in `SOURCES.md`. Session265 supplies the distinction between realizable Hamiltonian counts and mere arithmetic residues. We retain five concepts: literal puncture / signed affine word / table modulus / four-vertex tournament / path-pointed isomorphism class. The hostile is that the same directed four-vertex path has eight tournament completions. The least-used coordinate is the complete set of nonconsecutive chords, followed by exact integer height.

## 1. The two requested cuts have different topological types

Let C(x)=x/2 for even integers and 3x+1 for odd integers. On positive magnitudes the negative sheet is C_-(x)=x/2 on even x and 3x−1 on odd x: C(−x)=−C_-(x). Write C_+ for the usual positive map. Remove exactly the edges 5→16 and −2→−1, leaving their vertices present and with no outgoing edge.

**Puncture lemma.** For a function graph, deleting e=(a,f(a)) changes only forward trajectories that reach a; each such trajectory stops at its first visit to a. If e is off every directed cycle, all cycles survive. If e is on a directed cycle, exactly that cycle is destroyed, and its component becomes a rooted incoming graph ending at a. This is a literal directed-graph assertion, with no assumption that every component has a cycle.

*Proof.* Before the first visit to a all edges are unchanged. A cycle is deleted iff it used e. A function has at most one directed cycle in a weak component, since two forward orbits that merge never separate. In the cycle case every vertex of that weak component reaches the cycle: a finite undirected connecting path can be followed inductively using this eventual-meeting property. Removing the edge stops those trajectories. QED.

For the plus map, 5→16→8→4→2→1→4 proves that the edge from 5 is off-cycle. The known positive cycle 1→4→2→1 survives untouched. Its former basin separates into the part that first hits 5 and the part that never hits 5. For the negative sheet, −1→−2→−1 is a cycle and the requested edge destroys it. The known cycles through −5 and −17 do not use that edge and survive. This does not assert that these are all negative cycles, or that the known positive cycle is the only positive cycle.

On odd vertices, use T_s(x)=(3x+s)/2^{v2(3x+s)}. The plus cut deletes the accelerated edge 5→1 (exponent 4). The negative raw cut deletes the accelerated loop −1→−1: in positive magnitudes, the T_- loop 1→1 (exponent 1). These are a noncyclic entrance and a periodic loop, respectively. They are not symmetric cuts of two corresponding cycles.

The first few direct plus predecessors of odd target 1 are

    1, 5, 21, 85, ... , x_j=(4^j−1)/3, j≥1.

Removing the edge from 5 leaves 21→1 (exponent 6) intact. In contrast 7→11→17→13→5 now terminates at 5. Therefore this cut treats the two inherited missing tournament counts 7 and 21 differently. Neither deleting these arcs nor identifying their basins is an explanation of the theorem that 7 and 21 are unattainable Hamiltonian counts.

## 2. The table records the two holes exactly, but residues cannot

The inherited atlas assigns a positive source x on sign s the odd modulus m=4x+s. Odd sources use the central pair (x,3x+s) directly; even sources use the oriented halving decoder. Thus:

- 5→16 is the plus-sheet odd-source table with modulus 21.
- −2→−1 is the negative of the minus-sheet even-source edge 2→1, whose modulus is 7.

**This is the precise appearance of 21 and 7 requested by the two missing arcs.** They are table **moduli**, not vertex labels or Hamiltonian counts in this map. The pair of prescribed deleted edges has exactly the pair {7,21} as its signed table codes. The agreement with the inherited excluded Hamiltonian counts is exact arithmetic, but a tournament realization map is still missing. No nonattainment proof follows just by renaming these table codes as counts.

There is an exact accelerated table transport. If y=T_s(x), k=v2(3x+s), m=4x+s and m'=4y+s, then

    2^k (m'−s)=3m+s.

Indeed both sides equal 4(3x+s). The plus hole is m=21 to m'=5 with k=4. The negative magnitude loop is m=m'=3 with k=1; the raw cut itself has code 7 because its source is the intermediate even integer 2. Thus acceleration forgets the intermediate table code 7 unless its even-edge list is retained.

**Literal-hole obstruction.** No predicate depending only on an integer modulo any fixed q≥1 can recognize the single source 5. The integers 5 and 5+q have the same residue, but only one loses its edge. The same argument applies to source 2. Hence a finite congruence picture cannot implement either literal puncture exactly. A flag for exact source equality, or a retained height/quotient coordinate, is indispensable. The code atlas itself is exact because m is retained as an integer, not merely reduced modulo a second modulus.

Connection contract: source = signed central tables with parity and orientation; target = raw Collatz edges; preserved = source, target, and literal edge deletion; lost under acceleration = even intermediates and their table codes; lost under finite congruence = exact single-edge identity; repair = signed exponent plus full integer lift.

## 3. The exact reason modular ± pictures cannot classify positive cycles

For shortcut maps F_s(x)=x/2 if even and (3x+s)/2 if odd, fix a binary word ε_0,...,ε_(L−1), with r ones. Its formal composition is

    F_word,s(x)=(3^r x+sB)/2^L,
    B=Σ_{j: ε_j=1} 2^j 3^{Σ_{i>j} ε_i}.

For r>0, B>0. Its only periodic candidate is

    x_s=sB/(2^L−3^r).

**Signed-necklace criterion.** Such a word gives a genuine nonzero integer cycle exactly when D=2^L−3^r divides B; repeated words give repeated traversals, and cyclic rotations change the marked starting point. The integer candidate has positive magnitude on the plus sheet iff D>0 and on the minus sheet iff D<0. Integrality is the same on the two sheets; positivity is opposite.

*Proof.* Affine composition gives the formula. D is odd and nonzero. An integral fixed candidate automatically follows the prescribed parity itinerary: reducing the fixed equation modulo 2 forces the first prescribed parity. After its legal first division, repeat with the cyclically shifted word. Equivalently the final divisibility relation modulo 2^L forces each successive division. Negation conjugates F_+ and F_- and preserves binary parity. B>0 proves the sign assertion. QED.

The explicit known examples for the plus map on signed integers are:

| marked source | shortcut cycle | L | r | B | D |
|---|---|---:|---:|---:|---:|
| 1 | 1,2 | 2 | 1 | 1 | 1 |
| −1 | −1 | 1 | 1 | 1 | −1 |
| −5 | −5,−7,−10 | 3 | 2 | 5 | −1 |
| −17 | −17,−25,−37,−55,−82,−41,−61,−91,−136,−68,−34 | 11 | 7 | 2363 | −139 |

These are four witnesses, not an exhaustive classification. On every odd residue ring, negation conjugates the corresponding formal ± affine branch registers, so their modular closed-word information agrees. The positive/negative distinction lies in the **ordinary sign of D**, which the modular picture has forgotten. This locates an exact obstruction rather than turning “three known negative cycles versus one known positive cycle” into a theorem about tournaments.

### All-radix version

Let b≥2. For each digit d in {0,...,b−1}, choose integers a_d,c_d with gcd(a_d,b)=1 and c_d≡−a_d d (mod b), and use the integral branch f_d(x)=(a_d x+c_d)/b when x≡d (mod b). A digit word has

    f_word(x)=(A x+B)/b^L,
    A=Π a_(d_i),
    B=Σ_j b^j c_(d_j) Π_(i>j) a_(d_i).

If b^L≠A, its periodic candidate is B/(b^L−A). Integrality is sufficient for the itinerary, by the same successive congruence argument: all a_d are units modulo b. If b^L=A, the word is either inconsistent (B≠0) or an identity affine word (B=0), requiring separate treatment.

Negation changes digit d to −d mod b, and changes c to −c. In binary, each digit equals its negative modulo 2, so the two sheets keep the SAME parity word. In other bases the reflected digit word must be retained. This is a genuine radix-dependent distinction, not a decimal-only coincidence. The affine-word mechanism supplies a common algebraic language with polynomial-growth recurrences without claiming that their dynamics coincide.

## 4. What four-vertex tournaments actually say

There are four unlabelled tournament types on four vertices:

| type | sorted outdegrees | cyclic triangles | H | automorphism order | labeled copies |
|---|---|---:|---:|---:|---:|
| transitive | 0,1,2,3 | 0 | 1 | 1 | 24 |
| C3 followed by a sink | 0,2,2,2 | 1 | 3 | 3 | 8 |
| source followed by C3 | 1,1,1,3 | 1 | 3 | 3 | 8 |
| strong | 1,1,2,2 | 2 | 5 | 1 | 24 |

A quick classification: if there is no directed triangle the tournament is transitive. Otherwise fix a triangle; the fourth vertex either beats all three, loses to all three, or has mixed adjacency. The two mixed cases are isomorphic after relabeling/reversing the choice of displayed triangle. Representative path enumeration gives the table, and the executable checks all 64 labeled inputs. Thus H=1+2 c3 in this order. The middle two classes are converse; the outer two are self-converse. There are four directed isomorphism classes but three converse classes, with one strong and three nonstrong classes. These are different meanings of “three versus one.”

The residue classes 1,3,5,7 modulo 8 are not in bijection with the four tournament types by H: 3 appears twice and 7 not at all. Any correspondence needs the source/sink chirality sidecar at H=3. A Hamiltonian **cycle** convention would give different statistics and is not used here.

## 5. A real all-order bridge: complete an observed orbit path

A Collatz segment with n distinct successive vertices gives a directed path, not a tournament: most pairwise comparisons are absent. Retain that path and orient every missing chord arbitrarily. This yields 2^((n−1)(n−2)/2) tournament completions, without identifying the arithmetic vertices with their numerical magnitudes. Orienting by numerical order instead would always give a transitive tournament and discard the actual dynamical arrows whenever they decrease.

**Path-completion fiber theorem.** For an unlabelled n-vertex tournament T, exactly

    H(T)/|Aut(T)|

of the completions of the fixed path 0→1→...→n−1 are isomorphic to T. Consequently |Aut(T)| divides H(T), and

    Σ_[T] H(T)/|Aut(T)|=2^((n−1)(n−2)/2).

*Proof.* Every directed Hamiltonian path in T labels its vertices in order by 0,...,n−1, producing one completion. Two paths produce the same labeled completion iff an automorphism carries one ordered path to the other. The automorphism action on ordered Hamiltonian paths is free: fixing the ordered path fixes every vertex. Hence each fiber has size |Aut(T)|. QED.

At n=4 the eight completions split 1,1,1,5 among the four classes in the table. The last “5” is not a fifth tournament type: it is the number of path labelings of the unique strong class. A four-step orbit does not select one of these completions unless three chord choices are supplied.

For the broader multiplicity classification let C_n(h) count fixed-path completions with H=h, L_n(h) labeled tournaments with H=h, and A_n(h) unlabelled classes. Then

    L_n(h)=n! C_n(h)/h,
    A_n(h)=(1/h) Σ_(completion Q: H(Q)=h) |Aut(Q)|.

This gives an exact census route on only 2^((n−1)(n−2)/2) completions, a factor 2^(n−1) fewer edge assignments than all labeled tournaments. The first formula needs no canonicalization or automorphism calculation. The second needs the automorphism sidecar, or one can canonically deduplicate completions. No claim is made that this factor removes the general exponential complexity.

The map thus preserves the observed directed path and exposes a precise **Hamiltonian-size-biased fiber measure**. It loses arithmetic branch laws; those must remain attached separately. Conversely a tournament without an actual integer orbit embedding need not describe Collatz dynamics. This is the constructive bridge available now, rather than an unsupported bijection of cycle counts.

## 6. Focused frontier

1. Can a natural arithmetic comparison of nonadjacent Collatz orbit points select chords, rather than an arbitrary completion? State the comparison and tie policy first.
2. Can the exact table codes 7 and 21 of the two punctures be consumed by a tournament operation that explains their count exclusions? The equality of numbers is established; the realization map is OPEN.
3. Can cycle-necklace denominator divisibility be counted jointly with the sign/height coordinate rather than in modular registers alone?
4. Use fixed-path completion histograms for the main H-fiber classification. Their labels and automorphism weights are now explicit and cannot be mistaken for unlabelled multiplicities.

Reproduce only these new local claims with `python3 probe.py`. All formal all-order statements above have proofs; the finite receipt does not claim convergence or exhaustive integer-cycle classification.

## 7. Punctured path forests retain a richer, exact multiplicity interface

The previous construction extends from one observed path to a path with missing arcs. Let λ=(λ_1,...,λ_r) be a composition of n into positive lengths. Partition the fixed labels 0,...,n−1 into consecutive blocks of these lengths; prescribe only the forward consecutive edges **inside** each block. Keep no prescribed edge across a gap. There are n−r prescribed arcs, so the number of completions is

    2^(binom(n,2)−n+r).

For a tournament T define its ordered path-cover response

    C_λ(T)=Σ_(ordered partitions V=B_1⊔...⊔B_r, |B_i|=λ_i)
                    Π_i H(T[B_i]).

The blocks are distinguished and ordered even when their lengths agree. Each summand counts the choices of a Hamiltonian path on every block.

**Punctured-path fiber theorem.** The number of completions belonging to [T] is exactly C_λ(T)/|Aut(T)|. In particular the quotient is an integer and

    Σ_[T] C_λ(T)/|Aut(T)|=2^(binom(n,2)−n+r).

*Proof.* A choice counted by C_λ gives a complete vertex labeling: place the first selected block-path in the first λ_1 positions, then the second, and so forth. It produces precisely a completion of the prescribed forest. Two choices produce the same labeled tournament iff related by an automorphism. This action is free because every vertex occupies a specified position. QED.

For one missing edge in a four-vertex path, λ is (1,3), (2,2), or (3,1). In the same class order as §4, the completion multiplicities are

    λ=(1,3):  (4,2,2,8), total16;
    λ=(2,2):  (6,2,2,6), total16;
    λ=(3,1):  (4,2,2,8), total16.

For (1,3), C_λ=Σ_v H(T−v). Since every cyclic triple contributes H=3 instead of1, C_(1,3)=4+2c3=H(T)+3 in order4. Dividing by the class automorphism order yields the first row. For (2,2), every two-vertex tournament has H=1 and there are six ordered two-two partitions; hence C_(2,2)=6 independently of T, and dividing by the automorphism order gives the second row.

Thus **where** a path arc is deleted changes the tournament-completion measure. The source/sink converse distinction remains, but a gap at an endpoint and a gap in the middle weight the transitive and strong classes differently. This is a concrete connection from missing dynamical arcs to the repository's path-cover interface. It does not require claiming that the numerical Collatz vertices determine missing tournament chords.

### Binary no-carry law for every punctured profile

Accept the inherited tournament parity theorem H(U) odd for every tournament U. Every product in C_λ is odd, hence

    C_λ(T) ≡ multinomial(n;λ_1,...,λ_r) (mod2).

The path-completion theorem implies |Aut(T)| divides the odd integer H(T), so the automorphism order is odd. Division therefore preserves this parity:

    C_λ(T)/|Aut(T)| ≡ multinomial(n;λ_1,...,λ_r) (mod2).

By the binary multinomial Lucas criterion, the right side equals1 precisely when adding the λ_i in base2 has **no carries**. One elementary proof uses

    (z_1+...+z_r)^n = Π_(j: bit_j(n)=1)(z_1^(2^j)+...+z_r^(2^j))

in characteristic2. A specified monomial z_1^λ_1...z_r^λ_r appears once iff each binary digit of n is assigned to exactly one λ_i. This proves the no-carry criterion without a number-theoretic citation.

The striking consequence is class-independent parity of every completion fiber: either every class appears oddly or every class appears evenly, according solely to the composition. At n=2^k, every proper multi-block composition requires a carry, so every punctured-profile fiber is even. In particular all weights in the two order-four puncture rows are even. This is a real radix bridge: ordinary carry geometry determines the parity of path-cover multiplicity, while Collatz's signed cycle criterion additionally needs ordinary height.

For λ=(n), no carries are possible; all H(T)/Aut(T) are odd. At n≥3 their total 2^((n−1)(n−2)/2) is even, so the total number of unlabelled tournament classes is even. This is a consequence of the same completion measure, not a new exhaustive census.

The last parity statement does not extend automatically from base2 to an arbitrary prime: H is universally1 modulo2 but not modulo a general prime, and automorphism groups can have odd-prime factors. The general-radix affine word in §3 and the binary path-cover parity law have different scopes.

## 8. The multiplicity census is a weighted signed-edge atlas

The numerical hit in §2 admits a precise, if deliberately coarse, map. Every odd h≥3 has a unique expression

    h=4x+s,  x≥1, s∈{−1,+1}.

Choose s=+1 when h≡1 mod4, and s=−1 when h≡3 mod4. Send a tournament T with H(T)=h≥3 to the raw signed edge encoded by h: for s=+1 use x→C_+(x), and for s=−1 use −x→−C_-(x). H=1 is retained as a separate boundary symbol (formally x=0,s=+1, outside the nonzero-integer edge graph). This is a map from actual tournaments, not a claim that their vertices are Collatz integers.

If A_n(h) counts unlabelled n-vertex tournament classes with H=h, the exact weight of edge code h is

    w_n(h)=A_n(h).

Likewise cumulative weights are w_≤N(h)=Σ_(n≤N) A_n(h). Labeled weights or converse-orbit weights can be used instead, but must be named as different measures. The full fiber over an edge is exactly its H-isomorphism-class fiber; this is the classification requested, placed on the signed Collatz graph. At n=4, the negative edge −1→−2 (code3) has weight2 and the positive edge 1→4 (code5) has weight1; the transitive boundary symbol has weight1.

Using the inherited nonattainment of7 and21, the two requested edges have zero weight at **every size**. The claim that every other nonzero signed edge eventually has positive weight is exactly the assertion that every odd h other than7,21 occurs as a tournament Hamiltonian-path count. Thus the full-spectrum problem is equivalently an eventual-support problem on this specific punctured signed-edge atlas. This is an equivalence of realization problems, not a proof of either, and it has no implication that Collatz trajectories converge.

### Ordinal composition transports to an explicit edge-code multiplication

For ordinal sum T⊕U (all arcs from T to U), every Hamiltonian path is uniquely a path in T followed by one in U, so H(T⊕U)=H(T)H(U). In the signed code coordinates this becomes

    (x,s) ⊙ (y,t) = (4xy+t x+s y, st).

Indeed (4x+s)(4y+t)=4(4xy+t x+s y)+st. Include the boundary identity (0,+1). For x,y≥1 the resulting first coordinate is positive, even in the worst case s=t=−1. Therefore the signed-edge atlas inherits an associative commutative code product, and the map from tournaments under ordinal sum preserves that product at the level of H. Tournament ordinal sum itself is not commutative up to directed isomorphism: the code map loses the order of SCC factors. That lost order is exactly part of the multiplicity information the main session retains.

Example: a directed triangle has h=3, code (1,−1), corresponding to −1→−2. Two ordered triangle blocks have h=9 and code (2,+1), corresponding to 2→1, since (1,−1)⊙(1,−1)=(2,+1). This is an actual tournament construction inducing a specified algebraic operation on edge codes. It is **not** one step of the Collatz map on those code labels.

Connection contract: source = tournament isomorphism classes with H and SCC order; target = integer-labeled signed Collatz edges with multiplicity; map = h↦(x,s) above; preserved = exact H-fiber cardinality and ordinal-product law; destroyed = internal arcs, SCC order, graph order unless stored in w_n, and all claims of dynamical orbit adjacency; required sidecar = n, class identity/automorphisms, and SCC sequence. The forbidden code pair is exactly the requested puncture pair. This supplies a mathematically usable bridge while preserving the open support problem.

### The induced code dynamics and its smallest physical obstruction

Transport the raw Collatz source dynamics through h=4x+s. Its exact map on odd codes h≥3 is

    Φ(h)=4 C_s(x)+s
        = (h+1)/2,  if h≡1 (mod8),
        = 3h−2,     if h≡3 (mod8),
        = 3h+2,     if h≡5 (mod8),
        = (h−1)/2,  if h≡7 (mod8).

These formulas follow by substituting x=(h−s)/4 into the odd or even branch. The residue h mod4 preserves the signed chart. The familiar positive raw cycle becomes

    5→17→9→5,

while the negative magnitude two-cycle becomes

    3→7→3.

The prescribed cut at source code7 breaks the latter. The prescribed cut at source code21 removes the entrance

    21→65→33→17→9→5→17.

This code graph is literally conjugate to the signed raw graph, before puncturing; the statement is stronger than a coincidence of visible arrows. However the tournament map to code values is **not** a dynamical conjugacy. A triangle is an actual H=3 tournament, but Φ(3)=7, which is an inherited impossible tournament count. Hence no tournament operation defined on every tournament can implement Φ on Hamiltonian counts. Any proposed physical realization must restrict its domain, retain a failure state, or change its count law. This tiny obstruction remains even if every odd value except7,21 is eventually realized.
