# A commuting marked-module grammar and an exact gap barrier

**Status: PROVED (with a finite symbolic skeleton certificate for the explicit core).**
This is a positive two-operation extension of the selected-set seed grammar,
not a proof of Chan–Pak. The eight starting marked modules below have literal
commuting insertions and an exact three-coordinate state. Their odd-slope
spectrum is also classified; this fixed carrier alone is sparse.

## Inheritance, ownership, and live concepts

We use the selected-set quantities H, W, U and tail identity
H(tail_t)=H−W+2^t U from THM-5126 / Session258. The selected set here has
exactly two vertices. Session258's CORE_CLASSIFICATION proves every positive
unit seed has one marked strong core, with exterior source padding neutral and
exterior sink padding multiplying U by two. Its missing odd-slope size bound
is an existence question, not a claim that all representatives have bounded
order.

The five-vertex core and the isolated module ladders come from the concurrent
`../constructors/module_ladders.py`. The finite compressed-path method is the
concurrent `../constructors/marked_position.py`; we extend it to two marked
clones and three free intervals. The endpoint-swap lemma and the commuting
marked-module system below are the coverage lane's contribution. Root O/Z/J
inverse fibers, positive coordinates, and binary height were already established
in `session198/bilinear-descent/PROOF.md`, `session200/root-weighted-height/PROOF.md`,
and THM-4904; none is claimed new here.

Incoming work at commit `48ffb162415d`,
`05-knowledge/results/seed_recursion_births_20261111/PROOFS.md`, R1–R3,
proves simultaneous stability of60 inherited stable-port states, commuting
all-size inflations,322 old-port transport classifications, and complete head
grids for coefficients3·2^r+1 and21·2^r+1. We credit the common finite
run-tensor method and do **not** claim a first general commuting-inflation
principle. The distinction here is reselection of two vertices **inside** an
inflated module, classification of their intervening gap, the explicit boundary
defect δ, and its rearming/odd-head return obstruction. The incoming R2
statement explicitly leaves reselection and core changes open, so these
results are compatible.

Concept board: (i) collective target coverage, (ii) adjacent twin contraction,
(iii) a moving selected pair, (iv) the signed unit defect, (v) internal rather
than exterior dyadic padding, (vi) the remaining marked-gap coordinate.
The broad root-change route did not yield a new closure in this pass. Its
useful output was to redirect attention to an actual marked-state carrier.

## 1. The consecutive-twin lemma

Let x→y be vertices with identical orientations toward every other vertex.
They are consecutive transitive twins, so their contraction T/xy is a
well-defined tournament. Select O={x,y}. Then

    W_T(O)=2 H(T/xy),
    H(T)−W_T(O)=H(T)−2H(T/xy).                         (1)

To prove this, K_xy counts Hamilton paths of T containing the consecutive
edge x→y: concatenate the two covering paths. Contracting that edge gives a
bijection with Hamilton paths of T/xy. Also K_yx=K_xy. Indeed swap x and y
between the two paths in a cover ending at y and starting at x. Neither path
contains the other marked vertex, and all affected neighboring vertices are
external twins' neighbors, so every affected arc remains valid. This operation
is an involution onto covers ending at x and starting at y. Diagonal K terms
are zero because covering paths are disjoint. These facts prove (1).

The twins hypothesis is essential. A pair separated by other vertices of a
transitive module need not satisfy it; Section 3 retains that gap explicitly.

### A finite structural criterion for an entire twin-seed tower

Fix a tournament S with a vertex v and replace v by a transitive module of
size m≥1. Write H_m for its count. Every adjacent pair in that module gives
a unit seed at size m≥2 iff

    H_m=2H_(m−1)+1.

Thus the whole tower consists of unit seeds iff

    H_m=(H_1+1)2^(m−1)−1 for every m≥1.                 (2)

There is a finite combinatorial criterion for (2), rather than an infinite
sequence of tests. Let c_r count directed path skeletons through all exterior
vertices and r nonempty module runs. Consecutive module symbols are merged,
and their order along the path distinguishes the r run slots. If there are
k exterior vertices, r≤k+1 and

    H_m = Σ_(r=1)^(k+1) c_r r! S(m,r),                 (3)

where S(m,r) is a Stirling number of the second kind. Assign the m ordered
clones onto the r slots surjectively; each slot's assigned clones appear in
increasing order. This is a bijection with paths. In particular c_r≥0.

**Tower criterion:** (2) holds iff c_r=0 for r≥3 and c_1=2c_2−1.
Sufficiency uses S(m,1)=1 and 2!S(m,2)=2^m−2. For necessity, a positive
largest c_r with r≥3 would force growth of order r^m, incompatible with (2).
With only r=1,2 remaining, compare the constant and 2^m terms. This criterion
is independent of the explicit five-core below.

## 2. The concrete marked family

Take the tournament with zero-based outgoing-bit rows

    B=(24,1,3,6,14).

Choose v∈{1,2,3,4}. Replace v by the ordered transitive module

    P_a → x → G_c → y → S_b,       a,b,c≥0,

where P,G,S are unselected chains of the indicated sizes, and select exactly
x,y. Every clone has v's old orientations to the four external vertices.
The resulting tournament has a+b+c+6 vertices. It is strong: B is strong,
and replacing a vertex of a strong tournament by any nonempty tournament
preserves strong connectivity, using an external quotient cycle for travel
backward inside that module.

Use the convention 0^0=1. For every a,b,c≥0 and every listed v,

    H = 24·2^(a+b+c)−1,
    D := H−W
      = 0^c+2^c−1−12·2^(a+b)(3^c−2^(c+1)+1),          (4)
    U = C_v·2^a·3^b (2^(c+1)+3^(c+1)) − e_v·0^a,    (5)

with

    (C_1,e_1)=(12,0),  (C_2,e_2)=(18,1),
    (C_3,e_3)=(20,0),  (C_4,e_4)=(30,3).

The module count has the two-run signature c_1=11,c_2=6. Equation (1)
therefore proves D=1 for adjacent selected clones without calculating U.
Equations (4)–(5) supply the additional nonadjacent case and its exact boundary.

### Why the formulas are all-length identities

`marked_module.py` enumerates compressed path words with six fixed vertices
(the four exterior vertices and x,y) and three repeatable run symbols P,G,S.
A repeated symbol may not follow itself. The run symbols form an acyclic
ordered graph, so between successive fixed vertices there are at most three
run symbols. The symbolic enumeration is finite, independent of a,b,c.

For each word with r P-runs, s G-runs, t S-runs, expansion contributes
r!S(a,r)·s!S(c,s)·t!S(b,t). For two-path covers the run slots in both paths
are combined before surjective assignment; this enforces disjoint exhaustive
allocation of every clone. First-path emptiness is allowed only for U.
The two selected fixed vertices are never repeatable symbols.

The exact conversion

    r!S(n,r)=Σ_(j=0)^r (−1)^(r−j) binom(r,j) j^n

gives the following small exponential coefficient certificate. A tuple
(i,j,k):z denotes z·i^a·j^c·k^b:

    H: (1,1,1):−1, (2,2,2):24.
    D: (1,0,1):1, (1,1,1):−1, (1,2,1):1,
       (2,1,2):−12, (2,2,2):24, (2,3,2):−12.
    U: (2,2,3):2C_v, (2,3,3):3C_v,
       and (0,1,1):−e_v when e_v≠0.

The program checks this certificate by exact integer equality, not by fitting
a bounded list of numerical outputs. An additional 108 subset-DP controls
(a,c,b∈{0,1,2}, four v) are independent of the compressed-word calculation.
These are controls, not the justification for all-length extrapolation.

## 3. Exact gap classification

These marked modules are positive unit seeds **iff c=0 or c=1**. This statement
allows all prefix and suffix lengths and all four core positions.

For c=0,1 the final bracket in (4) vanishes and its first part equals one.
For c≥2 let A_c=3^c−2^(c+1)+1. Then A_2=2 and
A_(c+1)=3A_c+2^(c+1)−2≥2A_c, so A_c≥2^(c−1).
Consequently 12·2^(a+b)A_c>2^c−1, and (4) gives D<0.

Thus moving the pair is a real freedom, but the intervening gap has exactly
two admissible values. Forgetting it would predict false unit seeds.

## 4. A literal commuting two-operation grammar

Fix v and an admissible c∈{0,1}; write q=U and

    δ=e_v·0^a,   κ=C_v(2^(c+1)+3^(c+1)).

Then q+δ=κ·2^a·3^b. Two physical operations preserve the unit-seed property:

- L inserts one new unselected clone at the beginning of P.
- R inserts one new unselected clone at the end of S.

Their complete numerical state transitions are

    L: (q,δ) ↦ (2q+2δ, 0),
    R: (q,δ) ↦ (3q+2δ, δ).                            (6)

Both also send H to 2H+1; they retain v,c and the selected x,y. The two
insertions commute as marked tournaments up to the natural ordered isomorphism,
not merely as q-values: one adds the leftmost unselected clone, the other
the rightmost. Equation (6) likewise gives the same final state
(6q+6δ,0) in either order.

The offset coordinate δ is necessary. At v=2,c=0,a=b=0 one has q=89,δ=1.
R gives q=269, not 267, and L gives q=180, not 178. After L the boundary
defect is erased and subsequent moves really multiply q by two or three.
For v=1,3 the defect is always zero. This is a small genuine marked-state
semigroup, not an append-only scalar induction for arbitrary tournaments.

Keeping a+b fixed and sliding the pair rightward changes the pure interior
multiplier by 2/3. The underlying unmarked tournament is unchanged. That is
a concrete example of why the selected pair must be retained in a collective
coverage argument.

## 5. Strong internal dyadic padding and exact scope

At v=1,c=b=0,

    q=60·2^a, H=24·2^a−1, order=a+6.

All cores are strong, all are unit seeds, and all have the same odd slope15.
There are arbitrarily large strong representatives of this fixed slope,
without any exterior sink padding. The move L is an exact inverse-compressible
internal source of 2-adic content. This disproves a bound on **every** strong
representative's order in terms of its odd slope, but does not disprove a
bound on the **smallest** representative and does not settle Session258's
open seed-existence question.

For completeness the odd slopes of this entire eight-start grammar are
exactly the union, over b≥0, of the following twelve sequences:

    15·3^b, 39·3^b, 25·3^b, 65·3^b,
    45·3^b, 117·3^b, 75·3^b, 195·3^b,
    90·3^b−1, 234·3^b−1, 150·3^b−3, 390·3^b−3.

For v=1,3 this follows by stripping 2^(a+2). For v=2,4, a=0 yields the
last four odd expressions; a≥1 erases δ and stripping the remaining power
of two gives the middle four. In particular there are only O(log X) distinct
slopes ≤X. No universal ternary numeration theorem follows from these commuting
moves. A new operation must change the primitive carrier κ or regenerate
boundary defects in a way not already described by (6).

Each q yields actual tail counts 1+2^t q, t≥1. Individual small head points
and the Chan–Pak exclusions remain separate. The result is a positive grammar
extension, an exact missing coordinate, and a compression target for future
collective searches, not a resolution of full odd-count coverage.

## 6. A correct binary-return criterion, and its small hostile

Simply rearming δ=1 after L does **not** generate exact odd multipliers.
Every L output is even, and every R preserves the parity of q. Once an L
has occurred, no sequence of L,R and same-q rearmings can produce odd q.
This is an invariant of the proposed enlarged grammar, not just of our eight
starting modules.

There is nevertheless a genuine collective binary-coverage lemma. Let B≥1.
Assume:

1. There are physical unit seeds with exact multipliers 2p for every integer
   B≤p<2B (a finite bank).
2. Whenever an even multiplier q≥2B is produced, a physical same-q rebase
   equips its seed with δ=1 and a legal L move of (6). Its output remains
   an actual unit seed; merely assigning a formal δ label does not suffice.

Then every even integer q≥2B is a unit-seed multiplier. Indeed an arbitrary
seed can always double q by exterior sink padding. In the coordinate p=q/2,
this is p↦2p. Rearming and then L give q↦2q+2, hence p↦2p+1.
For p≥2B its parent floor(p/2) is smaller and at least B. Repeated binary
parent descent reaches the bank [B,2B), proving the claim by induction and
giving logarithmic construction depth beyond the finite bank.

Consequently **every odd slope** c has a seed under these two assumptions:
choose r≥0 with p=2^r c≥B, and use q=2^(r+1)c. The construction covers all
odd slopes, including 7 and21, which are not excluded slopes. It does not by
itself provide a seed of exact odd multiplier c, nor the first ray head
1+2c. In particular it is not yet a Chan–Pak proof.

If, in addition, every resulting q=2p≥2B admits a physical unit-seed
**deflation** to exact multiplier p, then every p≥B is a seed multiplier.
One tail vertex now gives every odd count N=1+2p≥2B+1. Actual tournament
witnesses for the remaining allowed odd N<2B+1 would finish Chan–Pak.
These three hypotheses separate binary reachability, physical rearming,
and the otherwise missing odd-head return.

The unrestricted hypotheses are false at small q, with inherited exact
hostiles rather than conjectural numerical holes. Session258's all-order
fixed-q classifier gives q=4 as a seed and excludes q=10. A physical rearming
of q=4 with legal L would create 2·4+2=10, contradiction. Likewise q=6 is
a seed but q=3 is not, so unrestricted physical halving is false. These
examples do not refute eventual rearming or eventual deflation.

On our explicit carrier the first boundary example q=89,δ=1 has L-output
q=180,δ=0. A successor-style half-return would ask for q=90. That multiplier
is absent from our exact eight-start image and also from the constructors'
finite two-core box (module lengths≤5,total order≤12). The latter is only a
finite-box observation; no all-tournament exclusion of q=90 is claimed.
Thus the next return must really change the marked core carrier. It cannot
be accomplished by invisible relabelling inside the fixed module family.

## Reproduction

    python3 marked_module.py
    python3 -O marked_module.py

Both write `marked_module_receipt.json`; no third-party dependencies are used.
The literal DP helper is imported from `../constructors/module_ladders.py`.
