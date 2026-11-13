# Short cluster returns and a measured paired-port code

**Status: PROVED elementary for the short-return and enumeration results;
INHERITED for the native-component classification; FINITE-EXACT for the
observation code.** The elementary mutation arguments were independently
reviewed in [the square-lane audit](../squares/AUDIT-cluster.md). No Chan–Pak
proof, new count coverage, or classification of
general quiver mutation classes is claimed. We classify the tournament locus
under ordinary matrix mutation and excursions of length at most three. The
code uses the declared 24 inherited port charts, not arbitrary coordinates.

## Inheritance and connection contract

Closest mechanisms: the source/sink remark of
[THM-4602, tournaments as chirotopes](../../../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md),
[Session290 exact profiles](../../session290_core_descent_20261113/descent/PROOF.md),
and the incoming
[Session291 paired-port theorem](../../session291_cluster_ports_20261113/ports/PROOF.md).
During this work, incoming commit e08896a818 promoted Session291 after
independent review. Its [cluster proof](../../session291_cluster_ports_20261113/cluster/PROOF.md)
already proves the source/sink and complete native-component classifications,
and the ambient pentagon. Propositions 1 and 2 below are therefore INHERITED
restatements, independently rederived and checked here. The enumeration
corollary, short-excursion classification and measured distance-five code
are the scoped extensions. Our code independently reads Session290's data.

Hostiles: mutating a directed triangle deletes an edge; the canonical paired
binary code has weight-two words. Corrected near miss: “12 base graphs” means
12 adjacency-plus-selection arrays, eight literal adjacency arrays, and 24
port charts, not twelve unmarked isomorphism classes. Least-used sidecar:
the full continuation profile rather than its terminal value.

Concept board: **SCC hull / return length / physical H / port marking /
modular observation / distance**. Anchor: lawful returns; niche: an exact
enumeration quotient; wildcard: error-correcting identification of ports.

| Source | Target/map | Preserved | Lost or required sidecar |
|---|---|---|---|
| Tournament | Skew matrix, +1 for i→j | Every edge | Mutation may leave tournaments |
| Native component | SCC hull and exterior singleton count | H and hull | Labelled exterior order |
| Twelve binary choices | Eleven observations mod 13 | All 4096 identities | Sum is not a single physical H |
| Paired 24 coordinates | Canonical paired binary code | Pairing, dimension 12 | Distance 2, not Golay distance 8 |

Matrix mutation is CITED from Fomin–Zelevinsky,
[*Cluster algebras I*](https://arxiv.org/abs/math/0104151).

## 1. Complete native component classification

Set b_ij=1 for i→j, b_ji=−b_ij, b_ii=0. Mutation at k is

    b'_ij = −b_ij                              if i=k or j=k,
    b'_ij = b_ij+(|b_ik|b_kj+b_ik|b_kj|)/2     otherwise.

A step is **native** when both matrices are tournaments. There is no sign
projection or arbitrary filling of zero entries.

**Proposition 1.** A mutation from a tournament is native iff its vertex is
a global source or sink.

Proof. If I and O are the in- and out-neighbors of k, each I×O entry becomes
b_ij+1, hence 0 or 2. All other nonincident entries stay unchanged. Thus
exactly I=∅ or O=∅ permits a tournament output. ∎

In the unique ordered strongly connected component (SCC) decomposition,
define the **SCC hull** B as the interval from the first to the last
nontrivial SCC, including singleton SCCs between them. If T is nontransitive,

    T = TT_a → B → TT_b,      t=a+b,

uniquely, where → is ordinal sum.

**Proposition 2.** Its labelled native component is a path on t+1 states.
At the endpoint where all exterior vertices precede B, write their order
v_1,…,v_t. The path states are

    v_(j+1) → … → v_t → B → v_1 → … → v_j,    0≤j≤t.

Unlabelled, the component still has t+1 distinct classes, represented by
TT_i→B→TT_(t−i). For transitive labelled tournaments it consists of the n
cyclic rotations of the total order: a cycle for n≥3, edge for n=2, singleton
for n=1. Unlabelled, transitive tournaments give one class.

Proof. No vertex inside a nontrivial SCC is a global source or sink.
Proposition 1 permits only the first or last exterior singleton to move
across the whole tournament. These moves increment or decrement j and stop
at the displayed endpoints. Internal singleton SCCs cannot move. SCC order
and uniqueness distinguish the unlabelled states. With no nontrivial SCC,
extreme moves simply rotate the total order. ∎

Every Hamiltonian path visits SCCs in order, so H(T)=∏H(SCC). Thus the whole
component has H(B). Native mutation cannot provide a count-changing operation.

**Corollary 3 (A000568 quotient).** Write a_n for ordinary n-vertex tournament
classes, a_0=a_1=a_2=1, and k_h for h-vertex SCC hulls. Then

    a_n = 1+Σ_(h=3)^n (n−h+1)k_h,
    k_n = a_n−2a_(n−1)+a_(n−2)                     (n≥3).

The number of native components on ordinary n-vertex classes is

    c_n = 1+Σ_(h=3)^n k_h = a_n−a_(n−1)+1          (n≥1).

This identity classifies a quotient; it is not a faster algorithm for unknown
a_n. The exact n=1,…,5 census gives a=(1,1,2,4,12), c=(1,1,2,3,9).

## 2. All returns of length at most three preserve H

An excursion leaves the tournament locus at its first step and returns only
at the end. Consecutive identical mutations cancel.

**Proposition 4.** No non-backtracking excursion has length two.

Proof. After a nonnative mutation at k every I×O entry is 0 or ±2. Every
l≠k belongs to some such pair. Mutating l negates its incident entries,
leaving a bad entry. Only k immediately undoes the first step. ∎

**Proposition 5.** For n≥3, the three-step excursions are exactly k,l,k when

    T = R→k→l     or     T = l→k→R,       |R|≥1.

The output is respectively l→k→R or R→k→l. Both counts equal H(R).

Proof. Restriction commutes with mutation when all mutated vertices are
retained. For a triple k,l,v, writing (a,b,c)=(b_kl,b_kv,b_lv), direct
substitution shows k,l,k returns a tournament exactly for
(a,b,c)=(1,−1,−1) or (−1,1,1). Applying this to every v forces the displayed
extreme transitive two-blocks. In those forms mutation leaves edges inside
R unchanged and gives the asserted output, establishing sufficiency.

For three distinct names k,l,m, the induced triple can return only when
(b_kl,b_km,b_lm)=(1,1,1) or (−1,−1,−1). Reverse all arrows to take the
positive case. If the first step leaves the tournament locus, some v→k.
Retain k,l,m,v and put e=b_lv, f=b_mv. Direct mutation gives:

| (e,f) | Final (b_kv,b_lv,b_mv) |
|---|---|
| (−1,−1) | (−5,−2,4) |
| (−1,1) | (−3,0,2) |
| (1,−1) | (−1,−2,2) |
| (1,1) | (1,0,0) |

Every row forbids a tournament. This four-vertex argument excludes distinct
names at all orders. Adjacent repeated names reduce to one step. ∎

**Corollary 6.** Any sequence of at most three mutations between tournaments
preserves H. A count-changing return must have length at least four and
pass through non-tournaments. This does not assert that length four works.

The checker exhausts all 1,099 labelled tournaments at orders 1 through 5,
5,405 one-step tests, and every non-backtracking three-step candidate. It
finds 12,48,320 three-step excursions at orders 3,4,5. H is independently
computed by subset DP and permutation enumeration. Matrix mutation is
checked against directed path addition/cancellation, including weighted
intermediate states. Four malformed-input controls remain active under -O.

## 3. A measured 13-ary port code

Read the 24 even-chart profiles U_j(m) from Session290. They are integer
linear combinations of b^m for b∈(1,2,3,4,6,8,9,12). For ε∈{0,1}^12 define
the baseline-subtracted formal aggregate

    W_ε(m)=Σ_i [U_(2i+ε_i)(m)−U_(2i)(m)]=Σ_i ε_i Δ_i(m),
    C(ε)=(W_ε(2),W_ε(3),…,W_ε(12)) mod 13.

All bases are nonzero mod 13, so the observations have period dividing 12.
The omitted m=1 phase is identically zero, because paired terminal marked
graphs coincide.

**Finite-exact theorem 7.** The code has 4096 words of length 11 over F_13,
minimum Hamming distance 5. Thus two erroneous symbols or four erased
symbols are uniquely recoverable, and mixed errors/erasures satisfy 2e+s<5.

Five is the minimum number of these periodic observation times that can
always identify ε. Exactly 225 five-time subsets work, including
{2,3,4,6,12}. No four-time subset works; the best has 4064 images, achieved
by {3,4,5,7}. Fewer than four cannot work since 13^3<4096.

The exact ordered distinct-codeword-pair distance histogram is:

| Distance | Ordered pairs |
|---:|---:|
|5|1,512|
|6|13,172|
|7|118,928|
|8|665,658|
|9|2,705,164|
|10|6,427,662|
|11|6,841,024|

The sum is 4096·4095. A distance-five witness is ε=2152, η=1431, as binary
masks whose bit i chooses chart 2i+1 rather than 2i. Their difference is

    (7,7,0,0,8,0,0,0,4,2,0) mod13   at m=2,…,12.

Certificate. Every difference of binary messages is a ternary vector in
{−1,0,1}^12; one with z zero coordinates has 2^z ordered binary preimages.
The script exhausts all 3^12−1 nonzero ternary differences and weights their
histogram accordingly. Independently it computes all 8,386,560 unordered
codeword-pair distances using packed symbols. Both histograms agree. All 330
four-time and 462 five-time subsets of the 11 active phases are tested.
Unique decoding follows from the elementary minimum-distance inequality. ∎

This code is nonlinear over F_13 because 4096 is not a power of 13. Its
notation is (11,4096,5)_13, not [11,12,5]_13. Only the full 11-symbol code
is assigned distance 5; a working five-symbol projection merely identifies
an uncorrupted message.

Each individual slope has physical tail H=1+2^t U_j(m), so its residue can
be recovered from H when the chart and t are known, since 2 is invertible
mod 13. W sums twelve separate observations. No tournament implementing
their addition is supplied. The result improves identification and caching
of this family, not its count-coverage theorem or the full metagraph census.

## 4. Relation and boundary to Leech

The canonical binary difference code of one-choice-per-pair indicators is

    C_pair={(x_1,x_1,…,x_12,x_12): x∈F_2^12}.

It has dimension 12, minimum distance 2, and weight enumerator (1+z²)^12.
Relabelling cannot change distance 2 into Golay's 8. Its elementary
Construction-A lattice is isometric to Z^24 and has norm-one vectors. These
points were already recovered in Session291 and are not new claims here.

The new observation code is a different, measured nonlinear F_13 code with
distance 5 and explicit decoding limits. It is a genuine structure attached
to the ports, but no Golay coupling or Leech isometry has been found.

A genuine distinct tournament-to-Golay map appears in
[THM-481, Paley gauge generates Gleason](../../../../01-canon/theorems/THM-481-paley-gauge-generates-gleason.md):
the bordered Paley tournament at 23 gives the extended quadratic-residue
code. A genuine Leech appearance that loses H is
[THM-4752, tournament scores in II_(25,1)](../../../../01-canon/theorems/THM-4752-tournament-scores-in-ii-25-1-backward-length-and-the-a24-deep-hole.md),
with recorded equal-score examples H=657 and 643. Neither map supplies the
missing physical addition or port-continuation operation here.

## Reproduction and limits

From the repository root:

```text
python 05-knowledge/results/triples_ports_sieve_20261112/ports/verify.py
python -O 05-knowledge/results/triples_ports_sieve_20261112/ports/verify.py --output receipt.optimized.json
```

The JSON receipts must be byte-identical. They contain the source profile
hash, full difference matrix, all working five-time sets, distance witness,
both distance controls, and mutation census. Checks use explicit exceptions.

OPEN: classify longer genuine returns; determine whether one changes H;
realize useful code couplings by physical tournament operations; turn
continuation-state identification into a target-preserving return.

The subsequent [four-step upgrade](FOUR_RETURNS.md) gives an explicit
finite-local argument and a depth-seven bounded scout; its review status
is recorded separately rather than inferred from the earlier audit.
