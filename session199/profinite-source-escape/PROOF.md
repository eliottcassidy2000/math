# Congruence-dense physical sources can escape every finite target

Session 199. **PROVED, architecture-specific:** the complete C2 spectrum is
dense in every congruence class compatible with oddness, already with one
fixed triangle source and a primitive unit-interface second source. A single
nonphysical integer profile is a simultaneous all-moduli limit of actual
primitive sources. Its formal C2 value is the genuine C2 gap 67.

This strengthens a finite-modulus obstruction to an all-moduli local/global
failure. It does **not** prove or refute universal coverage by C1/C2/C3.

## 1. The historical connection and what actually transfers

The archived Collatz result THM-4501, file
`01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md`,
constructs genuine positive integer orbit families with every fixed finite
prefix, while their infinite-horizon integer intersection can be empty. Its
last section also proves that a convergent Collatz basin can be dense in the
2-adic integers. The mechanism is exact affine word compilation followed by
congruence selection; it is not a convergence theorem for arbitrary inputs.

The related THM-4654,
`01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md`,
separates two load-bearing ingredients: real contraction puts rational states
in a bounded window, and 2-adic separation then makes finite memory possible.
Its correction explicitly rejects deriving finite memory from real
contraction alone. These are the surrounding themes used here, rather than
an analogy based merely on the words "inverse" or "automaton."

The source-to-target map is:

* Collatz finite parity-word realization becomes a finite modular source word.
* Agreement of arbitrarily long prefixes becomes agreement modulo all moduli.
* The lost datum is the ordinary integer height of the chosen witness.
* The required sidecar is a uniform bound on actual source cover counts.

The decisive test is not whether every residue occurs. It is whether a
specific nonphysical integer profile can pass **every** residue test. Section
3 answers yes, with a coherent sequence of actual witnesses and a quantitative
height escape. Section 5 states exactly the boundedness ingredient that
repairs compactness. No Collatz conjecture is assumed or concluded.

## 2. Physical C2 counts meet every odd congruence class

Use source response

    (D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)),

with exactly one Hamilton path starting at c, and the frozen physical source
and C2 theorems. The actual triangle source is b=(3,1,3,1). Substitution into
the N gluing form gives, for every actual nontrivial second source q,

    N(b,q)=5D-3E+11U-4d.                         (1)

The fixed-unit-interface modular theorem supplies actual strong sources

    q=(D,1,U,1)

with arbitrary independently prescribed odd residues of D,U modulo every
even modulus M. They can be chosen with arbitrarily large D. Every such q
is primitive because d=1. Therefore (1) becomes

    N(b,q)=5D+11U-7.                             (2)

**Theorem 1.** For every positive modulus m and every odd integer r, there
are infinitely many actual strong C2 tournaments, built with the same left
triangle source b, whose Hamilton-path counts are r modulo m. Both source
profiles can be required to be primitive.

**Proof.** Work modulo M=2m. Write D_0=1+2x and U_0=1+2y. Equation (2) is r
modulo M exactly when

    5x+11y=(r-9)/2 (mod m).

Since gcd(5,11)=1, choose integer x,y solving the displayed equation even
before reduction. The modular interface realizes their odd residues by one
actual q. The physical N construction then has C=2 and the required count.
The compiler gives arbitrarily large D, and (2) is at least 16D-7 because
U>=D. Thus its counts are unbounded and hence infinitely many are distinct.
Both b and q are primitive. QED.

This also gives simultaneous solvability for any finite collection of moduli,
by using their least common multiple. It is an all-order constructive theorem,
not extrapolation from a finite residue table.

## 3. One fixed phantom profile passes every modulus

Let

    p=(1,0,6,1).                                (3)

It is primitive and satisfies the elementary physical cone and parity
conditions. Nevertheless it is not an actual source. Indeed H=1 forces a
transitive tournament, and a source-marked transitive tournament has

    (D,E,U,d)=(1,0,2^k,1), k>=1,

or the singleton (1,1,1,1). One proof of H=1 implying transitivity is to
insert the remaining named vertices deterministically into each of the
three distinct paths of any directed triangle. Deleting the inserted names
recovers each old path. Thus any nontransitive tournament has H>=3.

**Theorem 2.** For every positive m there are actual primitive weak sources
q_m with E=0 and d=1 exactly and

    q_m=p (mod m) coordinatewise.               (4)

**Proof.** Apply the strong modular interface at M=2m to D=1,U=3 modulo M,
with E=d=1 exactly. Append one universal sink. The physical operation
O(D,E,U,d)=(D,0,2U,d) gives (4), in fact modulo 2m. QED.

Choosing m_j=24 j! gives a coherent sequence of residue approximations:
for every fixed modulus h, all sufficiently large q_(m_j) equal p modulo h.
Thus p lies in the profinite closure of the physical profiles, while p is
an ordinary integer point outside that set. The individual q_(m_j) need
not be initial segments of one source word, nor is that asserted.

Consequently even the infinite family of **all** modular membership tests
cannot recognize physical source profiles. This is stronger than the earlier
statement that no fixed finite list of tests suffices: here the same p passes
every modulus, not a different false positive for each finite list.

The escaped value has a concrete scalar consequence:

    N(b,p)=5+11*6-4=67.                          (5)

There is no physical C2 graph with H=67. The complete C2 inverse bounds both
source cover counts by floor((67-1)/6)=11. Its complete physical bank has 14
nontrivial profiles. Exhausting N and M on these profiles, and the exceptional
15, gives the support through 67 exactly

    {13,15,23,25,29,41,43,45,49,53,55}.          (6)

The independent verifier regenerates the finite bank using O/Z/J and tests
every ordered pair, without a graph-order cutoff. This is a finite-exact
calculation with a proved all-order cutoff. The separate C3 construction of
67 is already known; no unrestricted spectrum gap is claimed.

There is also a real-height lower bound on these particular witnesses.
Suppose m is divisible by 24 and q=(D,0,U,1) is any actual source satisfying
(4). If D=1 then U is a power of two, impossible since U=6 modulo 3. Hence
D>=m+1. A weak unit-deletion source has a strong or singleton initial core
followed by at least one sink. Its core has cover count at least D, so
U>=2D. Therefore

    U>=2m+2,
    N(b,q)=5D+11U-4>=27m+23.                    (7)

In particular the witnesses of (4) must escape to infinity in ordinary
height, although their coordinatewise modular limit is the fixed small
integer tuple (3). The large multiple of m in the scalar readout does not
disappear merely because the congruence error vanishes at every modulus.

## 4. A very sparse actual family can have the same congruence closure

**Theorem 3.** There is an explicitly constructible increasing sequence of
actual C2 counts h_i, all from the fixed-left-source construction in Theorem
1, that meets every odd-compatible residue class modulo every modulus
infinitely often and satisfies

    #{i:h_i<=X} <= floor(log_2 log_2 X), X>=2.   (8)

**Proof.** Enumerate all pairs (m,r), m>=1 and r an odd representative modulo
2m, with each pair repeated infinitely often. At stage i>=1 use Theorem 1's
effective compiler for that pair, increasing its neutral-loop repetitions
until

    h_i>h_(i-1),    h_i>=2^(2^i).

The compiler has unbounded actual D, so this finite search terminates. These
are literal physical C2 constructions, not freely assigned integer values.
The scheduled pair guarantees all the residue visits. If h_i<=X then
2^i<=log_2 X, proving (8). QED.

Thus even perfect all-moduli support together with actual graph recipes and
primitive inputs does not force positive density, let alone exact coverage.
The theorem is about a selected subfamily, not the density of the whole C2
spectrum. Its role is to expose precisely what a congruence-only finishing
argument would forget.

## 5. The exact bounded-height repair

Let S be the set of actual nontrivial source profiles, and let

    C2={15} union N(S,S) union M(S,S).

For a fixed integer n, the following are equivalent:

1. n belongs to C2.
2. There is a constant B such that, for every positive integer m, either
   15=n modulo m, or there are actual v,w in S with U_v,U_w<=B and
   N(v,w)=n modulo m or M(v,w)=n modulo m.

The forward direction uses one exact representation, or 15. For the reverse
direction, the bounded inventory S_B is finite: its other coordinates obey
1<=d<=D<=U<=B and 0<=E<=D. Thus the set of possible readouts in (2) is finite.
If n is absent from that finite set, choose m larger than every nonzero
absolute difference between n and those readouts. No readout is then n
modulo m, a contradiction. This proves the equivalence.

This is a compactness statement with a **uniform Archimedean bound**, not a
new source-existence claim. The modular interface supplies no such bound as
m grows; (7) demonstrates its necessary failure at a genuine C2 gap. The
source-word halving descent from Session 198 contracts an already given
representation. It does not give the missing uniform bound for a sequence
of congruence witnesses chosen independently at larger moduli.

For the target C1/C2/C3 union, the same finite-inventory argument applies
after adding the finitely many readout families and boundary seeds. A proof
that all odd target n outside {7,21} admit such uniformly bounded actual
witnesses would finish coverage. That uniform-height assertion remains OPEN.

## 6. Reproduction and scope

Run `python3 verify.py` and `python3 -O verify.py`; compare their JSON receipts.
The verifier independently implements the modular word compiler, both binary
evaluators, raw source moves, source-word reconstruction, C2 formulas, and
literal subset path/cycle counters. It checks all odd residue pairs through
modulus 16, scalar residue witnesses through modulus 32, the coherent phantom
examples, the complete C2 exclusion at 67, and small literal graph controls.
It does not use finite checks to infer all-moduli density or sparsity.

The frozen source and modular proof snapshots are under `dependencies/`.
The two historical dynamics snapshots give mechanism provenance only; all
new source statements are proved above without a Collatz or Mahler premise.
No literature-priority claim is made. The deliverable is a precise
local/global obstruction and its bounded-height repair, not the requested
universal odd-coverage proof.
