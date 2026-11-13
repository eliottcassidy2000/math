# Signed arithmetic, prime targets, and a physical obstruction to halving

2026-11-13, Codex. Statements below are individually **PROVED**, **CITED**,
or **OPEN**. Ordinary Collatz and Chan--Pak remain **OPEN**.

## 1. The target reduction is existential

The inherited [Session193 prime-atom induction, Theorems1--2](../session193_triangulated_proof_20261105/prime_induction/PROOF.md)
and [THM-4094, Hamiltonian matching deficit and two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
reduce the spectrum problem to the lanes

    odd primes p != 7;     7p with prime p>=5, p!=7,

plus the attained finite atoms49,63,343. Ordinal joins physically multiply
Hamilton counts. With the attained seeds1,3,5,23,49,63,67,343, it suffices
to construct a tournament with count N=3h+d whenever h>=3 is attained,
h is not7 or21, d is2 or4, and N lies in either infinite lane.

For an unseeded target N in those lanes the smaller odd parent is forced:
d=2 for N=5 modulo6 and d=4 for N=1 modulo6. The only forbidden-parent
cases that are themselves atoms are23 and67, already seeds. This is a
strong induction, not a conjecture about a single reusable gadget.

It permits replacing the entire parent witness. Thus a correct next
question is: **given this target N, which source representation of its
smaller count h admits the required exit?** It is unnecessarily stronger
to require both exits from every witness, the same marked root forever,
or a finite grammar closed under all outputs.

The conjecture's published formulation is
[Chan--Pak, section6.5, Conjecture6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
which credits the earlier bof/Royle formulation. Our constructions address
that spectrum question; they do not establish the universal premise above.

## 2. Exact coordinate maps, and their types

For odd H put r=(H+1)/2. Then

| count move | r move | physical information needed |
|---|---|---|
| 3H | 3r-1 | ordinal join with C3 always exists |
| 3H+2 | 3r | an actual selected ternary constructor |
| 3H+4 | 3r+1 | an actual selected ternary constructor |
| (H-1)/2 | r/2 | r even, plus an actual count realization |

The first three coordinate identities are inherited, including in
[clone-clock CLOCKS](../clone_clocks_negative_collatz_20261113/CLOCKS.md).
They are useful because they locate the signs exactly. They do not make
all three moves ordinary Collatz: r need not be odd, the middle map is
not its odd branch, and halving has not thereby acquired a graph operation.
For the standard maps and the distinction between arithmetic variants,
see [Lagarias, The 3x+1 Problem: An Overview, section1](https://arxiv.org/html/2111.02635v1).

There is already a genuine complementary minus-sign construction:
[the certified two-clock, CLOCKS sections3--4](../clone_clocks_negative_collatz_20261113/CLOCKS.md)
has H=1+60*3^j*2^k and a single-arc reversal changing (j,k) to(j+1,k-1),
for k>=2. This gives H'=(3H-1)/2. It uses two retained modules and an
explicit guard; neither is recovered from H alone.

## 3. A genuine plus branch, followed by a proved one-module obstruction

Take the new family [HIGHER_MODES, section1](HIGHER_MODES.md):

    H_m = 6*3^m + 6*2^m - 5,              m>=1,
    r_m = 3^(m+1) + 3*2^m - 2.

Every r_m is odd. The one-arc, one-new-vertex +4 replacement gives

    r_m -> 3r_m+1,
    H'_m = 18*3^m + 18*2^m - 11.

So this is literally the plus odd branch in the shifted coordinate,
on an infinite explicitly constructed family. The next arithmetic
halving would require

    Q_m = (H'_m-1)/2 = 9*3^m + 9*2^m - 6.              (1)

**Theorem (PROVED, one-module boundary).** No fixed finite rooted
tournament frame (A,r) satisfies H(A[r<-TT_m])=Q_m for every m>=1.
Nor can a fixed frame satisfy H(A[r<-TT_k])=Q_(m0+qk) for every
sufficiently large k, with fixed m0>=0 and q>=1. Thus an arithmetic
progression cannot rescue a single transitive-module clock by
reparametrizing its length.

**Proof.** In the power-basis expansion of a fixed frame, the greatest
nonzero base R is its greatest nonzero run sector, and its coefficient
is a_R, divisible by(R-1)! by freely permuting internal excursions.
Distinct exponential sequences are linearly independent. Equation(1)
would therefore force R=3 and a_3=9, which is odd, contradicting2|a_3.
On an arithmetic progression the greatest base is3^q and its coefficient
is9*3^m0, again odd. Since3^q>=3, the same contradiction applies. QED.

For the first member, the arithmetic path is

    r: 13 -> 40 -> 20 -> 10 -> 5,
    H: 25 -> 79 -> 39 -> 19 -> 9.

All five individual Hamilton counts are attained. The theorem does
**not** exclude39 or any other particular Q_m. It excludes one uniform
single-module realization of the whole halved family, even on an
arithmetic progression. Repeated halving of the scalar is not the same
as repeated use of the physical construction.

The follow-up [clone-parity theorem](CLONE_PARITY.md) proves that even
finitely many steadily growing transitive modules in a fixed frame
cannot realize the full halved family: their count stabilizes modulo4,
whereas Q_m alternates between1 and3. Residue-dependent changes of frame
or module type remain possible and need a separate constructor.

This gives a precise connection between the external problems: the
power-of-two division in Collatz can remove a factor forced by a free
symmetry action on tournament path decompositions. The extra module in
the known minus-sign construction supplies information absent from the
one-module attempt. It suggests an explicit two-module inverse problem,
not a claim that either conjecture implies the other.

## 4. A borrow identity, and the direction of an induction

The arithmetic identity

    3h+4 = 3(h+2)-2

allows one to exchange a +4 task for a companion source and a -2 exit.
The new h=15*2^m-1 family and the incoming source15*2^m+1 with a -2
child in [Session304 PHYSICAL_FORK](../session304_collatz_peeling_20261113/PHYSICAL_FORK.md)
make this identity concrete. Our one-arc construction supplies the
child directly. The intermediate companion is not asserted to arise
from a universal h->h+2 operation.

For an unseeded target N=3h+4, its companion h+2=(N+2)/3 is still smaller
than N. Therefore a strong-induction argument may use its attained count
when it avoids7 and21, provided an appropriate -2 constructor is proved
for some realization of that companion. This is another target-specific
option, not a replacement for the missing source-selection theorem.

The [Session304 signed descent](../session304_collatz_peeling_20261113/ARITHMETIC.md)
already proves short decreasing paths using both odd-part maps
oddpart(3n+1) and oddpart(3n-1), avoiding7 and21. It does not prove
ordinary plus-only Collatz. To turn that descending arithmetic tree into
a proof of count attainment, one needs **inverse physical lifts from
the smaller attained endpoint back to the target**, with every guard
supplied. A forward map from a hypothetical target to a known count is
insufficient. Equation(1) is a concrete hostile to forgetting this sidecar.

## 5. Connection ledger and next decisive tests

| source -> target | map and preserved predicate | lost data / required sidecar | decisive test |
|---|---|---|---|
| ordered set partitions -> path covers | onto bins yield K_j(TT_m)=j!S(m,j) | endpoint incidences / full skeleton | independent permutation-slots enumeration |
| excursion permutation symmetry -> clock parity | free S_(j-1) action gives factorial divisibility | scalar count hides the action / highest run sector | Q_m forces odd a_3 and is excluded uniformly |
| finite automata -> core grammar | directed rule graph preserves exact family identities | a path may exhaust resources / SCC and parameter increments | no cycle containing +4 in FINITE_PHASE |
| signed Collatz -> count recursion | r=(H+1)/2 preserves exact affine equality | parity, module positions and existence / explicit constructor | plus branch succeeds; first uniform one-module halving fails |
| prime factorization -> target selection | ordinal join preserves multiplication | chosen core and readiness / existential source witness | solve the one required digit at each atom |

The incoming [63 retained-port repairs](../sixty_three_bridge_20261113/README.md)
add a useful independent option: a neutral three-arc repair preserves a
root star, while a paired-module norm has a separate primitive-output
criterion. They do not remove the selection gap here. Their two-module
sidecar is precisely the kind of extra data that our finite-phase theorem
allows, and that a proposed Collatz lift should explicitly retain.

Three focused **OPEN** next tests follow: realize the halved mixed family
with residue-dependent frames or changing module types and a proved
transition; find an unbounded source
selector for the prime/7p targets rather than another finite frame bank;
or prove a companion-source -2 rule on a class whose entry follows from
the smaller-target induction hypothesis. No asserted connection to
unrelated number coincidences is used as a proof dependency.

## Final incoming boundary check

The newly received [Session305 synthesis](../session305_collatz_core_changes_20261113/README.md)
distinguishes its bounded-depth signed arithmetic compiler from the still
missing physical reverse-block compiler. Its selective fusion transports
an existing triangle factor by exchanging transitive and nontransitive
blocks with the same normalized two-cover clock. Such an exchange is
outside our fixed-transitive-block parity obstruction and is a concrete
candidate preparation to investigate; neither packet proves that it
repairs the halved mixed family. Its warning that a third visit can
distinguish equal two-cover clocks supports retaining our full run vectors.

The new [Session306 exact Collatz observer](../session306_relative_isolators_20261113/COLLATZ_LANGUAGE.md)
studies a different finite-state limitation, for prefix-slope recognition.
That language and our graph of exact frame identities are different
objects. Its finite-state boundary does not prove our construction
obstruction or conversely; integer registers and changing physical cores
lie outside the respective finite-state restrictions. These incoming
notes were read for compatibility and next questions, not used as proof
dependencies or treated as a resolution of ordinary Collatz.
