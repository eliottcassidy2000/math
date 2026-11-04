# The twelve-vertex module-converse residual must be prime

**PROVED** all-order inheritance and twin obstructions C1-C3.
**PROVED with FINITE-EXACT dependencies** the prime-core reduction C4.
**FINITE-EXACT** the explicitly prescribed twelve-vertex scout.
Universal two-port reconstruction remains **OPEN**. No complete census of
order twelve, exact collision, or literature-priority claim is made.

## Inheritance and the observer retained

The closest proved mechanism is [B3-B4, ambient module congruence
descent](../../two_port_seven_core_20261103/algebra/PROOF.md). A module M
whose internal reversal is invisible in an arbitrary ambient quotient
must satisfy the intrinsic necessary conditions

    S(M)=E(M) mod2,       C(M)=C(M)^t mod4.           (1)

The [complete order-ten sieve and independent
review](../../two_port_seven_core_20261103/residual/INDEPENDENT_REVIEW.md)
exclude their conjunction in all even orders through ten; the odd-module
current excludes every odd order>=3. Even self-converse cores are also
excluded at all sizes by endpoint parity. Twelve is therefore the first
remaining possible core size, not an established collision size.

The hostile here is a six-vertex tournament with symmetric C modulo4 but
unequal endpoint parity. It disproves the tempting shortcut that C modulo4
alone always separates a converse. The least-used sidecars are uniform
submodules, their contracted quotients, and exact differences between
endpoint counts at a pair of twins. The live concepts are module runs,
intrinsic modular residues, exact endpoint order, full named three-path
response, and finite versus universal enumeration.

All vertices remain named; bit v in row u means the actual arc u->v.
S_u and E_u count directed Hamilton paths starting and ending at u.
C_uv counts ordered spanning covers by two nonempty directed paths with
the first ending at u and the second starting at v. These include absent
arc gradient directions. A module is a nonempty set whose vertices have
identical incidences toward every exterior vertex. A prime tournament
has no proper module of size between two and n-1.

The full two-port function is the inherited bidegree-(2,2) polynomial

    F_T(a,c)=[V]{(1+e a)(1+(1-a)^t G c)(1+(1-c)^t s)
                  +(1+e c)((1-c)^t G a)(1+(1-a)^t s)},

where G is the named squarefree path matrix, e=1^tG and s=G1. Its
four-index three-path tensor is D_(i,j,k,l)=[V]e_i G_jk s_l; a singleton
middle path j=k is included and supports are disjoint. This packet uses
necessary information recovered from that full function, without asserting
that H,S,E,C alone suffice. Actual full-response controls are routed in
the surrounding packet; a modular survivor would still require those
stronger tests and would not itself be a collision.

## C1. Two different kinds of inheritance

Let T=Q[M], replacing a distinguished vertex z of Q by a nonempty module M.

**Theorem C1(a).** If C(T)=C(T)^t modulo4, then
C(M)=C(M)^t modulo4.

This is the unconditional transfer identity from B4:

    C_uv(T)-C_vu(T)
       = H(Q)(C_uv(M)-C_vu(M)) mod4,   u,v in M.

Its proof uses H(T)=sum_k c_k q_k(M), c_1=H(Q) odd, and the fact that
the directed derivative skew of q_k is divisible by2k!. Thus all higher
terms vanish modulo4 and the leading multiplier is a unit. No equality
between T and a modified ambient tournament is needed for this implication.

**Theorem C1(b).** If S(T)=E(T) modulo2, then S(Q)=E(Q) modulo2.

**Proof.** For each exterior vertex w, its fixed-start or fixed-end count
has a uniform-run expansion sum_k d_k q_k(M), with integer exterior
skeleton coefficients. For k>=2 the ordered-cover polynomial q_k is
divisible by k!, hence even, whereas q_1=H(M) is odd. The one-run
coefficients are S_w(Q) and E_w(Q). Consequently S_w(T)=S_w(Q) and
E_w(T)=E_w(Q) modulo2 for every w!=z. The hypothesis balances these
coordinates of Q. Since sum_w(S_w(Q)-E_w(Q))=0 exactly, its remaining
coordinate z is balanced too. QED.

These directions differ: the C-skew condition passes to the **module**,
while endpoint parity passes to the **quotient**. We do not assert that
both conditions automatically pass to each induced subtournament, or
that every modular survivor is prime at arbitrary order.

## C2. An exact ordered-cover sign recovers the orientation of twins

Suppose {a,b} is a module with a->b, and let Q be its contraction to z.
The exact endpoint identities are

    S_a(T)-S_b(T)=S_z(Q),
    E_b(T)-E_a(T)=E_z(Q).                            (2)

To prove the first, swap a,b in every Hamilton word starting at a in
which they are not consecutive. Uniform exterior incidences make this a
bijection to words starting at b, since b cannot be followed by a. The
unpaired words start with a,b; contracting this initial block gives
exactly the quotient words starting at z. The end identity is dual.

The stronger identities needed here concern ordered two-path covers. For
every exterior vertex w,

    C_bw(T)-C_aw(T)=C_zw(Q),
    C_wa(T)-C_wb(T)=C_wz(Q).                         (3)

For the first, swap the twins in covers counted by C_bw unless they are
consecutive. Every resulting cover ends its first path at a and starts
its second at w; all such covers occur. The unmatched covers have block
a,b at the end of their first path. Contraction gives exactly C_zw(Q).
The second identity is dual: the unmatched block a,b starts the second
path. Twins in different paths are part of the swap bijection. No
exterior arc changes inside a swapped cover, and no factorial is lost.

**Theorem C2.** Define the named contrast

    L_(a,b;w)(T)=C_bw-C_aw+C_wa-C_wb.

For twins a->b it satisfies

    L_(a,b;w)(T)=C_zw(Q)+C_wz(Q)=2w_zw(Q)>0.         (4)

Here positivity is Session176's [Positive cut circuit theorem, section1](../../session176_tournament_information_20261104/circuits/PROOF.md):
its K(u,v) is the present C_uv, and w_uv=(C_uv+C_vu)/2 is a positive integer.
For the same named twins oriented b->a in another tournament T', the
same contrast is strictly negative. The two quotient tournaments may
differ arbitrarily: each of their cut weights is positive. Consequently
two tournaments with equal ordered C cannot have the same named twin
pair oriented oppositely, even if other arcs also change.

The same obstruction survives a positive projective scale: full two-port
tables that are positively proportional have proportional C by the
inherited [Session27 linear tomography, section7](../../session27_repunits_lattices_20261022/session27-response/proof.md),
specifically C_uv=(B_uv+J_uv)/2 with B,J linear in the measured response
values, so they have the same sign of L. Its H and
one-port inputs are themselves two-port queries with universal source or
sink padding. Thus
opposite twins cannot have positively proportional full responses either.
This is a partial projective-faithfulness result; it does not assume that
response tables have primitive integer content. Their nontrivial content
was quantified in [B1-B2](../../two_port_seven_core_20261103/algebra/PROOF.md).

There is an exterior w whenever the tournament has at least three
vertices. At order two, the endpoint identities (2), with singleton
quotient counts one, separate the two orientations. Thus the complete
two-port response separates an opposite twin pair at every order.

In particular, if a module M contains twins, reversing all of M inside
any fixed larger quotient reverses a twin pair in both full tournaments.
Its external incidences remain uniform, so the sign obstruction applies
directly in that ambient tournament. No exact endpoint equality of M,
strong connectivity, or descent of the full response to M is assumed.

## C3. The directed-triangle quotient retains module endpoint parity

**Theorem C3.** Let T=M->a->b->M, with a,b singletons. For u in M,

    S_u(T)=S_u(M)+B_2(u),
    E_u(T)=E_u(M)+E_2(u),                            (5)

where B_2,E_2 count ordered covers of M by two nonempty paths, the first
fixed slot starting/ending at u. These follow by separating one M-run
from the two-run word with exterior segment a,b; permuting fixed slots
gives the stated end formula without a multiplicity factor.

The inherited squarefree parity identity e=Us, P=U-1 gives
e_uP=s_uP=e_u+s_u over F2. Extracting the full support therefore yields

    B_2(u)=E_2(u)=S_u(M)+E_u(M) mod2.

In particular (5) gives the swapped parities

    S_u(T)=E_u(M),   E_u(T)=S_u(M) mod2.             (6)

If T has balanced endpoint parity, then so does M. This is a special
quotient identity; C1(b) alone would not justify this direction.

## C4. A genuine twelve-vertex converse collision must have a prime core

The complete inherited representative bank through order eight is pinned
by SHA256 in both checkers. Independent residue evaluation here gives:

| Core order | Classes | C=C^t modulo4 classes |
|---:|---:|---:|
|1|1|1|
|2|1|1|
|3|2|0|
|4|4|0|
|5|12|0|
|6|56|1|
|7|456|0|
|8|6880|0|

The unique six-vertex class in the last column has code22693 and rows

    (10,20,33,22,37,11),
    H=45,
    S=(7,7,7,8,8,8), E=(8,8,8,7,7,7).

Its C is symmetric modulo4, but its endpoint parity differs in every
coordinate. This is a concrete **REFUTATION** of a C-mod4-only universal
converse obstruction. The strongest survivor is the joint condition (1).

The six-core has three directed twin pairs,5->0,3->1 and4->2; it is an
inflation C3[TT2,TT2,TT2]. Each pair's two vertices have identical
incidences to the other four. These twins remain twins under every further
uniform substitution into a larger tournament. Consequently C2 excludes
this six-core as a submodule of any hidden module-converse candidate,
at every ambient size, without a quotient enumeration.

There are48 endpoint-parity-balanced seven-vertex quotient classes.
Replacing each of their seven possible distinguished vertices by that
unique six-vertex core gives336 rooted expansions. Every one has
**unequal ambient endpoint parity**. This is a complete finite check of
all twelve-vertex candidates with a six-vertex module that could pass
C1: the module class is forced by C1(a), and the quotient class is forced
by C1(b). Considering every root covers all rooted quotient classes,
with harmless duplicate rooted orbits. Arbitrary renamings of the module
do not change either invariant under test. This336-case check corroborates
the structural exclusion; it is not load-bearing in C4 because the twin
sign already excludes the unique six-core in every ambient tournament.

The same reasoning gives an all-order small-module exclusion: a core
whose converse reversal is hidden in any ambient tournament has no
proper nontrivial submodule of any size2 through9. Odd sizes fail C1(a),
orders4/8 fail the complete residue table, and orders2/6 contain an
oppositely oriented twin pair detected by C2. More generally every
remaining proper nontrivial submodule must have even size at least10.

**Theorem C4.** Let M be a twelve-vertex module in an arbitrary ambient
tournament T=R[M]. If T and R[M^op] have equal full two-port response,
then the intrinsic core M is prime. There is no bound on the ambient size.

**Proof.** The inherited B3-B4 descent gives the intrinsic modular
conditions (1) for M. Suppose a proper nontrivial module N of M exists.
By C1(a), C(N) is symmetric modulo4. Odd |N|>=3 is impossible by the
all-order odd-current theorem. Orders4 and8 are excluded by the complete
table above. At order6, N must be the unique listed class, whose twins
are excluded by C2 directly in the actual ambient tournament. For |N|=10,
the quotient has order3 and balanced endpoint
parity by C1(b), so it is the directed triangle. C3 then gives balanced
endpoint parity inside N, contradicting the
[complete order-ten joint sieve](../../two_port_seven_core_20261103/residual/full10.json),
whose [independent review](../../two_port_seven_core_20261103/residual/INDEPENDENT_REVIEW.md)
records all9,733,056 classes,12,670 endpoint-parity candidates and zero
joint survivors. Finally |N|=2 is a twin pair in both ambient tournaments, reversed
by M->M^op. It is impossible by the exact ordered-cover sign theorem C2.
This exhausts2<=|N|<=11. QED.

A weaker statement applies to the intrinsic congruences alone: any
order-twelve core passing (1) can have no proper nontrivial module except
a two-vertex module. The exact contrast (4), read in the actual ambient
tournament, removes that last case for a hidden module reversal. This
distinction is essential: the proof does not assume that equality of the
entire response descends unchanged through an unknown exterior scale.
It uses only established modular descent and an actual ambient sign.

Universal module-converse separation at order twelve remains OPEN. C4
reduces a possible counterexample to a prime, even, nonselfconverse core;
it supplies neither a prime-core separation proof nor an exact collision.

## Exact twelve-vertex scout and retained hostile controls

The fixed scout tests311,872 specified labelled constructions, with
duplicates permitted. It is not a class census. Its complete named
families and deterministic streams are:

| Family | Tested | S=E modulo2 | Joint survivors |
|---|---:|---:|---:|
|All labelled six-vertex quotients, every vertex replaced by TT2|32,768|0|0|
|All labelled four-vertex quotients, modules TT2,TT2 and two order-four class representatives|1,024|0|0|
|All labelled three-vertex quotients, three order-four class representatives|512|0|0|
|Ordered pairs of all order-six class representatives|3,136|0|0|
|All births over TT11|2,048|0|0|
|All births over cyclic11|2,048|0|0|
|All births over Paley11|2,048|0|0|
|Specified fixed-involution self-converse12 stream|2,048|0|0|
|One specified edge flip from each paired source|2,048|1|0|
|Two specified edge flips from each paired source|2,048|2|0|
|Specified generic twelve-vertex stream|262,144|143|0|

The generic seed is0x0123456789abcdef, and the paired seed is
0x13579bdf2468ace1. Both use the exact uint64 xorshift13,7,17 operation
in the source, consuming one output per prescribed edge choice. Random
here names a deterministic generation rule, not a probabilistic proof or
a uniform class sample. The second flip may repeat the first; such cases
remain part of the declared stream. Final generator states are frozen.

All146 endpoint-parity candidates receive full integer H,S,E,C and a
separate residue-only modulo4 DP comparison. They all fail the second
condition. Additional regularly spaced and initial stream controls audit
the fast parity predicate against integer DP even when it rejects.

The sparse retained twelve-vertex hostile has rows

    (3664,3173,1185,839,1422,1113,2068,1899,615,3190,328,1468).

It has H=289025 and balanced endpoint parity, but24 unordered pairs have
nonzero C skew modulo4. It is index23714 in the prescribed generic stream;
this is the minimum such skew support among this scout's146 parity
candidates, not a global extremal claim. The complete integer sidecar is
saved in `normal.json`.

Another retained control, used by the parent full-three-path audit, is

    (3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002),

index1473 in the two-flip paired stream, with28 nonzero unordered skew
entries. Its full response and its converse must be compared using the
actual tensor and attachment battery; neither its endpoint parity balance
nor this finite scout is a collision claim.

Both retained twelve-vertex controls are prime. The module checker
exhausts all4,082 vertex subsets of sizes2 through11 separately for each
control and finds no module. The rows and subset counts are frozen in
the module receipts. These are actual prime-core controls, though both
already fail the ordered-cover congruence.

## Reproduction and validity

Run `python run.py` and `python -O run.py` for the full specified scout;
run `python module_check.py` and `python -O module_check.py` for the module
exclusions. Sources and the inherited compressed bank are SHA256 pinned
before and after execution. Compilers and temporary banks/executables
live outside the repository. All gates use explicit exceptions.

`scout.cpp` implements the exact endpoint recurrence both as a bit-parity
DP and an integer DP, plus a separately instantiated residue-only mod4
endpoint-and-cover DP. The latter shares the recurrence implementation
with the integer DP; it is an arithmetic cross-check, not an independent
algorithmic implementation.
At order12 every fixed-end count is at most11!, and every fixed C_uv
also is at most11!: concatenating its two paths with the distinguished
u,v boundary gives a permutation with u immediately before v. Counts and
products are accumulated in uint64, safely above all these bounds.

The module checker verifies75 exact twin endpoint identities and the
corresponding cross-cover identities and orientation signs (all labelled
quotients through order4),1,099 triangle parity identities (all labelled
module cores through order5), and all336 forced six-module candidates.
It also checks all8,164 proper nontrivial subsets of the two saved
twelve-vertex controls for primality.
Its normal and optimized receipts are `modules.json` and
`modules.optimized.json`; the scout receipts are `normal.json` and
`optimized.json`. Finite bank completeness is inherited and explicitly
separate from the all-order arguments C1-C3 and the bounded corollary C4.
