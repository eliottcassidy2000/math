# One-exit attachments: the missing penalty and two positive fixed gadgets

**PROVED IN NOTE:** the general one-exit attachment formula; necessary and
sufficient moment equations for a fixed gadget to produce a uniform return
from the incoming Q_k sources; all-order exclusions for multipliers 27 and
81 within this specified grammar; and two actual uniform constructions with
multipliers 375 and 377. **FINITE-EXACT, GLOBAL VIA INHERITED INVENTORY:**
the multiplier-243 gate has no solution at any order. **FINITE-EXACT:**
the complete inherited ordinary bank through order seven and all selected
cuts, plus the literal controls below. **OPEN:** an arbitrary nonzero-digit
return and Chan–Pak.

This note follows the successful native +2 construction in [PROOF.md](PROOF.md).
Its purpose is to retain the useful generalization of transitive inflation
without losing the interaction between the two boundaries of a path run.

## 1. The physical operation and its four source coordinates

Let X be a tournament with marked vertex x. Put

    h = H(X),
    F = sum_P H(X minus vertices(P)), over nonempty paths P ending at x,
    s = number of full X paths starting at x,
    d = H(X minus x), with H(empty)=1.

Let G be a nonempty disjoint tournament, and choose A contained in V(G).
Keep all internal arcs. Direct every comparison from G to X, except that
x points to every vertex of A. Thus x is the only possible exit from X
into G. Denote the resulting tournament by T(X,x;G,A).

For the gadget, write H=H(G) and define

    E = number of full G paths ending in A,
    U = number of ordered spanning two-path covers (R,S)
        with S nonempty and starting in A, while R may be empty,
    W = number of ordered spanning two-path covers (R,S)
        with R,S both nonempty, R ending in A and S starting in A.

The empty complementary path contributes one. Covers partition the vertices;
they impose no arc between the two paths.

**Theorem 1 (one-exit formula).** For every X,x,G,A as above,

    H(T(X,x;G,A)) = hH + FU - sE - dW.                  (1)

Proof. A full path has at most two G-runs, since every return from X to G
must use x. If it starts with its only G-run, the count is hH-sE: the
forbidden boundary is a G path ending in A followed by an X path starting
at x. If its only G-run follows x, the count is F times the number of full
G paths starting in A. If it has two G-runs, it has the form

    R → P → S → Q,

where R,S are nonempty G paths covering G, S starts in A, P ends at x,
and Q covers the remaining X vertices, possibly vacuously. The unrestricted
count is F times the corresponding nonempty two-cover count of G. The
boundary R→P fails precisely when P is the singleton x and R ends in A.
For each such G cover there are d choices of Q, so this subtracts dW.
Combining the two terms involving F gives FU. These cases are disjoint and
exhaustive. QED.

**Derivation boundary.** The subtraction dW is essential. Take X to be a
singleton, G the transitive three-vertex tournament, and A its first two
vertices. Then (h,F,s,d)=(1,1,1,1) and (H,U,E,W)=(1,6,0,2). Formula (1)
gives 5. Omitting dW would incorrectly predict the impossible count 7.
The uncorrected shortcut was caught before it entered a maintained result.

The [endpoint capacity theorem](../../session258_endpoint_seed_20261111/endpoints/PROOF.md),
Theorem 6, gives U=E+2V(A), W=2I(A), where I(A) is the internal selected
capacity and V(A) the selected row volume. Since the capacities are
nonnegative and I(A)≤V(A),

    U ≥ E+W.                                          (2)

This is the needed sidecar: selecting a large binary volume also selects
internal capacity, which incurs the double-boundary penalty.

## 2. Exact uniform comparison for the Q_k source family

The inherited [Q_k construction](../../berggren_target_return_20261111/reset/TARGET_RETURN.md)
has, for D=2^k and k≥1,

    h=297D+71,    F=297D+72,
    s=105D+24,    d=105D+25.                         (3)

Fix G and A independently of k. We ask for the actual target

    H(T(Q_k,x;G,A)) = λ(h+1)-1
                     = 297λD+72λ-1                 (4)

for every k≥1. This is a coefficient comparison, not a numerical fit.

**Theorem 2 (complete fixed-gadget gate).** Equation (4) holds for every
k≥1 if and only if there is an integer t≥0 such that

    H = 45t+1+E,
    U = λ-10t-1-E,
    W = 99t-E.                                     (5)

Proof. Equating the coefficients of D in (1) and (4) gives

    99(H+U-λ)=35(E+W).

Since 99 and 35 are coprime, H+U-λ=35t and E+W=99t for an integer
t≥0. Equating the constant coefficients then gives H+W=144t+1.
Solving these three equations gives (5), and reversing the calculation
proves sufficiency. QED.

By (2), every such gadget satisfies

    λ ≥ 109t+1+E.                                  (6)

Also H is odd, so E has the same parity as t. These are consequences of
the exact gate, not replacements for any of its equations.

**Corollary 3.** No fixed tournament G and selected set A in this one-exit
grammar realizes (4) for λ=27 or λ=81.

Proof. Inequality (6) forces t=0. Equations (5) give E=W=0, H=1 and
U=λ-1. A tournament with exactly one Hamilton path is transitive: a
nontransitive tournament contains a directed triangle, and vertex insertion
injects its three Hamilton paths into paths of the whole tournament.
In a transitive tournament, W=0 forces |A|≤1; every pair of selected
vertices contributes positively to W. For a singleton selected vertex,
U is a power of two. Neither 26 nor 80 is a power of two. QED.

The scope is fixed G,A, independently of k, with the specified one-exit
attachment. This does not exclude other physical sources of the target
count, larger grammars, k-dependent gadgets, or same-count rebases.

For λ=243 the possibilities are exactly the following two moment systems:

| t | H | U | W | parity of E |
|---|---|---|---|---|
| 1 | 46+E | 232-E | 99-E | odd |
| 2 | 91+E | 222-E | 198-E | even |

The t=0 case again fails because 242 is not a power of two; inequality
(6) excludes t≥3. These equations remain an explicit positive construction
target at all orders.

### A rooted three-coordinate constructor criterion

There is an equivalent inverse target that uses the existing rooted
tournament language directly. Adjoin one vertex z to G, with z→A and
G-A→z, and mark z. Write (h_z,e_z,u_z) for its Hamilton count, full
endpoint count at z, and weighted start-cover count at z. Then

    h_z=H+U-E-W,    e_z=H-E,    u_z=H+U.              (10)

The first identity is Theorem 1 with a singleton X. For the second, a
path ending at z must first traverse G and end outside A. For the third,
the distinguished path starting at z is either the singleton z, contributing
H, or begins z→A and then traverses a nonempty G path; its complementary
G path is unrestricted, giving U.

Consequently the fixed-gadget gate is equivalent to finding any rooted
tournament with the three exact coordinates

    (h_z,e_z,u_z)=(λ-64t,45t+1,λ+35t).               (11)

Conversely, delete z from any such rooted tournament and let A be its
out-neighborhood. Equations (10) recover (5), so it supplies the actual
one-exit constructor, with no extra deletion-count requirement. This is
a necessary and sufficient constructor criterion, not a loss of the
double-boundary sidecar: that information is retained in u_z-h_z=99t.

In particular λ243 asks for either

    (h_z,e_z,u_z)=(179,46,278) or (115,91,313).       (12)

The second alternative is impossible at every order. For any root z with
an out-neighbor, inject paths ending at z into paths not ending there:
remove z and insert it immediately before the first vertex it beats.
The preceding vertex, if any, beats z by minimality, so the new path is
valid. Moving z back to the end recovers its source, proving injectivity.
Thus e_z≤(h_z-1)/2. If z is a sink, then e_z=h_z and u_z=h_z. The target
(115,91,313) satisfies neither alternative. More generally every positive-t
solution of (11) must have a nonsink root and

    λ≥154t+3.                                      (13)

Therefore **the single remaining all-order λ243 target is (179,46,278)**.

The positive λ375 and λ377 gadgets correspond respectively to rooted
eight-vertex states (311,46,410) and (313,46,412). Thus existing rooted
compilers can attack the surviving target directly instead of enumerating
arbitrary cuts. Section 6 below decides that target globally using the
inherited exact multiplicity fibre.

For any source with F=h+1 and d=s+1, the whole attachment formula also
compresses to

    H(T)=(h-s)u_z+(s+1)h_z-e_z.                     (14)

Indeed, writing γ=H(G), recover U=u_z-γ, E=γ-e_z and
W=u_z-γ+e_z-h_z from (10), then substitute into (1). Every coefficient
of γ cancels. This is a universal three-coordinate constructor for that
marked source class, including the inherited strong and weak 1-ready
sources; the gadget's deletion count is genuinely unnecessary.

Thus changing u_z by two at fixed h_z,e_z changes the target by 2(h-s),
not by two. A small gap in the rooted coordinate must retain this actual
target multiplier when used in another operation.

## 3. Why a nonuniform transitive cut still fails

Order TT_m from first to last. Its selected volume is

    U(A)=sum_(i in A) 2^(m-i-1),

and E is one precisely when the last vertex is selected. But W is not
generally zero. For a pair of selected vertices i<j, its contribution
K_ij+K_ji equals 2 if j-i=1 and 2^(j-i-1) if j-i≥2. Hence W=0 if and
only if at most one vertex is selected.

The tempting binary cut for multiplier 27 uses TT_5 and A={0,1,3}, giving
(H,U,E,W)=(1,26,0,8). The actual all-k count is

    27(h+1)-1 - 8d = 7179D+1743,                    (7)

not 8019D+1943. At k=1 these values are 16101 and the intended 17981.
The literal control counts 16101. The failure comes from exactly eight
selected two-boundary covers, each costing H(Q_k-x), rather than an
unspecified limitation of binary encoding.

## 4. Two positive nontransitive fixed gadgets

The complete inherited ordinary bank of orders 1 through 7, with every
selected subset, contains the following two solutions of the positive-t
coefficient-cancellation gate. Outmasks use vertices 0 through 6.

| multiplier λ | G outmasks | A | (H,U,E,W) |
|---|---|---|---|
| 377 | (86,124,40,113,100,65,4) | {0,2,3} | (65,347,19,80) |
| 375 | (86,124,24,113,96,69,4) | {0,2,5} | (67,343,21,78) |

Both have t=1. Substituting the displayed moments into (1) and using (3)
proves, for every k≥1, actual strong tournament families with

    H(T_377,k)=111969·2^k+27143 =377(h+1)-1,
    H(T_375,k)=111375·2^k+26999 =375(h+1)-1.          (8)

Both G are strong, A is nonempty, and Q_k is strong. A path from Q_k can
reach x, cross to A, and reach every G vertex; G reaches Q_k-x directly.
These facts prove strongness of the combined graph. Each construction adds
seven vertices and preserves every old vertex and arc of Q_k.

The mechanism applies to the broader marked source class

    F=h+1,    d=s+1,    99s=35h-109.                 (9)

For the first gadget, (1) simplifies to 412h+267-99s=377h+376. For the
second it simplifies to 410h+265-99s=375h+374. Thus the shared endpoint
relation cancels the large selected-set penalties exactly.

There is also a physical sibling comparison. The gadgets differ in the
two arcs {2,4},{2,5}; their selected sets exchange 3 and 5, changing the
two comparisons {x,3},{x,5}. Under (9), these four arc changes alter H
by exactly 2(h+1). This is an actual pair of tournaments with a common
embedded X, not addition of independently realizable integers.

The literal k=1 counts are 251081 and 249749. In each graph the only
indegree-one vertex is G vertex 1, with predecessor G vertex 0. Its
weighted endpoint defect V-H is respectively 16960 and 12573. Thus
these particular rooted graphs are not a 1-/2-ready return. The positive
uniform counts and four-arc comparison stand independently of that next
interface condition. No ready-class iteration is asserted.

## 5. Finite universe, controls, and remaining gate

[one_exit_probe.py](one_exit_probe.py) reads the pinned inherited bank at
commit 63b86513f8be and expands converse-merged records using their ordinary
multiplicity field. The universe consists of 1,1,2,4,12,56,456 ordinary
representatives in orders 1 through 7; every subset A is tested. Bank
completeness is inherited, rather than re-enumerated here.

For each graph, endpoint subset DP computes exact full-path counts and
ordered two-cover moments. The script records all positive-t solutions
of (5), allowing every λ. Exactly the two gadgets in section 4 occur.
In particular no λ243 gadget occurs in this finite universe. This says
nothing about orders eight and above.

Independent literal whole-graph path counts test formula (1) for all rooted
sources and gadgets from the ordinary bank through order three and every
cut. Additional controls are the singleton hostile, the TT_5 failed binary
cut, and both actual 18-vertex Q_1 extensions. The receipt includes literal
rows, source-bank hash, complete universe size, and root defects.

Reproduce from the worktree root:

    python3 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/returns/one_exit_probe.py

The output is [one_exit_receipt.json](one_exit_receipt.json).

## 6. The prime fibre makes the multiplier-243 decision global

The all-order strong inventory in
[THM-4744, interval chains and the spectrum below 625](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md)
contains exactly 39 ordinary strong classes at H=179: 26 of order eight,
nine of order nine, and four of order ten. Its completeness is inherited
here. We use this existing mathematical object, rather than enumerate a
larger arbitrary gadget bank.

**Prime-prefix reduction.** Every tournament with H=179 has exactly one
nontrivial strong component S, with H(S)=179; its other strong components
are singleton vertices. This follows from multiplicativity across the
ordered strong components and primality of 179. Therefore it has the form
TT_a→S→TT_b. A root with endpoint count 46 cannot lie before the last
component, where its endpoint count would be zero, or be the final
singleton, where its endpoint count would be 179. Thus b=0 and the root
lies in S.

Removing the prefix TT_a preserves both E and U_start. Full paths traverse
the prefix in its unique order before the full S path, preserving E. In
a two-path cover whose distinguished path starts at the root, no prefix
vertex can follow that root. All prefix vertices therefore lie, uniquely
ordered, in the other path. Deleting them is a bijection with the covers
of S, including the case in which the other S path is empty. Hence the
rooted target (179,46,278) exists at any order if and only if it occurs at
one of the 329 roots of these 39 strong representatives.

The exact root check finds just two roots with E=46:

| order | inherited decimal code | root | (H,E,U_start) |
|---|---|---|---|
| 8 | 167805088 | 1 | (179,46,232) |
| 9 | 18254201090 | 1 | (179,46,276) |

No root has the required U_start=278. Consequently **there is no fixed
one-exit gadget for multiplier 243 at any order**, by the exact gate and
the inherited complete strong fibre. The closest displayed coordinate
miss is two in U_start, but (14) shows that importing it into Q_k loses
2(h-s)=384·2^k+94 paths. This is not a scalar +2 repair.

[rooted179_probe.py](rooted179_probe.py) reads the pinned inventory and its
literal code convention, checks exactly these 39 classes and 329 roots,
and records all endpoint vectors, both candidate graphs' literal rows,
the weighted-start counts, and the source-bank hash in
[rooted179_receipt.json](rooted179_receipt.json). Reproduce with

    python3 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/returns/rooted179_probe.py

This does not exclude multiplier-243 target counts, k-dependent gadgets,
or other source-changing grammars. The next constructive options are to
find a useful return from the positive 375/377 families, or use an operation
that preserves the target while changing relation (9) or the one-exit
interface. The native +2 theorem in PROOF.md is unaffected.
