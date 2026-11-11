# Session254: zero-count-cost padding resets and the mod27 curvature obstruction

**PROVED HERE, with FINITE-EXACT controls.** These are actual operations on
Hamilton-count fibres and complete modular substitution interfaces. They do
not give integer target coverage, a strong-source replacement, or Chan–Pak.
No literature-priority claim is made.

## 1. Inheritance and concept board

Session253's tournaments/PROOF.md supplies the factorial substitution jet,
its mod3 scalar/parity collapse, and the exact positive constructor
B_r(S)=C3[S ▷ TT_r,K1,K1], with H(B_r(S))=H(S)+2^r Q_S(2).
THM-3166 supplies Q_(A▷B)=Q_A Q_B; THM-3181 supplies
Q_S(-t)=(-1)^|S| Q_S(t). All dependencies are pinned in SOURCE_PINS.json.
They are used at their stated scope, not re-audited.

Closest mechanism: the modular jet and its legal parity-flipping padding.
Hostile: H alone cannot control cyclic substitution. Corrected near miss:
complete residue support is not bounded integer coverage. Least-used sidecar:
Taylor derivatives at the three residue centres 0,1,-1, before quotienting.
Board: positive count fibres; factorial jets; p-adic local Taylor data;
source-preserving versus source-changing operations; downstream integer cost.

The concurrent five-core reset (alternatives/RESET.md) replaces a port and
retains seven vertices at H141; its root-replacement variant is uniform but
still fails its425 source frontier. Our reset is orthogonal: it retains ALL
old vertices and arcs, preserves exact H, but only normalizes a modular
interface and becomes nonstrong. It cannot repair that integer425 failure.

## 2. Definitions and the evaluation-to-substitution bridge

Let p_S(k) count unordered spanning directed-path covers, set
F_S(k)=k!p_S(k), and Q_S(t)=sum_k p_S(k)(t)_k in Z[t]. Put

    H=Q(1), D=Q'(1), E=Q''(1)/2,
    K(S)=2HE-D²+HD (mod3).

The integer E is the second Hasse derivative, not division in a residue field.
We work with 3-unit sources: 3 does not divide H.

The complete jets are (F(1),...,F(5)) mod9 and
(F(1),...,F(8)) mod27, since the first factorial valuations reaching2 and3
occur at6 and9 respectively. By Session253 these jets are congruences for
ALL finite iterated tournament-substitution contexts.

Equality of the integer functions Q:Z -> Z/m is equivalent to equality of
the complete modular jet for m=9,27. In one direction
F(k)=Delta^k Q(0) recovers each jet coordinate. Conversely
Q(t)=sum_k F(k) binom(t,k); terms beyond the jet vanish, including at
negative integer t since binom(t,k) is integral. Thus evaluation equality
is exactly the right invariant, not merely a necessary test.

## 3. A small Taylor interface before normalization

For any integer polynomial f and integers c,u,

    f(c+3u) = f(c)+3u f'(c)                 (mod9),
    f(c+3u) = f(c)+3u f'(c)+9u² f^[2](c)  (mod27).

These follow by binomial expansion; discarded terms contain 3² or3³.
Use c in {0,1,-1}. Tournament parity reconstructs the -1 chart from1.
Since Q(0)=0, the arbitrary-source interface is:

- mod9: parity, H mod9, D mod3, Q'(0) mod3;
- mod27: parity, H mod27, D mod9, E mod3,
  Q'(0) mod9, Q^[2](0) mod3.

Parity forces the derivative at0 to vanish for even order, and the second
Hasse derivative at0 to vanish for odd order. These are sufficient tuples;
we do not claim all formal tuples are realizable by tournaments.

## 4. Physical zero-count-cost reset modulo9

Adjoin r universal sink vertices in transitive order. This is an actual
induced extension S_r=S ▷ TT_r, with

    Q_r(t)=t^r Q(t), H_r=H, D_r=D+rH.

All original vertices and arcs survive. For r>=2, Q_r(3u)=0 mod9.
If 3 does not divide H, and a desired order parity epsilon is specified,
there is a unique r in {2,...,7} satisfying

    D+rH=0 (mod3), |S|+r=epsilon (mod2).

This is the Chinese remainder theorem modulo6. For that padding,

    Q_r(x)=0                          if 3 divides x,
    Q_r(x)=H                          if x=1 (mod3),
    Q_r(x)=(-1)^epsilon H              if x=-1 (mod3),  modulo9.

**Theorem A.** Every 3-unit source has a constructive, exact-H-preserving
reset of its complete mod9 substitution jet with at most7 added vertices.
The normalized jet depends only on H mod9 and chosen parity. Thus any two
sources with equal H modulo9 can be made interchangeable modulo9 in EVERY
iterated substitution context, using separate prescribed paddings.

This is stronger than scalar congruence for a single quotient. It is weaker
than equality of integer responses, strongness, or integer target coverage.

## 5. Mod27 retains an invariant curvature

Under the same padding,

    E_r=E+rD+binom(r,2)H,
    2H E_r-D_r²+H D_r = 2HE-D²+HD.

The latter is an exact integer identity. Thus K is invariant under every
ordinal padding, not just the normalization range.

For r>=3, Q_r(3u)=0 mod27. For a 3-unit H and chosen epsilon there is a
unique r in {3,...,20} satisfying

    D+rH=0 (mod9), |S|+r=epsilon (mod2).

Write e=K/(2H) in F3, which is defined since H is a3-unit. Then E_r=e mod3.
For any x not divisible by3 choose s in {1,-1} with sx=1 mod3 and put
u=(sx-1)/3. The normalized evaluation formula is

    Q_r(x)=s^epsilon [H+9u² e] (mod27),
    Q_r(3u)=0 (mod27).                                  (1)

**Theorem B.** Every 3-unit source admits a zero-H-cost normalization of its
complete mod27 substitution jet using at most20 added vertices. For fixed
parity, two normalized jets agree if and only if their H mod27 and K mod3
agree. All3 values of K occur in actual tournaments.

Proof of sufficiency is (1) and the jet bridge. For necessity, Q(1) recovers
H; since H is a unit, Q(-1) recovers parity if not fixed. At x=4, (1) gives
Q(4)=H+9e mod27, recovering e and hence K. This also proves that K cannot be
silently discarded. QED.

More generally, among sources already padded by at least3 vertices, the
jet determines D mod9 and E mod3: the difference of evaluations at4,-2
is6D mod27, while their sum minus2H is18E mod27. Dividing the known factors
3 or9 in the integer congruence recovers the stated coordinates.
Consequently two sources with the same H mod27 and different K cannot have
equal mod27 jets after ANY paddings. Otherwise append3 more sinks to both;
jet equality is preserved by multiplication by t³, and the preceding
reconstruction contradicts their padding-invariant K values.

The boundary 3|H is essential to these normalizations: D+rH cannot vary
freely mod3, and need not be cancellable. No unit-free theorem is asserted.

## 6. An actual same-count obstruction and a quotient that sees it

Two explicit six-vertex tournaments, encoded by outgoing-neighbour bitmasks,
are

    A=(56,9,51,36,10,18), Q_A=t^6+14t^4+22t²,
    B=(10,4,57,50,3,19), Q_B=t^6+12t^4+24t².

Both have H37. Their derivatives are D_A106,D_B102, and their curvature
classes are K_A2,K_B0. Thus they occupy the SAME exact Hamilton-count fibre
but have different invariant modulo27 padding orbits. This is a concrete
failure of scalar-only source-fibre normalization, one prime-power above
the positive mod9 theorem.

For even normal parity, the prescribed minimal r>=3 are20 for A and6 for B.
Their mod27 normalized values at2 are19 and10 respectively. Therefore the
physical C3 two-ear, whose count is H+Q(2), distinguishes them by residues
2 and20 modulo27. This is a specific legal observer, not an existential
context-separation assertion. The graphs and polynomials are reproduced by
endpoint and anchored-set-partition recurrences in probe.py.

### 6.1 A positive single-arc surgery changes the curvature

The obstruction has a particularly economical repair. Reverse ONLY the arc
3->5 of A. The new tournament is

    A*=(56,9,51,4,10,26), Q_(A*)=t^6+12t^4+24t².

Thus H(A*)=H(A)=37 but K(A*)=0 differs from K(A)=2. This is a literal
count-neutral arc reversal: every vertex survives, and all14 other arcs
stay unchanged. It has zero vertex cost and zero Hamilton-count cost.
The profile change is the exact polynomial

    Q_(A*)-Q_A = 2t²(1-t²).

This vanishes at1 (preserving H), but its derivative and second derivative
supply precisely the change in the surviving curvature. Both graphs are
strong: both have directed cycle0->3->2->5->4->1->0, directly checked
in the displayed rows.

For ANY tournament R with 3 not dividing H(R), perform the same arc reversal
inside A▷R. The resulting A*▷R has exactly the same H=37H(R), the same
vertices, and changes curvature by1 modulo3, by the additive law below.
The same holds with R on the other side or further ordinal factors. This
is an infinite family of legal curvature-changing target-fibre operations,
not just equality of formal coefficients. We do not claim every source of
count37H(R) contains this patch or that the marked endpoint data survive.

## 7. Exact downstream cost and compositional curvature

Padding costs exactly r vertices and ZERO change in the target H. However,
a later two-ear reads

    H(C3[S_r,K1,K1])=H+2^r Q_S(2).

Thus the reset can be expensive after substitution. For a fixed normalized
parity/state and r>=a (a2 or3), all lawful r lie in a single progression of
step6 or18. Since Q_S(2)>0, the chosen representative in2..7 or3..20 uniquely
MINIMIZES that two-ear's actual output among this stable-padding family.
This is an exact costed constructor, not a global optimum among tournaments.
The bounds are H+2^7 Q_S(2) and H+2^20 Q_S(2), with final orders at most
|S|+9 and |S|+22. We retain Q_S(2) as the downstream-cost sidecar; the finite
modular state alone cannot recover it.

There is also an operation law. For 3-unit sources define

    C(S)=2E/H-(D/H)²+D/H in F3 = K(S)/H².

Since H²=1 in F3, C=K numerically. Product differentiation for
Q_(S▷T)=Q_S Q_T gives C(S▷T)=C(S)+C(T). Universal sink padding has C0.
So curvature is an additive ordinal-join charge, not an arbitrary label.
Attaching a source with nonzero charge changes H by multiplication; this
does NOT repair the charge within the same exact-H fibre for free.

## 8. Transfer ledger, limits, and next precise obligation

Source: actual tournaments with their path-colour polynomial. Target: finite
substitution jets. Map: evaluation/Taylor reduction at0,±1. Preserved:
Hamilton count under the normalization operation; every modular substitution
response after state identification. Lost: exact downstream output, order,
strongness, endpoint readiness. Sidecars: actual S, chosen r, Q_S(2), K.
Cheapest decisive test: the explicit H37 pair and its C3 two-ear observer.

The two prime-power levels explain why a first-order reset can succeed yet
fail at the next level: first logarithmic derivative can be gauged away,
second logarithmic derivative remains. This is an exact operation/interface
statement. It supplies no map from planar-JC derivatives or LRC phase to
physical tournament replacements, so such a cross-problem identification
would currently be an analogy only.

Next positive fibre question: extend the explicit H37 single-arc surgery to
a family usable by the five-core ternary constructor, retaining the actual
endpoints/readiness needed by its branch. Section6.1 supplies a genuine local
operation and its exact changed response; ordinal padding alone is now
conclusively insufficient. For the existing
H141 branch, H is not a3-unit, so this note's normalization cannot be applied
without a distinct nonunit analysis. This restriction is not hidden.

## 9. Reproduction and scope

Run `python3 probe.py > output.json` from this directory and compare with
`python3 -O probe.py`. The checker enumerates every labelled tournament
through order5 (1099 cases), plus500 seeded order6 samples. It computes
Hamilton counts on every subset by endpoint DP, then unordered covers by
anchoring the least vertex's block, and converts to integer Q coefficients.
It checks4052 normalizations against EVERY residue of9 or27 and compares
complete finite-difference jets across40 normalized states. All3 curvatures
occur; the first same-H unequal-curvature pair is the displayed H37 pair.
General proofs are algebraic; the sample is not an exhaustive order6 census.
No inherited large coverage computation or claimed conjecture was audited.
