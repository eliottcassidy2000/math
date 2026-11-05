# Compatible rooted-clock profiles inside affine symplectic rank four

**PROVED constructions and continuation-rank obstruction; FINITE-EXACT controls.
Independently reviewed.** Session 193, 2026-11-05.

This combines the actual Session 191 two-level physical fan, its multiple
incoming ports, and Session 192's affine symplectic substitution closure. It
supplies a reusable complete state for arbitrarily nested rooted clocks and
realizes every state of the two-level fan in rank four. It also identifies
precisely why this combination does not yet give the desired literal ternary
digit induction. The full odd Hamilton-path spectrum outside `{7,21}` remains
**OPEN**.

## 1. Inheritance and the active board

The closest proved mechanisms are:

- [Session 191 two-level cycle tree](../../session191_tournament_growth_20261104/ternary_fan/TWO_LEVEL_CYCLE_BUDGET.md): the physical two-level
  `3N+2 / 3N+4` fan, including its two-port root-change bridge.
- [Session 191 rerooted continuation](../../session191_tournament_growth_20261104/ternary_fan/REROOTED_CONTINUATION.md): the cross-cover
  correction needed when an exterior vertex has two incoming core ports.
- [Session 192 symplectic substitution](../../session192_surprising_approaches_20261105/symplectic_substitution/PROOF.md): the integral Heisenberg compiler
  for substitution, with no increase beyond affine symplectic rank four.
- THM-3121, `path-cover-walk-content-substitution-kernel`, and THM-3134,
  `tournament-endpoint-jet-and-c3-newton-profile-transform`: the full path-cover
  profile is the compositional sidecar, not just its scalar H coordinate.

The relevant hostile examples are Session 191's physical H13 ancestry and
unit-start obstructions, and Session 192's same-H singleton/transitive-edge
pair. None is evaded by giving vertices symplectic coordinates. The corrected
near miss is the assertion that geometric substitution closure itself closes
the arithmetic digit grammar. The least-used useful sidecar here is the
**whole marked path-cover polynomial**, including mixed end/start ports.

The live concept board is: literal roots; arbitrary clock modules; full
labelled-bin polynomials; rank-four Heisenberg geometry; cycle budget; context
rank. The anchor is compatible physical states. The orthogonal question is
whether those states admit fixed-dimensional linear compression. The wildcard
is a shift operator in the number of labelled path bins; it becomes exact
rather than an analogy.

No new literature-priority claim is made. The profile/substitution framework
is inherited. The compact marked-clock operator, its transported physical
fan, the exact finite rank-four seed geometries, and the continuation-rank
calculation are the contributions of this packet.

## 2. The complete state and marked sidecars

For a finite tournament T let

    P_T(d) = sum_j (d)_j pc_T(j),

where `pc_T(j)` counts spanning covers by j unordered nonempty directed paths
and `(d)_j=d(d-1)...(d-j+1)`. Equivalently, for every integer d>=0, P_T(d)
counts assignments of all vertices to d labelled bins, each bin carrying a
directed path, with empty bins allowed. In particular

    P_T(1)=H(T),
    P_T(2)=2H(T)+2pc_T(2),
    P_empty(d)=1,       P_T(0)=0 for nonempty T.

The complete polynomial P_T is equivalent to the complete path-cover profile,
by its expansion in falling factorials. It is not a new scalar compression.
For an ordered join, `P_(A->B)(d)=P_A(d)P_B(d)`: in each bin the A-path must
precede the B-path. For a transitive n-tournament, `P_TT_n(d)=d^n`.

Fix a nonempty core A and a nonempty set B of its named vertices. Retain the
following polynomials, evaluated at positive bin counts:

- L_u(d): A-covers in d labelled bins whose first bin is nonempty and
  **starts** at u.
- R_q(d): A-covers in d labelled bins whose first bin is nonempty and
  **ends** at q.
- K_qu(d), q!=u: A-covers in d labelled bins whose first bin is nonempty
  and ends at q and whose second bin is nonempty and starts at u.

The remaining bins may be empty. Only d>=2 is used for K. These are genuine
labelled covers; no unlabeled multiplicity or division is implicit.
Writing e_U(q),s_U(u) for endpoint path counts on an induced subset U gives
explicit finite formulas

    R_q(d) = sum_(q in U) e_U(q) P_(A-U)(d-1),
    L_u(d) = sum_(u in U) s_U(u) P_(A-U)(d-1),
    K_qu(d) = sum_(U,V disjoint; q in U,u in V)
                e_U(q) s_V(u) P_(A-(U union V))(d-2).

Thus every coefficient below can be computed exactly from the core without
rebuilding a substituted clock. At d=2, `R_q(2)=F_q(A)` is precisely the
Session 191 ordered two-path cover whose first path ends q and whose second
path may be empty. Also `K_qu(2)` is its mixed-port cross-cover constant.

## 3. The arbitrary-clock, multi-port transfer theorem

Let S be any nonempty tournament. Form W_(A,B)[S] on A, S and a new vertex w
with

    w -> S -> A,
    q -> w for q in B,       w -> q for q in A-B.

All other arcs are those of A and S. No transitivity or strongness hypothesis
on S is needed.

**Theorem 1 (complete reusable response).** For every integer d>=1,

    P_W(d) = d [
        (P_A(d) - sum_(q,u in B; q!=u) K_qu(d+1)) P_S(d)
        - (sum_(u in B) L_u(d)) P_S(d-1)
        + (sum_(q in B) R_q(d+1)) P_S(d+1)
    ].                                                     (1)

These are polynomial identities in d, since both sides are polynomials and
agree at every positive integer. In particular, putting `h=H(S)`,
`p=P_S(2)` and `C_B=sum_(q,u in B;q!=u)K_qu(2)` gives

    H(W) = (H(A)-C_B)h + (sum_(q in B)F_q(A))p,              (2)
    hc(W) = h sum_(q in B)e_q(A).                           (3)

Here hc counts Hamiltonian cycles modulo rotation. If A is strong then W is
strong, even when S is not.

**Proof of (1).** Choose the one bin containing w, in d ways. Every other
bin is an S-path followed by an A-path, either possibly empty.

If the w-bin contains no A vertex before w, then w must be its first vertex,
since no S vertex beats w. There are initially `P_S(d)P_A(d)` choices. The
invalid ones have an empty S-path in that bin and an A-path starting at some
u in B. They contribute exactly `P_S(d-1) sum L_u(d)`.

Otherwise, the w-bin is uniquely

    S_before, A_before, w, S_after, A_after,

with A_before nonempty and ending at some q in B. Its two S bins and the
other d-1 S bins give `P_S(d+1)` choices; the A choices give
`sum R_q(d+1)`. The only invalid boundary has S_after empty and A_after
nonempty starting at u in B. As q has already been used in A_before, q!=u.
Those choices contribute `P_S(d) sum K_qu(d+1)`.

The two cases are disjoint, exhaustive, and recover every selected path
cover uniquely. This proves (1). At d=1 the L term vanishes because S is
nonempty, yielding (2). A Hamiltonian cycle has a unique writing starting
at w: an S Hamilton path, then an A Hamilton path ending in B, then w.
This proves (3). Finally, A reaches a port in B and then w; w reaches every
S vertex, and S reaches A. Strongness follows. QED.

The negative terms are exact forbidden-boundary subtractions. They do not
claim that signed formal profiles are physically realizable. For every
physical input S, the operator outputs the actual nonnegative cover counts.
Equation (1) can be iterated with different cores and port sets: **the output
P_W is precisely the input state required by the next patch**. That is a
closure statement for the complete polynomial, not merely geometric closure.

## 4. The full two-level physical fan transports to arbitrary clocks

The finite data file `seed_data.py` records all fourteen Session 191 states
with their literal core adjacency rows and actual incoming port sets. In the
following table `(a,b,c)` means

    H(W)=a h+b p,       hc(W)=c h.

| State | a | b | c |
|---|---:|---:|---:|
| P | 9 | 11 | 1 |
| L | 29 | 33 | 3 |
| R | 31 | 33 | 3 |
| LL first birth | 51 | 57 | 3 |
| LL | 89 | 99 | 3 |
| LR first birth | 65 | 73 | 4 |
| LR old-arc flip | 27 | 27 | 1 |
| LR | 91 | 99 | 2 |
| RL two-port bridge | 9 | 67 | 6 |
| RL re-rooted | 31 | 34 | 3 |
| RL | 95 | 99 | 7 |
| RR old-arc flip | 25 | 27 | 2 |
| RR first birth | 49 | 51 | 2 |
| RR | 97 | 99 | 2 |

For each of `(P,L,R)`, `(L,LL,LR)`, `(R,RL,RR)`, the exact laws are

    H(left)=3H(parent)+2h,
    H(right)=3H(parent)+4h.                                (4)

Every original vertex, exterior incidence, newborn incidence, and named
arc reversal is preserved by replacing the old transitive clock with the
same S. The fourteen physical moves are listed in `EDGES` and independently
checked on labelled vertices. All states are strong. Each intermediate is
strictly below its own route's final target: the original coefficientwise
positive differences remain positive for `h,p>0`. Every state on these
specific routes has at most `7h` cycles, and the RL leaf has exactly `7h`.
No sharp lower bound for arbitrary realizations of these scaled counts is
asserted.

The two-port bridge genuinely uses the multi-port correction. Its core has
`H(A)=31`, `F_1=33`, `F_4=34` and `K_14(2)+K_41(2)=22`, so its count is
`9h+67p`, not `31h+67p`.

For S=TT_k, `(h,p)=(1,2^k)` and (4) recovers the original literal ternary
fan for every k>=1. For S=C3, `(h,p)=(3,12)`, so the first fan is

    159 -> 483,489,

with increments 6 and 12 after tripling, rather than 2 and 4.

**Exact rearming obstruction for this route.** The simultaneously literal
increments `+2,+4` in (4) require h=1. A tournament has H=1 if and only if
it is transitive. The forward implication follows because a nontransitive
tournament contains a triangle, which has three Hamilton paths, and adding
a vertex cannot decrease H: deterministically insert it at a legal position
in each old Hamilton word; deletion recovers that word, proving injectivity.
The reverse implication is immediate from the unique transitive order.
Consequently arbitrary nontransitive clocks cannot recursively rearm this
particular fan while retaining the literal digits. Replacing TT_k by more
expressive symplectic clocks gives homogeneous ternary laws, not the desired
universal arithmetic induction.

### 4.1 A permanent divisibility obstruction for retained clock modules

Put

    g(S)=gcd_(d>=1) P_S(d)=gcd_(j=1)^|S| j! pc_S(j).

The equality follows from `P_S(d)=sum_j binom(d,j) F_S(j)`, with
`F_S(j)=j!pc_S(j)`: binomial finite differences recover every F_S(j) from
integer values of P_S. It suffices to take values d=1..|S|.

**Proposition 2 (content cannot be removed).** For every rooted multi-port
clock operation in Theorem 1,

    g(S) divides g(W_(A,B)[S]).

Consequently g(S) divides every Hamilton count after any finite sequence of
these operations, with arbitrary later choices of cores and incoming ports.

**Proof.** Every coefficient in the right side of (1) is an integer at each
positive integer d. Every value P_S(d-1),P_S(d),P_S(d+1) is divisible by g(S),
including P_S(0)=0. Thus every positive-integer output value is divisible by
g(S), hence so is its gcd. Iterate this argument. QED.

In fact the same divisor is preserved by **every tournament substitution
context** with S as one entire factor, not just the rooted operators. The
inherited THM-3121 law for ordered profiles reads

    F_(Q[S_i])(d) = sum_c [x^c] W_Q(x)^d product_i F_(S_i)(c_i).

All displayed coefficients are integers, and each summand contains one
ordered-profile value of S. Hence g(S) divides every F of the output. This
extends the obstruction to changing the outer quotient or the other factors
while retaining S as an intact module; rank-four geometric closure does not
remove it.

In particular `g(C3)=3` and `g(C3[C3,C3,C3])=81`, exactly. A C3 clock therefore
cannot be rearmed to literal ternary children by *any* choice of nested
rooted-clock cores/ports: all reachable counts remain divisible by three,
whereas `3N+2` and `3N+4` do not. This excludes an entire module-preserving
strategy, rather than only the fourteen fixed fan seeds, and also applies
to arbitrary substitution contexts that retain the same clock module. It does not exclude
changing the retained module's internal arcs or using a different physical
operation that splits it. Clocks with g=1, such as the H13 witness, are not
excluded by this divisor argument.

There is also a target-level necessary condition: if H(S) is odd and an
intact-module substitution is to produce `3H(S)+2` or `3H(S)+4`, then g(S)=1.
Indeed g(S) divides H(S) and the proposed output, hence divides respectively
2 or 4, while it is odd.

No realization of H=3 has this property. The ordered strong-component
factorization of H forces an H3 tournament to be `TT_a -> C3 -> TT_b`:
Camion's cited theorem gives H(C)>=|C| for every nontrivial strong component
by cutting a Hamiltonian cycle, so the only possible component with H3 is
C3. Its full polynomial is

    d^(a+b) (d^3+2d).

Every integer value is divisible by 3, and its value at one is exactly 3,
so its profile content is exactly 3. Thus no intact-module substitution
context can send **any** H3 presentation to either H11 or H13. This is an
ancestry obstruction, not a spectrum gap. It does not apply when the old
tournament is used as a quotient and its individual vertices are substituted,
or when nonuniform exterior arcs split the old module. A global construction
must allow such a change of interface somewhere. Camion and the component
factorization are the same classical inputs used in Session 191's
`UNIT_START_BOUNDARY.md`; no new classical theorem is claimed here.

## 5. Every transported fan state remains in affine rank four

Use the alternating form

    omega(u,v)=u1 v2-u2 v1+u3 v4-u4 v3,

on points `(1,p,q,r)`. The fourteen exact integral coordinate lists in
`seed_geometry.json` realize the fourteen skeletons W_(A,B)[K1]. Every one
of their 476 pairings is checked as a nonzero integer with the prescribed
sign; the smallest absolute pairing is 139. These finite certificates do
not require a numerical optimizer to verify.

For illustration, the two-port bridge, whose old core has rows
`(14,52,40,34,45,1)` and ports `{1,4}`, has points, in the order core0..5,S,w,

    (1,631,63,15),  (1,651,78,80),  (1,198,2,10),
    (1,614,20,19),  (1,-164,-14,5), (1,37,-65,-21),
    (1,-514,-6,-7), (1,-482,-48,-19).

The bridge cannot be obtained merely by reusing the two-apex compiler: its
feedback vertex number is three, as is that of the re-rooted state and RL
leaf. Rank-four geometry is genuinely broader than two-apex geometry here.
The upper bound on rank, not a claim that these graphs all have minimum rank
four, is the property needed for substitution.

Now replace the skeleton's single clock vertex by any nonempty integral
common-affine-chart rank-four configuration S. The Session 192 compiler puts
factor points `(1,P,Q,R)` near outer points `(1,p_a,q_a,r_a)` as

    (1, M^2 p_a+M(r_a Q-q_a R)+P, Mq_a+Q, Mr_a+R).

For two different blocks the pairing is `M^2 A+M B+C`, while internal
pairings are exactly unchanged. With `D=max|B|`, `E=max|C|`, choose
`M=D+E+1`; since the integer outer pairing has `|A|>=1`, every cross sign is
preserved. Hence the output is precisely W_(A,B)[S].

**Corollary 3.** All fourteen transported states and their full physical
fan occur in the same affine symplectic rank-at-most-four class whenever S
does. Arbitrary nesting by these fourteen skeleton operators stays in that
class and has its complete count profile given by (1).

The nested operators are genuinely order-sensitive. Let D be the singleton
core operation of Section 6, and let E be the P-state wrapper from Section 4.
Starting with K1, the two nine-vertex constructions satisfy

    H(E(D(K1)))=159,       H(D(E(K1)))=675.

The second identity uses `P_(E(K1))(1)=31`, `P_(E(K1))(2)=644` and
`H(D(S))=P_S(1)+P_S(2)`. Both are checked by physical vertex DP and integral
geometric substitution. This is graph-level nested noncommutation, not a
commutator of symplectic coordinate changes (which preserve every arc).
Neither ordering is asserted to supply the missing literal digit closure.

The labelled graph transitions are retained, not the numerical coordinates
of old vertices: different states may use different coordinate realizations.
No continuous fixed-coordinate flip path, special rank-two chirotope, or
physical Lonely Runner interpretation is asserted.

## 6. A sharp obstruction to fixed-dimensional linear continuation states

Take the singleton core A={r}, B={r}. Then W is simply
`C3[S,K1,K1]`, and `P_A(d)=d`, `L_r(d)=R_r(d)=1`. Formula (1) becomes

    (D f)(d)=d^2 f(d)+d f(d+1)-d f(d-1).                   (5)

This is a reusable scalar-polynomial operator, not a map depending only on
H=f(1). Let U_0(n)=TT_n, n>=1, and recursively

    U_(t+1)(n)=C3[U_t(n),K1,K1].

Then `P_(U_t(n))=D^t(d^n)` and

    H(U_t(n)) = sum_(j=1)^(t+1) c_(t,j) j^n,              (6)
    c_(t,t+1)=t!.

Indeed, applying (5) shifts an evaluation index by at most one. Starting at
index one, reaching t+1 in exactly t steps requires every shift to increase
by one, whose coefficient is `1*2*...*t`. The first coefficient rows are
`(1)`, `(1,1)`, `(-1,5,2)`; in particular the n=1 path counts start 1,3,15.
All signs in this linear description refer to exact algebraic cancellation,
not negative physical path counts.

**Theorem 4 (infinite rational linear continuation rank).** For every m>=0,
let M_m have rows t=0..m, columns n=1..m+1, and entries H(U_t(n)). Then

    det M_m = (m+1)! [product_(k=1)^m k!]^2 > 0.           (7)

Therefore no finite-dimensional rational vector state, with a fixed linear
update for the C3 singleton-core patch and a linear H output, represents all
of these physical inputs and all their nested continuations.

**Proof.** The coefficient matrix of (6) is lower triangular, with diagonal
`0!,1!,...,m!`. The evaluation matrix has entries j^n for j,n=1..m+1; its
determinant is the Vandermonde product times `product j`, namely
`(m+1)! product_(k=1)^m k!`. Multiplying the two determinants gives (7).
A d-dimensional state representation factors every M_m through Q^d, forcing
rank at most d, contrary to (7) for m>=d. QED.

This lower bound allows the initial state encoding of S to be arbitrary;
only its **transition and observation are required to be linear**. It does
not rule out nonlinear compression, state-dependent arithmetic transitions,
a finite set of control types with unbounded counters, or a different proof.
It explains why the two-apex language's inherited 17-dimensional rational
continuation state cannot simply be reused after this nested closure.

These witnesses really remain low-dimensional geometrically and escape a
bounded-feedback language. Their orders and largest transitive induced
subsets satisfy

    |U_t(n)|=n+2t,       alpha(U_t(n))=n+t,
    tau(U_t(n))=t.                                        (8)

A transitive subset of a cyclic substitution cannot meet all three outer
blocks. Its maximum size is therefore `max(alpha(previous)+1,2)`, giving (8)
by induction. All t>=1 outputs are strong. For t>=2 an outer vertex dominates
an earlier directed triangle. A rank-two alternating determinant tournament
has transitive in- and out-neighborhoods, because each lies in an open
half-plane and is ordered by angle. Thus these U_t(n) require rank at least
four; the substitution compiler gives rank at most four. Their minimum
alternating sign-realization rank is exactly four for every n>=1,t>=2.

## 7. Verification, finite scope, and remaining target

From the workspace root run

    python3 verify.py
    python3 -O verify.py

The verifier uses only the standard library and explicit exceptions. The
coordinate source is literal integer data, not a floating point test.

- Exhaustive full-polynomial test: all 11 labelled tournaments of orders
  1..3 as cores, all their nonempty port sets, and all 11 labelled clocks of
  orders 1..3, for 693 full identities. Actual full profiles come from an
  independent permutation-and-cuts enumerator; H is also checked by subset DP.
- All fourteen seed configurations have every sign checked exactly.
- Seventy labelled geometric substitutions use clocks K1, TT2, C3, the
  two-apex H13 witness, and the nine-vertex `C3[C3,C3,C3]`.
- Forty-two actual fan outputs (the first three clocks) are independently
  counted by vertex-subset Hamilton and cycle DPs. Fourteen C3-clock outputs
  have their second full-profile coordinate checked separately.
- Seventy named physical move checks preserve exactly the stated old arcs.
- Fifteen nested examples, n=1..3 and t=0..4, are checked by vertex DP and
  exhaustive transitive-subset search. Nine context minors m=0..8 are checked
  by exact rational elimination, independently of the determinant formula.

These are control universes, not censuses of all affine rank-four tournaments.
The 693 tiny controls also verify the full ordered-profile gcd divisibility,
and the named source profile contents are independently evaluated exactly.
A further 319 arbitrary quotient contexts (all labelled quotients and clock
factors of orders 1..3, at every factor address, with other factors singletons)
verify full profile-content divisibility by independent permutation cuts.
The all-order conclusions follow from the displayed bijections, compiler,
and determinant proofs, not extrapolation of finite data.

Connection contract: the source is a labelled rooted-clock tournament with
its full marked core profile. The target is a complete output path-cover
polynomial and a rank-four labelled realization. The map is the multi-port
run decomposition (1) plus Heisenberg substitution. It preserves every arc,
every path-cover count, every named physical move, and the exact cycle
response. The scalar quotient H loses higher-bin evaluations; the new context
minors prove that unboundedly many linear response coordinates can be needed.
The geometric representation alone does not restore those lost coordinates.
The complete-profile gcd is an additional preserved arithmetic obstruction
when a clock module is kept intact.

The remaining global obligation is to choose physical core/port operators
and target-dependent state updates that achieve literal arbitrary admissible
odd H, including rearming both digit branches. Equation (1) supplies exact
compatible states, and the fourteen seeds supply actual low-rank physical
operators. The content proposition and equations (4) and (7) locate separate remaining obstacles:
a nontrivial clock-content divisor cannot be removed,
clock multiplication scales the digits, and nested continuation cannot be
compressed to a fixed-dimensional rational linear interface. None is a
new forbidden Hamilton-path value or a proof of the full `{7,21}` spectrum.
