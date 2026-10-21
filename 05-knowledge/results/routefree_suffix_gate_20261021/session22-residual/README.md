# Route-free terminal-zonotope suffix gate

2026-10-21. **PROVED:** the conditional all-future energy gate below, with no
assumption that an unfinished route ever reaches the frontier. **FINITE-EXACT:**
listed bounded-prefix certificates. **OPEN:** uniform 153N coverage for every
N>=4096, and Collatz entry/convergence. No theorem ID is reserved.

## Inheritance, change of coordinate, and honest target

The closest proved mechanism is the actual-clock source ledger in
`cross_scale_frontier_credits_20261021/check.py`, and the closest geometric
sidecar is the twelve-cut subset envelope in
`annulus_block_completion_20261020j/README.md`. The incoming
`finite_horizon_budget_20261021a/README.md` proves the short-horizon wedge
N>=2^(8d+12). The reconstructed Session21 remainder packet is historical
context only, not computational evidence for this packet.

The inherited moment proposal Y<=105N fails at N4096,d149: Y=433053>430080.
This is freshly reproduced in `moment-hostile.json` (not merely inherited).
Dropping actual
arrival times from a cancellation packet is also unsafe. Instead of trying to
bound an arbitrarily deep arrived fibre directly, this packet changes to the
**unfinished-source** coordinate at one finite cut. Its vectors have bounded
length even when their eventual destinations and clocks are wholly unknown.

Working board: finite-cut completed energy; unfinished charge zonotopes;
terminal linear supports; exact primitive clocks; source-owned residual
repartition. The first four produce a gate; the last prevents treating it as
new negative credit.

Use positive odd primitive sources s not congruent to5 mod8, the fixed
627-point stopped frontier V, and the native primitive map P. A step first
divides 3s+1 by 2^a, a=2 for s=1 mod8 and a=1 otherwise, then repeatedly
removes a B(x)=4x+1 factor while the odd intermediate is5 mod8. Its clock is
a+2j, not the fully odd Collatz clock. At scale N, old sources satisfy s<4N
and new sources 4N<=s<256N. A_N is the inherited two-coordinate charge;
A_(64N)=A_N on the old sources. Every charge is zero or in the angular
alphabet

    (1,0),(2,1),(1,1),(1,2),(0,1),(-1,1),
    (-1,0),(-2,-1),(-1,-1),(-1,-2),(0,-1),(1,-1).

Set Q(a,b)=a²-ab+b² and beta(v,w)=Q(v+w)-Q(v)-Q(w).
At a fixed cut h, W_u and H_u sum the arrived old and new signed charges at
terminal u, with the sign (-1)^t from the actual accumulated clock. Put
T_u=W_u+H_u and Delta0=sum_u[Q(T_u)-Q(W_u)]. Unfinished charges below are
**unsigned original source charges**, not guesses about future clock parity.
The zero charges may be omitted.

The physical target is Delta<=153N. If I=Delta-A, M=153N-A, and an owned
repartition has R'=I-J+K and J'=J-K, the equivalent residual target is

    R' <= M-J'.

The stronger R'<=M does not follow. K is still added to R'; no extracted
pair is spent again. Our proof bounds Delta directly and allocates no pair
credits at all.

## 1. Exact signed-zonotope scalar

For a finite indexed multiset F of charge vectors, define

    D(F)=max_{epsilon_s in {-1,1}} Q(sum_s epsilon_s v_s).

The indexing retains multiplicity. If c_i counts angular type i, let
Vsum=sum_i c_i v_i and C_i=sum_{j=i}^{i+5} c_(j mod12) v_(j mod12).
Then exactly

    D(F)=max_(0<=i<12) Q(2C_i-Vsum).                 (1)

Proof: at a maximizing signed sum z, flipping epsilon_s gives
Q(z-2 epsilon_s v_s)<=Q(z), hence

    epsilon_s beta(z,v_s)>=2Q(v_s)>0.

Thus all copies of each type have the same sign, determined by a strict
open half-plane through the origin. No populated direction is on its
boundary. A small perturbation through empty directions selects one of the
twelve consecutive six-direction half-planes. Its signed sum is2C_i-Vsum.
Conversely every displayed sum is an allowed sign choice. For F empty both
sides are zero. This is a signed-zonotope maximum, not the earlier arbitrary
subset maximum, although both use the same angular arrangement.

### The missing-destination lemma

Assign any subset of F to arbitrary terminals, with arbitrary signs, and let
z_u be the assigned sum at terminal u. Then

    sum_u Q(z_u) <= D(F).                          (2)

To prove it, combine the nonzero terminal sums one at a time. Flip the sign
of the entire next group if necessary to make its beta product with the
accumulated sum nonnegative. Each step increases Q by at least Q(z_u).
The combined sum therefore has Q at least sum_u Q(z_u), and is an allowed
signed sum of the assigned sources. Add the unassigned sources one at a time,
choosing a sign with nonnegative beta product. This can only increase Q and
produces an allowed signed sum using every source. Definition of D proves
(2). The group flips are a proof device, not modifications of physical clocks.
The missing destinations, signs, and absent sources have only been enlarged
in a safe upper bound. The argument works for every positive-semidefinite
quadratic form; the twelve-pattern evaluator specifically uses this alphabet.

## 2. All-future gate, including unfinished old sources

Let F_o,F_n be the old and new sources unfinished at cut h. Define

    D_o=D(F_o), D_n=D(F_n),
    L_n=sum_(s in F_n) max_(u in V) |beta(T_u,v_s)|,
    L_o=sum_(s in F_o) max_(u in V) |beta(H_u,v_s)|.

For every later horizon d>=h, including arbitrarily large finite horizons,

    Delta(d) <= Delta0 + D_n + L_n + L_o + 2 sqrt(D_o D_n).   (3)

No assertion that any unfinished source will ever arrive is required.
To see this, let z_o,u and z_n,u be the contributions of the unfinished old
and new sources which actually arrive between h and d. These are a subset
of the original sources, with some terminal assignment and actual clock
sign. Direct expansion, including cancellation of the old self-energy,
gives

    Delta(d)-Delta0
      =sum_u[Q(z_n,u)+beta(T_u,z_n,u)
             +beta(H_u,z_o,u)+beta(z_o,u,z_n,u)].

Equation(2) bounds the new quadratic term by D_n. Each linear term is
bounded source by source by L_n or L_o. Cauchy--Schwarz for the direct-sum
Q inner product bounds the mixed term by
2 sqrt(sum Q(z_o,u) sum Q(z_n,u))<=2 sqrt(D_o D_n).
This proves(3), including zeros and sources that never arrive. For integer
verification replace the last term by ceil_sqrt(4D_o D_n).

When all old sources have arrived, D_o=L_o=0, so the gate simplifies to

    Delta(d) <= Delta0+D_n+L_n.                    (4)

This is not valid if unfinished old sources are silently discarded. One
terminal with W=0,H=(1,0), one unfinished old charge(1,0), and no new charge
has Delta0=1 but future Delta=3. The full gate gives3 exactly.

A weaker but transparent all-old-completed corollary uses m unfinished
nonzero new sources and C=max_u Q(T_u):

    Delta(d) <= Delta0+3m²+2m sqrt(3C).            (5)

Consequently a uniform family satisfying Delta0<=(153-epsilon)N,
C<=cN and m<=eta sqrt(N), with
3eta²+2eta sqrt(3c)<=epsilon, has uniform all-future coverage.
This is a concrete square-root pending-count gate, not a proved estimate for
the native family. Formula(4) is substantially sharper in the finite probes.

## 3. Independent bounded-prefix certificates

`check.py` imports no historical producer. It authenticates frontier.json,
rechecks its grounded ranks, computes P using the actual first division
valuation, reconstructs charges from independent three-row coordinates, and
retains original source, mod15, mod35, terminal, full clock and depth. Each
new exploratory route walk has at most h steps. Previously proved finite
suffixes can be reused; unfinished source destinations and future signs are
never read by the gate. Cached suffix information beyond a source's cutoff
is discarded from that source's classification and bound. No unrestricted
route search or convergence hypothesis is used.

All finite cuts0..h are checked directly. Formula(3) covers every later cut
without a later replay. Hashes commit to every source row, not just an energy
aggregate. Mod15/mod35 are checked for their shared mod5 coordinate and joint
mod105 reconstruction, not treated as independent525 states. They constrain
the literal source receipt; the gate deliberately enlarges future terminal
and sign possibilities, never falsely claims to preserve unknown guards.

The default receipt contains:

| N | Cut h | Unfinished nonzero new | Delta0 | D_n | L_n | All-future upper | 153N |
|---:|---:|---:|---:|---:|---:|---:|---:|
|4096|105|386|421276|91876|104672|617824|626688|
|8192|112|615|793287|229957|188427|1211671|1253376|
|16385|121|961|1483533|537663|449028|2470224|2506905|
|32768|131|1320|2837313|1004644|921964|4763921|5013504|

Every old source is finished in these four certificates. These are individual
scales, not an interval of scales. N16385 lies just beyond the inherited
4096..16384 finite range. Native N7,d8 remains hostile:
Delta=1084>1071; the upper1415 also correctly refuses certification.

Independent small controls exhaust all signed choices for every multiset of
at most4 angular vectors; they also exhaust absence/terminal/sign choices
for a mixed old/new four-source test. A separate hostile records why pending
old sources cannot be ignored. Standard and python -O immutable replays use
explicit checks rather than removable asserts.

Reproduce from any checkout, writing only to an explicit output:

    python3 check.py --repo-root /path/to/math --output receipt.json
    python3 -O check.py --repo-root /path/to/math --verify receipt.json
    python3 check.py --repo-root /path/to/math --cases 4096:149 --verify moment-hostile.json

## 4. What this changes and does not change

Source/target/map: actual stopped-prefix ledger -> pending indexed charge
multiset -> signed zonotope plus terminal supports. The preserved target is
an upper bound on the physical energy increment at every future horizon.
Discarded future destinations, eventual convergence, arrival order and clock
parity become adversarial choices, so no sidecar is secretly assumed. The
retained source sidecar restores multiplicity, actual past clock, frontier
and the mod15/mod35 guards.

The structural advance is an all-future stopping certificate from a bounded
prefix, not merely more completed orbit enumeration. It can stop with many
unfinished sources, including a hypothetical nonterminating one. The
remaining uniform obligation is to prove that every intended N admits a cut
whose explicit D/L gate lies below153N, together with prefix coverage to that
cut. A uniform bound on pending counts and terminal supports would suffice;
neither is proved here. The incoming short-horizon wedge and this tail gate
must not be declared overlapping uniformly without that missing estimate.
