# Switched Collatz returns: exact recharge gates and a common height certificate

**Status, 2026-10-19:** PROVED elementary conditional descent theorems and
source-wise obstructions. The transition graph must be justified on the
actual source. No universal source coverage or Collatz convergence is claimed.

## 1. What is new, and what is inherited

The [local-return-rank packet](collatz_return_rank_obstruction_20261019.md)
proves exact valuation budgets for one word and exhibits unbounded recharge
when its rational center changes. The [fixed-translation packet](fixed_translation_descent_20261019.md)
provides smaller-source receipts outside a linear core, conditional on an
isolated merge actually occurring. Neither statement controls arbitrary
switches of return word for an unresolved integer source.

This note supplies three sharper statements.

1. Every nontrivial change of word center has an exact, bounded-rank
   resonance gate, and that gate admits arbitrarily large lawful recharge.
   Merely taking all the finitely many word-centered valuations together
   still cannot remove this obstruction.
2. A finite graph of lawful affine blocks has a common weighted real-height
   rank exactly when every directed cycle has contracting multiplier. The
   certificate permits individually expanding blocks and gives an explicit
   logarithmic deadline to a finite core or an uncertified exit.
3. An actual expanding adjacent-pair return exists. A following reset gives
   a genuine smaller-source receipt on an explicitly identified cylinder,
   but is not forced by the expanding word. Arbitrarily long repetitions
   of this expanding word obstruct any finite contracting graph that
   purports to include all lawful finite return histories.

The graph-directed criterion is a self-contained application of finite
multiplicative difference inequalities, not a claim of novelty for general
affine-system stability. The source guards and explicit Collatz return are
the mathematical interface established here. No residue census is used.

Throughout, T(x)=x/2 for even x and T(x)=(3x+1)/2 for odd x. A nonempty
branch word i has length n_i, positive odd count a_i, map

    F_i(x) = lambda_i x + b_i = (3^{a_i}x+C_i)/2^{n_i},
    p_i = C_i/(2^{n_i}-3^{a_i}).

Both 2- and 3-adic valuations of rational differences are well defined;
p_i is integral at both primes. Every asserted forward or inverse path
retains its parity cylinder, integrality, and positive-source guards.

## 2. The exact resonance gate, and why more fixed centers do not fix it

Put r_i(x)=v_2(x-p_i), excluding x=p_i when necessary. At a lawful i step,
for any next center p_j,

    F_i(x)-p_j = lambda_i(x-p_i) + (p_i-p_j).       (2.1)

If p_i=p_j, its rank simply drops by n_i. If p_i differs from p_j, put
s_ij=v_2(p_i-p_j), a finite nonnegative integer. Then

    r_j(F_i(x)) = min(r_i(x)-n_i, s_ij)
                    whenever r_i(x) differs from n_i+s_ij.    (2.2)

Only the exact shell r_i(x)=n_i+s_ij can recharge the new rank beyond s_ij.
This is the ultrametric equality for summands of unequal valuation.

**Universal recharge theorem.** For every two words with distinct centers
and every M>s_ij, there are arbitrarily large positive integers x such that
i is lawful, r_i(x)=n_i+s_ij, and r_j(F_i(x))=M. If M>=n_j, the actual next
word is j, and exactly floor(M/n_j) complete j blocks can be followed.

**Proof.** Choose a positive integer z whose difference from p_j has exact
2-adic valuation M, and additionally require

    2^{n_i}z = C_i  (mod 3^{a_i}).

The Chinese remainder theorem combines these independent conditions. Set
x=(2^{n_i}z-C_i)/3^{a_i}; sufficiently large lifts make x positive. The
integer inverse of the i branch follows that branch. Since M>s_ij,
v_2(z-p_i)=s_ij, hence r_i(x)=n_i+s_ij. The subsequent j budget is the exact
one-word theorem. All forward iterates of a positive integer are positive. ∎

The dual statement at prime 3 also holds. For a legal inverse i step y→z,
replace n_i by a_i and r by v_3(·-p). Arbitrarily large recharge is obtained
by prescribing z near p_j modulo 3^{M+1}, while placing z in i's independent
parity cylinder modulo 2^{n_i}, and setting y=F_i(z). If a prescribed finite
inverse j run needs a positivity guard, choose the CRT lift sufficiently
large; every finite inverse affine prefix then remains positive.

There is a stronger finite-center obstruction. Let S be any finite set of
2-adically integral rational centers containing the p_i. If there are at
least two distinct p_i, some inverse image q=F_i^{-1}(p_j), with p_j in S,
lies outside S. Otherwise every injective map F_i^{-1} permutes S. Each
q in S then lies on a finite F_i orbit. Since lambda_i is positive and not
one, F_i^k(q)=q implies q=p_i. This would force S={p_i} for every i, contrary
to the distinct-center hypothesis.

For such a q outside S, choose positive integers x arbitrarily close to q
2-adically. At sufficiently high precision all initial values

    (v_2(x-p))_{p in S} = (v_2(q-p))_{p in S}

are the same fixed finite vector, while v_2(F_i(x)-p_j) is arbitrarily large.
The i cylinder is automatic once v_2(x-q)>=n_i: F_i(q)=p_j is 2-adically
integral, and branch integrality determines that cylinder. Thus, in
particular, no positively weighted sum or maximum of these finitely many
valuations is a universally decreasing switching rank. A different,
nonmonotone encoding is not excluded by this argument. For the concrete
F,H,J bank below, F^{-1}(p_H)=-9/431 lies outside all three centers:
all three initial ranks can remain exactly (3,1,0), while the next vector
is (1,M,0) for arbitrarily large M. Adding J's valuation to the old two
therefore does not repair their rank.

## 3. A common weighted height that does survive switching

Let G be a finite directed graph. Every edge e:v→w carries an actual
positive Collatz branch F_e(x)=A_e x+B_e, with rational A_e>0, B_e>=0.
The edge may have an arbitrary exact source guard. A source path must obey
that guard and the graph, rather than merely choose a convenient edge.

**Cycle-contraction theorem.** The following are equivalent:

- every directed cycle has product of its edge multipliers strictly below 1;
- there exist rational q in (0,1) and rational weights omega_v>0 such that
  A_e omega_v <= q omega_w for every edge e:v→w.

**Proof.** Multiplication around a cycle proves necessity. Conversely,
there are finitely many simple cycles. Choose rational q<1 so close to 1
that every such cycle has product of A_e/q below 1. Define omega_v to be
the maximum product of A_e/q over all simple directed paths ending at v,
including the empty path of value 1. Any nonsimple path can have a cycle
removed, increasing its product. Appending an edge therefore proves the
required inequality. All weights are rational. An acyclic graph has no
cycle condition and the same construction works for any rational q in
(0,1). ∎

Given these weights, choose a rational K>=0 such that, for every edge,

    K >= B_e/(omega_w-A_e omega_v),

and put h_v=K omega_v. Denominators are strictly positive. More efficient
rational caps h_v can be supplied directly, provided

    F_e(h_v) <= h_w.                              (3.1)

For a source above its vertex cap, use

    R(v,x)=(x-h_v)/omega_v.

Each admitted edge then satisfies

    R(w,F_e(x)) <= q R(v,x).                       (3.2)

This is a common rank even when an individual A_e exceeds 1. Let

    delta = min_v (floor(h_v)+1-h_v)/omega_v > 0.

For an integer source above its cap, R>=delta. Thus any prefix of k edges
whose endpoints all stay above their caps obeys

    delta <= q^k R(v_0,x_0).                       (3.3)

It has at most floor(log(R(v_0,x_0)/delta)/log(1/q)) edges. Including a
possible final edge into the core adds at most one edge to this deadline;
multiplication by the largest edge length bounds total Terras time. The
inequality can also be checked using exact rational powers, without logs.
A path must enter the finite core 0<x<=h_v or leave the guarded graph.
An exit is a missing certificate, not a claim of convergence.

This controls every valuation recharge along an admitted path. Endpoint
height is at most

    max_v h_v + (max_v omega_v) R(v_0,x_0).

For a word center C_i/D_i distinct from that integer endpoint,
v_2(x-p_i)<=log_2|D_i x-C_i|. Hence every finite local rank is O(log(x_0+1)),
and the number of switches before core/exit is O(log(x_0+1)), with constants
depending only on the certificate. One can enlarge caps to include any
positive integral centers, so that infinite valuations occur only in the
finite core. The sum of positive rank jumps is consequently
O((log(x_0+1))^2); the bounds concern a fixed supplied source, unlike the
varying-source CRT examples of section 2.

**Conditional coverage corollary.** If every positive source can be assigned
a vertex, every source above its cap has a guarded outgoing edge, and all
positive integers in the finite core have actual receipts to 1, then all
positive sources converge. These are explicit additional premises. The
height theorem proves termination after those premises; it does not prove
the guard-coverage premise.

## 4. A genuinely expanding oriented return of adjacent pairs

In addition to the previously proved returns

    F(x)=(3x+1)/4,             x=1 (mod 4),
    H(x)=(81x+101)/512,        x=43 (mod 512),

consider

    J(x)=(4782969x+10623730)/4194304,
    x=1817598 (mod 4194304).                       (4.1)

The exact lower and upper words on this cylinder are respectively

    0111111111001101101000,
    1111111111011000000011.

Each has length 22 and odd count 14. Their affine offsets are 10623730
and 10035065, whose difference is 4194304-4782969. Consequently the two
endpoints remain one apart after the return. Checking the words at the
stated residue fixes their full cylinder; the displayed affine identity
then proves the adjacent return for every positive lift.

Here

    lambda_J=4782969/4194304 > 1,
    p_J=-2124746/117733 < 0.

Thus J strictly increases every positive source. The center separations
are v_2(p_F-p_H)=1 and v_2(p_F-p_J)=v_2(p_H-p_J)=0. Arbitrarily many J
repetitions have positive integer representatives, but a fixed positive
source has only floor(v_2(x-p_J)/22) repetitions.

The following reset is genuinely descending:

    H(J(x)) = (387420489x+1284146834)/2147483648,
    p_HJ=1284146834/1760063159 < 1.                (4.2)

Therefore H(J(x))<x for every positive integer x>1 in its lawful cylinder.
The source guard is exact, not inferred from the contraction:

    x=2119941118 (mod 2147483648).                 (4.3)

Indeed write x=1817598+4194304t. Then
J(x)=2072698+4782969t, which is 43 modulo 512 exactly when t=505 modulo
512. At the smallest source in (4.3), the actual receipt is

    2119941118 →[22 steps] 2417472043
               →[9 steps] 382451632 < 2119941118.

J alone does not force H: x=1817598 maps to 2072698, outside H's cylinder.
Nor does it force any reset before another J: x=4555502500862 follows JJ.
These are exact positive controls against silently asserting source entry.

## 5. A rational two-vertex certificate, and its exact boundary

Use vertices A,B with these edges:

- F and H from A to A;
- J from A to B;
- F from B to B;
- H from B to A.

Retain each edge's actual word cylinder. In particular, a B source in J's
cylinder has no edge here. Set

    q=3/4,   (omega_A,omega_B)=(1,2),
    (h_A,h_B)=(1,4).

Every weighted slope inequality and every cap inequality (3.1) holds
exactly. The only expansive edge J has

    lambda_J/2 < 3/4,    J(1)=15406699/4194304 < 4,

and H resets B with

    2 lambda_H < 3/4,    H(4)=425/512 < 1.

F's cap images are F(1)=1 and F(4)=13/4<4. All positive core sources are
among 1,2,3,4 and have literal receipts to 1. Here delta=1/2, so a path
starting at A above the core has at most

    floor(log(2(x_0-1))/log(4/3))

edges that remain above the core. This graph permits arbitrary repeated F
blocks and any number of J expansions separated by H resets. Their local
valuation ranks may recharge, but the common height cannot.

The graph's restriction is essential. It is not a source-covering graph.
For example its J→H cylinder is only one residue class among 512 lifts of
the J cylinder. Lawful J→J paths are deliberately not certified by it.

There is a precise finite-state barrier to repairing this by adding enough
states while retaining only contracting cycles. For every L>=1, the
parity cylinder J^L has a positive integer representative. Any finite graph
that represents every one of these finite J-only paths contains a J-only
path longer than its number of vertices. Repeated vertices then give a
cycle with multiplier lambda_J^ell>1. Thus no finite cycle-contracting graph
can represent all lawful finite J histories. Additional unbounded
source information, such as the exact valuation budget, or a different
source-sensitive termination proof is indispensable.

The progress is therefore a sharp conditional mechanism and a sharp
limitation: common weighted height controls arbitrary recharge on a
source-guarded graph; the newly exhibited lawful expansive return proves
why source coverage cannot be smuggled into that finite graph.

## 6. Exact controls

The companion [script](../experiments/collatz_switching_height_certificate_20261019.py)
performs **2,701 explicit checks**. It reconstructs the return words, coefficients and full cylinder residues;
checks the resonance identity and CRT recharge; verifies every rational
height inequality; and replays lawful mixed paths on the two-vertex graph.
It checks the explicit JH descent and JJ exclusion, and constructs long
J-only positive prefixes to test the finite-state boundary. These are
focused theorem controls, not a new census or an assumed global tail.
The [normal receipt](collatz_switching_height_certificate_20261019.json)
and [optimized receipt](collatz_switching_height_certificate_20261019.optimized.json)
are required to agree byte-for-byte.

```sh
python3 05-knowledge/experiments/collatz_switching_height_certificate_20261019.py > /tmp/switch-height.json
python3 -O 05-knowledge/experiments/collatz_switching_height_certificate_20261019.py > /tmp/switch-height-O.json
cmp /tmp/switch-height.json /tmp/switch-height-O.json
cmp /tmp/switch-height.json 05-knowledge/results/collatz_switching_height_certificate_20261019.json
```
