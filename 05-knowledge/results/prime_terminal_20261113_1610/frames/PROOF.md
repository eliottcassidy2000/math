# Exact ternary targets through changing module charts

**PROVED scoped response constructions and obstructions; FINITE-EXACT
order-nine classification; Chan--Pak remains OPEN.** A two-vertex shift in
the growing module supplies an actual infinite +2 family even though its
unshifted response is impossible. The +4 branch is still a one-shot success.

## Inheritance and the retained target

The incoming [shared-source frames](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md)
already produce both 3h+2 and 3h+4 on h=15*2^m-1. Their complete matching-guarded
frame catalogue through order eight is the input here, not a new census.
The strongest hostile is their arbitrary-module identity: the two digits
become 2H(M),4H(M). The corrected near miss is to treat an impossible fixed
module response as an impossible target ray. The least-used coordinate is
the division between fixed core vertices and the growing module.

Write h=H(M), p=P_M(2), and K_j(M) for ordered covers by j nonempty paths.
Thus K_1=h and K_2=p-2h. A rooted frame with at most two visits has response

    H(F[r<-M])=a_1 h+a_2 K_2(M)=a p+b h,
    a=a_2, b=a_1-2a_2.                                  (1)

The [two-run cone](../absorption/TWO_RUN_CONE.md), C1--C4, proves the exact
matching criterion and, when r is cyclic,

    a+2 <= a_1 <= 3a,       -a+2 <= b <= a.               (2)

For a transitive m-vertex module, h=1,p=2^m. All graph rows below encode
out-neighbours as bit masks, with zero-based module ports.

## F1. Constant digits in a fixed two-run comparison

Fix two two-run frames A,B. Suppose their difference on the same nonempty
module M is

    H(B[M])-3H(A[M])=u h+v p=delta,   delta in {2,4}.       (3)

If v is nonzero, the order n of M satisfies

    n+1 <= (|u|+delta)/|v|.                              (4)

Indeed cutting each Hamilton path at any of its n-1 positions injects
(path,cut) into ordered two-path covers. Hence K_2>=(n-1)h and p>=(n+1)h.
Equation (3), divided by h, gives the stated upper bound for p/h.

If v=0, the positive odd integer h divides delta, so h=1. A tournament
has H=1 exactly when it is transitive. Consequently a fixed finite bank of
two-run comparisons cannot supply constant digits on nontransitive modules
of unbounded order. This is a corollary of the inherited cover calculus and
path-cut injection, not a global obstruction to changing-core constructions.
The all-rank retained-module obstruction C5 in the cone packet is stronger
on the specific Q-power family.

## F2. Both digits by a real nontransitive cancellation

Let Q=(6,12,8,1), so h=5,p=32. The following complete run vectors are exact:

| Frame rows | Port | (a_1,a_2) | Count after substituting Q |
| --- | ---: | ---: | ---: |
| (2,4,25,19,3) | 1 | (9,7) | 199 |
| (2,12,73,81,103,79,3) | 5 | (45,17) | 599 |
| (2,12,9,17,103,79,15) | 5 | (41,18) | 601 |

All higher run coefficients vanish. The two differences are

    26h-4p=2,       20h-3p=4.                            (5)

Thus the same actual H199 source count has both prime targets599 and601.
These are actual tournament substitutions, not formal Bezout identities.
Neither identity stays constant along Q^r or arbitrary substitutions.

## F3. A shifted all-size +2 construction

Use the strong frames

    S=(16,5,9,35,46,7),                      port0,
    F=(258,4,265,51,39,135,191,287,122),      port1.

Their complete run vectors are respectively (37,27) and (599,324).
Therefore for every integer k>=3,

    h_k=H(S[TT_k])=27*2^k-17,
    H(F[TT_(k-2)])=324*2^(k-2)-49=3h_k+2.              (6)

The source and target orders are k+5 and k+6. The first instances are
199->599,415->1247,847->2543. Each is a full Hamiltonian count. The
construction selects a different core and uses two fewer vertices in the
uniform module. It does not assert an induced one-vertex extension.
The companion [rearming packet](../repairs/REARMING.md) records the explicit
marked embedding test: none of the eighteen order-nine target presentations
found here contains the required marked eight-vertex source presentation.

The source formula can also be recovered by a cyclic triangle of modules:
for C3[D,Q,1], direct run enumeration gives

    H=(P_D(2)-H(D))*32+(3H(D)-P_D(2))*5
     =27P_D(2)-17H(D).

Taking D=TT_k gives the same scalar ray; the explicit S frame above is the
literal certificate used by the checker. No new general triangle formula
is claimed beyond the existing cover calculus.

## F4. Why the unshifted chart really fails

There is no rooted tournament frame whose uniform transitive-module
response is either

    81*2^m-49,       81*2^m-47                           (7)

for all positive m, or for all sufficiently large m.

First, the nonnegative ordered-run expansion and
K_j(TT_m)=j!S(m,j) imply that a largest nonzero run coefficient at j>=3
has positive leading exponential j^m. It cannot agree with an affine
2^m response on an infinite tail. Thus (7) forces at most two visits,
with a_2=81 and a_1=113 or115.

Second, factor the frame into its ordered strong components. The root
component carries every repeat of the root; all other components multiply
every run coefficient by their H. Since gcd(81,49)=gcd(81,47)=1, all those
other factors have H=1 and hence are singletons. Delete that padding.
The root is cyclic because a_2>0, leaving a strong frame of H113 or115.

The exact strong minimum in Arthur Busch,
[Theorem1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf),
is f(3j)=3*5^(j-1), f(3j+1)=5^j, f(3j+2)=9*5^(j-1).
In particular f(10)=125, and the minimum is increasing. The remaining
strong frame has at most nine vertices. The inherited complete order-eight
catalogue excludes both profiles through eight. The exact extension
calculation below excludes order nine. This closes the all-order argument
for (7).

This does not exclude the shifted chart in (6): its actual frame response
has amplitude324 and base count599. Nor does it exclude higher-order
frames for the shifts (162,-49), (162,-47), or (324,-47). Their reported
absence is only through order nine.

As another application of (2), the proposed source profile
(a_1,a_2)=(95,11), obtained in the H9011 inverse-module experiment, is
impossible in every order because95>33. Its impossibility is structural;
the bounded census is merely a corroborating control. H9011 itself is
attained, as [MODULES](../arithmetic/MODULES.md) constructs explicitly.

## F5. Exact search universe and hostile controls

The input is all6880 ordinary order-eight tournaments in the pinned
`multicore_metagraph_20261113/rooted_growth.json.gz`. Enumerate every one-
vertex attachment to representatives whose H<=601. The other97 parents
cannot occur by induced-subtournament monotonicity of H, inherited from
the positive odd-cycle formula, THM-002 in
[THM-002-ocf.md](../../../../01-canon/theorems/THM-002-ocf.md).
Every order-nine tournament with H<=601 is represented, possibly repeatedly,
by this extension method. We retain strong results and each matching-guarded root.

The exact universe has6783 retained parents,1736448 extensions, and2478
strong presentations with H in{95,113,115,275,277,599,601}. Of the5372
guarded rooted presentations, eighteen have a desired profile, all
(a_1,a_2)=(599,324). The 80 presentations with H95,113,115 and all534
guarded roots therein are independently recounted in Python using a
backward endpoint recurrence, separate from the C++ forward recurrence.
All eighteen positive presentations are also recounted literally.

For the +4 continuation, a further bounded experiment reverses every set
of at most four arcs in the ten-vertex H601 witness of F2. It tests164221
graphs;141879 have at least one two-run root,257 retain H601, and none
has a two-run root with a_2=324. This is a local stopping reason only.
It does not imply that the missing shifted +4 frame does not exist.

The source bank and programs are pinned in the receipts. No isomorphism
count is inferred from presentation counts. The independent skeleton
enumerator checks every run coefficient of the named positive frames,
so the all-k identities do not depend on interpolation from small k.

Reproduce from the worktree root:

    python -B 05-knowledge/results/prime_terminal_20261113_1610/frames/search.py
    python -B 05-knowledge/results/prime_terminal_20261113_1610/frames/check.py
    python -B -O 05-knowledge/results/prime_terminal_20261113_1610/frames/check.py
    python -B 05-knowledge/results/prime_terminal_20261113_1610/frames/neutral_search.py

## The next constructive question

The shift in (6) is a genuine escape from a proved unshifted obstruction.
The precise residual target is an actual frame with response324*2^m-47,
or another chart for81*2^k-47. A successful +4 chart would provide both
children on a third dyadic source family, but would still need a proof
that source selection covers every critical prime or seven-prime target.
The module-preservation obstruction explains why changing only the
exterior of a fixed nonunit module is sometimes insufficient: the split
between module and core is itself a resource to change.
