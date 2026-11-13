# A complete intrinsic inverse and the first ready hole on the 405 ray

**PROVED reduction / FINITE-EXACT all-order consequence.**
There is no tournament of Hamilton count

    N = 53,084,231 = 405·2^17+71

with an indegree-one root r satisfying
H(T−r)−s_T(r)=1. This is the first such failure on the specified ray:
the targets for k=1,...,16 have explicit unit-ready sources.
**N itself is attained. Chan–Pak remains OPEN.**

The exclusion is over every tournament order. A new intrinsic path-cover
bound reduces this marked question to a complete finite core search.
Two independent algorithms then exhaust that search. Neither the 881
fixed-block obstruction nor a guessed core-order cutoff is used.

## 1. Inheritance and concurrent overlap

The complete strong unique-start grammar is inherited from
[FUNNELS F1–F3](../polarized_path_calculus_20261103/FUNNELS.md).
Its six-state compiler is inherited from
[TARGET_RETURN §3](../berggren_target_return_20261111/reset/TARGET_RETURN.md).
The binary-shell recurrence and section inverse were recovered in
[the H13 grammar](../h13_paired_descent_20261111/square_graceful/GRAMMAR.md)
and the earlier recursive-seed switching inverse.

During this session, incoming commit9b8aa35188 supplied the independently
developed [arbitrary-tail target interface](../changing_core_returns_20261113/returns/PROOF.md).
That note already proves the observer formula below, supplies 29 of the
first 30 ray targets, and excludes k17 over all unique-start cores through
order17. We credit that overlap explicitly. Our additional result is the
intrinsic complete core bound and its exhaustive all-order consequence;
the formula and square-root-scale inversion are not claimed as new.

## 2. Exact target observer with a physical source address

For a strong rooted unique-start core C, or the singleton, retain

    (h,e,s,p,f,t) = (P(1),E(1),S(2),P(2),E(2),S(3)),  G=s+p.

P(j) counts ordered j-bin path covers, allowing empty bins. E and S
require the first bin to be nonempty and end or start at the root.
Encode a binary tail word low bit first by x≥0; x=0 is the empty word.
Its nonempty words end in1. The literal funnel adds that transitive tail
after C and a new root which beats only the old root inside C and beats
exactly the zero-bit tail vertices. Append a final vertex beaten only
by this new root. The result is a unit-ready strong tournament, except
that singleton C plus empty word is outside the strong constructor.

The full Hamilton count, including every path sector, is

    Q_C(x) = h g(x)+Gx+e(x mod2)+G−e−f−1,                 (1)
    g(0)=g(1)=1,
    g(2a+1)=g(a),  g(2a)=g(a)+a+(a mod2).

Indeed the tail coefficients are L2=x+1, B1=g(x), D1=1−(x mod2).
These are the inherited compiler identities. They retain an actual graph
constructor, not merely a numerical representation.

The inverse tests both binary sections. For a g(x)+b x+c(x mod2)+d=N,

    even: (a,b,c,d) → (a,2b+a,a,d),
    odd:  (a,b,c,d) → (a,2b,0,d+b+c).                   (2)

The bounds 1≤g(x)≤x+1 give the necessary integer interval

    max(0,ceil((N−d−a−c)/(a+b))) ≤ x ≤ floor((N−d−a)/b). (3)

Empty intervals are rejected; intervals of width at most3 are directly
evaluated. Binary sections partition all nonnegative addresses. This
does not assume monotonicity of g, which has large downward jumps.
For fixed coefficients the search has square-root-scale branch count,
with logarithmic-cost integer evaluation. That computational bound is
different from the previous fixed-seed construction's sparse target image.

## 3. An intrinsic complete core bound

**Lemma.** Every core C occurring in (1) at target N satisfies

    P_C(2) ≤ 2N,       |C| ≤ floor(log2(2N)).             (4)

Proof. Since g(x)≥1 and the other x-dependent terms are nonnegative,
Q_C(x)≥Q_C(0). Moreover

    Q_C(0)−(p−f)=h−e+s−1 ≥0.

Here e≤h and s≥1. Swapping the two bins injects the covers whose first
bin ends at the root into a disjoint subset of all two-bin covers, so
f≤p/2. Thus Q_C(x)≥p/2. Every one of the 2^|C| assignments of vertices
to two bins has at least one pair of directed Hamilton paths, including
empty-bin conventions. Hence p≥2^|C|, proving (4).

For a nonempty outer word, x≥1 and Q_C(x)−Q_C(0)≥G, giving the optional
stronger bound Q_C(x)≥3p/2+1. The code uses the uniform weaker bound.

The complete unique-start decomposition supplies every strong unit-ready
source: deleting the ready root leaves a strong unique-start parent,
whose last funnel has one core C and one such word. The inherited weak
classification adds only a transitive suffix at unchanged H and defect.
Thus weak tournaments create no extra target values. No Busch order
bound or extrapolation from a finite ready-source bank is needed here.

## 4. Safe pruning at every recursive core level

All ancestors of a retained core also satisfy p≤2N. One can see this
directly from vertex insertion: for each old two-bin cover and each
choice of a bin, deterministically insert the new vertex into that bin's
path. Deleting it recovers both the old cover and the chosen bin. Hence
adjoining any vertex multiplies P(2) by at least2.

The implementation additionally uses a stronger tail-length bound. The
six-state compiler gives

    p' = 2[L3 t+B2 p−D2 f].

The independently proved cover injections give t=S(3)≥p/2 and f≤p/2.
For a legal nonempty word of length l,

    2B2−D2≥1,       L3≥2^l+3^(l−1),
    p'≥(L3+1)p≥(2^l+3^(l−1)+1)p.                    (5)

For the empty word p'≥2p. The first inequality follows by choosing
the last bit in the before bin, giving B2≥2^(l−1), while
D2≤2^l−1. The final1 contributes 3^(l−1) to L3 in addition to its
2^l term. These are positive inequalities on actual covers.

The search keeps all grammar states below the bound, not only distinct
scalar observers. Source addresses and filtered-parent indices remain
available. Rejecting a partial state above the bound cannot hide a later
accepted core, by monotonicity. The order cap in (4) is absolute.

## 5. The exact exclusion and its independent path

For N=53,084,231, the bounds are p≤106,168,462 and |C|≤26.
The complete retained-core counts at orders1,3,...,19 are

    1,1,3,8,21,55,144,377,987,2584,6765,17711,46368,
    121393,317811,832040,2108788,2184969.

Orders20 through26 have no retained cores. There are **5,640,026**
retained cores in total. The producer applies the exact binary inverse
to every core and finds no target address.

The independent [audit](ports/COMPLETE_INVERSE_AUDIT.md) uses a separate
depth-first address traversal and directly evaluates **every** allowed
outer address, rather than importing the producer's inverse or core bank.
It tests **74,845,805** addresses, matches every retained-core layer,
and obtains zero hits. Its separate unpruned small-core controls agree.
Positive881 returns exactly the two inherited ready addresses; the known
1997 hostile is independently excluded with287 retained cores.

This promotes the incoming bounded-order miss at k17 to an all-order
marked obstruction. It refutes an all-k unit-ready405 ray in this exact
indegree-one sense. It does not refute a changed selected interface,
nonunit defects, a direct larger surgery, or Chan–Pak.

## 6. Positive target-directed sources and the honest remaining gap

`target_inverse_receipt.json` retains a source for every k=1,...,44
except17 over a complete core bank through12, with no tail-length cutoff.
The adjacent source certificate records full chronological words and
outgoing graph rows, and checks the compiler and small literal controls.
These are finite constructions, not an all-k selector theorem. The first
16 positive witnesses together with §5 establish the first failure.

The target immediately after the hole, 3N+2=159,252,695, has a ready
singleton-core source with outer address45,157,534. Thus the missing
interface at N is not a missing next target. [AT_TARGET](AT_TARGET.md)
records a physical selective-port bridge and nonunit-defect exits.

Reproduce from the repository root by compiling the two C++ files and
running the commands in REVIEW.md. The proof arguments establish the
infinite scope; receipts certify the stated finite reductions. No Lean
formalization or new global theorem identifier is claimed.
