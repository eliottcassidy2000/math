# Positive path-cover operator for triangle substitution

Status: proved exact formulas and a finite modular certificate. The constructions are strong tournaments, not a classification of their Hamilton-path spectrum. No claim of novelty is made.

## The full bin profile

For a finite tournament X and integer r>=0, let F_r(X) count ordered r-tuples of vertex-disjoint directed paths whose union is V(X), allowing empty paths. Thus F_0(X)=0 for nonempty X and F_r(empty)=1. Put h(X)=F_1(X). For nonempty X, S158 uses F_2(X)=C_2(X)+2h(X).

Define G(X) by substituting X and singleton vertices a,b into a directed triangle X→a→b→X. Then G(X) is strongly connected for nonempty X, and

    F_r(G(X)) = r F_(r+1)(X) + r^2 F_r(X) - r F_(r-1)(X),  r>=1.   (1)

This is valid for empty X too, under the above conventions. Here is a direct proof, avoiding any experimental extrapolation. In an r-slot path cover of G(X), consider the slots containing a,b. If they are in the same slot, there are r choices. That slot contains an induced G(S), S subset V(X), and the rest of X is split among the other r-1 slots. For nonempty S, S158's triangle bijection gives h(G(S))=h(S)+F_2(S). For S empty, h(G(empty))=1 instead of 2. Summing produces r[F_(r+1)(X)+F_r(X)-F_(r-1)(X)]. If the two singleton vertices occupy distinct slots, there are r(r-1) choices. In each such slot, the singleton is forced to the appropriate end of the X-path, yielding r(r-1)F_r(X). Add both cases.

For any ordinal sum X⊕Y, every slot path cuts uniquely into an X-path followed by a Y-path, hence F_r(X⊕Y)=F_r(X)F_r(Y). Also F_r(T_n)=r^n. These operations provide an exact arithmetic calculus on full profiles, with triangle substitution supplying strong outputs. The r=1 case recovers h(G(X))=h(X)+F_2(X).

## A subtraction-free species operator

Let p_k(X) denote the number of unordered k-block path covers, p_k=C_k/k!, with p_0(empty)=1 and p_0(X)=0 for nonempty X. Let

    P_X(z) = Σ_(k>=0) p_k(X) z^k,   F_r(X)=Σ_k p_k(X) (r)_k,

where (r)_k=r(r-1)...(r-k+1). Formula (1) is equivalent to

    P_G = (z+z²)P_X'' + (3z+2z²)P_X' + (z+z²)P_X.   (2)

Equivalently, p_j (r)_j contributes

    j(j-1)(r)_(j-1) + j(j+2)(r)_j
      + (2j+1)(r)_(j+1) + (r)_(j+2),

and consequently

    p'_k = k(k+1)p_(k+1) + k(k+2)p_k + (2k-1)p_(k-1) + p_(k-2),

where out-of-range terms are zero. At k=0 the negative-looking coefficient multiplies p_(-1)=0. This transfer is entrywise nonnegative. It encodes singleton insertion and collision as a finite-band positive operator on the carrier's path-cover profile. Differential operator (2), for example, maps z to z³+3z²+3z, the three cover counts of a directed triangle.

## Exact tropical refinement

Mark any original directed arcs of X, and weight each path cover by t to the number of marked arcs it uses. Assign both new singleton vertices and all new arcs weight 1. Define p_k(X;t) as the resulting polynomial in N[t]. The same proof establishes (2) coefficientwise in t: restrictions to X preserve all used marked arcs, while all new edges are unmarked. It can be repeated for any nesting of G and arbitrary separately marked initial blocks.

Let ell_k be the minimum exponent of t with nonzero coefficient in p_k, with min(empty)=infinity. Because all four terms in the recurrence have nonnegative coefficients, exact min-plus propagation is

    ell'_k = min(ell_(k+1) if k(k+1)>0,
                 ell_k if k(k+2)>0,
                 ell_(k-1) if 2k-1>0,
                 ell_(k-2)),

omitting terms with nonexistent indices or zero multiplier. Likewise, maximum exponents propagate by max. Thus these tropical states compute exact optimal numbers of marked original arcs among k-path covers of the literal enlarged tournament. No use is made of cancellation or of a matchgate/Pfaffian solvability claim. One must not substitute p-adic valuations for ell: positive integer sums can cancel modulo a prime and increase p-adic valuation.

## Iterated strong family and a modular invariant

Starting with X_0 a singleton, define X_(m+1)=G(X_m). Each X_m for m>=1 is strong and has 2m+1 vertices. Iterating (2) gives

    m              0  1   2    3     4      5       6
    h(X_m)         1  3  15  123  1479  23955  496479.

For every m>=1, v_3(h(X_m))=1. This infinite assertion has a finite, exact certificate: modulo 9 the vector (p_1,...,p_8) is closed under the recurrence, since the only omitted influence is 8*9*p_9=0 mod 9. Starting with (1,0,...,0), the values are:

    m  p1 p2 p3 p4 p5 p6 p7 p8  (mod 9)
    0   1  0  0  0  0  0  0  0
    1   3  3  1  0  0  0  0  0
    2   6  3  6  1  1  0  0  0
    3   6  6  6  8  5  3  1  0
    4   3  3  6  7  1  6  8  0
    5   6  6  3  8  5  3  7  0
    6   3  3  6  4  7  6  8  0
    7   6  6  3  2  8  3  4  0
    8   3  3  6  1  4  6  2  0
    9   6  6  3  5  2  3  1  0
    10  3  3  6  7  1  6  5  0
    11  6  6  3  8  5  3  7  0

Rows 5 and 11 coincide, so the deterministic state repeats with period dividing six. As p1 is 3 or 6 throughout the transient and cycle, the valuation follows. In fact h(X_m)≡3 mod 9 for even m>=4 and ≡6 mod 9 for odd m>=3, while h(X_1)≡3 and h(X_2)≡6. This assertion pertains to this family only.

## Reproduction

`python3 check.py` and `python3 -O check.py` enumerate all 75 labeled tournaments on one through four vertices, construct the literal G(X), independently compute induced-subgraph Hamilton counts by subset/endpoint dynamic programming, and compare actual slot assignments against (1) and the positive transfer for r=1,2,3. The script checks the first seven iterated counts and the displayed mod-9 cycle using explicit conditions that remain active under `-O`. The general identities follow from the proofs above; the finite modular table is itself a complete induction certificate.
