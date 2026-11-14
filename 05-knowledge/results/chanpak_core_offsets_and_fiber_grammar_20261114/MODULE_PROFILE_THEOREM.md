# Exact modular Hamiltonian fibers and finite continuation states

Status: PROVED consequences of the cited path-cover/gate/modular identities; FINITE-EXACT independent census through order six. Local research note. No claim of literature priority.

## Definitions and known inputs

For a nonempty tournament T on n vertices, pc_T(k) is the number of spanning covers by k **unordered**, nonempty, vertex-disjoint directed paths, including singleton paths. Put F_T(k)=k! pc_T(k), and

    P_T(q) = sum_k pc_T(k) (q)_k,

where (q)_k is the falling factorial. In particular P_T(1)=H(T). These are ordinary path counts, not tournament-class counts. The power form of the known odd-cycle formula is

    P_T(q)=q^n+sum_{j=1}^{floor((n-1)/2)} w_j(T) q^(n-2j).

Here w_j is the sum of 2^|K| over collections K of pairwise vertex-disjoint directed odd cycles with sum_(C in K)(|C|-1)=2j. Thus w_j are nonnegative even integers. The empty collection supplies the leading term, and there is no constant term.

Sources at the inspected repository head 8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90:
- [THM-3121, exact quotient walk/path-cover kernel](https://github.com/eliottcassidy2000/math/blob/8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
- [THM-4900, known gate operator and size-weighted odd-cycle formula, with earlier credits](https://github.com/eliottcassidy2000/math/blob/8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90/01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md).
- [THM-4803, canonical modular census and its corrected automorphism boundary](https://github.com/eliottcassidy2000/math/blob/8f3e642bd29de5bfbb7b0bc1fa558fbfe18b6b90/01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md).

## 1. Coarsest unmarked module-substitution state

**Theorem.** For tournaments T,U of the same order n, the following are equivalent:
1. P_T=P_U, equivalently their full path-cover profiles agree.
2. Every one-hole tournament module-substitution context has the same H after inserting T or U.
3. H(Gamma^r T)=H(Gamma^r U) for r=0,...,d-1, where d=floor((n-1)/2) and Gamma(X)=C3[X,1,1]. For n=1,2 the profile is unique and no probes are needed.

**Proof.** The exact kernel in section 2 proves (1) implies (2); (2) implies (3) because nested gates are legitimate module contexts. The known gate operator is

    (LP)(q)=q(P(q+1)+qP(q)-P(q-1)),  P_(Gamma T)=LP_T.

Expanding (L^rP)(1) uses samples P(1),...,P(r+1). Terms reaching zero vanish, and reaching r+1 requires all positive shifts; its coefficient is r!. Thus the gate responses recover P(1),...,P(d) triangularly. The polynomial (P(q)-q^n)/q^(n-2d) has degree at most d-1 in q^2; interpolation at the distinct squares 1^2,...,d^2 recovers it. This proves (3) implies (1).

An independent family of d probes is H(T) and H(C3[T,TT_k,TT_k]), 1<=k<=d-1. The latter is triangular in F_T(1),...,F_T(k+1), with last coefficient (k!)^2. Indeed, the C3 walk kernel allows run counts only within distance one; the maximal run count k+1 in T forces exactly k in each transitive block. This also recovers P(1),...,P(d).

This is coarsest observational equivalence, not a claim about the fewest bits in an encoding or the fewest probes among all possible contexts. It applies to whole-module attachments, NOT arbitrary attachments to individual vertices, rooted/endpoint states, automorphisms, or metagraph arc reversals.

**Smallest fixed-order scalar-H hostile.** Index pairs (i,j), i<j, lexicographically from bit zero, with a set bit reversing i->j. The order-five masks 69 and 76 have

    P_69=q^5+8q^3+6q,   pc_69=(15,39,33,10,1),
    P_76=q^5+10q^3+4q,  pc_76=(15,45,35,10,1).

Both have H=15, but H(Gamma T) is 123 and 135. No such failure is possible at fixed n<=4, because d<=1. This disproves scalar-H substitution closure at its first possible order.

## 2. Exact class histogram recursion, not merely profile realization

Let a_(n,p) be the number of tournament **isomorphism classes** of order n with P=p, and s_(n,p) the analogous strongly connected count. Write a_(0,1)=1. Every tournament has a unique ordered SCC decomposition, and

    P_(A->B)(q)=P_A(q) P_B(q).

For integer q this follows by assigning each path component injectively to q ordered bins: a bin is a directed path, possibly empty. On an ordinal join, the A and B contents of each bin concatenate uniquely, so the choices factor. Polynomial equality follows. Therefore

    a_(n,p)=sum_(k=1..n) sum_(u v=p) s_(k,u) a_(n-k,v).       (1)

For each isomorphism type of prime quotient Q of order r>=3, let G=Aut(Q). A strongly connected tournament is either a singleton or has a unique prime quotient with nonempty arbitrary tournament children, up to the action of G. Here prime means no nontrivial proper module; it includes C3.

The exact substitution profile Phi_Q is supplied by THM-3121. With A_Q its directed adjacency matrix and X=diag(x_i), put

    W_Q(x)=1^T X (I-A_Q X)^(-1) 1.

For child profiles p_i and F_i(c)=c! pc_i(c),

    pc_(Phi_Q(p_i))(d)
       = (1/d!) sum_c [x^c] W_Q(x)^d prod_i F_i(c_i).        (2)

All child contents c_i are positive. The factor d! divides only after actual path components have been decorated: arbitrary bare quotient words need not have a free permutation action. Formula (2) then converts from pc to P in the falling-factorial basis.

**Profile-refined Burnside theorem.** For sigma in G, assign to every cycle C of sigma a size/profile pair (n_C,p_C), with total sum_C |C| n_C=n. Expand p_C at all quotient vertices in that cycle. Add the weight

    prod_(cycles C of sigma) a_(n_C,p_C)

to the resulting profile Phi_Q, and average these histograms over sigma in G. Sum over prime Q and add the singleton to obtain s_(n,p). Together with (1), this is a triangular exact recursion for every Hamiltonian-count fiber:

    number of order-n isoclasses with H=h = sum_(p:p(1)=h) a_(n,p).   (3)

**Why the weights count classes correctly.** A sigma-fixed decoration requires the SAME ACTUAL CHILD ISOTYPE around each sigma-cycle. If c is the number of isotypes with profile p, a cycle contributes c choices, NOT c^|C|. Distinct cycles choose independently, possibly choosing the same type. Grouping choices by profile loses no multiplicity because (2) depends only on profiles. Thus child profiles need not determine child automorphism groups. Computing an individual composite automorphism group still requires more information, as corrected THM-4803 emphasizes.

The quotient data must retain enough information to compute Phi_Q for the expanded cycle assignments. An unmarked cycle-length census alone is not justified for this nonlinear refinement. Keeping Q, its actual automorphisms, and their vertex-cycle partitions is sufficient. This note does not claim to have proved a smallest counterexample to cycle-length-only compression.

For Q=C3 the formula becomes

    ( sum_(p1,p2,p3) a_p1 a_p2 a_p3 [Phi_C3(p1,p2,p3)]
      + 2 sum_p a_p [Phi_C3(p,p,p)] ) / 3,

with sizes retained. In particular, c distinct child isotypes sharing one profile give (c^3+2c)/3 rotation-orbits of three children, not one orbit or c^3/3.

## 3. Finite continuation states at fixed H

Put h=H(T), m=(h-1)/2. All h are positive odd and sum_j w_j/2=m.

An elementary cutoff is w_j=0 when j>m. A contributing packing with lengths 2l_i+1 gives h>=prod_i H(T[V_i])>=prod_i(2l_i+1)>=1+2j: opening each chosen cycle gives its length many Hamiltonian paths, and the positive odd-cycle formula includes all unions of independently chosen internal packings.

There is a sharper cutoff using the Busch bound for a strongly connected component C of order t:

    H(C)>=5^((t-1)/3).

Each directed cycle lies in one SCC. A disjoint nonempty cycle packing within a t-vertex SCC has total excess sum(|cycle|-1)<=t-1. Since H(T) is the product of its SCC counts, every nonzero w_j satisfies

    5^(2j) <= h^3.                                           (4)

This is an exact integer test, with no logarithmic rounding. Thus e is the largest j<=d satisfying (4); only e gate/collar responses, including already known H, are needed. The normalized weight vector has at most

    binom(m+e-1,e-1)

possible values when e>=1, independently of n except for the cap d. At h=1 the unique profile is q^n. The bound counts POSSIBLE PROFILES, usually including unrealizable ones; it does not bound the number of classes by the same expression.

Sharpness of (4): by [Busch’s attained bound at orders 1 modulo 3](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), an order-seven strong extremal has H=25; by Camion’s theorem it has a spanning 7-cycle. The ordinal sum of r copies has h=25^r and a packing of r such cycles at j=3r, giving equality. This uses the supplied size-seven extremal and Busch theorem, not a new proof of their existence/bound in this note. The elementary cutoff above is independent of these inputs.

## 4. Exact computation and falsification guards

Run `python3 module_profile_verify.py` (also valid under `python3 -O`). The verifier independently generates isotypes by all vertex relabelings through n=6. It obtains pc from all permutations and cut sets, checks H by a separate subset DP, and compares the ENTIRE profile-class and strong-profile histograms with (1)-(2) and Burnside. This is not a comparison of scalar totals only.

Results through n=6:

    n:             1  2  3  4   5   6
    A000568:       1  1  2  4  12  56
    strong:        1  0  1  1   6  35
    profiles:      1  1  2  3   8  24
    prime classes: 0  0  1  0   3  15

There are 94 additional independent direct substitution controls. The mask69/76 counterexample is checked by the subset DP on the actual seven-vertex gates. The C3-only subclass recurrence gives counts at n=0,...,9:

    1,1,1,2,4,9,24,65,182,528.

These are NOT A000568 beyond order four. They are checked against a separate scalar Burnside/SCC recurrence. Every logical check uses explicit exceptions and therefore survives optimized Python. The recorded first complete run took approximately four seconds in this environment; this is a replay duration, not an asymptotic speedup benchmark.

## 5. What this does and does not solve

The state theorem gives a sharp qualitative boundary: full P is exactly the unmarked modular continuation information; H alone first fails at order five. Profile aggregation can share computation between many isotypes, and (4) bounds the normalized state space for a fixed target h. The positive kernel also gives the useful pruning bound H(Q[S_i])>=H(Q) prod_i H(S_i), by using single-run quotient Hamiltonian paths.

The prime-quotient census remains an input. This recursion does not turn arbitrary tournament counting into a polynomial-time algorithm, enumerate new prime quotients for free, classify individual isotypes from P, or prove that any proposed h is attained. It does not identify isotype counts with metagraph lines, H values, or extremal H counts. Toward Chan-Pak, it supplies a finite, rigorously closed continuation-state grammar at each target; the missing step is proving realizability/coverage of the desired states or targets, not collecting additional finite examples.
