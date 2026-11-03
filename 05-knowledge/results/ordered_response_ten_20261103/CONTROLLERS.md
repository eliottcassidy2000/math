# A terminating arithmetic program on actual metagraph routes

Status: **PROVED scoped**. The program is an explicit terminating factorization normalizer, not a Collatz convergence argument or a universal machine. Its strict height suggests a design method: choose semantics first, then orient value-preserving relations by a well-founded rank.

## C1. Four FRACTRAN fractions with a complete trace law

Encode a factorization a^alpha c^gamma d^delta e^epsilon by

    N=2^alpha 3^gamma 5^delta 7^epsilon,
    V(N)=3^alpha 49^gamma 343^delta 63^epsilon.

Use the ordered fractions

    48/49, 20/21, 36/35, 27/25.

They implement e^2->ca^4, ec->da^2, ed->c^2a^2 and d^2->c^3 respectively. A fraction applies exactly when its denominator divides N. Every application preserves V and strictly decreases `2v_7(N)+v_5(N)`. Other prime factors are inert. Hence every positive integer input terminates, in at most that many steps, with no exceptional-input qualification.

The entire run is explicit. Put u=floor(epsilon/2) and c'=gamma+u. Number rules1,...,4.

| Condition | Trace |
|---|---|
| epsilon even | 1^u 4^floor(delta/2) |
| epsilon odd, c'>0 | 1^u 2 4^floor((delta+1)/2) |
| epsilon odd, c'=0, delta>0 | 3 4^floor((delta-1)/2) |
| epsilon=1, gamma=delta=0 | empty |

The third case forces u=0. Therefore the set of possible traces is exactly

    1* (empty | 2) 4*  union  3 4*.

Every word in this language is attained by the corresponding exponents. There are2n+2 traces of length n>=1 and one empty trace; the generating function is `(1+2z-z^2)/(1-z)^2`. This is a language count, not a count of integer inputs.

## C2. Thirteen feasible snapshots, ten complement-compatible ones

The guard tuple records divisibility by(49,21,35,25). For any proposed enabled set let L be the lcm of its denominators, with lcm(empty)=1. The tuple is feasible iff no disabled denominator divides L. Necessity follows from L|N; sufficiency uses N=L.

There are13 feasible tuples. Exactly these10 have feasible bitwise complements:

    0000 0001 0011 0100 0111
    1000 1011 1100 1110 1111.

Their one-bit adjacency graph has12 edges. Quotienting by bitwise complement gives K_{2,3}, on representatives0000,0001,0011,0100,0111, with bipartition{0000,0011} and the other three. This ten-vertex graph is not the Petersen graph: its degrees already differ. The older ten-address Kneser/Petersen calculation is useful for its insistence on addressed incidence, not a claimed isomorphism.

Realize a tuple as a marked six-vertex tournament with Hamilton path0->1->...->5. Store its four bits on chords(i,i+2), i=0,...,3; orient every other non-path chord by bit0. There are ten non-path pairs in total. Native path complement reverses all these pairs, sending this exact marked source to the bitwise-complement source. The chart is physically closed under native complement precisely on the ten feasible tuples above.

This finite closure does not close arithmetic evolution. Inputs25 and625 share guard0001; rule4 sends them to27 and675, with different next guards. More strongly,5^(2k) needs k applications of rule4. Thus the guard is a snapshot; exact unbounded counters and the remaining clock cannot be inferred from it.

## C3. Actual black-line compilation

Use the inherited seven-leaf cycle hub from `session160_ordered_programs_20261103/metagraph/NOTE.md`. Its W has forward pairs

    01,03,04,05,12,23,26,34,36,45,46,56

and all other pairs reversed. It contains C=(0,1,2,3,4,5,6,0). Let T_g reverse cycle edge(g,g+1), let P_g be C cut at that edge, and let Z reverse every non-cycle edge of W. The literal native operation K_P reverses all arcs except the path arcs. Then

    K_(P_g)(T_g)=Z,  K_(P_h)(Z)=T_h.

Hence every ordered leaf transition has a two-native-move implementation, including a self-transition by going out and back. Distinct converse-normalized score sequences certify eight different merged nodes; each chosen spoke is black. The actual H values are103 at Z and129,37,75,59,73,111,49 at the leaves.

For the arithmetic program keep (N,current leaf,target leaf,phase). Initially take leaf6; rule i targets leaf i-1. Each arithmetic step uses the paths P_current and P_target, updates N by its fraction, and costs exactly two native moves. All seven physical vertex names are retained. The intermediate Z must retain the target and phase. The hub's Hamilton count is not V(N); this is a product controller with explicit external arithmetic state.

The same compiler implements all70 transitions of the minimal seven-state primitive selector in [BANK.md](BANK.md), one transition per source letter. The actual ordered bank word must still be retained if it will be used as a tournament source. Its refined thirteen-state parity selector needs more control information; it is not silently identified with the seven hub leaves.

`controllers.py` constructs all charts and all49 leaf routes, independently applies native operations as pair reversals, checks all70 selector transitions and replays supplied integer executions with their exact external state and clock. No arbitrary graph edge is assumed to be a native metagraph line.
