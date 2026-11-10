# Two exact even-difference compilers and their prime boundary

**PROVED scoped constructions; all-even difference existence and infinitude are INHERITED, not new.** The new quantitative application is a small one-vertex extension with both counts at most5N, plus an endpoint-only universal word comparison. No arbitrary count or prime pair is selected.

## Inheritance and correction of novelty

`session199/physical-interpolation/PROOF.md`, §4, already proves that every even difference has infinitely many distinct strong single-arc pairs. `05-knowledge/results/session221_computation_to_construction_20261107/physical/PROOF.md` gives explicit same-order strong single-arc pairs on5+ceil(log2 N) vertices, with counts12N-9-f(N-1),14N-9-f(N-1). THM-5000-coprime-interface-composition-consecutive-odd-windows.md also supplies arbitrarily long strong odd-count intervals. Those are stronger than a bare gap-existence assertion, and are expressly credited here.

THM-5001-every-positive-endpoint-quotient-has-a-small-strong-realizer.md supplies every positive endpoint count and its exact binary selector. The word maps below are inherited from Session199 and Session237. The present comparison gives a low absolute-count bound for a different physical relation (one-vertex extension), and isolates exactly why a universal endpoint-based count-gap compiler can be arithmetically useless for prime pairs.

## 1. The actual O/Z/J operations

For rooted source state (D,E,U,d), O appends a sink and retains the root. Z appends a vertex beating only the old root and retains the root. J is the same graph addition but marks the new vertex. Their states are

    O: (D,0,2U,d),
    Z: (D+U-E,d,2U-d,d),
    J: (D+U-E,D-E,D+U,D).                               (1)

For example E_J=D-E because a Hamilton path ending at the new vertex is exactly an old path not ending at the old root. For U_J, either the new root is its entire distinguished path (D covers), or removing its initial edge to the old root leaves one of the U old distinguished covers. The other formulas follow from the singleton graft law in PROOF.md and the sink path split.

Every Hamilton path splits at the root into a distinguished suffix and the complementary prefix. This injection gives U>=D. Starting at the edge0->1 rooted at0, the root reaches every vertex throughout any O/Z word. After a final Z, every old nonroot vertex reaches the new vertex, which reaches the root, so the output is strong.

## 2. A one-vertex extension of prescribed size, below5N

**Theorem.** For every N>=1 there exist strong tournaments L_N and R_N such that R_N minus one designated vertex is isomorphic to L_N, and

    H(R_N)-H(L_N)=2N,
    2N+1 <= H(L_N) < H(R_N) <= 5N,
    |R_N| <= ceil(log2 N)+4.                              (2)

For N>=2 the orders are exactly ceil(log2 N)+3 and ceil(log2 N)+4. For N=1 they are3 and4.

**Construction and proof.** At the fixed-root edge the state is(1,0,2,1). Along O/Z, d=1 and the response U changes by

    O: U->2U,             Z: U->2U-1.

Every U>=2 has a unique finite word from2: invert an even U by U/2 and an odd U>2 by(U+1)/2. The integer strictly decreases until2. The word length is ceil(log2 U)-1; its final parity selects its final letter, so the inverse is unambiguous.

For N>=2 build the source with response N, then append O. Write its unchanged count as D=g(N). It now has U=2N,E=0,d=1, and1<=D<=N by the injection U>=D before O. At N=1 use the original edge, D=1. Compare the two literal suffixes Z and OZ. Formula(1) gives

    H(L_N)=D+2N,           H(R_N)=D+4N.                    (3)

This proves the value bounds and gap. Both final Z moves make strong outputs by the reachability argument. Deleting the O vertex added immediately before the final Z in R_N leaves exactly the graph L_N: all arcs on the old source and the final bridge agree. This proves the one-vertex assertion, not merely an equality of counts. The source's word length and the two final suffix lengths give the stated orders. ∎

**Sharp equality within this chosen family.** For N>=2, H(R_N)=5N iff N=2^k+1 for an integer k>=1; equality also holds at N=1. No global extremality is claimed.

To prove this, track delta=U-D and E in the unique inverse word before the final O. Initially(delta,E)=(1,0). O gives delta'=U+delta>=2,E'=0. Z gives delta'=delta+E-1,E'=1. The first Z from the initial edge makes delta0; successive Zs retain0. Once an O has occurred, its delta is at least2; the next Z can reduce it by only1 and later Zs cannot reduce it again. Thus D=U occurs exactly at nonempty all-Z words. After k Zs, U=2^k+1. Equation(3) proves the claim. At N=3 the pair is9,15; at N=5 it is15,25.

This does not produce infinitely many one-vertex extensions for each fixed gap, nor a +2 extension at an arbitrary prescribed count. At N=1 it gives precisely the allowed3-to5 extension. The forbidden5-to7 and the inherited no-+2-extension theorem remain intact.

## 3. A fixed five-letter endpoint compiler

For any physical rooted state, multiply the elementary integer maps(1). One obtains

    H(A.JOJOZ)=9D-3E+9U,
    H(A.JOOOZ)=9D-E+9U.                                  (4)

Thus these two universal fixed words differ by exactly2E. The equality is an identity of linear forms, not a fit to a graph census. The verifier evaluates the four basis vectors to check the coefficient rows and then checks the graph operations separately.

If the starting source is strong, both outputs are strong. Each J creates a root reaching all vertices: it first reaches the old root, which already reached all. Each O retains this reachability, and each final Z makes the graph strong. Both outputs have exactly five more vertices than the same source. They need not differ at one edge or be induced deletions of one another.

To select E=N, take the inherited endpoint selector. For any

    n >= max(4,ceil(log2(N+1))+2),

let the base be TT_n with its extreme arc reversed. Its endpoint counts at interior vertices1,...,n-2 are1,2,...,2^(n-3). Adjoin root x, with incoming neighbors the interior vertices selected by the binary digits of N, and with x beating all remaining vertices. Then E_x=N. The base is strong, and x has both incoming and outgoing neighbors, so the source is strong. The resulting pair has common order n+6; at minimum this is max(10,ceil(log2(N+1))+8).

Let n grow without changing the chosen binary digits. The endpoint count remains N, but D>=2^(n-2)+1 by the path-insertion injection. Since U>=D, the lower count in(4) is at least18D-3N and is unbounded. This recovers infinitely many distinct strong pairs for every even difference through a different endpoint-preserving construction; the global existence statement was already established in Session199.

**Decisive prime obstruction.** The lower count in(4) is divisible by3 and exceeds3 for every such strong source. Therefore this entire universal endpoint compiler supplies **no prime pair at all**. If3 divides N, both counts are divisible by3. This is stronger than saying that primality has not been proved: arithmetic excludes it identically. The offset alone carries no Polignac consequence.

## 4. Exact controls

From this directory:

    python3 -I -B verify_gaps.py > gaps_replay.json
    cmp gaps_replay.json GAPS.json
    python3 -I -B -O verify_gaps.py > gaps_optimized.json
    cmp gaps_optimized.json GAPS.json

The standard-library verifier imports only the committed sibling verify.py. It checks all one-letter state transfers on every labelled tournament of orders1 through4 with every root; actual graph pairs for N=1,...,64 in(2); the integer inverse and size bound for N=1,...,10000; and nine literal endpoint-pair examples including forbidden full counts N=7,21 used lawfully as endpoint coordinates. Those endpoint values do not assert that7 or21 is a full count.

Its finite word search also finds no pair with coefficient difference(0,2,0,0) when both O/Z/J words have length at most4. This is only a finite comparison in the specified alphabet, not a universal minimality theorem for all graph operations.

The constructor preserves a chosen even difference and explicit graph relation, while losing arbitrary absolute placement and prime status. The one-vertex and endpoint-only lanes have different advantages; neither supersedes the inherited single-arc construction or resolves Chan–Pak or Polignac.
