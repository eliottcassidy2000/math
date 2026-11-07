# Paired-black leaves: exact rigid reduction and the remaining center alphabet

Status: PROVED structural reductions; FINITE-EXACT new order-nine exclusion.
The strong deficit-two center alphabet remains OPEN. In particular this note
does not classify all bridges, prove that every bridge is pendant, or assert
that S4 is the only strong deficit-two center at every order.

## 1. Inheritance, objects and the lost coordinate

The paired-black quotient Q divides every literal black-line multiplicity
by two, pairing lines by tiling reflection. Loops have degree two. For a
converse-merged tournament node, let M be its tiling-fibre mass and B its
number of reflection-fixed (blue) tilings. The proved degree identity is

    deg_Q = (M-B)/2.

Here SC means self-converse, not strongly connected; "strong" always means
strongly connected. An ordinary tournament has fibre mass
f(T)=H(T)/|Aut(T)|, an odd positive integer. Its Hamiltonian paths and its
automorphism group factor over its unique, linearly ordered strong
components. A converse-nonself-converse (NS) node has M=2f and B=0; an SC
node has M=f and B odd.

Dependencies and controls:

- The degree identity and constructive lower-bound families are in
  `../merged_metagraph_black_quotient_bridges_20261020/README.md`,
  sections 2–3; its complete census is through order eight.
- The all-order ordinary rigid classification is inherited from section 3
  of `../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md`:
  strong path-rigid atoms are exactly singleton and directed C3. It uses
  the cited tournament Hamiltonian-bypass theorem, not an extrapolation.
- The SCC center law is THM-2454,
  `01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md`, section 2.
  Its product formulas, not just its purity consequence, are required here.
- The corrected near miss is regular C5: f=B=3 although it is not rigid.
  Therefore neither nonrigidity nor nonpurity alone measures black deficit.
- The exact SC frame catalog
  `../odd_blue_spectrum_hollow_cubes_20261020/pairing_n9.json` tests a small-deficit
  question without constructing the full 2^28 fixed-path cube.

The map from a tournament to its ordered SCC word retains each component's
fibre and blue counts. This preserves the full
paired-black degree through the product formula, but discards the actual
neighbor and whether a nonleaf edge is a bridge. The necessary sidecar is
the literal complement incidence data; counts cannot restore its wiring.

## 2. All NS leaves, with no missing higher-order rigid atoms

**Theorem.** At every order n>=3, the NS Q-leaves are exactly the unordered
reversal pairs of nonpalindromic words in letters 1 and 3, interpreted as
ordered stacks of singleton and C3 components.

Proof. For an NS node, deg_Q=f. Thus it is a leaf exactly when f=1.
Since component fibre masses are positive integers and multiply, each
component is rigid. The all-order rigid classification makes every
component singleton or C3. Conversely every such word has f=1.
The SCC decomposition is canonical; each of its letters is self-converse,
and converse reverses the word. Thus the node is NS exactly when the word
is not a palindrome. The two orientations are exactly the reversal pair.
A degree-one vertex cannot carry a loop. QED.

Let a_0=1, a_j=0 for j<0, and a_n=a_(n-1)+a_(n-3).
Let p_n count palindromes in these letters. Then

    p_(2r)=a_r,              p_(2r+1)=a_r+a_(r-1),
    N_n=(a_n-p_n)/2.

This promotes Family I in the earlier bridge packet from a constructed
subfamily to the complete NS leaf classification.

## 3. All mixed leaves reduce to a single strong center

Write any SC tournament in its canonical SCC form as

    A -> Z -> A^op,

where A is an ordered list of strong components and Z is either a single
strong self-converse component or empty. Put

    t = product of the component fibre masses in A,
    h = f(Z),    b = B(Z).

For the empty center take h=b=1. The center law says

    M=t^2 h,       B=t b,       deg_Q=t(th-b)/2.                 (1)

The integer t is odd. If deg_Q=1, equation (1) gives t(th-b)=2,
so t=1 and h-b=2. Therefore A consists entirely of rigid atoms,
hence singleton/C3 letters. The center cannot be empty because an empty
center has h-b=0. Conversely every strong SC center with h-b=2, wrapped
in any rigid word A and its converse, gives a mixed Q-leaf. All these
presentations are unique by ordered SCC decomposition.

**Exact classification.** Every mixed leaf is uniquely

    A -> Z -> A^op,

where A is a singleton/C3 word and Z is a strong SC tournament with
exactly two nonblue Hamiltonian-path orbits, equivalently f(Z)-B(Z)=2.

The word "mixed" is automatic: b>=1 for a self-converse tournament, and
the deficit is positive. The reduction is an iff, but identifying every
possible Z remains an independent classification problem.

### A sharp low-degree extension

Since h>=b>=1,

    deg_Q = t((t-1)h+(h-b))/2 >= t(t-1)/2.                    (2)

Thus degree one OR degree two forces t=1. Outer nonrigidity can first
occur at degree three. Equality there with t>1 forces

    t=3,       h=b=1.

This is realized by C5 -> C5^op with empty center:
M=9, B=3, deg_Q=3. It is a decisive hostile boundary to extending
"small degree forces rigid outer components" beyond degree two.

## 4. Exact generating series in the unknown center alphabet

Let d_s count strong SC isomorphism classes of order s with f-B=2,
and let D(z)=sum_(s>=1) d_s z^s. This reduction does not require
finiteness, rationality or an eventual recurrence for D.

The complete leaf counting series is

    Leaf(z) =
      1/2 [ 1/(1-z-z^3) - (1+z+z^3)/(1-z^2-z^6) ]
      + D(z)/(1-z^2-z^6).                                   (3)

The first summand counts all NS leaves modulo reversal; the second counts
all mixed leaves by their unique half-word and central atom. Coefficients
below order three are zero, so the formal inclusion of empty words causes
no graph-convention issue.

The unique strong S4 has f=5 and B=3, so d_4=1. The full old census
through n=8 and the new SC9 test give d_s=0 for s<=9, s!=4.
Thus D(z)=z^4+O(z^10), a FINITE-EXACT initial segment, not an all-order
identity. The exact leaf counts for n=3,...,9 are

    0, 2, 1, 3, 3, 6, 7.

The positive S4 control can be labelled with out-neighbor bitsets
(2,4,9,3), that is arcs 0->1, 1->2, 2->3, 3->0, 2->0, 3->1.
Direct permutation enumeration gives H=5 and |Aut|=1. In the earlier
explorer coordinates its fibre is {1,3,5,6,7}, with blue subset {1,6,7}.

## 5. Cheap exact order-nine exclusion

The complete input consists of all 2,752 self-converse order-nine classes,
represented by normalized anti-involution frames. Exactly 2,661 are strong.
This input was independently enumerated from 65,536 frame masks and from
all 65,536 blue tilings in the odd-blue packet; it is not a sample of the
2^28 full fixed-path tilings.

For each class, subset DP computes the exact Hamiltonian-path count H.
The equality f-B=2 would require

    |Aut| = H/(B+2).

Filtering is safe and one-sided:

1. Required order must be an integer: 184 candidates.
2. It must be odd and divide 9!: 37 remain.
3. It must be divisible by the already known anti-involution centralizer
   order: all 37 remain.
4. Aut preserves the degree partition and every equitable refinement.
   Hence its order divides the product of the refined color-cell
   factorials: only five remain.

Only these five classes require automorphism backtracking. They have
(frame mask,H,B,required Aut,actual Aut) equal to

    (32231,15,3,3,5),
    (64526,3,1,1,3),
    (64559,9,1,3,9),
    (64605,27,1,9,27),
    (64607,9,1,3,9).

All five are nonstrong pure-blue classes, and none has the required order.
In particular every strong order-nine class is eliminated before any
automorphism search. The result excludes all SC leaves at n=9, not merely
strong ones.

Reproduction, from this published packet directory:

    python3 check_deficit_two.py --source ../odd_blue_spectrum_hollow_cubes_20261020 \
      --output verification.json > verification.out

The script records every H value, candidate order, input SHA-256 and
filter count in `verification.json`. It uses exact integers only, and
does not enumerate 9! labelings for each class or the full tiling cube.
Controls directly enumerate all vertex permutations through n=5,
cross-check the subset DP, and verify S4 versus pure-blue regular C5.
All reconstructed frame triangle counts are checked against the input.

The completed independent audit in `audit/AUDIT.md` reconstructs all
anti-involution frames, verifies that the supplied classes partition
the entire frame universe, independently recomputes every blue count,
and recomputes H with a backwards subset recurrence. It verifies the
five surviving automorphism groups by direct permutations within degree
cells, and checks every producer row at n=3,5,7,9. All checks pass.
The filters are necessary conditions only; failure of equality is not
replaced by a heuristic group estimate. A full order-nine adjacency
census is unnecessary.

## 6. Bridge honesty and research boundary

Every Q-leaf edge is a bridge, but two leaves can share one edge: at n=4,
the NS leaf and the S4 mixed leaf are its two ends. For the constructed
singleton/C3 and S4 families, the earlier triangle-current separation
proof shows distinct bridges for n>=5, preserving the all-order bound

    #bridges(Q_n) >= (a_n-p_n)/2 + [n even] a_((n-4)/2).

This note does not upgrade that bridge lower bound to equality. Even a
complete leaf classification would not exclude bridges joining two
nonleaf regions. Nor does the finite-exact SC9 leaf test compute all
order-nine bridges or prove active black connectivity.

The additional all-order cycle-cut obstruction is stated and audited in
`CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md`. It rules out all odd orders
n>=7 not divisible by three, and all orders divisible by eight.
Combined with the complete finite checks, any still-unknown center has
odd order divisible by three and at least 15, or even order at least 10
with v2(n)<=2. These are the congruence classes left by this arithmetic
obstruction; no finite cutoff is proved here. The obstruction does not
assert that purity alone classifies near-purity.
