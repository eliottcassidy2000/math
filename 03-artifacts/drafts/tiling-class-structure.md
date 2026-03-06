# Tiling Class Structure Analysis

**Instance:** kind-pasteur-2026-03-05-S10
**Status:** Computationally confirmed at n=3,...,6

---

## Summary of Confirmed Structural Facts

### F1: class_size = H(T) / |Aut(T)|

For every isomorphism class at every n tested (n=3,...,6), the number of tilings
in the class equals H(T)/|Aut(T)| where T is the class representative.

When |Aut(T)| = 1 (the generic case), class_size = H(T).

This is an orbit-counting identity: each tiling corresponds to a Ham path of T
with a specific vertex labeling; automorphisms identify |Aut(T)| labelings that
give the same tiling.

### F2: Transitive class always has size 1

Class #0 (the all-zeros tiling) gives the transitive tournament T_n with
scores (n-1, n-2, ..., 0). It always has H(T) = 1 (unique Ham path),
|Aut| = 1, and is self-converse.

### F3: Full tiling class H(T) sequence (CORRECTED)

The all-ones tiling (every arc flipped from transitive) belongs to a class
whose size = H(T) (since |Aut| = 1 for n >= 4). The H(T) values are:

- n=3: H=3 (class size 1, |Aut|=3)
- n=4: H=5 = 2^2+1
- n=5: H=9 = 2^3+1
- n=6: H=17 = 2^4+1
- n=7: H=31 (NOT 2^5+1=33)
- n=8: H=57 (NOT 2^6+1=65)

The conjecture H(T) = 2^(n-2)+1 holds for n=4,5,6 but FAILS at n=7,8.

The actual sequence 3, 5, 9, 17, 31, 57, 105, 193, 355, 653 satisfies the
**Tribonacci recurrence**: H(n) = H(n-1) + H(n-2) + H(n-3) for n >= 6,
with initial values H(3)=3, H(4)=5, H(5)=9.

This is OEIS A000213 (Tribonacci numbers starting 1,1,1), specifically
H(n) = A000213(n). The growth rate converges to the Tribonacci constant
~1.83929, the real root of x^3 - x^2 - x - 1 = 0.

Verified computationally through n=12.

### F4: External blue pair connects transitive to full

At every n, the flip (bit inversion) of the transitive tiling gives the full
tiling. Both are grid-symmetric (sigma-fixed), so the connecting line is blue.
This is the unique blue line involving class #0.

### F5: Self-paired class hierarchy

A class is "self-paired" if flipping some member lands back in the same class.
- n=4: 1 BLUE self (#2, the full class itself)
- n=5: 2 BLACK self (#8 size=11, #10 size=13)
- n=6: 2 BLUE self (#46 size=41, #54 size=15) + 6 BLACK self (3 pairs)

Key: the self-paired classes at n+1 are CHILDREN of self-paired classes at n.
- n=4 blueself #2 -> n=5 blackself #8, #10 (both have parent #2)
- n=5 blackself #8 -> n=6 blueself #54 (parent #8)
- n=5 blackself #10 -> n=6 blackself #51 (parent #10)

### F6: Black self classes come in transpose pairs

At n=6, all 6 black-self classes form 3 transpose (RED) pairs:
{#14,#18}, {#29,#40}, {#51,#52}. Each pair has equal sizes.
This follows from: if T has a self-flip, T^op does too.

### F7: Vertex-orbit child count conjecture

The sum of vertex orbits over all classes at n predicts the number of classes
at n+1:
- n=3: sum=4, classes at n=4: 4 (match)
- n=4: sum=12, classes at n=5: 12 (match)
- n=5: sum=48, classes at n=6: 56 (FAILS)

The conjecture breaks at n=5->6 because some vertex removals from different
classes at n=5 give isomorphic tournaments at n=4, causing "collisions" in the
parent map that are not accounted for by the simple orbit count.

### F8: "Perpendicular" geometry

The user's observation: blueself/blackself pairs are "perpendicular" to the
external blue pair (transitive <-> full).

Interpretation: The external blue pair represents the "radial" axis in tiling
space (from no flips to all flips). Self-paired classes represent a transverse
structure: the flip stays within the same isomorphism class, meaning the
tournament has a non-trivial relationship with its "anti-tournament" (arc
complement relative to the fixed Ham path).

### F9: H(T) = H(T^op)

Confirmed by transpose pairs always having equal H values. This follows from
the path-reversal bijection Ham(T) <-> Ham(T^op).

### F10: Score sequences don't uniquely determine classes

Multiple classes can share the same score sequence. At n=5: classes #4 and #6
both have scores (3,3,2,1,1); classes #8,#9,#10 all share (3,2,2,2,1).

---

## Connection to Tex Paper Concepts

### self-converse = grid-symmetric under transpose

A tournament T is self-converse (T ~ T^op) iff the class is its own transpose
target (not in a RED pair). The grid transpose map (x,y) -> (n-y+1, n-x+1)
corresponds to T -> T^op.

### BlackSelf(n) from Open Problem 7

The tex defines "black-self" for a tournament with T ~ T^op, |Aut|>1, and
specific parity conditions on Fix(beta). The computation finds:
- n=5: 2 black-self classes (#8, #10), but both have |Aut|=1
- n=6: 6 black-self classes, all with |Aut|=1

The tex's definition requires |Aut|>1, so these are DIFFERENT from the tex's
BlackSelf concept. The HTML's "black" refers to flip-pair symmetry (not both
grid-symmetric), while the tex's "black-self" refers to the beta-involution
parity structure.

### Grid symmetry = sigma-fixed tilings

The count of grid-symmetric tilings = |Fix(sigma)| = 2^{floor((n-1)^2/4)},
which is distributed across classes. Grid-symmetric tilings correspond to
self-converse tournaments with a specific anti-automorphism structure.

---

## Open Questions

1. ~~Why does the full class have size 2^(n-2)+1?~~ RESOLVED: H(T_full) follows
   the Tribonacci recurrence H(n) = H(n-1)+H(n-2)+H(n-3) (OEIS A000213).
   The 2^(n-2)+1 pattern was a coincidence for n=4,5,6.
   Ham paths of T_full biject to "run decompositions" of {0,...,n-1}: ordered
   sequences of consecutive intervals where max(I_k) >= min(I_{k+1})+2.
   Proof sketch: interval containing max is always first, min always last;
   conditioning on first interval size gives Tribonacci via auxiliary sequence g(n).
2. What determines whether a self-paired class is blue vs black?
   PARTIAL: Blue = both tilings in the self-flip pair are grid-symmetric.
   Black = at least one is not. Deeper structural fact:
   - BLUE self => class is SELF-CONVERSE (T ~ T^op). Always true.
   - BLACK self at n=5: both classes are self-converse (not in transpose pairs).
   - BLACK self at n=6: all 6 form 3 TRANSPOSE PAIRS, so NOT self-converse.
   The transition from n=5 to n=6 is: at n=5, black-self classes are individually
   self-converse; at n=6, they lose self-converse status and pair up.
3. Why does the vertex-orbit conjecture fail at n=5->6?
4. Is there a formula for the number of self-paired classes at each n?
   Data: n=4: 1, n=5: 2, n=6: 8. Sequence: 1, 2, 8, ...
5. Do the self-paired classes have special OCF/Claim A properties?
   All self-paired classes are SELF-CONVERSE (T ~ T^op). OCF verified for all.
   All self-paired H values are ODD (5, 11, 13, 15, 25, 41, 43, 45).

---

## F11: Perpendicular Geometry (confirmed)

Every self-flip pair {t, flip(t)} has its Hamming-distance midpoint exactly at
m/2 (the center of the hypercube {0,1}^m). This is necessarily true since
flip exchanges 0s and 1s, so bit-count(t) + bit-count(flip(t)) = m.

The "perpendicular" structure: the transitive<->full axis is the main diagonal
of the hypercube (from (0,...,0) to (1,...,1)). Self-paired classes sit in the
hyperplane perpendicular to this diagonal at the center. Their members are
NOT concentrated at the center — they span a range of bit counts — but their
self-flip pairs always straddle the center symmetrically.

## F12: Tribonacci Recurrence for Full Tiling (PROVED computationally)

H(T_full_n) = H(T_full_{n-1}) + H(T_full_{n-2}) + H(T_full_{n-3}) for n >= 6,
with H(3)=3, H(4)=5, H(5)=9. This is OEIS A000213.

The full tiling tournament has the adjacency rule: v_i beats v_j iff j=i+1
(path edge) or i > j+1 (backward arc). Its Ham paths biject to "run
decompositions": ordered sequences of consecutive intervals covering
{0,...,n-1} where max(I_k) >= min(I_{k+1}) + 2.

Structural facts used in the proof:
- The interval containing max(S) must be FIRST (or the jump-down condition
  from any previous interval to it fails, since max(S) > all other elements)
- The interval containing min(S) must be LAST (same argument, reversed)
- f(n) = g(n) + g(n-1) where g(n) counts decomps with first interval size >= 2
- g(n) = g(n-1) + g(n-2) + g(n-3) (Tribonacci)
- Therefore f(n) = f(n-1) + f(n-2) + f(n-3)

Verified computationally through n=12.

---

## Full Proof: H(T_full) satisfies the Tribonacci recurrence

**Theorem.** Let T_n denote the "full tiling" tournament on n vertices (the
tournament obtained by reversing every arc of the transitive tournament).
Then H(T_n), the number of Hamiltonian paths in T_n, satisfies

    H(n) = H(n-1) + H(n-2) + H(n-3)    for all n >= 3,

with initial values H(0) = H(1) = H(2) = 1 (equivalently, H(3) = 3, H(4) = 5,
H(5) = 9). This is OEIS A000213.

### Step 1: Bijection with run decompositions

A **run decomposition** of the set S = {0, 1, ..., n-1} is an ordered sequence
of intervals (I_1, I_2, ..., I_m) such that:

  (R1) Each I_j = [a_j, b_j] is a non-empty interval of consecutive integers.
  (R2) The intervals are pairwise disjoint and their union is S.
  (R3) I_1 contains max(S) = n-1 (the maximum element is always first).
  (R4) I_m contains min(S) = 0 (the minimum element is always last).
  (R5) **Gap condition:** For each consecutive pair, max(I_j) >= min(I_{j+1}) + 2.

The Hamiltonian paths of T_n biject to run decompositions of {0, ..., n-1}.
(Each Hamiltonian path decomposes the vertex ordering into maximal ascending
runs; the arc structure of T_n forces max(S) to come first and min(S) last,
and consecutive runs must satisfy the gap condition because the arc from the
end of one run to the start of the next must go "backward" by at least 2.)

Let f(n) denote the number of run decompositions of {0, ..., n-1},
with f(0) = 1 (the empty decomposition of the empty set).

### Step 2: Conditioning on the first interval

Since I_1 must contain n-1 and consist of consecutive integers, it has the
form [n-k, n-1] for some k in {1, 2, ..., n}. After removing I_1, the
remaining set is {0, ..., n-k-1}, and the gap condition requires

    max(I_1) = n-1 >= min(I_2) + 2,   i.e.,   min(I_2) <= n-3.

**Case k >= 2:** The next interval I_2 must contain max({0,...,n-k-1}) = n-k-1,
so min(I_2) <= n-k-1. Since k >= 2, we have n-k-1 <= n-3, so the gap
condition is automatically satisfied. The remaining decomposition of
{0, ..., n-k-1} is therefore **unconstrained**, contributing f(n-k).

Summing over k = 2, 3, ..., n:

    [k >= 2 contribution] = f(n-2) + f(n-3) + ... + f(0) = S(n-2),

where S(m) := f(0) + f(1) + ... + f(m).

**Case k = 1:** The first interval is the singleton {n-1}. The remaining set
is {0, ..., n-2}, and the gap condition forces min(I_2) <= n-3. Since I_2
must contain n-2 (the maximum of the remaining set), I_2 has the form
[a, n-2] with a <= n-3, so I_2 has length j := n-1-a >= 2.

After removing I_2, the remaining set is {0, ..., n-2-j}. The gap condition
from I_2 requires min(I_3) <= max(I_2) - 2 = n-4. Since the next interval
contains max({0,...,n-2-j}) = n-2-j and j >= 2, we have n-2-j <= n-4,
so the gap is again automatically satisfied. The remaining decomposition of
{0, ..., n-2-j} is unconstrained, contributing f(n-1-j).

Summing over j = 2, 3, ..., n-1 (where j = n-1 means I_2 = {0,...,n-2}):

    [k = 1 contribution] = f(n-3) + f(n-4) + ... + f(0) = S(n-3).

### Step 3: The intermediate formula

Combining both cases:

    f(n) = S(n-2) + S(n-3)
         = [f(n-2) + S(n-3)] + S(n-3)
         = f(n-2) + 2 S(n-3).                              ... (*)

This holds for all n >= 3.

### Step 4: Telescoping to get the Tribonacci recurrence

Applying (*) at n and n-1:

    f(n)   = f(n-2) + 2 S(n-3)
    f(n-1) = f(n-3) + 2 S(n-4)

Subtracting:

    f(n) - f(n-1) = [f(n-2) - f(n-3)] + 2 [S(n-3) - S(n-4)]
                  = [f(n-2) - f(n-3)] + 2 f(n-3)
                  = f(n-2) + f(n-3).

Therefore:

    f(n) = f(n-1) + f(n-2) + f(n-3).                       QED

### Step 5: Base cases and the sequence

Direct enumeration:
- f(0) = 1 (empty decomposition).
- f(1) = 1 (the single interval {0}).
- f(2) = 1 (only {0,1}; the singleton {1} then {0} fails the gap: 1 < 0+2).
- f(3) = f(2)+f(1)+f(0) = 1+1+1 = 3.

The resulting sequence is:

| n   | 0 | 1 | 2 | 3 | 4 | 5  | 6  | 7  | 8  | 9   | 10  | 11  | 12  |
|-----|---|---|---|---|---|----|----|----|----|-----|-----|-----|-----|
| f(n)| 1 | 1 | 1 | 3 | 5 | 9  | 17 | 31 | 57 | 105 | 193 | 355 | 653 |

This matches H(T_n) for all n tested computationally (n = 3 through 12).
The growth rate converges to the Tribonacci constant tau ~ 1.83929, the
real root of x^3 - x^2 - x - 1 = 0.

### Remark on the coincidence with 2^(n-2)+1

For n = 4, 5, 6, the Tribonacci values 5, 9, 17 happen to equal 2^(n-2)+1.
This is a numerical coincidence: the Tribonacci sequence grows as ~tau^n
(exponential with base ~1.839), while 2^(n-2)+1 grows as ~2^n. They diverge
starting at n = 7, where f(7) = 31 but 2^5+1 = 33.
