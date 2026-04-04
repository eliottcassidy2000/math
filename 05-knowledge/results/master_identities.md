# Master Identities from 2^m = Σ H/|Aut|
## opus-2026-04-03-S27

All verified at n=3,4,5,6. Let [T] range over isomorphism classes.

### The Five Master Equations

**[ME1] Σ H(T)/|Aut(T)| = 2^m**
  The tiling equation: tiling counts sum to hypercube size.
  m = C(n-1,2). Each term = number of base-path-fixed labelings in class [T].

**[ME2] Σ n!/|Aut(T)| = 2^C(n,2)**  
  Standard orbit-stabilizer: labeled tournament counts sum to total.
  Ratio ME2/ME1 = 2^C(n,2)/2^m = 2^(n-1) = ... wait:
  C(n,2) - C(n-1,2) = n-1. So 2^C(n,2) / 2^m = 2^(n-1).
  This means: each iso class has 2^(n-1) times more labeled tournaments than tilings.
  But 2^(n-1) ≠ n!/H in general. So this ratio 2^(n-1) is an AVERAGE.

**[ME3] Σ 1/|Aut(T)| = 2^C(n,2)/n!**
  From ME2 divided by n!. This is just Burnside restated.
  Equals A000568(n)? No: A000568 = number of classes = Σ 1.
  But Σ 1/|Aut| ≠ A000568 (different sum).

**[ME4] E[H | uniform over labeled tournaments] = n!/2^(n-1)**
  EXACT at all n tested! The average Hamiltonian path count over all
  labeled tournaments (uniform random) equals n!/2^(n-1).
  This is because each directed arc is present with probability 1/2,
  and each of n! permutations is a Hamiltonian path with probability 1/2^(n-1)
  (each of n-1 arcs must be present).

**[ME5] Σ gs_tilings = 2^((m + floor((n-1)/2))/2)**
  Grid-symmetric tilings from hypotenuse fixed-tile counting.
  Combined with ME1 restricted to SC classes: Σ_{SC} H/|Aut| = SC tiling sum.

### Derived Identities

**From ME1 + ME4:**
  Σ H²/|Aut| = (Σ H×size) = sum over all tilings of H(tiling's class)
  n=3: 4, n=4: 32, n=5: 632, n=6: 29696

**From ME1 + THM-281:**
  Σ_{SC} H/|Aut| is a sum of an EVEN number of ODD terms = EVEN
  Σ_{NS} H/|Aut| is a sum of EVEN terms = EVEN
  Both halves of the decomposition are even.

**From |Aut| distribution:**
  At n=6: 73% of classes have |Aut|=1 (trivial automorphism group).
  For these: H = tiling count (all Hamiltonian paths create distinct base-path labelings).
  Only |Aut| ∈ {1, 3, 5, 9} appear at n=6. All odd!

### |Aut| IS ALWAYS ODD for tournaments

At n=3..6, every tournament automorphism group has ODD order.
|Aut| values: {1, 3, 5, 9} only (all odd).

This follows from: tournament automorphism groups have odd order (Alspach 1970).
A permutation σ with even order has a 2-cycle, which would require both arcs
between two vertices — impossible in a tournament.

### Consequence: H = |Aut| × (odd number)

Since |Aut| is always odd and H/|Aut| = tiling count = integer:
  H = |Aut| × tiling_count
Both factors are odd (|Aut| odd by Alspach, tiling_count odd for SC by THM-281).
So for SC classes: H = (odd) × (odd) = odd. ✓ (H is always odd for all tournaments.)
For NS classes: tiling_count is not necessarily odd, but H = |Aut| × tiling_count
with |Aut| odd, so H has the same parity as tiling_count.
Since H is always odd (known), tiling_count is always odd... but we showed NS
paired sizes are even! Contradiction?

No: NS paired size = 2 × tiling_count_of_one_partner. Each partner has 
tiling_count = H/|Aut|. H is odd, |Aut| is odd, so H/|Aut| is a ratio of
two odd numbers. It CAN be non-integer if |Aut| doesn't divide H? 
No — we proved tiling_count IS an integer. So |Aut| | H, and since both
are odd, H/|Aut| is odd/odd. But odd/odd can be odd or can be any integer...

Actually: we showed NS paired sizes (= 2 × single_class_size) are even.
single_class_size = H/|Aut| = (odd)/(odd). This is odd iff the division is exact
and the quotient is odd. The data shows NS class sizes include both odd and even
values (e.g., at n=6: NS sizes 1,3,5,9,11,13,...,43 — all odd!).

So NS single-class sizes ARE also always odd? Let me check...

YES! At n=6 all 44 NS class sizes are: 1,1,3,5,5,9,9,9,11,11,13,15,19,23,23,25,29,31,33,37,37,43 (each ×2 for the pair)
These are ALL ODD!

So ALL class sizes are odd, not just SC ones!

This follows from: H is odd, |Aut| is odd, so H/|Aut| = (odd)/(odd dividing it).
An odd number divided by an odd divisor gives... well, 15/3 = 5 (odd), 15/5 = 3 (odd).
In general: if a is odd and b|a with b odd, then a/b is odd.
Proof: if a/b were even, then a = b × (even) = (odd)(even) = even. Contradiction since a is odd.

THEOREM: ALL tournament class tiling counts are odd.
(Not just SC ones! THM-281 is actually a special case of a STRONGER theorem.)
