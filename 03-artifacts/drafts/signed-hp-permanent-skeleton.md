# The Signed HP Permanent and the Skeletal Structure

**Instance:** opus-2026-03-07-S35c10
**Status:** NEW DISCOVERY — multiple theorems proved

## The Signed HP Permanent S(T)

**Definition.** For a tournament T with signed adjacency B[i][j] = 2A[i][j] - 1 (so B = +1 if i->j, -1 if j->i), the *signed HP permanent* is:

S(T) = sum_{all permutations P} prod_{i=0}^{n-2} B[P_i][P_{i+1}]

This sums over ALL n! permutations (not just Hamiltonian paths), with each permutation weighted by the product of signed adjacency entries along its path.

## Main Results

### THM-A: S(T) = 0 at even n (PROVED)

**Proof.** B is skew-symmetric: B[i][j] = -B[j][i]. The reversal P^rev = (P_{n-1}, ..., P_0) gives:

prod B[P^rev] = prod_{i} B[P_{n-1-i}][P_{n-2-i}] = prod_{j} B[P_{j+1}][P_j] = prod (-B[P_j][P_{j+1}]) = (-1)^{n-1} prod B[P_j][P_{j+1}]

At even n: (-1)^{n-1} = -1, so each pair {P, P^rev} contributes 0 to S(T). Since P != P^rev for all permutations of n >= 2 distinct elements, S(T) = 0. QED.

**Consequence:** At even n, the constant term c_0 of the W-polynomial vanishes universally.

### THM-B: H = c_0 + 3*t3 at n=5 (PROVED)

The W-polynomial W(r) = tr(M(r)) has only even powers of r:
  W(r) = c_0 + c_2*r^2 + c_4*r^4

where:
- c_4 = n! = 120 (universal)
- c_2 = (n-2)! * (2*t3 - C(n,3)/2) = 6*(2*t3 - 5) (THM-054)
- c_0 = S(T)/2^{n-1} = S(T)/16

At r = 1/2: H = c_0 + c_2/4 + c_4/16 = S/16 + 3*(2*t3-5)/2 + 120/16 = S/16 + 3*t3.

Verified for all 12 tournament classes at n=5.

### THM-C: GS flip always crosses the S=0 boundary at odd n=5 (PROVED)

For every GS tiling at n=5, exactly one of {T, flip(T)} has S(T) = 0.

Verified: 16/16 GS tilings cross the boundary.

This is EQUIVALENT to THM-060's bipartiteness: since S=0 iff t3 is odd (at n=5), and THM-060 proves t3 parity flips under GS flip at odd n, the S=0 boundary crossing follows.

### THM-D: S(T) divisibility (VERIFIED)

S(T) is always divisible by 2^4 = 16 at n=5 and n=7.

At n=5: S/16 in {-1, 0, 1, 3} (always integer).
At n=7: S/16 takes values {-1, -5, -9, 0, 3, 15, ...} (integer, but S/64 is not always integer).

## Connection to Skeleton Structure

The skeleton of isomorphism classes under GS flip has two "sides" at odd n:

**Side A** (S=0, odd t3): H determined entirely by cycle counts. H = 3*t3 at n=5.
**Side B** (S≠0, even t3): H has extra contribution from S(T)/2^{n-1}.

Every GS flip crosses between sides. This is the "perpendicularity" of the skeleton.

The unpaired (self-flip) classes at n=5 are classes 6 (H=11) and 8 (H=13), both on Side B (S≠0). No class on Side A is self-flip at n=5.

## Position Parity Connection

The diagonal M[a,a] = sum_P (-1)^{pos(a,P)} (THM-053). The M-diagonal vector determines (with vertex labeling) the class up to score-equivalent isomorphism.

Classes on Side A (S=0) have M-diagonal vectors with sum = H = 3*t3 but M-diag sorted values in {[3,1,1,1,-3], [3,3,1,1,1]} — always containing negative entries.

Classes on Side B (S≠0) have M-diagonal vectors that can be all-positive (e.g., class 7: [3,3,3,3,3]) or mixed. The regular tournament (class 10, H=15) is on Side A (S=0).

## Cross-Scale Pattern

| n | #classes | skeleton edges | bipartite? | S=0 characterization |
|---|----------|---------------|------------|---------------------|
| 3 | 2        | 1             | YES        | S=0 iff t3=1 (cyclic) |
| 5 | 12       | 35            | YES        | S=0 iff t3 odd |
| 7 | 456      | ~10000        | YES        | S=0 NOT determined by t3 parity alone |

The simplification "S=0 iff t3 odd" holds at n=3,5 but breaks at n=7, where S(T) depends on deeper structural invariants (t5, t7, bc33, etc.).

## The W-Polynomial Hierarchy

At general odd n, the W-polynomial has floor((n-1)/2)+1 coefficients:
- c_{n-1} = n! (universal, top coefficient)
- c_{n-3} = (n-2)! * (2*t3 - C(n,3)/2) (depends on t3, THM-054)
- c_{n-5} = depends on t3 and t5 (THM-055)
- ...
- c_0 = S(T)/2^{n-1} (the "residual" invariant)

Each coefficient from top down peels off one more cycle-count dependence. The bottom coefficient c_0 is the MOST refined invariant — it distinguishes tournaments that agree on all cycle counts.

At even n: c_0 = 0 universally, so one fewer degree of freedom.

## Open Questions

1. **Prove S(T) divisibility.** Why is S(T) always divisible by 2^4 at n=5,7? Is the minimum v_2(S) equal to floor(n/2) + 1?

2. **Algebraic formula for S(T).** S(T) = sum_P prod B_e is like a "path permanent" of the skew-symmetric matrix B. Is there a closed form in terms of other tournament invariants?

3. **Cross-scale embedding.** How does the n=5 skeleton embed inside the n=7 skeleton via vertex deletion? Each n=7 class has sub-tournaments at n=5 — how do their S-values relate?

4. **Self-flip class characterization.** At n=5, self-flip classes have S≠0. Is this true at all odd n?
