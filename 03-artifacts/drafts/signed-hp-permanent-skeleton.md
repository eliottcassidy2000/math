# The Signed HP Permanent and the Skeletal Structure

**Instance:** opus-2026-03-07-S35c10, S35c11
**Status:** NEW DISCOVERY — multiple theorems proved, major corrections at n=7

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

### THM-C: GS flip always crosses the S=0 boundary at n=5 (PROVED)

For every GS tiling at n=5, exactly one of {T, flip(T)} has S(T) = 0.

Verified: 16/16 GS tilings cross the boundary.

This is EQUIVALENT to THM-060's bipartiteness: since S=0 iff t3 is odd (at n=5), and THM-060 proves t3 parity flips under GS flip at odd n, the S=0 boundary crossing follows.

### THM-D: S(T) divisibility (VERIFIED)

S(T) is always divisible by 2^4 = 16 at n=5 and n=7.

At n=5: S/16 in {-1, 0, 1, 3} (always integer).
At n=7: S/16 values ≡ 3 mod 4 always. S/64 has fractional part exactly 3/4.

### THM-E: S(T) NEVER ZERO at n=7 (VERIFIED — MAJOR CORRECTION)

**S(T) ≡ 48 mod 64 for ALL tournaments at n=7.** (20000/20000 checked, 0 zeros found.)

This is a STARK contrast to n=5 where S=0 for exactly half (6/12) of the isomorphism classes.

Congruence pattern across odd n:
- n=3: S ∈ {-2, 6}. S ≡ 2 mod 4. c_0 ∈ {-1/2, 3/2}. S NEVER ZERO.
- n=5: S ∈ {-16, 0, 16, 48}. S ≡ 0 mod 16. c_0 ∈ {-1, 0, 1, 3}. S=0 OCCURS.
- n=7: S ≡ 48 mod 64 ALWAYS. c_0 has fractional part 3/4. S NEVER ZERO.

So: S=0 is the exception (only at n=5), not the rule!

### THM-F: class_size = H / |Aut(T)| (PROVED by orbit-stabilizer)

In the backbone tiling encoding, each isomorphism class has size = H(T) / |Aut(T)|.

**Proof.** Each directed HP of T gives a unique labeling with the backbone path. Labelings related by Aut(T) give the same tiling class. By orbit-stabilizer, #distinct backbone-compatible labelings = H / |Aut(T)|.

Verified at n=3,4,5,6 for all classes.

**Corollary:** |Aut(T)| always divides H(T). Since tournament automorphism groups have odd order and H is odd, H/|Aut(T)| is always a positive odd integer.

### THM-G: perm(B) = 0 for skew-symmetric B at odd n (VERIFIED)

The standard permanent perm(B) = sum_P prod B[i][P(i)] vanishes at odd n=5 for ALL tournaments. (At even n, perm(B) ≠ 0 in general.)

This means S(T) (the "path permanent") captures information that the standard permanent cannot.

## Connection to Skeleton Structure — REVISED

The earlier claim that the skeleton has "Side A (S=0)" and "Side B (S≠0)" is **specific to n=5**.

At n=7, ALL classes have S≠0 (S ≡ 48 mod 64). The bipartition (THM-060) is by t3 parity, not S=0/S≠0. The perpendicularity at n=7 manifests differently than at n=5.

### S(T) under GS flip at n=7

S + S_flip is NOT constant at n=7 (unlike the clean S=0/S≠0 partition at n=5). However:
- S/16 ≡ 3 mod 4 for BOTH T and flip(T) (universal congruence)
- (S + S_flip)/32 is always ODD

The self-flip question at n=7 remains open (no self-flip classes found in sampling of 443/456 classes).

## Position Parity Connection

The diagonal M[a,a] = sum_P (-1)^{pos(a,P)} (THM-053). The M-diagonal vector determines (with vertex labeling) the class up to score-equivalent isomorphism.

At n=5, classes on Side A (S=0) have M-diagonal vectors with sum = H = 3*t3 containing negative entries. Classes on Side B (S≠0) can have all-positive M-diagonals. The regular tournament (class 10, H=15) is on Side A (S=0).

## Cross-Scale Pattern

| n | #classes | S=0 classes | S congruence | c_0 type |
|---|----------|------------|--------------|----------|
| 3 | 2        | 0          | S ≡ 2 mod 4 | half-integer |
| 4 | 4        | all (THM-A)| S = 0        | zero |
| 5 | 12       | 6          | S ≡ 0 mod 16 | integer |
| 6 | 56       | all (THM-A)| S = 0        | zero |
| 7 | 456      | 0          | S ≡ 48 mod 64| quarter-integer (3/4) |

The pattern of S=0 occurrence is: NEVER at n=3, ALWAYS at n=4, SOMETIMES at n=5, ALWAYS at n=6, NEVER at n=7. The "sometimes" at n=5 is the unique case that creates the S=0/S≠0 bipartition.

## The W-Polynomial Hierarchy

At general odd n, the W-polynomial has floor((n-1)/2)+1 coefficients:
- c_{n-1} = n! (universal, top coefficient)
- c_{n-3} = (n-2)! * (2*t3 - C(n,3)/2) (depends on t3, THM-054)
- c_{n-5} = depends on t3 and t5 (THM-055)
- ...
- c_0 = S(T)/2^{n-1} (the "residual" invariant)

At n=7: c_0 is NOT determined by {t3, t5, t7, bc33} alone (regression max error ~6). S(T) encodes genuinely deeper structural information.

At even n: c_0 = 0 universally, so one fewer degree of freedom.

## S(T) NOT related to Pfaffian or standard permanent

For skew-symmetric B:
- det(B) = 0 at odd n (universal)
- perm(B) = 0 at odd n (verified)
- S(T) ≠ 0 at odd n in general

S(T) is a "path permanent" (product along Hamiltonian paths) rather than a "cycle permanent" (standard permanent = product over cycle covers). These are fundamentally different objects. The "path permanent" S(T) captures information invisible to both det and perm.

## Open Questions

1. **Prove S(T) ≡ 48 mod 64 at n=7.** Why this specific congruence? Predict the congruence at n=9.

2. **Prove S(T) divisibility by 16.** Why is v_2(S) ≥ 4 at n ≥ 5? (At n=3, v_2(S) = 1.)

3. **What determines S(T) at n=7?** It's not a linear function of {t3, t5, t7, bc33}. What additional invariants are needed?

4. **n=9 prediction.** Will S(T) = 0 occur at n=9? The pattern (never/sometimes/never) suggests it depends on n mod 4 or n mod 8.

5. **Cross-scale embedding.** How does the n=5 skeleton embed inside n=7? The S=0 structure is completely different.

6. **Self-flip classes at n=7.** Do they exist? None found in 443/456 classes sampled.
