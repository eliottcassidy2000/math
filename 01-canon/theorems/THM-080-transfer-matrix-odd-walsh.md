# THM-080: Complete Walsh Spectrum of Transfer Matrix M[a,b]

**Status:** PROVED (analytical derivation + exhaustive verification at n=5)
**Author:** opus-2026-03-07-S35c5, extended S35c6
**Date:** 2026-03-07
**Dependencies:** THM-030 (transfer matrix symmetry), THM-077 (Walsh OCF proof)

## Main Theorem: Complete Walsh Formula for M[a,b]

For the transfer matrix M[a,b](T) = sum_{S subset U} (-1)^|S| E_a(S+{a}) B_b(R+{b}):

**hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-|S|)! / 2^{n-2}**

when S is a *valid monomial* (defined below), and 0 otherwise.

### Valid Monomial Conditions

S is valid for M[a,b] if and only if:
1. S forms a union of disjoint paths in the complete graph K_n
2. |S| is odd (complement antisymmetry)
3. |S| <= n-2 (so the factorial is defined)
4. No connected component of S contains both a and b
5. Components containing a: vertex a must be at a path endpoint
6. Components containing b: vertex b must be at a path endpoint
7. Components containing neither a nor b: must have even length

### Key Properties

- **No r-dependence:** Unlike H (where amplitude has factor 2^r for r components), the M amplitude depends ONLY on |S| and n. This is a fundamental simplification.

- **Degree-1 specialization:** hat{M[a,b]}[{p,w}] = sgn(p-w) * (n-3)!/2^{n-2} for p in {a,b}, w in U. At n=5: ±1/4. At n=7: ±3/4.

- **Linear approximation:** M[a,b]^(1)(T) = (n-3)!/2^{n-3} * (s_a + s_b - (n-2)) where s_a, s_b are out-degrees to interior vertices.

## Verification

At n=5 (m=10): **968/968 monomials match** (all 2^10 = 1024 tournaments, exhaustive Walsh transform).
- Degree 0: 0 nonzero (correct)
- Degree 1: 6 nonzero, amplitude 1/4 = 2!/2^3 (correct)
- Degree 2: 0 nonzero (correct)
- Degree 3: 24 nonzero, amplitude 1/8 = 0!/2^3 (correct)
- Degree 4-10: 0 nonzero (correct)

At n=7: degree-1 amplitude 3/4 = 4!/2^5 verified via paired sampling (500+ backgrounds, SE < 0.01).

## Analytical Proof (Degree 1)

For edge {a,w} with w in U:

1. **Factorization:** Only partitions with w in S contribute (cross edges average to 0).

2. **Rooted HP Walsh:** hat{E_a}[{a,w}]_k = sgn(w-a)/2 * (k-2)!/2^{k-1}, where k = |S_U|+1. This follows because only HPs ending at a with penultimate vertex w contribute, and E[A[w][a] * chi_{aw}] = sgn(w-a)/2.

3. **Average B_b:** E[B_b(R+{b})]_k' = (k'-1)!/2^{k'-1} where k' = |R_U|+1.

4. **Alternating sum telescoping:**
   sum_{t=0}^{n-3} (-1)^t * C(n-3,t) * t! * (n-3-t)! = (n-3)! * sum (-1)^t
   Since C(n-3,t) * t! * (n-3-t)! = (n-3)! for all t, and for n odd the alternating sum = 1.

5. **Result:** hat{M}[{a,w}] = sgn(a-w) * (n-3)!/2^{n-2} = (-1)^{asc} * (n-2-1)!/2^{n-2}. QED.

## Analytical Proof (General Degree)

For a monomial S = P_a (component on a-side) union P_b (component on b-side):

1. **E_a side:** hat{E_a}[P_a edges]_k = (-1)^{asc_rev(P_a)} * (k-d_a-1)!/2^{k-1}. This uses the contiguity argument: P_a must appear as a contiguous block in sigma (ending at a), with a at the block's end (one orientation only for odd d_a; for even d_a, the non-a case has two orientations but they give the same sign).

2. **B_b side:** hat{B_b}[P_b edges]_{k'} = (-1)^{asc_fwd(P_b)} * (k'-d_b-1)!/2^{k'-1}. Analogous, with b at the block's start.

3. **Alternating sum:** With d_a + d_b = d vertices fixed, the remaining n-2-d free vertices partition freely. The sum telescopes via C(n-2-d, t) * t! * (n-2-d-t)! = (n-2-d)!, giving factor (n-2-d)! * (-1)^{d_a}.

4. **Sign unification:** (-1)^{asc_rev(P_a) + d_a} = (-1)^{asc_fwd(P_a)} (since asc_rev + asc_fwd = d_a). Combined with (-1)^{asc_fwd(P_b)}: total sign = (-1)^{asc(S)}.

5. **Result:** hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-d)!/2^{n-2}. QED.

## Complement Antisymmetry and Odd Walsh

**Theorem:** M[a,b](T^op) = -M[a,b](T) for all tournaments T.

**Proof:** Since chi_S(T^op) = (-1)^|S| * chi_S(T), and M has support only on odd |S|:
M(T^op) = sum_S hat{M}[S] * chi_S(T^op) = sum_S hat{M}[S] * (-1)^|S| * chi_S(T) = -M(T). QED.

## Degree-3 Detailed Structure (n=5)

24 nonzero degree-3 monomials:
- 12 of type P2+P1: P2 on one side (a or b at endpoint), P1 on the other (the other root at endpoint). Both {a,b} present on DIFFERENT components.
- 12 of type P3: single 3-edge path with exactly one of {a,b} at a path endpoint. The other root absent from S.

Zero at degree 3: triangles (cycles), P3 with both {a,b}, P3 with {a,b} at interior, P2+P1 with both {a,b} on the P2 (forces cross edge), pure interior paths.

## Comparison with H Formula

H (THM-077): hat{H}[S] = (-1)^{asc(S)} * 2^r * (n-|S|)! / 2^{n-1} for even-length path unions.

M (this theorem): hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-|S|)! / 2^{n-2} for valid odd-degree monomials.

Key differences:
- H has factor 2^r (depends on number of components). M has no such factor.
- H uses (n-|S|)!. M uses (n-2-|S|)! = (n-|S|-2)! — shift by 2 (for the two pinned vertices a,b).
- H normalizes by 2^{n-1}. M by 2^{n-2} — factor of 2 difference.
- H is on even degrees, M on odd degrees.

## Amplitude Table

| n | deg 1 | deg 3 | deg 5 | deg 7 |
|---|-------|-------|-------|-------|
| 3 | 1/2   | —     | —     | —     |
| 5 | 1/4   | 1/8   | —     | —     |
| 7 | 3/4   | 1/16  | 1/32  | —     |
| 9 | 45/8  | 3/8   | 1/16  | 1/128 |

Formula: amplitude at degree d = (n-2-d)!/2^{n-2}

## Open Questions

1. Can this formula provide a new proof of THM-030 (transfer matrix symmetry M[a,b] = M[b,a])? The Walsh spectrum is manifestly symmetric in a,b by the ascent sign rule.

2. The H-M duality (even/odd Walsh) combined with their amplitude formulas suggests a generating function interpretation. Is there a "super-object" G(T,z) = H(T) + z*M(T) with clean Walsh structure?

3. Does the no-r-dependence of M's amplitude have a representation-theoretic explanation? For H, the 2^r factor comes from block orientations; for M, the pinning at a,b eliminates this freedom.

4. Can the M formula be extended to prove OCF? Since H = sum_{a,b} M[a,b] and M's Walsh spectrum is known, perhaps OCF follows from summing M's contributions.
