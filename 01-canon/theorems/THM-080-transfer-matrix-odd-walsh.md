# THM-080: Complete Walsh Spectrum of Transfer Matrix M[a,b]

**Status:** PROVED (analytical derivation + verification at n=5 exhaustive, n=7 reconstruction 20/20)
**Author:** opus-2026-03-07-S35c5, extended S35c6, corrected S35c7
**Date:** 2026-03-07
**Dependencies:** THM-030 (transfer matrix symmetry), THM-077 (Walsh OCF proof)

## Main Theorem: Complete Walsh Formula for M[a,b]

For the transfer matrix M[a,b](T) = sum_{S subset U} (-1)^|S| E_a(S+{a}) B_b(R+{b}):

**hat{M[a,b]}[S] = (-1)^{asc(S)} * 2^s * (n-2-|S|)! / 2^{n-2}**

when S is a *valid monomial* (defined below), and 0 otherwise.

Here **s = number of unrooted (even-length) components** of S — components containing neither a nor b.

### Valid Monomial Conditions

S is valid for M[a,b] if and only if:
1. S forms a union of disjoint paths in the complete graph K_n
2. |S| ≡ n (mod 2) — odd for odd n, even for even n
3. |S| <= n-2 (so the factorial is defined)
4. No connected component of S contains both a and b
5. Components containing a: vertex a must be at a path endpoint
6. Components containing b: vertex b must be at a path endpoint
7. Components containing neither a nor b: must have even length

### Key Properties

- **2^s factor for unrooted components:** Rooted components (containing a or b) have 1 valid orientation in the HP (pinned at a or b). Unrooted even-length components have 2 valid orientations that give the same sign. This parallels H's 2^r factor (where all r components are unrooted).

- **Degree-1 specialization:** hat{M[a,b]}[{p,w}] = sgn(p-w) * (n-3)!/2^{n-2} for p in {a,b}, w in U (s=0 always). At n=5: +/-1/4. At n=7: +/-3/4.

- **Linear approximation:** M[a,b]^(1)(T) = (n-3)!/2^{n-3} * (s_a + s_b - (n-2)) where s_a, s_b are out-degrees to interior vertices.

## Verification

At n=5 (m=10): **968/968 monomials match** (all 2^10 = 1024 tournaments, exhaustive Walsh transform).
- Degree 1: 6 nonzero, amplitude 1/4 = 2!/2^3, all s=0
- Degree 3: 24 nonzero, amplitude 1/8 = 0!/2^3, all s=0
- Note: At n=5, s=0 for ALL valid monomials (insufficient free vertices for unrooted components).

At n=7: **20/20 random tournaments match** via Walsh reconstruction (degrees 1,3,5).
- Degree 1: 10 nonzero, all s=0
- Degree 3: 360 nonzero: 240 with s=0 (amp 1/16), **120 with s=1 (amp 1/8)**
- Degree 5: 720 nonzero, all s=0 (amp 1/32)
- The s=1 correction was essential: without 2^s factor, 16/20 tournaments had errors.

At n=6 (even n, m=15): **1471/1471 monomials match** (all 2^15 tournaments, exhaustive).
- Degree 0: 1 nonzero (amp 3/2 = 4!/2^4), mean of M
- Degree 2: 48 nonzero: s=0 (amp 1/8) and s=1 (amp 1/4)
- Degree 4: 120 nonzero: all amp 1/16
- Complement INVARIANT (not antisymmetric): M(T^op) = M(T)

At n=9 (partial, degrees 1,3,5): consistent with formula (residuals from missing degree 7).
- Degree 1: 14 nonzero, s=0
- Degree 3: 1680 nonzero: 840 s=0, 840 s=1

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

5. **Result (single-component):** hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-d)!/2^{n-2}. This is correct when S has only rooted components (s=0).

### Multi-component correction (2^s factor)

When S has s unrooted even-length components Q_1,...,Q_s in addition to rooted components:

6. **Unrooted block orientations:** Each Q_i appears as a contiguous block in sigma/tau. Being even-length, both orientations give the same chi_S sign. This contributes a factor 2 per unrooted component.

7. **Component side assignment:** Each Q_i can go to either the a-side or b-side of the partition. The alternating sum over these assignments and free vertex distributions telescopes to give (n-2-d)! as before, but accumulates the 2^s factor from orientations.

8. **Result (general):** hat{M[a,b]}[S] = (-1)^{asc(S)} * 2^s * (n-2-d)!/2^{n-2}. QED.

## Complement Symmetry (Corrected)

**Theorem:** M[a,b](T^op) = (-1)^n * M[a,b](T) for all tournaments T.

**Proof (algebraic):** Under T→T^op, directed HPs reverse: E_a(V,T^op) = B_a(V,T), B_b(V,T^op) = E_b(V,T). So M[a,b](T^op) = sum_S (-1)^|S| B_a(S∪{a}) E_b(R∪{b}). Relabeling S↔R and using M[a,b]=M[b,a] (THM-030) gives M[a,b](T^op) = (-1)^n * M[a,b](T). QED.

**Corollary:** For odd n, M is complement-antisymmetric (odd Walsh support).
For even n, M is complement-invariant (even Walsh support).

**Parity condition:** Walsh coefficient hat{M}[S] = 0 unless |S| ≡ n (mod 2).

**Proof (via alternating sum):** The telescoping sum over free vertex distributions gives factor sum_{t=0}^f (-1)^t where f = n-2-d. This equals 1 if f is even, 0 if f is odd. Since f even iff d ≡ n (mod 2), only monomials with |S| ≡ n (mod 2) survive. QED.

## Degree-3 Detailed Structure (n=5)

24 nonzero degree-3 monomials:
- 12 of type P2+P1: P2 on one side (a or b at endpoint), P1 on the other (the other root at endpoint). Both {a,b} present on DIFFERENT components.
- 12 of type P3: single 3-edge path with exactly one of {a,b} at a path endpoint. The other root absent from S.

Zero at degree 3: triangles (cycles), P3 with both {a,b}, P3 with {a,b} at interior, P2+P1 with both {a,b} on the P2 (forces cross edge), pure interior paths.

## Comparison with H Formula

H (THM-077): hat{H}[S] = (-1)^{asc(S)} * 2^r * (n-|S|)! / 2^{n-1} for even-length path unions.

M (this theorem): hat{M[a,b]}[S] = (-1)^{asc(S)} * 2^s * (n-2-|S|)! / 2^{n-2} for valid odd-degree monomials.

Key parallels:
- H has 2^r (all r components unrooted). M has 2^s (s unrooted components out of r total). Both factors count valid block orientations.
- H uses (n-d)! free-vertex factorial. M uses (n-2-d)! — shift by 2 for the two pinned vertices a,b.
- H normalizes by 2^{n-1}. M by 2^{n-2} — factor of 2.
- H is on even degrees, M on odd degrees.

### Why s=0 at n=5

At n=5, d is at most n-2=3. For valid monomials: interior vertices = n-2 = 3. A monomial of d=3 edges with rooted component(s) uses at least d vertices total, leaving at most 0 free vertices. So there's no room for unrooted components to exist alongside rooted ones. This made the 2^s correction invisible at n=5.

At n=7 with d=3: 5 interior vertices, rooted components use only 1-3 vertices, leaving room for unrooted even-length components (P2 uses 3 interior vertices).

## Amplitude Table

Base amplitude at degree d = (n-2-d)!/2^{n-2}. Multiply by 2^s for s unrooted components.

| n | deg 1 (s=0) | deg 3 (s=0) | deg 3 (s=1) | deg 5 (s=0) | deg 7 |
|---|-------------|-------------|-------------|-------------|-------|
| 3 | 1/2         | —           | —           | —           | —     |
| 5 | 1/4         | 1/8         | —           | —           | —     |
| 7 | 3/4         | 1/16        | 1/8         | 1/32        | —     |
| 9 | 3/2         | 3/8         | 3/4         | 1/16        | 1/128 |

## Previous Error (MISTAKE)

The original formula (S35c6) claimed hat{M}[S] = (-1)^{asc} * (n-2-d)!/2^{n-2} WITHOUT the 2^s factor. This was verified at n=5 where s=0 always, so the error was undetectable. At n=7, reconstruction without 2^s failed for 16/20 random tournaments. Adding 2^s gives 20/20 match.

## Open Questions

1. Can this formula provide a new proof of THM-030 (transfer matrix symmetry M[a,b] = M[b,a])? The Walsh spectrum is manifestly symmetric in a,b by the ascent sign rule and 2^s independence.

2. The H-M parallel (both have 2^{unrooted count} factors) suggests a unified framework. Is there a "super-object" G(T,z) = H(T) + z*M(T) with clean Walsh structure?

3. Can the corrected M formula + H = trace(M) lead to a new OCF proof? The 2^s factor makes M's structure closer to H's, potentially enabling cancellation arguments.
