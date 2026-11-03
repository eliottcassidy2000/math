# Independent synthesis audit

Status: **PASS for the norm theorem and cross-domain scope; strengthen the proposed next target.**

Read NOTE.md and check.py. The quadratic identity and inert-prime argument are correct for arbitrary nonempty A. Its valuation is exactly twice min(v_q(F_A(1)),v_q(F_A(2))) for the stated primes. The modulus-8 consequence correctly uses inherited tournament Hamilton-path oddness. The code only evaluates positive norms, so its valuation loop terminates. The cross-domain discussion explicitly distinguishes complex roots, decimal residues, and tournament Hamilton counts; the norm analogy claims no predicate-preserving map between them.

The individual lattice, polynomial-root and external record claims are attributed to separate packets; this audit does not independently reprove those packets or verify external records.

Independent controls in ../spectrum/audit_norm7.py cover all 75 labelled tournaments of orders 1–4 using literal wrapper Hamilton DP, plus all 2401 pairs modulo49. Normal and optimized certificates agree.

## Additional theorem: all transitive-feedback mirrored cyclic wrappers exclude 7 and 21

Let T_A(r) have nonempty A, its dual A*, and a transitive feedback tournament Z_r, oriented A→A*→Z_r→A, with r≥1. Then H(T_A(r)) is never 7 or21.

For r=1 this follows from the norm obstruction already proved in NOTE.md.

For r=2, the positive skeleton expansion is

    H=5q1²+6q1q2+7q2²+4q2q3+2q3²,

where qk=k!p_A(k), with nonexistent coefficients zero. If |A|≥2, every tournament has a Hamilton path, so q1≥1. Splitting one such path at any arc gives a two-path cover, so q2≥2. Every summand is nonnegative, and H≥5+12+28=45.

Hamilton-path counts are nondecreasing under induced tournament extensions. To prove this without an unverified monotonicity assumption, take a Hamilton path v1,...,vn and adjoin vertex x. If x→v1, prepend x. Otherwise, if vn→x append it. In the remaining case v1→x and x→vn, so there is a first index i with x→vi; then v_(i−1)→x→vi, allowing insertion. Use a deterministic choice. The map from old paths to extended paths is injective because deleting x recovers the original order. Iteration proves monotonicity under any number of adjoined vertices.

For r≥2 the wrapper contains an induced T_A(2), so |A|≥2 implies H≥45 throughout.

If |A|=1, name the two old vertices a,l. Vertex a beats only l. If a is not last in a Hamilton path it must be immediately followed by l. Such paths are exactly

    Z_S (in increasing order), a,l, Z_complement (in increasing order),

for any subset S of the r feedback vertices, giving 2^r paths. If a is last, the unique path is l followed by all feedback vertices in increasing order, then a. Thus

    H(T_A(r))=2^r+1.

Neither 7−1=6 nor21−1=20 is a power of two. This completes all orders and all r≥1. ∎

The argument excludes this entire construction class, not all tournaments. It supplies no global proof of the reported missing spectrum values. Consequently NOTE.md §5 should not propose trying other seeds or merely increasing transitive feedback size within the mirrored cyclic family. The possible escape by changing the feedback interior is addressed and ruled out in the stronger section below.

`audit_wrappers.py` independently checks the singleton formula by literal Hamilton DP for r=1,...,12 and the r=2 bound and positive expression for every labelled tournament A of orders1–4. The general theorem is proved above, not inferred from these tests.

## Stronger closure: arbitrary feedback interiors also fail

The same exclusion holds when Z is **any nonempty tournament**, not necessarily transitive, with uniform cyclic orientation A→A*→Z→A.

If |A|≥2 and |Z|≥2, selecting any two feedback vertices gives the induced two-feedback wrapper, so the preceding ≥45 bound applies. If |Z|=1, the norm theorem applies.

It remains to consider |A|=1, writing the two old vertices as a,l and m=|Z|. If a is not last, it must be immediately followed by l. The two segments of Z before and after this adjacent pair form two labelled path slots, so these paths number F_Z(2). If a is last, l must be first, and the middle is a Hamilton path of Z. Therefore

    H = F_Z(2)+H(Z) = 3H(Z)+2p_Z(2).                      (4)

For every subset S of V(Z), each of Z[S] and Z[V(Z)\S] has a Hamilton path, counting an empty path as one. Their products count two-slot covers, so F_Z(2)≥2^m and H≥2^m+1. Hence m≥5 gives H≥33.

For m=1,2 the only tournaments are transitive and H=3,5. For m=3 the transitive and cyclic cases give9 and15. For m=4, if Z is transitive, F_Z(2)=16 and H(Z)=1, giving17. Otherwise Z contains a cyclic triangle; induced-extension monotonicity gives H(Z)≥3. Let c≥1 count its cyclic triangles. Every two-path cover has block sizes1+3 or2+2. The former contribute Σ_v H(Z−v)=4+2c, because a triple has one Hamilton path if transitive and three if cyclic. The latter contribute exactly3 (the three unordered pairs-of-pairs). Thus p_Z(2)=7+2c≥9 and (4) is at least27. This excludes7 and21 in every case. ∎

This proof uses no unproved cyclic-triangle formula for H on four vertices and no computational classification premise. Additional literal tests enumerate all75 labelled feedback tournaments through order4 and confirm (4); observed values by feedback order are {3}, {5}, {9,15}, {17,27,37}.

The revised next target must therefore change more than the feedback interior: break the mirror relation A* or the uniform cyclic intermodule orientation, or use a genuinely different architecture. The norm theorem remains informative because it gives infinitely many prime-valuation restrictions, although exclusion of the two specific targets now holds for the larger class by the combined argument.
