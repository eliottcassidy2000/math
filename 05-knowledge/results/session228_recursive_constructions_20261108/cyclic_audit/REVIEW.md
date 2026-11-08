# Independent audit: binary cyclic physical inverse

**PASS, with the U convention made explicit before freezing.** Audited producer
PROOF.md SHA-256: 5b84f5de851da40ca0de4659132649964442c864e5cd5450f9cd6eef533cfce0.
Repository inheritance checked at 1a178a44f13f9602a3be9e84714be91cbfc7d7a6.

## All-order conclusions

1. The full colour-polynomial law is valid. In the distinguished bin containing
u, the only forbidden boundary patterns are an A-part immediately before u
without the B-part, and a B-part after u without the A-part. Inclusion-exclusion
has precisely the intersection counted in Section 1. Other bins necessarily
have A followed by B. The bijection proves the integer evaluations; polynomiality
proves the polynomial identity.
2. The six possible maximal block-run words give D=3ab+as+br+rs.
Only ABu ends at u, and deletion is the order join, hence E=d=ab.
The U counting argument gives (a+r)(2b+s)+ab, equal to D+br.
3. Crucially U permits its first path to be empty. Ordinary nonempty ordered
two-path counts F(2) do not. In C3 the latter root-start statistic is 2 whereas
U=3. This was flagged during audit; the frozen proof now states the convention
explicitly, consistent with the parent's Session225–227 convention. There is
no correction to the formula itself.
4. The inverse is an exact iff for existence of a representative with the
prescribed numerical state. Indeed t=br and z=(a+r)s force all proposed divisor
candidates. Conversely d,t,z reconstruct D=3d+t+z and U=D+t. Positivity excludes
zero denominators. gcd(d,0)=d handles the singleton boundary. This is not a
recognition theorem for a supplied graph.
5. Physical pair membership really is finite: F_T(2)>=2^n-2 follows by assigning
Hamilton paths to every ordered nontrivial bipartition. Thus an individual pair
(H,r) has order at most floor(log2(r+2)). The sharper target cap used by the
producer follows from D>=(2^m-1)(2^n-1)+2>=2^(m+n-1)+1. The extra inequality
has difference 2(2^(m-1)-1)(2^(n-1)-1), including all equality boundaries.
No size cap for arbitrary tournaments of bounded H is implied.
6. The scalar 13 obstruction is valid for this construction only. Its bound
forces m+n<=4, whose complete possibilities give 3,5,9,11,15. The independent
bank actually contains (H,F2)=(13,78), so unrestricted 13 realization is not in
doubt. No all-odd realization theorem follows.

## Attribution and scope

THM-3121-path-cover-walk-content-substitution-kernel.md already supplies the
complete substitution mechanism. THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md
already supplies complete cyclic profile transfer and the endpoint-jet sidecar.
The producer correctly calls its formulas explicit specializations, rather
than a new general calculus. The finite physical numerical inverse is the
specific useful outcome. The Session198/THM4882 source-grammar scope should not
be silently imposed here: these factors are arbitrary nonempty tournaments,
and their root-start count is ab, not necessarily one. Current proof does not
make that mistake. No literature-priority assertion is warranted or made.

## Independent computation

check.py uses literal permutations, counts backward adjacency positions, and
counts permissible cuts. It imports no producer algorithms. It checks:

- All 1,099 labelled factor tournaments through order five and their F2 floor.
- All 505 ordered labelled factor pairs with each order at most four and summed
  order at most six, by literal composite permutations (at most seven vertices).
- Full rooted state, root-start count ab, and order cap for all 505 composites.
- Full Q identity evaluated at 1 through the composite order, using
  Q(k)=sum_d binom(k,d)F(d); the identity at 0 is automatic. Both sides have
  degree at most composite order, so these checks cover the complete polynomial
  of every finite case, not merely the first three values.
- Inverse soundness and completeness for all 225 combinations of the 15 physical
  factor pairs through order five; singleton and formal false-positive controls.
- Complete scalar image through 64: 3,5,9,11,15,17,23,27,33,37,47,51. The target cap
  proves factors through order five suffice. It agrees with the producer's
  independently computed image. Producer's through-128 census is not independently
  recomputed in this audit; its code and finite scope were inspected.

The hostile formal state (7,1,7,1) passes sign and divisor tests but would require
(H(B),F_B(2))=(1,4), impossible: the floor forces order at most two, with F2 0 or 2.
This is an all-order illustration of why integrality alone is insufficient.

Normal and optimized independent transcripts are byte-identical. No Python
assert statements are used for verification. Reproduce from this directory:

    python3 check.py > normal.out
    python3 -O check.py > optimized.out
    cmp normal.out optimized.out

The script reads the frozen producer_CHECK.json next to its own source only for
its final image comparison. All mathematical counts and witness computations
preceding that comparison are independent. The packet is portable. Input and
output byte hashes are recorded in SHA256SUMS.

## Final physical-hostile addition

Independently verified the tournament with outgoing masks (20,1,2,7,14), rooted
at vertex 1. Literal permutations and cuts give state (9,3,9,3); directed
reachability from every vertex verifies strong connectivity. Here t=z=0, so
any cyclic factorization would force r=s=0, hence both factors singleton and
d=1, contradicting d=3. Thus the target is itself a physically realizable strong
rooted state but lies outside this construction. Arithmetic factor candidates
(3,0),(1,0) or their swap fail physical membership. This strengthens, and does
not alter, the earlier inverse statement. The final producer script also adds
literal reconstruction of its 47 inverse states through 128; inspected without
claiming an independent recomputation of the full order-six bank.
