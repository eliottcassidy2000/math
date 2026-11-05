# Independent review: odd-prime radix constraints at fixed physical ports

**PASS for the clarified frozen statement.** The prime spectral split, factorial divisibility of the amplitude, port lower bound and offset parity all follow, including an arbitrary finite initial transient. No physical p=5 or p=7 affine clock is inferred.

## Frozen revision and clarification

Reviewed in `session197-hybrid-constructions/fixed-port-prime-radix/`:

- `PROOF.md`: **b9ebccedba223520982c6cc66abbe048848b632d4243cae8cd485fc1a797a25c**
- `graphs.py`: 3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c
- `verify.py`: 0d572f8b6bfba28f1f04f2badbb62a5791bf5c13a19a5ea0c91a980327a1e6dd
- `VALIDATION.json`: 23587c9d4da79bed1940762ec4733c7972daf33f19ea45f8cc93b5e0ad40ae38

During review, the author clarified the hypothesis to say explicitly that the **count sequence is not eventually constant**. This is necessary: a constant count h also satisfies h_next=p*h-(p-1)h, even though the affine map is nonconstant. For example, three universal-sink births per period keep a singleton seed's count equal to1 and would violate the direction-count conclusion if constant tails were admitted.

The earlier wording, “a nonconstant recurrence,” was ambiguous on this point. The proof already used the intended nonconstant-tail hypothesis. This PASS is bound to the new hash above, not the preclarification proof hash115a968a.... The old fixed-mixed foundation remains unchanged; its guide-level statement should likewise carry the explicit intended hypothesis.

The source verifier's normal and Python-optimized receipts are identical to the frozen receipt. The code and receipt hashes did not change with the clarification.

## Finite-interface and spectral review

The retained cover condition c-P<=1 is closed under newborn deletion. For internal deletion the component count rises by one, but at least one selected-port endpoint is restored. Endpoint and singleton deletions do not increase c-P either. The finite shape bound c<=2k+1 follows from P<=2k, and insertion/deletion are inverse with the correct multiplicities for repeated nonport component types.

For right births, the potential2k-P+(nonport starts) cannot decrease. The only equality case is appending at a nonport end, which leaves the shape unchanged. Its diagonal multiplicity d is at most k+1, because c-P<=1 implies the number of nonport ends is at most the number of port starts plus one. The diagonal is independent of the mask. Thus a right word of length r has bases d^r; the analogous left statement is valid.

Opposite-side births commute while fixing every original vertex and port. Their newborn-to-newborn arc is always left-to-right. The finite cover matrices therefore commute as well, and simultaneous triangularization puts a mixed period's eigenvalues among d^r e^l. A nonconstant eventual affine sequence forces p to occur as an eigenvalue. Since p is prime, at least one nonempty side has length one, and p<=k+1. The absent-side case is handled without a0^0 convention.

## Mixed factorial amplitude and the transient

The free action on the j-1 nonempty internal fragments is valid: every fragment has a nonempty, disjoint labelled vertex set, and each still connects through the same repeated root after any permutation. Hence (j-1)! divides a_j.

When both directions occur, reversing arcs if needed makes the left subword a single birth. Its n newborns form an actual transitive module after commuting the opposite-side births. Every internal fragment ends at a different selected old port, so the uniform bound j<=k+1 holds even though the right-side source grows.

For a fixed number q of left clones, the remaining right continuation starts from a fixed physical graph and has the same uniform finite-interface annihilator. Triangular Stirling inversion therefore expresses each a_j(n) as an exponential-polynomial sequence with bases d^r, after a finite zero-eigenvalue transient. Jordan blocks contribute only polynomial factors. This is a uniform recurrence argument, not an assumption of bounded graph order.

The positivity bound a_j(n)=O((p/j)^n) is valid on a tail, since j!S(n,j) is asymptotic to j^n for fixed j. For j>p it forces eventual zero by integrality. For j=p it gives a bounded sequence; among the one-sided bases only1 can remain, and boundedness removes any nonconstant polynomial factor. Thus a_p(n) is eventually a constant integer divisible by(p-1)!.

Expanding the Stirling factor gives bases i*d^r. Primality leaves exactly the two cases in the proof. The i=p,d=1 contribution is the eventual a_p constant. The i=1,d^r=p case requires r=1; the growth bound excludes this base from every a_j with j>=2. Its sole possible contribution is a_1(n), which is itself an ordinary fixed-root clone sequence with one fixed left vertex. Positivity of that clone expansion excludes run indices above p, so its p coefficient is also divisible by(p-1)!.

Consequently A=a or a+b is an integer with the required factorial divisor. These coefficients are formed using n from the original seed, not a reindexed tail. No unnoticed division by p^N is introduced when the transient is removed. Nonconstant positive count growth forces A>0. Odd Hamilton counts then make beta an odd integer, and c=-(p-1)beta gives the exact claimed residue. In particular zero offset is impossible.

## Independent controls

`check_prime_radix.cpp` is a separate implementation. Its optimized and debug outputs agree. It checks:

- all5,405 labelled rooted tournaments through order5 by direct repeated-root subset DP for factorial divisibility;
-41,055 physical opposite-side commutations on all labelled seeds through order4, all port subsets and all two masks, checking every arc after swapping only the newborn labels;
- six cyclic three-block examples with positive top run indices2 through7 and exact top coefficients(1,4,36,576,14400,518400), equal to((j-1)!)^2. These test nonzero factorial cases at prime indices5 and7 without claiming a pure affine clock;
- both separate branches of the mixed amplitude proof.

For the latter controls, L(mask0) followed by R(mask7) on the negative ternary seed has a=0,b=180. Reversing that seed and using one L(mask7) followed by two R(mask0) births has a=180,b=0. The first construction has a growing source prefix; the second has a growing sink suffix. Direct graph counts and the conditional left-module run coefficients verify nine finite stages. No example with both a and b nonzero is needed by the argument or claimed by these checks.

A separate constant-tail control verifies why the explicit hypothesis is necessary. The author's two strong one-sided clocks and their parity conventions also reproduce exactly under source replay.

Independent checker SHA256:18d2450d4da29ed4042a39db95e2f6d2625c0c36efe9788f7216f8f2c56c88b1.
Independent receipt SHA256:46f756c77aca74fe596b9a81b729257a1567661122f1909f10362525216ac6e9.

The all-prime, all-k result rests on the finite-interface proof and the positive root-run decomposition. The finite controls do not claim existence at new radices or coverage by adaptive, port-changing or graph-rebuilding languages.
