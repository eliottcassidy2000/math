# Universal modular clone invariance detects primes and Carmichael numbers

**Status: PROVED elementary theorem below, using the inherited substitution law. Finite controls are separately scoped. Chan–Pak remains OPEN.**

## Inheritance and concept board

The closest proved mechanism is THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`, equation (4): tournament substitution requires full path-cover data, not H alone. THM-3134's endpoint jet is the neglected sidecar. The hostile is a cyclic triangle, whose transitive root clone gives 2^m+1, rather than a scalar multiple of its original value 3. The corrected near miss is a base-2 Fermat pseudoprime being mistaken for a universal substitution congruence. Session258's same-base endpoint budgets and Session259's insistence on a genuine operation map are inherited design constraints, not claims proved again here.

Live board: actual tournament substitution; positive Stirling skeletons; finite-field Frobenius; strong alternating-root tests; exact versus congruence realization. Anchor: Chan–Pak constructor interfaces. Niche: Carmichael detection by physical substitutions. Wildcard: prime-power failure as a p-run path-cover defect.

## 1. The actual operation and inherited expansion

Fix a rooted tournament (T,v). Replace v by a transitive tournament TT_m, giving all its vertices the old external orientation of v. Write h_v(m) for the Hamiltonian-path count, with m>=1. Split a Hamiltonian path into maximal runs inside the clone module. Collapse each such run to the symbol v. If there are r runs, the resulting skeleton visits every outside vertex exactly once and v exactly r times, with no adjacent v symbols. Let a_r(T,v) be the number of such directed skeletons. Then

    h_v(m) = sum_{r>=1} a_r(T,v) F_r(m),
    F_r(m) = r! S(m,r)
           = sum_{j=0}^r (-1)^(r-j) binom(r,j) j^m.             (1)

Here S(m,r) is the Stirling number of the second kind. The coefficients are nonnegative integers, a_1=H(T), and a_r=0 for r>floor((|T|+1)/2). Two successive occurrences of v require at least two outside vertices between them, since v->x->v is impossible. Formula (1) is a specialization of THM-3121, not a new substitution theorem. Its bijection assigns the m ordered clone vertices surjectively to the r successive run positions; within each run their order is forced.

An immediate positive survivor: if v lies in a nontrivial strong component, a_2>=1, hence

    h_v(m) >= H(T) + 2^m - 2.                               (2)

One proof that a_2>=1 avoids trying to extend a triangle skeleton: duplicate v transitively. The duplicated tournament contains an induced copy of T, so OCF counts all its old packings. A directed triangle through v gives a second triangle using its duplicate, hence the count strictly increases. Since F_1(2)=1, F_2(2)=2, and all other F_r(2)=0, (1) forces a_2>=1. Every vertex of a nontrivial strong tournament belongs to a triangle: otherwise no edge leads from its out-neighbourhood to its in-neighbourhood, contradicting strong connectivity. For a general tournament apply this inside its strong component. Conversely, if v is a singleton strong component, substituting TT_m only replaces that singleton by consecutive singleton components, so h_v(m)=H(T). Thus a_2=0 iff v is a singleton strong component. In particular no nonconstant arithmetic progression can arise from increasing a fixed transitive clone block.

## 2. Main theorem: an exact universal congruence classification

**Theorem.** For an integer m>=2, the following are equivalent.

1. For every finite rooted tournament (T,v),

       H(T[v <- TT_m]) == H(T) (mod m).                     (3)

2. The same congruence holds for every rooted **strong** tournament.
3. m is squarefree and p-1 divides m-1 for every prime divisor p of m.
4. m is prime or a Carmichael number, where a Carmichael number is a composite m for which j^m==j (mod m) for every integer j.

The familiar equivalence of (3) and (4) in this list is included in the proof; no unverified classification is imported. The theorem is a new use of the inherited skeleton law, not a claim of external priority.

### Sufficiency

If m is squarefree and each p-1 divides m-1, Fermat's theorem gives j^m==j modulo every prime divisor p, including j divisible by p. Chinese remaindering gives the congruence modulo m. Substituting in (1) gives

    F_1(m)==1 (mod m),     F_r(m)==0 (mod m), r>=2,

since an r-th finite difference annihilates the polynomial j when r>=2. Therefore (3) holds for every rooted tournament.

### Genuine strong test objects

For r>=2 construct T_r on v and 2r-2 outside vertices ordered transitively as

    a_1 < b_1 < a_2 < b_2 < ... < a_(r-1) < b_(r-1).

Let v beat every a_i and lose to every b_i. The tournament is strong: v reaches every outside vertex through a_1, and every outside vertex reaches v through b_(r-1). There are at most r occurrences of v in any skeleton. At exactly r occurrences, v must occupy both ends, with each gap consisting of a pair a_i,b_j satisfying i<=j. Covering every a_i and b_i forces this matching to be i=j: summing i<=j over the matching leaves equal totals. The matched pairs may occur in any order. Therefore

    a_r(T_r,v) = (r-1)!,   a_s(T_r,v)=0 for s>r.             (4)

This supplies a triangular family of physically realizable tests, not an arbitrary choice of formal coefficients.

### Necessity: recover Frobenius modulo p

Assume (3) holds for every rooted strong tournament. Fix p^e dividing m. Testing T_2,...,T_p and using (4), induction on r gives

    F_r(m)==0 (mod p^e),              2<=r<=p,              (5)

because the leading coefficient (r-1)! is a unit modulo p^e when r<=p. The a_1 terms cancel as F_1(m)=1.

Newton inversion says, for every nonnegative integer j,

    j^m = sum_{r=0}^j binom(j,r) F_r(m).

Taking 1<=j<=p-1 in (5) gives j^m==j (mod p). A generator of the cyclic group of nonzero residues then forces p-1 | m-1. For p=2 this divisibility is automatic.

### Necessity: the p-run defect excludes square factors

Suppose e>=2. Stirling recurrence gives

    S(m,p)==S(m-1,p-1) (mod p).

As p-1 divides m-1, the finite-difference formula and binom(p-1,j)==(-1)^j (mod p) give

    (p-1)! S(m-1,p-1)
      == sum_{j=1}^{p-1} j^(m-1)
      == -1 (mod p).

For p=2 the same formula is valid (minus one equals one). Wilson's theorem therefore gives S(m-1,p-1)==1 (mod p). Consequently

    v_p(F_p(m)) = v_p(p! S(m,p)) = 1,

contradicting (5). Hence m is squarefree. This proves (2)=>(3), while (1)=>(2) is immediate. Applying the same argument to universal scalar congruences, or standard Fermat/Wilson reasoning above, proves equivalence with the stated Carmichael definition.

**Finite constructive witness bound.** If m is neither prime nor Carmichael, some T_r with 2<=r<=P(m), where P(m) is the largest prime divisor of m, violates (3). Thus a counterexample can always be chosen strong of order at most 2P(m)-1, and its coefficient vector can be evaluated modulo m without constructing the much larger clone graph. More precisely, use any prime p for which squarefreeness or p-1 | m-1 fails; the argument supplies r<=p.

## 3. Sharp positive and hostile controls

C3 only detects

    h_v(m)-h_v(1)=2^m-2.

Thus m=341 passes this first test, since it is a base-2 Fermat pseudoprime. It fails the next actual strong test T_3, with h_v(341)-h_v(1)==330 (mod 341). The first skeleton vectors are

    T_2: (3,1),
    T_3: (9,10,2).

The exact receipt records all tested vectors rather than assuming their values. All prime m and every Carmichael m, including 561,1105,1729, pass **all** rooted tournament tests, not merely the small tests used in the receipt.

The user sequence 20,104,464,650,1952 belongs to sigma(m)=2m+2, not to this classification; its members fail already the C3 test. The theorem gives no implication from quasiperfection sigma(m)=2m+1. Matching an aliquot offset and a Hamilton-path offset without an operation is explicitly not a bridge.

## 4. Multimodule survivor and the exact-construction boundary

If several vertices are replaced by TT_m modules with a common prime/Carmichael size m, the Hamiltonian count remains congruent to the original H modulo m. Apply the universal one-module theorem successively to an actual representative of each still-unreplaced vertex. More generally independent module sizes m_i give the inherited multivariate Stirling expansion, and reducing modulo a common d works whenever j^(m_i)==j modulo d for every j and every i.

For one fixed root, if m==1 (mod p-1), then h_v(m)==H(T) modulo p. Thus any attempt to use such clone sizes to hit an exact target N must retain the original residue modulo p. When p divides H(T), all those actual outputs remain divisible by p. This is a useful branch sieve, not an obstruction to arbitrary target-dependent bases or other operations.

**Connection ledger.** Source: integers with universal Frobenius congruence. Target: actual transitive-module substitutions. Map: m maps to replacing a vertex by TT_m. Preserved predicate: H modulo m, uniformly over every base. Destroyed information: exact H and all endpoint/path-cover coordinates. Necessary sidecar: skeleton coefficients and, for endpoint constructors, their own ordered-cover profiles. Cheapest separating test: T_3 at m=341. Missing Chan–Pak step: congruence invariance does not supply a target, an inverse operation, a seed offset of one, or an interval. This lane contributes a precise realizability filter and an all-order arithmetic bridge, not spectrum completeness.

## 5. Reproduction and scope

Run `python3 05-knowledge/results/session260_quasiperfect_chanpak_20261111/bridges/check.py` from the repository root. The receipt contains literal subset-DP checks for all labelled bases through order four, every root and clone sizes 1..4; additional alternating-root controls; and a finite-difference/Korselt comparison for 2<=m<=600. The theorem above is all-order and does not depend on extending this census. No navigation, canon, or other lane files are modified.
