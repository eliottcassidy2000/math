# A common two-bin divisor bounds prime-target core depth

Status: **PROVED, independently audited in the stated scope.** This companion is separate from the frozen entropy proof. It gives a target-sensitive exact reduction of an unbounded-depth physical grammar, not an occupancy theorem.

## 1. Physical operations and inherited transfer

For a rooted tournament use the nonnegative retained state

    (h,x,y,d)=(H(T),R_r(T)-H(T-r),H(T-r)-end_r(T),H(T-r)).

For an actual non-singleton rooted core C write H_C=H(C-r), beta_C=P_(C-r)(2), and

    lambda_C=L_C-H_C, alpha_C=R_C-H_C,
    sigma_C=H_C-start_r(C), eta_C=H_C-end_r(C).

The audited root-amalgamation law in `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`, Section 3, gives the literal graph product T*C:

    h' = H_C h + lambda_C x + sigma_C y + (H(C)-H_C)d,
    x' = beta_C x + alpha_C d,
    y' = eta_C d,
    d' = H_C d.                                        (1)

Identify the roots, orient every old nonroot vertex toward every new nonroot vertex, and retain the shared root. This specifies every arc.

The two singleton-block operations have exact laws

    O(h,x,y,d)=(h,2x+d,d,d),
    Z(h,x,y,d)=(h+x+y,2x,0,d).                          (2)

O appends a universal sink; Z appends a vertex beating only the retained root. Begin at the rooted singleton, state (1,0,0,1). Arbitrarily many O/Z letters are allowed before, between and after core operations.

## 2. Divisibility filtration

**Theorem.** Let c>=2 be an integer dividing both H_C and beta_C for every core in a specified bank. The bank may be infinite. After m core operations and arbitrary intervening O/Z letters,

    c^m divides d,
    c^(m-1) divides x and y            (m>=1),
    c^(m-2) divides h                  (m>=2).         (3)

**Proof.** At m=0 no divisibility of h,x,y is required, and d=1. At m=1 the last two assertions, when applicable, are vacuous; the d assertion follows from (1). At the next core, d is multiplied by a multiple of c. The term beta_C x is divisible by c^m, and alpha_C d is divisible by c^m, so the new x is divisible by c^m. The new y is a multiple of the old d and therefore of c^m. For m>=2, H_C h is divisible by c^(m-1), lambda_C x and sigma_C y are divisible by c^(m-1), and the remaining term is divisible by c^m. Hence the next h is divisible by c^(m-1). At the boundary m=1 -> 2 its required exponent is zero, so it is automatic.

Each operation in (2) preserves all current divisibility assertions: the only additions are from x,y into the lower-grade coordinate h and from d into x or y. Thus arbitrary binary intervals do not change the induction. QED.

The proof is an integer filtration on a physically returning interface. It does not require choosing profile coordinates independently or silently deleting internal binary letters.

## 3. Exact consequences for target searches

Deleting the root leaves the ordered join of the core deletions and all singleton blocks. Thus d is the product of the H_C. In particular d>=c^m. Hamilton paths of T-r inject into those of T by deterministic insertion, so h>=d>=c^m.

If m>=3, (3) gives the proper nontrivial divisor c^(m-2) of h. Therefore a prime output can use at most two cores. More generally, for any prescribed positive target n with gcd(n,c)=1, at most two cores are possible. The exact depth condition is

    c^(m-2) divides n                         (m>=2). (4)

This is necessary, not sufficient. For c prime it says m<=2+v_c(n).

For a bank whose core deletions consist of one directed triangle plus any number of ordered singleton components, H_C=3 and beta_C=12*2^k, where k is the number of singleton deletion vertices. Thus c=3 applies. Every odd prime other than 3, and every target 7p for prime p>=5, needs at most two such cores. Arbitrary binary spacers remain allowed and unbounded.

This is a useful positive simplification of the inverse problem: for those prime-atom targets, arbitrary-depth search in this bank reduces exactly to zero, one or two core insertions with at most three intervening binary segments. It does not assert that the bank represents every prime atom. For a finite bank the discrete core choices are finite, while each binary segment still has an unbounded actual-word/index parameter. An extra dyadic-length parameter may be needed to encode a middle binary segment; unlike the initial singleton prefix, its leading Z letters are not generally neutral.

## 4. The fixed triangle-block family has only five prime values

For the narrower family with no binary spacers, each rooted core is a triangle plus the shared root. Its isomorphism type is determined by j=0,1,2,3, the number of triangle vertices beaten by the root. Independent cyclic relabellings of the blocks preserve all cross arcs. The exact core data are

    H_C=3, beta_C=12,
    H(C)=(3,5,5,3)[j], lambda_C=9-3j,
    alpha_C=3j, sigma_C=3-j, eta_C=j.                 (5)

These data are directly counted in the inherited root-amalgamation packet and checked by literal graphs here. Formula (1) and the prime-depth theorem reduce every prime output to m<=2. At m=1 the values are 3 and 5; the sixteen j-pairs at m=2 have value set

    {9,15,31,41,45,61,75,99}.

Consequently the exact prime palette of the entire arbitrary-length fixed-triangle shared-root family is

    {3,5,31,41,61}.                                  (6)

Arbitrary ordered joins cannot add prime values, because every joined factor has a positive integer count. This statement deliberately excludes internal singleton spacers: they preserve the two-core bound but can produce additional prime values with zero, one or two triangle cores.

## 5. Verification and outcome

`verify_content.py` checks (5) against literal adjacency-DP root cores; checks all fixed-triangle incidence words through three blocks against the inherited physical matrix law; checks the divisibility filtration on a declared complete short alphabet universe; and checks the prime palette through the theoretically complete m<=2 bound. It separately provides a literal H=13 graph using one triangle core and a binary prefix, illustrating why singleton spacers cannot be silently excluded from the richer grammar.

This adds no universal occupancy theorem and gives no unconditional threshold near 23. It isolates a bounded core-depth arithmetic problem for the prime and 7p lanes, while the entropy companion rules out one tempting no-spacer reconstruction class altogether.
