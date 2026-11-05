# The exact universal Hamilton-count content of an intact tournament module

Session197. **PROVED all-order divisibility and an explicit finite physical test bank.** For a nonempty tournament T, the gcd of its Hamilton counts in every whole-module tournament context is

    G(T)=gcd_(1<=r<=|T|) (r-1)! r! c_r(T),                (1)

where c_r is the unordered r-path-cover count. A physical eight-vertex example has ordinary bin-value content9 but universal content27. The refinement affects valuations; it does not supply an additional prime/7p obstruction beyond the ordinary content and the input count bound.

## 1. Exact run coefficients and their free permutation factor

Let Q be any rooted tournament with root q. Substitute T as an intact module: every exterior vertex sees all of T uniformly as it saw q. Put O=N_Q^+(q), I=N_Q^-(q). Every output Hamilton path has r>=1 maximal T-runs. These are an ordered r-path cover of T and contribute

    f_r(T)=r!c_r(T).

The exterior consists of an optional initial path ending in I, an optional final path starting in O, and r-1 nonempty internal paths, each starting in O and ending in I. All these exterior paths are vertex-disjoint and together cover Q-q. The two outer paths may be empty; internal paths cannot be empty because the T-runs are maximal.

Let a_(Q,r) count these exterior arrangements with their internal paths ordered. Then cutting and reconstructing the Hamilton word are inverse, so

    H(Q[q:T])=sum_r a_(Q,r) f_r(T).                      (2)

The actual internal paths are distinct nonempty objects on disjoint vertex sets. Permuting their r-1 positions freely preserves validity: every one still starts in O and ends in I, so all intervening module-to-exterior and exterior-to-module arcs are present. This is a free action of the symmetric group on r-1 letters, leaving the two outer paths fixed. Consequently

    a_(Q,r)=(r-1)! b_(Q,r), b_(Q,r) a nonnegative integer. (3)

No fragment is identified merely because it has the same endpoint type. The argument concerns literal paths and all exterior vertices; it does not assume an arbitrary formal operator is realizable. Equation(2) is the one-module maximal-run specialization of the inherited S195 cover calculus. The internal-fragment factorial is proved here explicitly.

Since f_r(T)=0 for r>|T|, equations(2)-(3) show that the right side of (1) divides every context count.

## 2. A unit-triangular bank of actual contexts

For every r>=1 construct Q_r with root q and r-1 out-neighbors O_1,...,O_(r-1) and r-1 in-neighbors I_1,...,I_(r-1). Orient both sides transitively. Across the sides use

    O_i -> I_i,
    I_j -> O_i whenever i!=j.

The incidences at q are q->O_i and I_i->q. For r=1 this is just the identity context.

There are at most r maximal input runs: every internal exterior path uses at least one O and one I vertex, so there can be at most r-1 such paths. At equality, all2r-2 exterior vertices are consumed by these internal paths, the outer paths are empty, and every internal path must be the matching edge O_i,I_i. Those r-1 distinct paths can be ordered in exactly(r-1)! ways. Therefore

    a_(Q_r,j)=0 for j>r,
    a_(Q_r,r)=(r-1)!.

By (3), the matrix b_(Q_r,j), 1<=r,j<=n, is integer lower triangular with diagonal1. Write w_j=(j-1)!j!c_j(T). Then

    H(Q_r[q:T])=sum_(j<=r) b_(Q_r,j) w_j.                (4)

A lower unit-triangular integer matrix has an integer inverse by back substitution. Thus the gcd of the n physical values in (4) is exactly gcd(w_1,...,w_n), proving equality in (1). These n contexts have at most2n-2 exterior vertices, so the test graphs have at most3n-2 vertices. Section3 improves this to only ceil(n/2) contexts using the actual parity of the input polynomial.

The gcd is not asserted to be an attained count. In general it is only an integer linear combination of context counts. The proof also does not claim every multiple of it is attainable.

## 3. Central-factorial coordinates and a shorter physical bank

Let g(T)=gcd_(d>=1) B_T(d) be the ordinary bin-value content. Put d=ceil(n/2). The integer same-parity polynomial B_T has a unique expansion

    B_T(x)=sum_(k=0)^(d-1) u_k q_k(x) for odd n,
    B_T(x)=sum_(k=0)^(d-1) u_k r_k(x) for even n,

where q_k=x product_(j=1)^k(x^2-j^2) and r_k=x^2 product_(j=1)^k(x^2-j^2). Every u_k is an integer: these monic integer polynomials are unit triangular on the monomials of the relevant parity. Nonnegativity is unnecessary for the following content argument. Set

    gamma_k=(2k+1)! for q_k,
    gamma_k=(k+1)(2k+1)! for r_k.

Then the exact contents are

    g(T)=gcd_k gamma_k u_k,
    G(T)=gcd_k k! gamma_k u_k.                           (5)

To prove this, use the polynomial identities

    q_k/gamma_k = binom(x+k,2k+1),
    r_k/gamma_k = binom(x+k+1,2k+2)+binom(x+k,2k+2).

The second follows by adding the adjacent binomial coefficients, or directly by their falling products. Vandermonde's identity expands these normalized polynomials in binom(x,r) with nonnegative integer coefficients N_(r,k). They vanish for r<k+1 and equal1 at r=k+1. Explicitly the coefficients are binom(k,2k+1-r) in the odd case and binom(k+1,2k+2-r)+binom(k,2k+2-r) in the even case, with out-of-range binomial coefficients zero.

The Newton coordinates f_r of B_T therefore satisfy f_r=sum_k N_(r,k) gamma_k u_k. Their first d rows are integer lower unit triangular on gamma_k u_k, and all later rows are integer combinations. This proves the first formula. After multiplying row r by(r-1)!, rewrite on coordinates a_k=k! gamma_k u_k. Each coefficient is ((r-1)!/k!)N_(r,k), an integer whenever nonzero because r>=k+1, again with diagonal1. This proves the second formula and shows

    G(T)=gcd_(1<=r<=d) (r-1)! f_r(T).

Consequently only the first d matching contexts Q_1,...,Q_d from Section2 are needed. They have at most2d-2 exterior vertices, so the physical test graphs have at most n+2d-2 vertices. This is at most2n-1 for odd n and2n-2 for even n.

A sharp uniform valuation bound is

    g(T) divides G(T), and G(T)^2 divides g(T)^3.         (6)

For a prime p, let a=v_p(g), and choose a nonzero gamma_k u_k attaining this minimum. Both definitions of gamma_k are divisible by(k!)^2. Since u_k is an integer, v_p(k!)<=floor(a/2). The second formula in (5) gives

    a <= v_p(G) <= a+floor(a/2).

This includes equality of prime supports and equivalence of divisibility by p^2. The physical example below attains a=2 and v_3(G)=3, so the exponent3/2 cannot be reduced as a uniform bound. The argument applies to any nonzero integer polynomial of one parity with zero constant term, using its Newton-coordinate gcds; physicality is needed for their interpretation as context counts.

## 4. Comparison with ordinary bin content

Define

    g(T)=gcd_(d>=1) B_T(d).

Since B_T(d)=sum_r binom(d,r) f_r(T), evaluating at d=1,...,n is an integer unit-triangular change of coordinates. Therefore

    g(T)=gcd_r f_r(T),
    g(T) divides G(T).                                  (7)

The two contents have exactly the same prime divisors. Indeed if p divides every (r-1)!f_r, then for r<p the additional factorial is a p-unit, so p divides f_r. For r>=p, f_r=r!c_r is automatically divisible by p. Thus p divides every f_r. The converse follows from (7).

For an independent verification of the low-valuation consequence of (6), use the actual parity of bin polynomials:

    p^2 divides G(T) if and only if p^2 divides g(T)      (8)

for every prime p. Only the forward direction needs proof. If p^2 divides G, then for1<=r<=p, (r-1)! is a p-unit, so p^2 divides f_r, with f_r=0 past the input order. Hence B_T(1),...,B_T(p) vanish modulo p^2. Write P=B_T; it has integer coefficients, P(0)=0 and P(-x)=epsilon P(x), epsilon=(-1)^|T|.

For1<=a<=p-1, Taylor expansion modulo p^2 gives

    P(p-a)=P(-a)+pP'(-a) (mod p^2).

Both P(p-a) and P(-a)=epsilon P(a) vanish, so p divides P'(-a), equivalently P'(a). Likewise P(p)=pP'(0) modulo p^2 forces p|P'(0). Now for every integer x=a+tp,

    P(x)=P(a)+tpP'(a)=0 (mod p^2).

Thus p^2 divides g. Integer coefficients, not merely integer-valued coefficients, justify these Taylor congruences. This handles p=2 as well, though physical H(T) is odd and both contents are odd.

Consequently strict valuation growth from g to G can first occur when g already has that prime to exponent at least2. It is not justified to infer a new7p obstruction from a hypothetical content7-to-content49 jump.

More precisely, for a target N which is prime or7 times a prime, the two necessary conditions g(T)|N and H(T)<=N already imply G(T)|N. If N is squarefree this follows from prime support and (8). The only nonsquare case is N=49. Here g is1,7 or49; the first two cases are controlled by the same statements, while g=49 and G|H(T)<=49 force G=49. Thus (1) is a sharper exact invariant, but not a new filter for those atoms once the existing input-count bound is enforced.

## 5. A physical strict-valuation example

Take the actual tournament

    T=[30,28,136,192,172,15,55,99].

Its unordered cover vector is

    (c_1,...,c_8)=(459,2889,5310,4125,1560,300,28,1),

and its bin polynomial is

    B_T(x)=x^8+34x^6+214x^4+210x^2.

The ordered-cover vector is

    (459,5778,31860,99000,187200,216000,141120,40320),

whose gcd is9. Multiplying coordinate r by(r-1)! instead gives

    (459,5778,63720,594000,4492800,25920000,101606400,203212800),

whose gcd is27. In particular the ordinary content does not determine the exact universal content, even for a physical positive bin polynomial.

The values B_T(1),...,B_T(4) are459,6696,50571,262944. The last has exactly two factors of3. The identity and triangle contexts have actual Hamilton counts459 and7155, whose gcd is27. These two small physical contexts already certify that the universal divisor27 is sharp, while (1) proves it divides every other context count.

It is not an attained count: every context graph contains T as an induced subtournament and has H at least459. To see monotonicity directly, extend each Hamilton word of T by inserting the additional vertices one at a time by a deterministic Hamilton-path insertion rule. Deleting those added vertices recovers the original word, so this is an injection.

## 6. Composition and exact scope

If S is obtained from T by one whole-module substitution, then G(T) divides G(S). Any later whole-module context applied to S is still a whole-module context applied to T, because the old T remains an intact module in the composed template. Taking gcds proves the claim.

This divisibility is not asserted under arbitrary nonuniform births, deleting vertices, or graph rebuilding. The theorem computes a gcd over all contexts; it is not a construction of a prescribed multiple and does not close the unrestricted inverse problem.

## 7. Reproduction and inherited sources

Run `python3 verify_universal_module_content.py` and its `-O` form. The checker counts exterior initial/internal/final fragments directly, checks the free factorial on all rooted labelled contexts through order5, verifies physical clone counts, builds the perfect-matching test bank and its unit diagonal, checks the exact central-factorial/Newton transformations and valuation bound, compares the shortened-bank gcd descriptions on small actual inputs, and directly counts the eight-vertex example and its triangle wrapper.

The maximal-run framework is the pinned S195 proof https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md . The frozen `PROOF.md`, SHA256 `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`, records the bin/parity conventions. The computational helper is `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`. The proof of the factorial and finite bank is included here and does not depend on an unproved content-equivalence conjecture.
