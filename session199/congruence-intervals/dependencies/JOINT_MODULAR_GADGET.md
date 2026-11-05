# Joint modular universality of physical one-apex gadgets

Session 197. **All-order constructive residue coverage, not exact-target coverage.** Every parity-compatible Hamilton-path residue modulo every positive integer occurs in infinitely many strong tournaments with exactly one Hamilton cycle and a one-vertex feedback set. For powers of two, a witness uses only quadratically many vertices in the modulus bit length.

## 1. Physical input and binary recurrence

Use the reviewed integer gadget compiler in `session197-lowcycle/ROOTED.md`, Section4, and `ARITHMETIC_INVERSE.md` (SHA256 `4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7`). For every L>=1 it gives a root b and a transitive base x_0,...,x_(s-1), where s=bit_length(L−1). The root beats x_i exactly when bit i of L−1 is zero. Its Hamilton-path count is g(L), its root-ending two-cover count is L, and exactly one Hamilton path ends at b. Its final base vertex beats b when s>0.

Put f(n)=g(n+1). The exact recurrence is

    f(0)=1,
    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod2).                         (1)

It follows immediately that f(n) is odd. Reading a binary word from left to right updates (n,f) by (1). Leading zero bits are harmless when starting at n=0.

For an even modulus M, the state (n mod M,f mod M) determines both transitions, including the parity of n. All statements below concern actual nonnegative integers and actual physical gadgets; residues are used only to design their binary words.

## 2. A modular loop and complete joint residues

**Theorem1.** For every even M>=2, every u modulo M and every odd z modulo M, there are arbitrarily large integers L with

    L=u (mod M),   g(L)=z (mod M).                 (2)

Write M=2^k q, with k>=1 and q odd. Choose an integer a>=2 with 2^a=1 mod q. For example a=2 phi(q) works, with a=2 when q=1. Append the word

    W=1^a 0^k

to any binary state having n=0 mod M. After the ones, n equals y=2^a−1 modulo M, and f has not changed. The k zeros add

    Delta=y(1+2+...+2^(k−1))+1
         =(2^a−1)(2^k−1)+1.                       (3)

The extra1 occurs only at the first zero, since y is odd. The new n is y2^k=0 mod M. Thus W is a loop on n mod M, with f increment Delta.

Modulo q, Delta=1. If k>=2, Delta=2 mod4 because a>=2. If k=1, Delta is even and is still coprime to q. Consequently

    gcd(Delta,M)=2.                               (4)

Now set r=u−1 mod M in {0,...,M−1}, and let v be its fixed-width binary word of width ceil(log2 M). Starting at n=0 mod M, appending v changes f modulo M by f(r)−1, independently of its initial f. Repeating W j times and then appending v therefore gives

    n=r (mod M),   f(n)=f(r)+j Delta (mod M).

By (4), choose j modulo M/2 to obtain any prescribed odd z. This also handles M=2, where the reduced modulus is1. Adding any positive multiple of M/2 to j preserves both residues and makes the binary word nonempty, with a leading1. Its length grows, so the resulting L=n+1 are arbitrarily large. This proves (2).

This is an effective integer algorithm. Taking a=2 phi(q), a<=2q<=M, and choosing M/2<=j<M gives a word of length at most

    (M−1)(M+k)+ceil(log2 M),

hence O(M^2) vertices in the physical gadget. No efficient factorization claim is needed for this existence or construction bound.

## 3. Short dyadic witnesses

**Theorem2.** When M=2^k, a witness in (2) can be constructed with at most 4k^2+2k+1 bits in n=L−1. Thus its gadget has at most 4k^2+2k+2 vertices.

Begin with the neutral word 10^k. It has n=0 mod M and changes f by M, so the state remains (0,1) modulo M while the integer word is nonempty.

The loop 1^a0^k now changes f by

    2−2^a = −2(2^(a−1)−1) (mod M).                (5)

Let r=u−1 mod M, and choose 0<=j<M/2 with

    f(r)−2j=z (mod M).

Write j in binary. A bit 2^i with i=0 is supplied by one a=2 loop. For i>=1, use one a=i+1 loop and one a=2 loop, since

    2^i=(2^i−1)+1.

At most 2(k−1) loops suffice; each has length at most 2k when k>=2. Finish with the k-bit word for r. For k=1, j=0 and no adjustment loops are needed. The stated generous bound follows in all cases. Every update is an exact binary-word operation, not an approximate carry normalization.

Appending additional copies of the neutral loop before these adjustments gives arbitrarily large witnesses in the same residue class, although the short bound is asserted only for the displayed initial construction.

## 4. Actual strong one-apex count residues

Given a nonempty-base gadget G(L), add a vertex a with all arcs a->G except b->a. The resulting tournament T(L) has

    H(T(L))=g(L)+L−(L mod2).                        (6)

This is the inherited singleton-prefix seam formula; it can also be counted by whether the path uses b->a. Deleting b leaves the transitive order a,x_0,...,x_(s-1). Every Hamilton cycle must use b->a, and then must follow that entire increasing order back to b. Thus T(L) is strong, has exactly one Hamilton cycle, and has feedback vertex number1. The gadget G(L) itself may be weak; this conclusion concerns the displayed singleton-prefix output.

Set u=1 in Theorem1. Then (6) is congruent to z modulo M. Since H(T(L))>=L, arbitrarily large L produce unbounded, hence infinitely many distinct, counts.

**Corollary3.** For every positive modulus m, every residue compatible with an odd integer is attained by infinitely many such strong one-apex tournaments. If m is even these are its odd residues; if m is odd these are all residues. Apply Theorem1 with M=2m, choosing an odd lift of the desired residue. Modulus1 is included.

For modulus 2^k the short construction gives T(L) at most 4k^2+2k+3 vertices. This is a residue witness bound, not a bound for realizing a prescribed integer exactly.

## 5. Every fixed physical prefix has full dyadic residue image

Let (h,e,R) be any fixed state of the reviewed rooted grammar. Its next-gadget count is

    F(L)=h g(L)+R L−e(L mod2).

Here h is odd and R=e mod2. The latter follows, for example, because the L1 seam gives the odd tournament count h−e+R. Fix any residue u modulo M=2^k and any odd target residue N. Since h is invertible modulo M, set

    z=h^(−1)(N−R u+e(u mod2)) (mod M).

This z is odd. Theorem2 constructs a physical L having the required joint residues, and hence F(L)=N mod M. The same holds at arbitrarily large L. Thus no dyadic residue class of odd targets can exclude the final-gadget image of even one fixed prefix.

This statement does not extend without qualification to odd moduli for every prefix: a prefix with common odd content can force every future count to remain divisible by that content.

## 6. Scope and verification

This does not prove that a fixed prefix image is cofinite, that the rooted grammar covers every sufficiently large odd integer, or that infinitely many of the constructed counts are prime or seven times a prime. A set may meet every congruence class allowed by parity and still have many exact-value gaps. The result shows specifically that a fixed-modulus obstruction cannot explain the strong one-apex spectrum's gaps, and that dyadic congruence obstructions cannot explain any fixed physical prefix's gaps.

The standalone verifier implements both word constructors, checks joint residues using the original g recurrence and a separate left-to-right recurrence, tests arbitrary starting modular loop states, and reconstructs small actual adjacency graphs. Independent subset path DP checks those literal graphs and their singleton-prefix counts. Finite controls corroborate the constructive all-order proof; they are not the source of universality.

Run `python3 verify.py` and `python3 -O verify.py`. Neither run writes files. Both receipts must agree.
