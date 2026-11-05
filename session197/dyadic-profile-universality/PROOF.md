# Independent dyadic control of a physical one-apex profile

Session 197. **PROVED, constructive residue statement; not exact-target coverage.** The integer gadget compiler can prescribe the Hamilton count and both rooted two-cover responses simultaneously modulo any power of two. The witnesses have order quadratic in the modulus bit length. Consequently actual strong start-one sources with both endpoint counts and root-deletion count equal to one have independently arbitrary odd Hamilton-count and starting-two-cover residues.

## 1. Physical conventions and inherited laws

Use ROOTED.md Section 4 (SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e). For L>=1, put n=L-1 and s=bit_length(n). The actual gadget G(L) has a transitive base x_0,...,x_(s-1) and root b. The arc x_i->b occurs exactly when bit i of n is one; otherwise b->x_i. When s>0 its last base vertex points to b.

Its Hamilton count is g(L), its ending-two-cover response is L, and exactly one Hamilton path ends at b. Its starting-two-cover response is Q(L), while the two-bin count of the root deletion is B(L)=2^s. Empty path slots are allowed as in ROOTED.md. Write f(n)=g(n+1). The inherited physical recurrence is

    f(0)=1,
    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod 2).                         (1)

In particular every f(n) is odd. The previously reviewed modular construction is MODULAR_GADGET_UNIVERSALITY.md, SHA256 74df1989d75ff743177937c108ac36dce133412e48701e4d412297d8c47b7d63, published at commit 48e10c6a4773de23f9ca19ef120b17012a0d7812. The new ingredient here is independent control of Q by the leading bits, while the loop controls f and the final suffix controls L. We reprove the needed loop law.

For a fixed-width binary word w, let revc(w) be its reversed complement, interpreted as a nonnegative integer. The exact gadget response is

    Q(L)=1+revc(bin_s(n)).                          (2)

Indeed a rooted first path consists of b followed by an increasing chosen base subset. If its first base vertex is x_i, then b->x_i is necessary and the later chosen vertices are arbitrary. The complementary vertices form the second increasing path. Thus Q=1+sum_(i:bit_i(n)=0) 2^(s-1-i), proving (2). The width is the actual minimal width s; arbitrary leading-zero padding would change the physical response.

## 2. Simultaneous joint-profile theorem

**Theorem 1.** Let k>=1 and M=2^k. For every residue u modulo M and odd residues z,q modulo M, there are arbitrarily large integers L such that

    L=u,   g(L)=z,   Q(L)=q,   B(L)=0  (mod M).     (3)

One witness can be chosen with at most 4k^2+3k bits in n=L-1, hence G(L) has at most 4k^2+3k+1 vertices. All witnesses below have a nonempty base and their binary words have a leading one.

**Proof.** Choose the unique k-bit word P such that

    revc(P)=(q-1) mod M.

Since q is odd, the last bit of revc(P) is zero; therefore the first bit of P is one. Begin the binary word for n with P, followed by k zeros. This fixes its first k bits permanently and puts n=0 modulo M. Let f_* be the exact current value obtained from (1).

For any a>=2, appending 1^a0^k from a state with n=0 modulo M returns n to zero modulo M. The ones leave f unchanged. The k zeros increase it by

    (2^a-1)(2^k-1)+1 = 2-2^a (mod M).             (4)

The extra one arises at the first zero, whose preceding integer is odd. All later preceding integers are even. Thus this loop subtracts 2(2^(a-1)-1) from f modulo M.

Set r=(u-1) mod M, and reserve its k-bit word as the final suffix. From n=0 modulo M, that suffix changes f by f(r)-1 modulo M. This follows directly from (1), since M is even and therefore both the previous residue and its parity determine each update.

Choose 0<=j<M/2 with

    2j=f_*+f(r)-1-z (mod M).                       (5)

The right-hand side is even, so this is well-defined, including k=1 where the reduced modulus is one. Write j in binary. A bit 2^0 is supplied by one a=2 loop. A bit 2^i with i>=1 is supplied by an a=i+1 loop followed by an a=2 loop, since (2^i-1)+1=2^i. These loops subtract exactly 2j modulo M.

Finally append the reserved suffix. The resulting integer has L=u and g(L)=z modulo M. Its first k bits remain P. In (2), the low k bits of the reversed complement depend exactly on these first k bits, so Q(L)=1+revc(P)=q modulo M. The full word has length at least 3k, so B(L)=0 modulo M.

For k>=2 there are at most 2(k-1) correction loops, each of length at most 2k. The initial prefix/reset and final suffix total 3k bits. This gives the stated generous bound 4k^2+3k; for k=1 there are no correction loops and exactly three bits suffice.

To obtain arbitrarily large witnesses, insert any number of copies of 10^k after the initial reset. From n=0 modulo M this returns n to zero and increases f by M, hence changes none of its residues. It also preserves the first k bits. The integer length grows strictly, completing the proof.

## 3. Actual strong start-one sources with independent residues

Here restrict to L>=2, as holds for every witness of Theorem 1. Add a vertex a to G(L), with every arc a->G(L) except b->a, and mark b. Call this T(L). The excluded case L=1 would give only a transitive edge, with root-ending count zero and no Hamilton cycle. The inherited one-seam identity, or a direct split at the sole reverse arc, gives

    H(T(L)) = g(L)+L-(L mod 2),
    e_b(T(L)) = s_b(T(L)) = H(T(L)-b) = 1,
    R_b(T(L)) = Q(L)+B(L).                         (6)

The equality for root deletion is literal: a followed by the base is transitive. Every Hamilton path starting at b must next use b->a and then traverse that transitive deletion. Every path ending at b must traverse a and the base in increasing order before b. The rooted two-cover formula is the singleton-prefix case of ROOTED.md equation (1).

Every Hamilton cycle must use b->a, the sole incoming arc of a, and then follow a,x_0,...,x_(s-1),b. Thus T(L) is strong, has exactly one Hamilton cycle and has feedback vertex number one. These assertions concern T(L), not the possibly weak gadget G(L).

**Corollary 2.** For every k>=1 and every odd residues H_0,R_0 modulo 2^k, there are arbitrarily large physical rooted sources T with

    H(T)=H_0, R_root(T)=R_0 (mod 2^k),
    s_root(T)=e_root(T)=H(T-root)=1.

They can be chosen strong with exactly one Hamilton cycle and feedback vertex number one. A witness has at most 4k^2+3k+2 vertices. Apply Theorem 1 with u=1, z=H_0 and q=R_0. The Hamilton counts are unbounded because H(T(L))>=L, so this is not merely an assertion about arbitrarily many padded presentations of a bounded count.

Appending t>=1 universal sinks after T and retaining the root gives the further physical response

    (H,e,R,d)=(H(T),0,2^t R_root(T),1).             (7)

Each appended sink preserves H and d, kills the root-ending count, and doubles R by its two possible path slots. The root-start count remains one. This enlarged source need not be strong; it is an actual start-one source. Thus, for any fixed t, every odd H residue and every R residue of exact 2-adic valuation t can be prescribed modulo 2^k when k>t, with O(k^2+t) vertices. No claim of independent arbitrary four-coordinate control is made.

## 4. Scope and reproducibility

The proof prescribes congruences, not exact target values. It does not prove cofinite coverage, prime infinitude, or unrestricted realizability of formal profiles. In particular the exact endpoint and deletion coordinates in Corollary 2 are fixed at one; the statement does not assert that arbitrary values of those coordinates can be prescribed. The leading-bit construction is tied to powers of two and is not asserted for general odd moduli.

The standalone verifier uses both the halving recurrence and a left-to-right binary evaluator, checks Q by a separate physical recurrence, and tests arbitrarily large neutral-loop insertions. It reconstructs literal small gadgets and strong outputs and counts their Hamilton paths, endpoint counts and rooted two-covers by subset dynamic programming. Larger bit-length controls use exact integers, not floating point. Run `python3 verify.py` and `python3 -O verify.py`; neither writes files. Their receipts must agree. The finite controls check the formulas and physical interface; the all-order conclusion is proved above.
