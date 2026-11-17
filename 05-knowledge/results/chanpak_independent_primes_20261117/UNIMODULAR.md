# A unimodular inverse for distinct prime targets

Take the seven-vertex exterior R₀ with outgoing masks (20,57,50,21,32,73,31). Replace its last vertex by any nonempty tournament B, uniformly at its contacts, to get R. Append q₀ beating precisely old exterior vertices with mask 42, and one new return vertex d beating only q₀. The old R is induced, with no orientation changes. Both quotients are strong; the substitution preserves strength.

Writing h=H(B), P=P₂(B), the exact reads of the two tournaments are

    H(R)=44P−13h = p,
    H(R+q₀+d)=7(71P−21h) = 7q.

In each quotient the variable slot has only one entrant, so a full path can meet it in at most two runs, and both H counts are linear reads of (P,h). Literal singleton and transitive-pair substitution determines the two coefficients in each read. [The replay](replay_pairs.py) also checks a nontransitive substitution.

**Target inversion.** The two-coefficient map has determinant −1 in the order (P,h), so the independently supplied p,q prescribe a unique profile:

    P=21p−13q,     h=71p−44q.

No lattice-index obstruction remains at the level of integer pairs. Physical profiles are much more restrictive. For instance for both p,q prime >3, B must have even order and p+q≡0 mod 3: modulo 3, P≡−q and h≡2p+q, while P≡(−1)^|B|h by the odd-pack expansion. For |B| even nonzero P is divisible by four, forcing p≡q mod 4. With |B|=n at least two, the split bound forces

    ((n+1)·71−21)p ≤ ((n+1)·44−13)q,     q < 71p/44.

Condition P≥2^n also bounds the order of any compatible source. The inverse therefore supplies an exact test at absolute p,q, not a universal pairing claim. At p=163, q=263 it returns (P,h)=(4,1), and the target inequalities force n=2, so they select the transitive pair. This gives a literal strong 163-path tournament inside a strong 1841-path tournament (7·263).

A nontransitive source B=(52,1,10,35,14,22) has (h,P)=(37,376), giving the distinct primes p=16063 and q=25919 and the induced sevenfold count 181433. For distinct primes the determinant also forces gcd(h,P)=gcd(p,q)=1. Every automorphism of B acts freely on each ordered path cover (fixing its slot paths would fix every position and thus every vertex), so its group order divides both coordinates. Hence any such source selected by distinct primes is rigid.

A useful negative control is (1747,2819), two certified primes. They pass the stated parity and slope bounds (with a possible even order) and invert to (P,h)=(40,1). But h=1 forces a transitive tournament, whose P is a power of two. So the exact pair is not physical for this compiler despite clearing the integer lattice and basic bounds.

This is not a missing-spectrum control: the separate strong tournaments with outgoing masks (220,289,250,434,162,321,410,290,21) and (222,148,568,1186,1000,387,1710,4,207,427,951) have counts 1747 and 19733=7·2819, respectively. Literal subset DP checks both. So even at a pair whose two targets are both realizable independently, the unique shared profile can fail to exist.

For even source order at least four one can tighten p≡q mod4 to p≡q mod8. In the odd-pack expansion of P₂, a term with c odd cycles and s singletons has weight 2^(2c+s). For even order ≥4, no-cycle terms are divisible by 16, one-cycle terms have s≥1, and two-or-more-cycle terms are divisible by 16. Thus P₂ is divisible by 8 and P₂/8 modulo two records exactly the parity of directed cycles of length |B|−1. The order-two pair (163,263) is the small exceptional case.
