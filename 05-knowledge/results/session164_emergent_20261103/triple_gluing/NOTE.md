# Literal endpoint gluing, a sign interface, and Gaussian carries

**Status:** proved compositional synthesis with exact finite controls. The endpoint/Pythagorean codec, Gaussian factorization, rational-circle group law, Cassini identity, and two-valued quotient construction have prior sources. This packet's contribution is their explicit composition with a literal tournament gluing, its exact information-loss ledger, and the consequences below. No literature-priority claim is made. It does not construct the missing strong ternary spectrum branches.

The central result is a useful distinction: forgetting global conjugation is harmless for **one** primitive Pythagorean triple, but forgetting each input's conjugation independently before composition loses a necessary relative bit. Shared-prime normalization can then erase further bits, with exactly computable multiplicities.

## 1. A literal two-terminal gluing

For an integer r>=1, let E(r) be the endpoint-reversal tournament with vertices `0,1,...,r,r+1`: every arc points forward in this order except `r+1 -> 0`. Its two terminals are 0 and r+1; its r interior vertices retain their order. Write I for its increasing Hamilton path and N(r) for all other Hamilton paths. The inherited decoder is

    P_A = (A increasing, r+1, 0, A^c increasing),  A subset {1,...,r}.

Every nonincreasing path must use the sole backward edge, so this is a bijection and `H(E(r))=1+2^r`.

Given E(r) and E(s), identify their 0-terminals and their final terminals, retain every old arc, and orient all first-factor interior vertices toward all second-factor interior vertices. This actual tournament, denoted `E(r) diamond E(s)`, is E(r+s) in the concatenated interior order. It is strong: its vertices lie on the directed Hamilton cycle `0,1,...,r+s,r+s+1,0`.

**Theorem 1 (path multiplication).** There is an explicit bijection

    N(r) x N(s) -> N(r+s),
    (P_A,P_B) -> P_(A disjoint-union B),

and therefore

    H(E(r) diamond E(s)) - 1
      = (H(E(r))-1)(H(E(s))-1).

Here B is shifted to the second interior block. More concretely, concatenate the first path's pre-terminal interior run with the second's pre-terminal run, traverse the single shared backward edge, and concatenate the two post-terminal runs. Each cross-block arc used is one of the specified new arcs. The inverse restricts the output address to the two interior blocks. This proves the bijection, not merely the count identity. Repeated gluing is associative up to the canonical ordered relabelling. The exceptional path is excluded on both sides; it is not another Gaussian factor or a degenerate triple.

The inherited weighted path signature also composes:

    W_C(z) = product_(p^e || C) (z_p^e + z_p^(-e)),
    W_(glued slots)(z) = W_C(z) W_D(z).

When C and D share a prime, the gluing still keeps two distinct interior slots. Identifying their arithmetic variable is a later decoding step, not identification of those vertices.

## 2. The one-bit interface theorem

Let C,D>1 have only primes congruent to 1 modulo 4 and have disjoint prime supports. Fix a Gaussian prime pi_p of norm p for each prime. An oriented address chooses either pi_p^e or its conjugate at each slot. Its product alpha has norm C or D. Squaring alpha and taking absolute coordinates yields the parity-ordered primitive triple

    PPT(alpha) = (|Re(alpha^2)|, |Im(alpha^2)|, |alpha|^2).

Complementing all address bits conjugates alpha. The inherited codec identifies precisely these complementary address pairs, giving `#PPT(C)=2^(omega(C)-1)`.

**Theorem 2 (relative conjugation).** Restriction of the output address gives a two-to-one map

    PPT(CD) -> PPT(C) x PPT(D).

For any chosen input lifts alpha,beta, the two output classes are

    PPT(alpha beta),  PPT(alpha conjugate(beta)).

They are distinct. Choosing one relative-conjugation bit gives a bijection between the input pair plus that bit and the output fiber.

**Proof.** The four oriented lifts of the input classes are `(alpha,beta)`, `(conj(alpha),beta)`, `(alpha,conj(beta))`, and `(conj(alpha),conj(beta))`. Simultaneous conjugation identifies them in pairs. No further equality is possible: uniqueness of Gaussian prime allocations and disjoint nonempty supports would otherwise equate an input address with its own complement. The literal path gluing in Theorem 1 realizes the two surviving classes. Equivalently,

    #PPT(CD) = 2 #PPT(C) #PPT(D).

The factor 2 measures the missing interface, rather than a new unaccounted family of paths. QED.

Example: the unique norm-5 and norm-13 input triples are `(3,4,5)` and `(5,12,13)`. The lifts `2+i` and `3+2i` give

    (2+i)(3+2i) = 4+7i       -> (33,56,65),
    (2+i)(3-2i) = 8-i        -> (63,16,65).

Thus no operation on the unmarked input triple classes can agree with Gaussian multiplication for **every** independent choice of lifts. A deterministic convention is possible, but it chooses the missing bit; it does not show the bit was unnecessary.

For k inputs with pairwise disjoint nonempty prime supports, the same argument gives exactly `2^(k-1)` distinct outputs over a fixed input tuple. Consequently k-1 binary choices are both sufficient and necessary to select one output without loss. They can be stored as relative signs along any tree on the k inputs. If signs are specified on a larger connected graph, they are consistent precisely when their product around every cycle is +1: assign a root sign, propagate along a spanning tree, and check the remaining edges. This is an exact sign-interface statement. It is not the physical observer port `(h,x,y)` used by the tournament-return packets.

## 3. Shared primes: normalization has an exact carry

Let v be a finitely supported integer vector indexed by split primes `p=1 mod 4`. Define

    alpha(v) = product_p pi_p^(max(v_p,0)) conjugate(pi_p)^(max(-v_p,0)),
    N(v) = product_p p^|v_p|.

The two integer coordinates of alpha(v) are coprime: a rational prime dividing both would require both members of a conjugate prime pair. For v!=0, they give a nondegenerate primitive triple. For v=0, alpha(v)=1 and the convention is the degenerate `(1,0,1)`.

**Theorem 3 (carry and reassociation).** Set

    g(v,w) = product_p p^((|v_p|+|w_p|-|v_p+w_p|)/2).

Then

    alpha(v) alpha(w) = g(v,w) alpha(v+w),
    gcd(Re(alpha(v)alpha(w)), Im(alpha(v)alpha(w))) = g(v,w),
    N(v+w) = N(v) N(w) / g(v,w)^2.

Moreover,

    g(u,v) g(u+v,w) = g(v,w) g(u,v+w).

**Proof.** At one prime, signs that agree add exponents and produce no rational factor. Opposite signs pair exactly `min(|v_p|,|w_p|)` copies of pi_p with conjugate(pi_p), each pair producing the rational integer p. The remaining exponent is v_p+w_p and is entirely on one conjugate side. This proves the product formula and exact coordinate gcd. Norms give the square loss. The carry exponent on either side of the final identity telescopes to

    (|u_p|+|v_p|+|w_p|-|u_p+v_p+w_p|)/2.

QED. Primitive normalization is therefore associative when the oriented vector is retained, and its accumulated scale is independent of parenthesization.

After forgetting the global sign, the well-defined operation is the **unordered multiset**

    [v] star [w] = [[v+w], [v-w]],  where [v]={v,-v}.

It is two-valued, not an ordinary group law. Its associativity follows by expanding either parenthesization to the four classes `[u+v+w]`, `[u+v-w]`, `[u-v+w]`, `[u-v-w]`, counting multiplicities. Its unit yields two copies of the input; the zero class must be retained for cancellation. This is the standard quotient-by-negation mechanism, expressed in the marked prime lattice.

### Exact fibers after a cancellation

Fix exponent profiles e_p for C>1 and f_p for D>1. At a common prime, choose aligned signs or opposed signs. The resulting absolute exponent is respectively `e_p+f_p` or `|e_p-f_p|`. These choices determine a unique output norm, since the two absolute exponents differ.

Suppose the selected branch cancels z primes completely, meaning `e_p=f_p>0` with opposed signs. If the output norm is greater than one, every primitive output triple has **exactly `2^z` input address-pairs modulo simultaneous complement**. If the output norm is one, the sole degenerate output has `2^(z-1)` such pairs.

Indeed every nonzero output coordinate reconstructs its input signs uniquely. A zero coordinate admits the two opposed input sign pairs; these are independent. A nonzero output class has two oriented lifts, whose common global factor two cancels the simultaneous-complement quotient. At zero there is only one output lift, so the final division by two remains. This proves the exceptional formula too.

At equal shared exponent e, the local weighted signature is

    (z^e+z^(-e))^2 = z^(2e) + 2 + z^(-2e).

The coefficient 2 is precisely the sign information erased by complete cancellation. Calling the three monomials three states does not make them the equally weighted ternary spectrum branches `3h,3h+2,3h+4`.

## 4. The divisor-balance shapes become actual composition branches

Use the established definitions: F counts proper nontrivial divisors, S counts the squarefree ones among them, and U counts proper prime divisors. For N≥2, the classification `F=S+U` has exactly the shapes p, p^3 and p^2qr, with distinct primes in the last shape. Impose the additional split-prime condition only to use the Pythagorean decoder.

1. **The two prime-power shapes.** Compose the unique norm-p and norm-p^2 triple classes. Their oriented valuations are ±1 and ±2 at p. The two normalized branches have exponents ±3 and ±1, so their norms are exactly p^3 and p. Both outputs satisfy the divisor balance. At p=5 the outputs are `(117,44,125)` and `(3,4,5)`.
2. **The three-prime shape.** Compose the norm-pq and norm-pr fibers. The shared p-slot has exponent one in each input. Alignment produces p^2qr with no rational carry; opposition produces qr with carry p. There are eight input address-pairs modulo simultaneous complement. Four aligned pairs give the four distinct primitive p^2qr triples bijectively. Four opposed pairs give the two qr triples, each twice. Thus the normalized output satisfies divisor balance **if and only if the shared p signs align**. The other branch is squarefree composite and has F=S but U=2.

For `(p,q,r)=(5,13,17)`, the balanced branch gives

    (2163,5084,5525), (4557,3124,5525),
    (5133,2044,5525), (5427,1036,5525),

each once; the cancelled branch gives `(21,220,221)` and `(171,140,221)`, each twice. The source tournament is the literal `E(2) diamond E(2)=E(4)`, whose 16 nonexceptional paths become eight complementary pairs. The normal forms have address supports three and two respectively. This is a map of path sets with explicit fibers, not an assertion that contracting vertices in E(4) preserves Hamilton paths.

This gives a compositional interpretation of the previously classified shapes, not a new proof that no other exponent profile satisfies the divisor equation. The balance is not preserved by arbitrary multiplication or normalized composition.

## 5. A Fibonacci bit selects a whole positive square-sum family

Write F_0=0,F_1=1 and

    alpha_n = F_(n+1) + i F_n.

For every n>=1, Cassini and the doubling identity give the exact Gaussian equality

    alpha_n conjugate(alpha_(n+1))
      = F_(2n+2) + i (-1)^(n+1).                         (1)

The real part is `F_(n+1)(F_n+F_(n+2))=F_(2n+2)`; the imaginary part is `F_n F_(n+2)-F_(n+1)^2=(-1)^(n+1)`. Both identities follow by induction from the Fibonacci recurrence (or from powers of its 2-by-2 companion matrix).

Take n=3j+2, j>=0. Adjacent Fibonacci numbers are coprime, and `F_m` is even exactly when 3 divides m. Thus both alpha_n and alpha_(n+1) have coprime coordinates of opposite parity. Their odd norms are

    C=F_(2n+1), D=F_(2n+3),

using `F_m^2+F_(m+1)^2=F_(2m+1)`. These norms are coprime: `gcd(F_a,F_b)=F_gcd(a,b)` and `gcd(2n+1,2n+3)=1`. Alternatively apply the Fibonacci Euclidean algorithm directly. All their primes split modulo 4 because they are primitive sums of two squares. Hence Theorem 2 applies to the actual endpoint path fibers.

**Theorem 4 (Cassini branch).** The relative-conjugate branch in (1) selects the primitive gap-two triple

    (K^2-1, 2K, K^2+1),  K=F_(2n+2), n=3j+2.

Here K is even and K>=8. In particular it selects one explicit primitive triple in the norm fiber `CD=K^2+1`, not just the correct number of triples. Choosing arithmetic representatives alpha_n and alpha_(n+1) fixes which interface bit is meant.

This selected branch always continues to a **positive** arithmetic-progression square-sum triple. Put a=K^2-1, b=2K, c=K^2+1 and set

    L = 2(c^2-4ab, c^2, c^2+4ab).

Its pairwise sums are respectively `4(a-b)^2`, `4c^2`, `4(a+b)^2`; its entries form an arithmetic progression. Positivity follows from

    c^2-4ab = K^3(K-8)+2K^2+8K+1 > 0  for K>=8.

The first instances are

    n=2: C,D=5,13;   PPT=(63,16,65);       L=(386,8450,16514),
    n=5: C,D=89,233; PPT=(20735,288,20737).

The other n=2 bit gives `(33,56,65)` and a negative first AP label, so the sign interface already distinguishes an admissible positive continuation from a rejected one. The assertion about positivity is specific to the selected Cassini branch; no claim is made that all other branches fail. These are square-sum triangles, not Hamilton paths in an entire consecutive-label square-sum graph.

## 6. What does not transfer to the ternary guard

The rational Collatz triple chart retains the signed even leg through `(s,t)` with odd coprime positive s,t. In the canonical positive-even-leg chart of a primitive triple with Euclid parameters n>m>0, `s=n+m` and `t=n-m`. The prior ternary normalization reaches the integer stratum exactly when t is a power of 3; with an arbitrary prime-to-3 factor in t, that factor remains.

The present multiplication does not preserve this guard. With pi_5=2+i, the canonical charts are

    norm 5:   (s,t)=(3,1),
    norm 125: (s,t)=(13,9),
    norm 625: (s,t)=(31,17).

The first two have the ternary denominator guard. Their same-orientation primitive product has norm625 and denominator17, so loses it. The carry is one in this example; cancellation is not the reason for failure. The first Cassini branch `(63,16,65)` also has canonical `(s,t)=(9,7)`, outside that guard. Thus Fibonacci positivity, Gaussian compositional closure, divisor balance, and ternary source closure remain genuinely different predicates.

A physical ternary-return implementation would still have to supply actual successor tournaments with the required observer state. Neither the arithmetic sign interface nor the literal endpoint gluing yields the missing affine offsets or their closure.

## 7. Sources and metanalysis

- [Endpoint/Pythagorean fibers, 2026-10-19](https://github.com/eliottcassidy2000/math/blob/f18780daea94f262eb1b7d6fd174015e2657ea86/05-knowledge/results/endpoint_pythagorean_fibers_20261019.md): inherited path decoder, marked Gaussian bijection, weighted signature, and AP square-sum continuation. The new step here is to put a literal composition before the quotient, exposing its relative-sign loss and repeated-prime fibers.
- [S163 triple codec](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/triple_codec/NOTE.md) and [strands](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strands/NOTE.md): supply the divisor-balance and Fibonacci strands. Here the Cassini product turns two separate fibers into a selected all-size output family.
- [Divisor balance](https://github.com/eliottcassidy2000/math/blob/f18780daea94f262eb1b7d6fd174015e2657ea86/05-knowledge/results/arithmetic_braids_20260917_divisors.md) and [ternary bridge](https://github.com/eliottcassidy2000/math/blob/f18780daea94f262eb1b7d6fd174015e2657ea86/05-knowledge/results/ternary_bridge_20260925.md): the precise definitions and guard boundary are inherited, not inferred from matching small numbers.
- [Amnon Yekutieli, *Pythagorean Triples, Complex Numbers, Abelian Groups and Prime Numbers*](https://arxiv.org/abs/2101.12166): prior rational-circle group structure and Gaussian-prime organization. The signed-valuation arithmetic is classical in substance.
- [Buchstaber and Dragovic, *Two-Valued Groups, Kummer Varieties, and Integrable Billiards*](https://amj.math.stonybrook.edu/pdf-Springer-final/018-0085.pdf): prior two-valued quotient framework; their elementary sum/absolute-difference example is the rank-one version of the sign quotient used here. No claim is made that the endpoint tournament itself is a Kummer variety.

The research lesson is operational: compare **quotient-then-compose** with **compose-then-quotient**. In this family their discrepancy is exactly measurable. A successful reusable source interface must retain every choice that distinguishes future accepted outputs, not merely enough information to decode the current count.

## 8. Reproduction

Run `python3 check.py > check.json` and `python3 -O check.py > check_optimized.json`, then compare the two JSON outputs. Both modes currently pass 88,873 explicit checks with identical output. All checks use explicit exceptions, not removable assertions.

Coverage: literal permutation enumeration through E(6) (eight vertices), independent subset DP through E(8) (ten vertices); all 900 address pairs for r,s in 1..4 under actual terminal gluing; 15,625 signed-valuation pairs over primes 5,13,17 with exponents -2..2; all 19,683 triples with exponents -1..1 for the carry identity and two-valued associativity; direct Euclid enumeration of all resulting fixed-norm fibers in five overlap examples, including total cancellation; and 34 Cassini families n=2,5,...,101 with positive labels and all square-sum identities checked exactly. The all-size statements rest on the proofs, not the finite ranges.
