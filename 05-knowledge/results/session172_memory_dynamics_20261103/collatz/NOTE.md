# Collatz itineraries as contextual tournament memory

**Status: PROVED finite realization and observer obstruction; exact descent certificates.** This packet connects S171's Hamilton read heads to Collatz parity dynamics. The arithmetic facts about parity cylinders and density-one finite stopping are classical (Terras/Everett), not claimed new. The contribution here is their explicit tournament realization, a rigorous memory/continuation contract, and a certificate interface distinguishing itinerary information from height descent. No Collatz convergence theorem is claimed.

## 1. The arithmetic state and its finite response

Use the normalized map on positive integers

    T(n)=n/2 if n is even; T(n)=(3n+1)/2 if n is odd.

Write b_j(n)=T^j(n) mod 2. For a word w=(b_0,...,b_(k-1)), put s_j=sum_(i<j)b_i and define

    a_0=0,  a_(j+1)=3^(b_j)a_j+b_j 2^j.

Equivalently,

    a_k=sum_(i=0)^(k-1) b_i 2^i 3^(s_k-s_(i+1)).

**Theorem 1 (finite itinerary interface).** For every k>=1:

1. Every binary k-word occurs as the first k parities of exactly one residue class r(w) modulo 2^k.
2. On that class, T^k(n)=(3^(s_k)n+a_k)/2^k.
3. Its residue is r(w)=-a_k (3^(s_k))^(-1) mod 2^k.
4. Two positive integers have identical first k parities iff they are congruent modulo 2^k.

**Proof.** The affine identity follows by induction from the branch formula. Suppose the first j parities are fixed. Replacing n by n+2^j preserves these parities and changes T^j(n) by the odd integer 3^(s_j). Thus the two lifts of a given class modulo 2^j give opposite next parities. Starting with the sole class modulo 1 proves the bijection by induction. Divisibility of the affine numerator gives (3), since 3 is invertible modulo 2^k. This also proves (4). All classes have positive representatives. ∎

In particular these are **all** binary words, not a restricted Reed–Muller family. A code restriction cannot be imposed on raw Collatz parities without losing integers. Error correction must be added as a separate injective encoding.

## 2. A literal constant-profile Collatz tape

Use S171's five-vertex strong seeds, specified by outgoing-neighbor masks,

    B_0=(16,1,3,7,14), B_1=(20,1,2,7,14).

They differ only at arc {0,2}, are nonisomorphic, and have the same labelled-path-slot polynomial

    F(u)=u^5+6u^3+2u,  H(B_b)=9.

Adding a vertex that beats only old vertex 0 gives Hamilton counts 17 and 19 respectively. These inherited finite identities are checked independently by the accompanying program.

For any binary word w define the ordinal-sum tournament

    G_w=B_(b_0) ⊕ ... ⊕ B_(b_(k-1)).

All earlier blocks beat later blocks. To read block j add one vertex within that block, using the seed attachment, beaten by all earlier blocks and beating every later block.

**Theorem 2 (finite Collatz tape).** The 2^k residue classes modulo 2^k have representations by 2^k pairwise nonisomorphic tournaments G_w on 5k vertices such that:

- Every complete labelled-slot profile is F(u)^k and every Hamilton count is 9^k.
- Labelled arc distance between two tapes equals the Hamming distance between their parity words.
- The addressed one-vertex context has Hamilton count (17+2b_j)9^(k-1).
- The k local read heads recover r(w), the exact affine map T^k on the entire residue cylinder, and every prefix's descent threshold.
- Any adaptive scheme using only these binary local read heads requires k queries in the worst case to identify an arbitrary residue class.

**Proof.** Path-slot covers multiply under ordinal sums, because each slot visits blocks in their unique order. The seeds are strong, so the ordered strong-component types recover w up to graph isomorphism. Exactly one labelled arc changes at each differing bit. The read count follows from replacing one factor 9 by 17 or 19. Theorem 1 then gives the arithmetic reconstruction. The binary decision-tree lower bound is 2^q>=2^k. ∎

This is an exact finite realization of Collatz **observations**, not an autonomous Collatz graph evolution. The original bulk profile forgets every one of the k residue bits. The controlled contexts recover them from the graph, not from the old profile.

### Robust and rich-output variants

An injective binary encoder E:{0,1}^k→{0,1}^L of minimum distance delta gives, by the same construction, 5L-vertex tapes of minimum labelled arc distance delta. Unique decoding is guaranteed for 2e+s<delta with e arc errors and s erased pairs. This encodes, rather than restricts, the Collatz words. No claim is made that arc errors cause local output-bit errors.

For a contrasting bandwidth tradeoff, repeat block b_j exactly 2^j times. There are N=2^k-1 blocks. Add a read vertex in every block in one experiment. Its single integer response is

    17^(N-z)19^z,  z=sum_j 2^j b_j.

The distinct prime factors recover z and hence the entire itinerary. This uses N additional vertices and exponentially many storage blocks. It is not a counterexample to the k-query lower bound for the one-vertex local interface.

## 3. Why the tape cannot evolve autonomously at fixed capacity

**Theorem 3 (no finite safe parity quotient).** Suppose X is a set, q:N_+→X, A:X→X, and c:X→{0,1} satisfy

    q(T(n))=A(q(n)),  c(q(n))=n mod 2

for every positive n. Then q is injective. In particular X cannot be finite.

**Proof.** Equal q-states have equal observations after every iterate of A, hence identical parity prefixes of every length. Theorem 1 implies their difference is divisible by 2^k for every k, so they are equal. ∎

This is the nonlinear infinite-state version of S171's safe-forgetting principle: the present parity observer has two fibers; closure under every allowed future continuation separates every positive integer.

The general b-adic contextual observer packet developed alongside this note extends this conclusion to affine full-branch systems. Here the following finite witness identifies precisely what a tournament tape lacks.

For a sharper finite witness, n and n+2^k have the same k-bit tape, but their images have different k-bit tapes. Their original parities first differ at index k, since after k common branches their iterates differ by the odd number 3^(s_k). Thus no deterministic update on the 2^k tape states can implement n↦T(n) while retaining a fresh length-k prediction. One must supply the missing bit, retain more of n, or reduce the remaining prediction horizon.

A graded interface does work: discard the first bit of a k-tape to obtain the (k-1)-tape of T(n). This forgets one bit and consumes one time step. A fixed-length shift register with an externally supplied new parity bit is also possible, but the external stream carries precisely the missing information. Neither is a finite autonomous proof of convergence.

The obstruction is interface-specific. It does not rule out a finite proof, a symbolic rewrite system, or a finite collection of rules acting on unbounded integers.

## 4. Exact descent certificates, not just readable histories

For a prefix of length j write D_j=2^j-3^(s_j). The affine identity gives

    T^j(n)<n  iff  D_j n>a_j.

Since a_j>=0, a prefix with D_j<=0 never certifies strict descent. If D_j>0, it certifies descent exactly when

    n>floor(a_j/D_j).

The pair (s_j,a_j), with its length j, is therefore a sufficient arithmetic interface for that endpoint's height comparison. The parity weight s_j alone determines the multiplicative drift but omits the additive cocycle a_j. For example words (0,1) and (1,0) have the same weight 1 but a_2=2 and 1 respectively. A simultaneous read of all j unweighted blocks returns 17^(j-s_j)19^(s_j), so it reads the slope and forgets this order-sensitive additive term.

**Proposition (the additive cocycle retains the missing order).** At fixed length k, the exact affine maps associated with distinct parity words are all different. At fixed weight s, their additive terms a_k are all different, giving exactly binom(k,s) distinct offsets. Consequently, under uniform residues modulo 2^k, the simultaneous unweighted Hamilton read leaves conditional entropy log_2 binom(k,s) after response s, and average hidden entropy k-H(Binomial(k,1/2)).

**Proof.** If two words have the same s and a, Theorem 1 gives the same residue and hence the same word. Equality of their affine maps first forces the same slope, thus the same s, and then the same a. The uniform source is uniform on the binom(k,s) words of given weight. ∎

This is a sharper loss statement than saying that a slope statistic may miss an offset: within every fixed-weight class, **all** the unresolved itinerary information is carried by the additive cocycle. An exact affine branch interface is a faithful record of the full k-word, even though the formula uses only two integers. Counting integer registers is not counting bits.

For an entire k-cylinder w, define

    h(w)=min_{1<=j<=k, D_j>0} floor(a_j/D_j),

with h(w)=infinity when the index set is empty. Then, for every positive integer n in that cylinder,

    some T^j(n)<n with 1<=j<=k  iff  n>h(w).

This is a fully exact threshold, including the finite exceptional small representatives. The least positive representative is r_+=r(w) when r(w)>0 and r_+=2^k otherwise. If h is finite, the number of cylinder representatives failing to descend by time k is

    max(0, floor((h-r_+)/2^k)+1).

These formulas provide a machine-checkable cover of all large representatives in each contracting cylinder and explicitly list its exceptions. They do not silently convert coefficient contraction into pointwise descent.

**Corollary (what would suffice for convergence).** If every n>1 has some finite prefix with D_j n>a_j, then every positive orbit reaches 1.

**Proof.** Strong induction on n: descend to a smaller positive integer and apply its induction hypothesis. The base n=1 is already at 1. ∎

Thus the missing theorem is universal descent, not universal observability of finite histories.

## 5. Exact density certificates and their quantifier boundary

Let S_k be the binary k-words satisfying 3^(s_j)>2^j for every 1<=j<=k. Equality is impossible for j>=1. By the preceding section, the positive integers having no strict descent in their first k iterates differ by only a finite set from the union of the |S_k| residue classes corresponding to S_k. Their natural density is therefore exactly

    d_k=|S_k|/2^k.

The survivors are counted without orbit simulation by the dynamic program

    v_(0,0)=1,
    v_(j,s)=1_[3^s>2^j] (v_(j-1,s)+v_(j-1,s-1)),
    |S_k|=sum_s v_(k,s),

with out-of-range entries zero. Every condition is an integer comparison, not a floating-point logarithm decision.

Since survival implies s_k>k alpha, where alpha=log(2)/log(3)>1/2,

    d_k <= Pr[Binomial(k,1/2)>=alpha k]
        <= exp(-k D(alpha || 1/2)),

where D(p||1/2)=p log(2p)+(1-p)log(2(1-p)). For completeness the last inequality follows from Markov's inequality applied to exp(t S), giving exp(-t alpha k)((1+exp(t))/2)^k, minimized at exp(t)=alpha/(1-alpha). Consequently d_k→0 and the integers with finite strict stopping time have natural density one: the exceptional infinite-stopping set is contained in each finite-horizon failure set, so its upper density is at most every d_k. This is a classical stopping-time conclusion, included here as a proved control and not as a new Collatz advance.

The coefficient-first-contraction words form a prefix-free set. Under uniform binary-prefix measure their total cylinder mass is 1, because the surviving mass d_k tends to zero. Thus a sequential Hamilton reader, supplied with progressively longer tapes or fresh stored bits, looking for coefficient contraction terminates almost surely under that measure, but its worst-case depth and required tape capacity are unbounded. This statement is about binary cylinder probability, not a uniform probability distribution on positive integers.

**Theorem 4 (no bounded descent horizon).** For every K>=1 there are infinitely many positive integers with T^j(n)>n for all 1<=j<=K.

**Proof.** Take n=q2^K-1 with q>=1, choosing n>1. Its first K parities are all one, and

    T^j(n)=3^j q2^(K-j)-1 > q2^K-1=n.

∎

In particular no finite collection of bounded-length parity cylinders, each justified only by descent at an inspected prefix, can cover all n>1. A successful approach of that type needs unbounded depth or an additional induction/composition argument. Density one, even with exponentially small finite-horizon failure density, does not eliminate a sparse exceptional set. Nor does reading every available bit prove that the next unread bit must be favorable.

## 6. Concrete next research target

Use the exact cylinder threshold h(w) to grow a prefix-free descent certificate tree on the domain n>1, with n=1 treated separately as the induction base. A node may close only when all representatives in this domain lie above its verified threshold, or when its finitely many exceptions in this domain have separately verified descent. The representative 1 is not required to have strict descent. Store unresolved nodes and exact affine cocycles, not merely parity weights. The natural mathematical target is a proven rule that maps every unresolved node to previously certified smaller integers or to a well-founded symbolic rank. Such a rule would add a genuinely new convergence mechanism. The tournament tape supplies authenticated finite histories and a precise capacity limit; it supplies no such ranking theorem by itself.

## Sources, reproduction, and provenance

- S171 storage/read-head packet: [S171 information](../../session171_knots_information_20261103/information/NOTE.md), routed by the [published S171 index](https://github.com/eliottcassidy2000/math/blob/1ad21a4dcf30ff54428b3e8290db5b00483dbe5b/00-navigation/S171.md).
- S171 contextual framework and synthesis: [S171 framework](../../session171_knots_information_20261103/framework/NOTE.md) and [S171 synthesis](../../session171_knots_information_20261103/synthesis/NOTE.md).
- Riho Terras, [A stopping time problem on the positive integers](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/30/3/101028/a-stopping-time-problem-on-the-positive-integers), Acta Arithmetica 30 (1976), 241–252, DOI 10.4064/aa-30-3-241-252. Original research source for the historical stopping-time attribution.
- Jeffrey Lagarias, [Behavior of the stopping time function](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/node4.html), author-hosted exposition of the parity bijection and Terras's coefficient/actual stopping-time distinction. Consulted to check provenance; all needed identities are proved above.

Run `python3 check.py > certificate.json`, repeat with `python3 -O check.py > certificate.optimized.json`, and compare. The checker validates seed profiles and 17/19 readouts by enumeration; all parity cylinders through depth 12, their affine formulas and exact thresholds on multiple positive representatives; literal two-block tournament contexts; nonautonomous tape witnesses; exact survivor counts; and rich-output injectivity. It imports no repository code and uses explicit exceptions rather than assert statements.

Validation snapshot: 8,190 cylinders and 49,140 positive representatives passed. The exact depth-256 survivor count is 37933744616918826336506027646167154052909131283110018145640477926649259 out of 2^256; its density is approximately 3.2760221243761753e-7.

SHA-256:

    check.py: 4f367bbe7c9fe66c96ef786bfc9821affa4c22bc52982cf2d9e01c7de0adf212
    certificate.json: b621a57cebc76ad69f2fcd3f5b7cfb85f4681716951f7176a2616556cdf4da5c
