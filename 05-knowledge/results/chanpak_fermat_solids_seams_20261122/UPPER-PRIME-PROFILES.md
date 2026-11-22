# The two prime fibers of the exponent-61 upper target

**PROVED and FINITE-EXACT.** The upper integer has the certified prime factorization
\[
  E_{61}+1=2,658,455,991,569,831,744,654,692,615,953,842,177
          = 2,432,582,681\cdot1,092,853,292,237,112,554,142,488,617.
\]
Write these factors as \(s<q\). A concrete strong tournament of order 34 has count \(s\). The remaining scalar question for this factorization is exactly whether \(q\) is realizable. If it is, this joins two concrete prime fibers to realize \(E_{61}+1\); one of the two cores can always be the specified 34-vertex witness. The strong upper realization required by the physical seam remains a separate question. [Strong lower rail](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md) · [Upper rigidity in every order](README.md)

## The small factor has a constructive address

Take the inherited one-apex/one-apex join \(F_1(b,e)\) at \((b,e)=(49,41466580)\). Its height is
\[
             g(b-1)g(e-1)+be-[b\text{ odd}][e\text{ odd}]
             =47\cdot8525963+49\cdot41466580=s.
\]
Here \(g(0)=g(1)=1\), \(g(2u+1)=g(u)\), and \(g(2u)=g(u)+u+[u\text{ odd}]\) for \(u>0\). Explicitly, the rooted core for address \(e\) uses \(v=e-1\) and \(k=\operatorname{bitlength}(v)\). The spine has arcs \(j\to i\) for \(j>i\); the root beats precisely the spine positions where bit \(i\) of \(v\) is one. Use the converse of the \(b\)-core on the left and the \(e\)-core on the right. Direct all arcs from right to left except for the root-to-root arc from left to right. The replay checks strong reachability in both directions in these literal rows and evaluates the compressed path count. [Join constructor and path formula](../chanpak_prime_frontier_selector_20261118/join_loads.py)

## Rigidity and orders of the two cores

The sizes of any strong realizers belong to the following windows. The minima and maxima in these windows are bounds.

| Prime fiber | Possible strong order | Certified realization |
|---|---:|---|
| \(s\) | \(16\) through \(41\) | \(34\) vertices |
| \(q\) | \(32\) through \(117\) | Open here |

For the lower bounds, if a permutation forms a Hamilton path, then each consecutive block of positions forms one as well. For block widths \(k_i\) summing to \(n\), assign the vertex sets to the blocks and bound the internal orders by their census maxima \(M(k_i)\). Thus \(M(n)\le n!\prod_i M(k_i)/k_i!\). The partition \(15=10+5\) gives \(M(15)\le709233525<s\); the partition \(31=10+10+10+1\) gives \(M(31)\le671677290585888037883932500<q\). The maximum is nondecreasing by ordinal singleton padding. For the upper bounds Busch's strong tournament inequality gives \(5^{n-1}\le H^3\). Direct integer comparisons give 41 and 117. [Audited maxima through order ten](../first_creator_generation_20261108/audit/AUDIT.md) · [Busch's lower bound](https://doi.org/10.37236/1141)

In a tournament of prime count, deleting the singleton strong components leaves precisely one nontrivial strong component, since the counts of components multiply. Its automorphism group acts freely on the Hamilton paths. A nontrivial group would therefore have order equal to the count and contain a cycle of that prime length, impossible within the displayed windows. Ordering of strong components fixes each padding vertex. At the full product, every realization with two nontrivial strong components has one of each prime type in either order. This is equivalent to the existence of a tournament at the larger prime because the explicit smaller witness is available. [LEM-003, free action on paths](../../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md)

## Primality and the remaining sparse test

For each large integer \(N\) in the certificate, the replay factors \(N-1\) completely and provides \(a\) with \(a^{N-1}=1\bmod N\) and \(\gcd(a^{(N-1)/r}-1,N)=1\) for every prime \(r\mid N-1\). For any prime divisor of \(N\), the order of \(a\) modulo that divisor must then be \(N-1\), forcing the divisor to equal \(N\). The only intermediate certification above eleven million uses \(N=1468888833652033002879689\). The witnesses for \(s\), this intermediate, and \(q\) are 3, 3 and 5 respectively. All leaves are proven by trial division through the integer square root. The script includes the factorization trees and uses integer arithmetic throughout.

**FINITE-EXACT.** At the larger prime, no representation with transitive deletion of one apex has at most eight incident wins or at most eight losses on that apex, without an order bound. After stripping redundant leading zeros and trailing ones, let the top mark be at \(x\). Its contribution to the last zero alone gives \(H\ge1+2^{x-1}\) for \(x\ge2\). All marks together give \(H<2^{x+1}\) by summing their individual powers of two. Thus \(2^{89}<q<2^{90}\) leaves just \(x=89,90\). For each top and budget the previously derived [monotone sparse inverse](README.md#sparse-asymmetric-contacts-at-the-absolute-exponent-61-target) exhausts every possible descending selection and tests the exact smallest mark. All these searches are empty. Reversal and complementation treat a small number of losses as well.

From this directory run `python3 -B check_prime_factors.py` for the factorization, certificates, order bounds and explicit small core. Run `c++ -O3 -std=c++17 sparse_upper.cpp -o /tmp/chanpak_sparse_prime && /tmp/chanpak_sparse_prime --large-factor` for the large prime's sparse test.
