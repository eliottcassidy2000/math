# Exact finite congruence states for whole tournament modules

Session197. All input tournaments are nonempty. A context is a finite tournament with one distinguished vertex, replaced by the entire input as an intact module. The results concern Hamilton-path counts modulo a fixed integer, not exact integer spectra or prime infinitude.

## 1. A parity-free universal quotient

Let c_r(T) be the unordered r-path-cover count and B_T(x)=sum_r r! c_r(T) binom(x,r). Define

    w_r(T)=(r-1)! r! c_r(T), r>=1.

For M>=2 put R=min{r>=1: M divides (r-1)!r!}. This is finite. Because these factorial products form a divisibility chain, w_r=0 modulo M for every r>=R, whatever the input order.

**Theorem.** Two actual nonempty tournaments have the same Hamilton count modulo M in every whole-module context if and only if their vectors (w_1,...,w_(R-1)) agree modulo M. Every fixed context induces a computable integer linear map on these finite vectors. Consequently every finite context alphabet has an exact finite modular automaton, valid under every continuation.

Here and below a vector space of formal residues is only an upper bound on the set of physically realized states.

Proof. Cut an output Hamilton word at its maximal input runs. Its r input runs form an ordered r-path cover. The exterior consists of two optional outer paths and r-1 nonempty internal paths, each starting at an out-neighbor of the distinguished vertex and ending at an in-neighbor. Permuting those distinct internal paths acts freely. Thus every context Q has

    H(Q[T])=sum_r b_(Q,r) w_r(T), b_(Q,r) nonnegative integers.     (1)

Only finitely many coefficients are nonzero. This proves sufficiency and the finite cutoff.

For necessity use the physical matching contexts Q_j. They have root q, out-neighbors O_1,...,O_(j-1), in-neighbors I_1,...,I_(j-1), transitive orientations within each side, matching arcs O_i->I_i, and all other cross arcs I_l->O_i. At most j input runs are possible; equality forces the j-1 matching edges as the internal paths. Therefore b_(Q_j,r)=0 for r>j and b_(Q_j,j)=1. The first R-1 responses form a lower unit-triangular integer transform of w_1,...,w_(R-1), so they recover those coordinates modulo M.

For closure, fix any context C. Integer triangular inversion gives each output coordinate w_i(C[T]) as an integer linear combination of H(Q_j[C[T]]), j<=i. Each composed graph is still one whole-module context of T. Applying (1) to each composition proves

    w_i(C[T])=sum_r A_(C,i,r) w_r(T), A_(C,i,r) integers.           (2)

Each row has finite support. Reducing and discarding r>=R is exact, even when C changes parity. Coefficients can be computed by finite enumeration of the exterior arrangements in the composed contexts, followed by integer triangular inversion. The transition can have negative integer coefficients; only the original scalar coefficients in (1) are asserted nonnegative.

There are at most

    product_(1<=r<R) M/gcd(M,(r-1)!r!)

possible physical vectors. Breadth-first search from any specified physical seed, with transitions (2), gives exactly all residues reachable by words in a finite alphabet, and each recorded word constructs an actual tournament. If the empty word is allowed, a shortest witness has length at most the number of reachable states minus one. To require a positive word, include a separate initial flag; a bound equal to the unflagged number of reachable states suffices. This is a decision procedure for modular targets, not integer targets. Strongness needs an additional hypothesis on the letters; the two letters in Section4 always produce strong graphs.

The zero vector means M divides the universal content G(T)=gcd_r w_r(T). It is absorbing under every context. No claim that arbitrary formal nonzero vectors are physical is needed.

## 2. Smaller parity-tagged coordinates

For odd/even input order respectively, use

    q_k=x product_(j=1)^k(x^2-j^2),
    r_k=x^2 product_(j=1)^k(x^2-j^2).

Write B_T=sum_k u_k q_k or sum_k u_k r_k. These monic integer bases give integer u_k. Put

    gamma_k=(2k+1)!                  (odd),
    gamma_k=(k+1)(2k+1)!             (even),
    p_k=k! gamma_k, a_k=p_k u_k.

The reviewed universal-content theorem proves an integer unit-triangular transformation between the first d weighted Newton coordinates and (a_0,...,a_(d-1)), with all later weighted Newton coordinates integer combinations. This also follows directly by Vandermonde expansion of

    q_k/gamma_k=binom(x+k,2k+1),
    r_k/gamma_k=binom(x+k+1,2k+2)+binom(x+k,2k+2).

Set K_epsilon(M)=min{k>=0:M divides p_k}. Since p_(k+1)/p_k is integral, all a_k for k>=K vanish modulo M. Within one parity, equality of the retained a coordinates is equivalent to equality in every whole-module context. One may recover a_k instead from the first k+1 repeated-triangle responses: that response matrix is lower unit triangular on a.

The same composition proof as in Section1 gives an integer matrix on weighted central coordinates for every context, including parity-changing contexts. To compute a column, expand q_k or r_k as an integer combination of positive monomials x^s=B_(Tr_s); compute the output profiles on those actual transitive inputs; form the signed sum; multiply output central coefficient j by its p'_j and divide by input p_k. Division is exact over the integers by the preceding theorem. Never divide a nonunit modulo M.

A parity tag gives a correct finite machine with at most

    sum_epsilon product_(0<=k<K_epsilon) M/gcd(M,p_k)

states, but may distinguish equivalent states of opposite parity. In particular the two zero states are equivalent. We claim minimality only within a fixed parity, except for the exact global counts in Section4.

For an odd prime p, K_odd(p)=K_even(p)=(p-1)/2. At p^2 these values are p and p-1 respectively: before k=p-1 the factorial (2k+1)! has at most one p-factor; at k=p-1 the even extra factor k+1 contributes the second, whereas the odd basis needs k=p. More generally K_epsilon(p^e)<=p ceil(e/3), since at k=pc, k! and (2k+1)! together contain at least 3c factors p. These are sufficient finite bounds, not optimality claims for every exponent.

## 3. Uniform triangle recurrences independent of input order

Let D denote the physical wrapper W(T)=C3[T,K1,K1]. Its exact polynomial action is

    Df=x^2 f+x(f(x+1)-f(x-1)).

The weighted-coordinate matrix has diagonal alpha_k=k^2+6k+3 for odd input and k^2+6k+5 for even input, upper diagonal1, and lower diagonal

    beta_k=2k^2(2k+1)          (odd),
    beta_k=2k(k+1)(2k+1)      (even).

Truncation to K_epsilon(M) is exact modulo M. Thus Cayley-Hamilton gives a uniform order-K recurrence for h_j=H(W^j(T)) modulo M, for every input order. No restriction on its initial profile or Hankel signs is used. In the unweighted path interpretation, any term reaching height K before returning to readout0 traverses every lowering edge1,...,K, whose product is p_K; that is another proof of the cutoff.

For M=7 the odd matrix is

    [3 1 0; 6 10 1; 0 40 19],

with characteristic x^3-32x^2+231x-336, congruent to x^2(x-4). Hence

    h_(j+3)=4h_(j+2) (mod7), j>=0.

The even matrix is [5 1 0;12 12 1;0 60 21], with characteristic x^3-38x^2+345x-708, congruent to (x+1)(x^2-4x+6). The quadratic has nonsquare discriminant6 in F7 and its roots have norm6=-1, so their eighth powers are -1. Therefore this characteristic polynomial divides x^16-1; the even matrix has sixteenth power I. Every even-input sequence has period dividing16 from j=0.

For T=Tr2 the sixteen residues are

    1,5,2,4,6,5,0,3,0,1,6,2,2,1,1,3.

Thus 7 divides H(W^j(Tr2)) exactly when j is6 or8 modulo16. The first two such values are composite:

    h6=2399173=7*89*3851,
    h8=2296349461=7*241*367*3709.

Divisibility alone does not establish a 7p witness or an infinitude statement.

## 4. Exact physical quotients modulo5 and7, and short steering

Use the second strong wrapper E(T)=C3[T,Tr2,K1], whose full bin operator is

    Ef=x[(x^2-2)f+(1-2x)f(x-1)+(1+2x)f(x+1)].

D preserves parity and adds2 vertices; E reverses parity and adds3. Both are strong for every nonempty input, because their three nonempty blocks form a directed cycle.

At modulus5 the weighted central matrices, with columns as input coordinates, are

    D_odd=[3 1;1 0], D_even=[0 1;2 2],
    E_odd->even=[0 3;2 2], E_even->odd=[1 3;3 3].

At modulus7 they are

    D_odd=[3 1 0;6 3 1;0 5 5],
    D_even=[5 1 0;5 5 1;0 4 0],
    E_odd->even=[5 3 0;5 4 5;0 1 3],
    E_even->odd=[4 3 0;1 2 5;2 6 2].

The reproducible complete finite checks give:

- Mod5, the orbit of K1 is all48 nonzero tagged states; maximum shortest-word length7. Both letters are bijections on these states. Every nonzero state reaches every scalar residue using a positive word of at most4 letters and at most10 added vertices, simultaneously.
- Mod7, the orbit of K1 is all686 tagged states; maximum shortest-word length14. Every one of the684 nonzero states can return to the K1 state, while zero is absorbing. Every nonzero state reaches every scalar residue using a positive word of at most5 letters and at most13 added vertices, simultaneously.

These are finite exact certificates using the displayed matrices and exhaustive graph search, not extrapolations to other primes. Any input with p not dividing G(T), p=5 or7, can therefore be turned into a strong graph with any prescribed scalar residue by the stated bounded number of physical wrappers. If p divides G(T), every context count is0 modulo p. Mod5 this particular D/E grammar cannot take a nonzero behavior class to zero, although it can make the scalar H divisible by5; mod7 it can. This distinguishes scalar divisibility from universal divisibility.

There are exactly49 and685 global whole-module behavior classes modulo5 and7 respectively, including inputs of both parities and all orders. To prove the upper bound for any odd prime p, note that w_1,...,w_(p-1) are equivalent to B_T as a polynomial function on F_p: factorials below p are units, B_T(0)=0, and the Newton evaluation matrix on0,...,p-1 is unit triangular. Odd input gives an odd function and even input an even function, each of dimension(p-1)/2 with value0 at0. Their intersection is the zero function because p is odd. Thus at most

    2 p^((p-1)/2)-1

global classes are possible. Within either parity the central coordinate description is injective, so the finite orbits above realize every nonzero class. Mod7 they realize zero too. Mod5 use the cyclic regular five-vertex tournament: its order5 rotation extends to every whole-module context and acts freely on Hamilton words, proving universal divisibility by5. Hence the upper bound is attained in both cases. Every class has a strong representative: nonempty D/E words realize the nonzero classes (a nonempty return word replaces the empty K1 presentation), and the cyclic zero representative is strong.

## 5. Verification, inheritance, and scope

Run python3 verify.py and python3 -O verify.py; outputs are deterministic and do not edit this packet. The checker derives all matrices by exact polynomial operations, checks actual profiles for small arbitrary contexts, verifies finite weighted-coordinate closure, regenerates both complete state graphs and bounded steering words, and checks the uniform recurrences. It includes composite-modulus and parity-changing controls, rather than relying only on prime examples.

Inherited ingredients are explicitly separated: the free internal-fragment factorial and matching-bank triangularity are from UNIVERSAL_MODULE_CONTENT.md, SHA256 4f86f7ada5eff79e61bb2ef4b8139b76df8efedf6996a71dae327bc66b24d9cb; full bin linearity and the D/E operators are from UNIVERSAL_CONTEXTS.md, SHA256 8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51, and PROOF.md, SHA256 b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b, in the path-polynomial packet. The central triangle-response bank is from our reviewed triangle-content/PROOF.md, SHA256 242a28612c17288da56ea785341e5189f9b03eea3e70fec16b4110d7fc88fcde, published at https://github.com/eliottcassidy2000/math/blob/0bd8d9c49c7ee5fca8650d3aee966942c6f5f34d/session197/triangle-content/PROOF.md . The unchanged helper core.py is that packet's verifier, SHA256 bbd76a9897f30d10bd846e39468ab644098d19d7ed9f36c90241f6bcb357fbcc.

The present additions are exact contextual congruence, composition-closed finite machines, uniform modular truncation, and the two certified physical quotient/steering examples. We make no literature-priority claim and no unrestricted prime/7p coverage claim. Whole-module substitutions are essential: arbitrary nonuniform births, deletion, and port-dependent operations are outside this quotient theorem.
