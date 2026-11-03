# Two-bin path-cover multiplication and a 7-divisible strong family

Status: proved exact constructions, with finite controls. The construction gives infinitely many strong counts divisible by 7 and particular exact counts 7p with prime p. It does not prove infinitely many primes of the displayed shape, nor does it characterize the spectrum. No claim of novelty is made.

## The multiplicative state

Let h(T) be the number of directed Hamilton paths of a nonempty tournament T. Define f(T) to be the number of ordered pairs of vertex-disjoint directed paths covering T, allowing an empty path in either slot (but not both, since T is nonempty). Then

    f(T) = C_2(T) + 2 h(T).

The two extra cases place a Hamilton path in one slot and the empty path in the other.

For nonempty tournaments X and Y let X⊕Y be their ordinal sum, with every vertex of X pointing to every vertex of Y. Restricting a path in X⊕Y to X and Y shows that it consists of an X-path followed by a Y-path, either part possibly empty. Conversely any such concatenation is a path. Applying this separately in two labeled slots gives bijections

    h(X⊕Y) = h(X)h(Y),       f(X⊕Y) = f(X)f(Y).

For the transitive tournament T_n, h(T_n)=1 and f(T_n)=2^n: assign each labeled vertex independently to either slot, where the internal order is forced. The same identities may be applied with n=0 by defining h(T_0)=f(T_0)=1.

Now let G(B) be the directed 3-cycle substitution with B in one block and a singleton in each other block. By the cyclic path-run bijection in NOTE.md,

    h(G(B)) = 3 h(B) + C_2(B) = h(B) + f(B).

Every G(B) is strong, since all its blocks are nonempty. Combining the formulas, for every nonempty X and integer n>=0,

    h(G(X⊕T_n)) = h(X) + 2^n f(X).                 (1)

Thus (h,f) is a two-coordinate exact state under ordinal sum, and adjoining a transitive block is an exact scaling operation on f alone. This is stronger than knowledge of h alone.

## An explicit 7-divisible family

Let X have vertices 0,1,2,3, every arc i→j for i<j, except 3→0. The five Hamilton paths are

    0123, 1230, 1302, 2301, 3012.

Thus h(X)=5. For f(X), 2 assignments put all vertices in one slot, giving 2·5; the six splits with two vertices in each slot contribute 6; and the four induced triples have Hamilton-path counts 1,3,3,1, contributing 2·8 across the (1,3) and (3,1) splits. Hence f(X)=10+6+16=32. Equivalently C_2(X)=22.

By (1), each n>=0 produces a strong tournament on n+6 vertices with exactly

    5+2^(n+5)

Hamilton paths. At n=0 this is 37, a prime. Since the order of 2 modulo 7 is 3, the count is divisible by 7 exactly when n≡2 mod 3. Setting n=3t+2 gives the infinite family of distinct strong tournament counts

    7 p_t,       p_t = (2^(3t+7)+5)/7,       t>=0.

Two prime-cofactor instances are

    t=0: n=2,  h=133=7·19   (8 vertices);
    t=2: n=8,  h=8197=7·1171  (14 vertices).

Here 19 and 1171 are prime by trial division through their square roots. The intervening cofactor p_1=147 is composite, illustrating why the parameter does not establish an infinite prime family. More precisely, p_t≡5+2t mod 7, so 7 divides p_t exactly when t≡1 mod 7; this follows from 2^(3t+7)=128·8^t and the binomial theorem modulo 49.

For a state-preserving but generally non-strong enlargement, ordinal-summing any transitive T_m before or after an existing witness leaves h unchanged. In particular 133 has witnesses of every order at least 8, with no claim that these enlarged witnesses are strong.

## Reproduction

`python3 extension_check.py` and `python3 -O extension_check.py` independently enumerate X, verify the two-slot multiplicative law on small tournaments, count literal G(X⊕T_n) by subset/endpoint DP for n=0,1,2,4,8, and check both prime cofactors by trial division. The argument above proves the infinite family; computation is a control.
