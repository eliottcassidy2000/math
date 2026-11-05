# Palindromic three-port clocks and every Mersenne target

Session 197. **Constructive all-order family, independently reviewed.** Every integer 2^k-1 with k>=4 is the Hamilton-path count of an explicit strong feedback-one tournament on exactly k+2 vertices. This order is minimal within the strong feedback-one class. The construction is a repeated physical left/right word with three fixed ports, not transitive cloning of one vertex.

## 1. Inherited physical integer gadget

ROOTED.md Section 4, SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e, constructs an actual gadget from L>=1. With n=L-1, s=bit_length(n), its transitive base x_0,...,x_(s-1) has x_i->b exactly when bit i of n is one. Put f(n)=g(n+1). The physical recurrence is

    f(0)=1,
    f(2x+1)=f(x),
    f(2x)=f(x)+x+(x mod2).                         (1)

For L>=2 add a first base vertex a, with a->x_i for every i and b->a. The resulting T(L) is the singleton-prefix one-seam output. It has

    H(T(L))=g(L)+L-(L mod2)=n+f(n)+(n mod2).        (2)

Deleting b leaves the transitive base a,x_0,...,x_(s-1). Its unique Hamilton cycle is that increasing base order followed by b and back to a. Thus T(L) is strong, has exactly one Hamilton cycle and feedback vertex number one. The graph has s+2 vertices. The L=1 transitive edge is excluded throughout.

For a binary word c of length ell, let rev(c) denote its reversed word as an integer. When its first and last bits are one, both c and rev(c) are odd and have length ell. Its rooted starting-two-cover response satisfies Q(c+1)=2^ell-rev(c), by the direct reversed-complement formula; this also follows from ROOTED.md's response recurrence.

## 2. An exact two-block formula

Fix any positive odd integer c, with ell=bit_length(c), and integers p>=1, q>=0. For t>=0 form the canonical binary word

    w_t = 1 0^(p+t) bin(c) 0^(q+t),
    n_t = value(w_t),  L_t=n_t+1.                 (3)

The word has length p+q+ell+1+2t. Define

    A=2^(p+q+ell+1),
    B=2^(q+1)c-2^p rev(c),
    delta=c-f(c)-1.

**Theorem 1.** These actual strong feedback-one graphs satisfy, at every t>=0,

    H(T(L_t)) = A 4^t + B 2^t - delta.            (4)

**Proof.** After reading 1 0^(p+t), the integer is x=2^(p+t), and f(x)=x+1 because p+t>=1. Appending bin(c) gives

    y=2^(p+t+ell)+c,
    f(y)=(2^ell-rev(c))2^(p+t)+f(c).              (5)

For completeness, appending any fixed suffix v to an even prefix x adds revc(v)x+f(value(v))-1 to f. This follows by summing the prefix contribution at each zero of the suffix; its first parity contribution vanishes because x is even. Equation(5) is that identity with f(x)=x+1 and 1+revc(c)=2^ell-rev(c).

If q+t>=1, the following zeros add y(2^(q+t)-1)+1 to f: y is odd, so precisely the first zero contributes the parity correction. Substitute this into(2) with n_t=y2^(q+t), which is even. Collecting the terms yields(4).

If q=t=0, no such zeros occur and n_t=y is odd. Equation(2) is instead y+f(y)+1. Substitution of(5) gives exactly the same expression(4). Thus the zero-length boundary is included without treating a nonexistent zero as a real transition.

## 3. A literal fixed three-port period

In T(L_t), retain the three vertices a, b, and v=x_(s-1), the final base vertex. Repeat these two births:

1. Add a vertex u which loses precisely to a and b among all existing vertices; it beats every other existing vertex.
2. Add a vertex z which beats precisely v among all existing vertices; it loses to every other existing vertex, including the new u.

The first birth inserts a new zero-bit vertex immediately after a in the transitive base. The second inserts a new zero-bit vertex immediately before v. The apex b beats both new vertices. Therefore the new graph is exactly T(L_(t+1)) after relabelling the base in its increasing order. These births retain the same three physical ports throughout: no port is replaced by a fresh vertex.

In the usual orientation language this is one left birth losing only to {a,b}, followed by one right birth beating only {v}. The process enlarges two different transitive modules simultaneously. In particular it is not a one-root transitive-clone clock, and no amplitude bound for that narrower class applies to it.

## 4. Pure clocks are exactly the palindromic balanced cases

By(4), an eventual affine recurrence H_(t+1)=4H_t+c_0 can hold only when B=0: its difference is -2B2^t+3delta. Conversely B=0 gives the exact recurrence from t=0.

Since c and rev(c) are odd, B=0 implies p=q+1 by comparing powers of two, and then c=rev(c). Consequently:

**Corollary 2.** Within the complete language(3) with p>=1, q>=0, the pure base-four clocks are exactly

    p=q+1 and bin(c) is a palindrome.

Their exact count and returning rule are

    H_t=2^(ell+2)4^(q+t)-delta,
    H_(t+1)=4H_t+3delta.                           (6)

This is an iff statement for the displayed language, not a classification of all mixed-port clocks.

The possible offsets have a finite exact inverse. From(1), f(c)<=c+1 for all c>=0, and for odd c>=1,

    1<=f(c)<= (c+1)/2.                            (7)

One can prove the first inequality inductively, using f(2x+1)=f(x) and the sharper odd-input bound before treating an even input whose half is odd. Thus delta=-1 occurs only at c=1; every c>=3 gives a positive odd delta. For prescribed positive delta, any possible c satisfies

    delta+2 <= c <= 2delta+3.                     (8)

Testing the finitely many odd palindromes in that interval is a complete inverse for the constant offset within(6). In particular delta=1 permits exactly c=3 (binary11) and c=5 (binary101). Their respective amplitudes with q=0 are16 and32. The opposite unit constant, delta=-1, has c=1 and amplitude8. Adding q merely advances the phase by multiplying the amplitude by4^q.

## 5. All Mersenne targets, with a sharp feedback-one order

**Theorem 3.** Every target 2^k-1 with k>=4 has an actual strong feedback-one realization on k+2 vertices. Explicitly:

- For k=2t+4, t>=0, use c=3, p=1, q=0 and

      L=2^(k-1)+3*2^t+1.

- For k=2t+5, t>=0, use c=5, p=1, q=0 and

      L=2^(k-1)+5*2^t+1.

The word width is exactly k, and(6) gives the required count. Both rays return by the same physical three-port word as

    H -> 4H+3.

No primality hypothesis is needed. Thus every Mersenne prime with exponent at least5 has a graph from the displayed formula, whether or not infinitely many Mersenne primes exist.

**Minimal order in the feedback-one class.** Every strong tournament whose deletion of one vertex is transitive has the form T(L) after choosing that vertex b. Indeed its transitive deletion has a first vertex a and a last vertex v. Strongness forces b->a and v->b; all remaining incidences with b are arbitrary. This is exactly the gadget construction above. Conversely the displayed construction is strong and feedback-one.

We use g(L)<=L, which is(7)'s first bound. If L is even, the gadget recurrence gives the stronger g(L)=g(L/2)<=L/2. A graph in this class with at most k+1 vertices has binary width at most k-1, hence L<=2^(k-1). If L is odd then L<=2^(k-1)-1, and

    H=g(L)+L-1 <=2^k-3.

If L is even then

    H=g(L)+L <=3*2^(k-2) <2^k-1       (k>=3).

So no smaller strong feedback-one graph realizes the target. The constructed order k+2 is optimal in this class. This is not an unrestricted minimum-order claim: H15 already occurs in the regular five-tournament, whereas this feedback-one construction uses six vertices.

## 6. The unique seven-times-prime member

Among the targets of Theorem3, exactly one has the form 7p with p prime:

    2^9-1=511=7*73.

Indeed divisibility by7 forces k=3j. For j>=2,

    2^(3j)-1=(2^j-1)(2^(2j)+2^j+1).

Both factors exceed one, so a product of7 and a prime forces one factor to be7. The first equals7 exactly at j=3, giving the other factor73. The second equals7 only at j=1, outside k>=4. This includes p=7: no resulting target49 occurs. Trial division certifies73 prime.

The verifier additionally supplies literal graphs for H31, H127, H511 and H8191 and checks primality of31,127,73,8191 by complete trial division. These are examples, not a claim of infinitely many primes in a fixed ray.

## 7. Verification and scope

The standalone verifier evaluates the physical graph family, both binary recurrences, every fixed-port incidence under the actual two births, the complete offset inverse over a stated finite range, and literal Hamilton counts for small graphs. The all-order formula and iff cancellation criterion are proved symbolically above, not inferred from interpolation. All exact graphs and prime certificates in the witness file are regenerated by the verifier. Run normally and with `-O`; receipts and witnesses must agree.

The new contribution is the explicit two-block cancellation, its fixed-port interpretation and the sharp feedback-one Mersenne realization. The scalar gadget recurrence and singleton-prefix transfer are inherited. The broader unrestricted {7,21} spectrum remains open; this packet covers a precise infinite exact-target family, not every prime or7p atom.
