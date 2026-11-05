# Self-converse low-cycle extremizers from palindromic seams

Session 197. **Explicit all-order constructions, exact extremal values in five residue classes, and a construction-family optimum in the sixth.** This does not settle the unrestricted Hamiltonian-path spectrum. The product formula specializes the inherited one-seam calculus, and the lower bound is Busch's published theorem. No literature-priority claim is made.

## 1. The literal graph and its exact product

Choose a nonempty list of integers `d_1,...,d_k>=2`. Put `s_0=0`, `s_i=d_1+...+d_i` and `n=1+s_k`. On vertices `0,...,n-1`, orient every arc forward in the natural order except

    s_i -> s_(i-1),       1<=i<=k.                       (1)

Call the resulting tournament T(d_1,...,d_k), with root s_k. Then

    H(T)=product_i (1+2^(d_i-1)),
    C(T)=1,   s_root(T)=1,   R_root(T)=H(T).             (2)

In particular every such graph is strong.

Here H counts directed Hamiltonian vertex words, C counts directed Hamiltonian cycles modulo rotation, and R counts ordered spanning path covers (P,Q), with P possibly empty and Q nonempty and starting at the root.

**Proof.** Begin with the singleton rooted at 0, whose H, R, and rooted starting-path count are all one. For each gap d, append a transitive tournament G on d new vertices, mark its sink b, and reverse only the cross arc b->a to the previous root a. This is exactly (1).

For this marked transitive gadget the constants in `ROOTED.md` are

    g=1, t=0, E=1, Q=1, L=B=2^(d-1), H(G-b)=1.

The L covers correspond to choosing an arbitrary subset of the d-1 nonroot vertices to precede b; the remaining vertices form the other increasing path. Q=1 because a path starting at the sink consists only of that sink. B assigns each of the d-1 transitive vertices independently to one of two ordered slots.

The inherited one-seam laws consequently give

    H'=H+2^(d-1)R,       R'=H+2^(d-1)R.

Thus R=H is maintained and each step multiplies H by 1+2^(d-1). The new root-start count is the old one times H(G-b), hence remains one. Every spanning cycle uses the sole reverse cross arc b->a; cutting there gives uniquely an old Hamilton path starting at a and the increasing G-path ending at b. Therefore there is exactly one Hamiltonian cycle. Total order is at least three because d>=2, so the two-singleton exception does not arise. The actual cycle proves strongness. QED.

The required one-seam identities and conventions are proved or explicitly sourced in `ROOTED.md`; this paragraph records their direct specialization rather than claiming a new general transfer law.

## 2. Reflection gives self-converseness

If the gap list is a palindrome, the vertex reflection

    rho(v)=n-1-v

is an isomorphism from T to its converse. Indeed a forward ordinary pair is sent to a reversed ordinary pair. Palindromicity gives `n-1-s_i=s_(k-i)`, so every exceptional pair `{s_(i-1),s_i}` is sent to the exceptional pair `{s_(k-i),s_(k-i+1)}`. Its orientation is also reversed. This checks every arc, not just the Hamiltonian cycle.

Hence every palindromic gap list supplies a self-converse strong tournament with exactly one Hamiltonian cycle and the product (2).

We do not assert that palindromicity is necessary for an arbitrary presentation to be self-converse, or that these graphs exhaust self-converse low-cycle tournaments.

## 3. Sharp self-converse minima whenever six does not divide the order

The external lower bound is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3, Theorem 1, [DOI](https://doi.org/10.37236/1141), [official primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). Its exact residue-class bound, for strong order n>=3, is

    L(n)=3*5^((n-3)/3)        if n=0 mod 3,
         5^((n-1)/3)          if n=1 mod 3,
         9*5^((n-5)/3)        if n=2 mod 3.             (3)

The primary theorem was checked in the preceding low-cycle packet and its independent audit. It is imported here, not inferred from tests.

For every n>=3 with 6 not dividing n, a palindromic gap list attains (3):

- `n=3k+1`, k>=1: use k gaps equal to 3. The count is 5^k.
- `n=3k+5`, k>=0: use `2, 3 repeated k times, 2`. The count is 9*5^k.
- `n=6t+3`, t>=0: use `3 repeated t times, 2, 3 repeated t times`. The count is 3*5^(2t).

These cases cover exactly the five nonzero residue classes modulo six. Each graph is self-converse, strong, and has C=1. Since (3) applies to every strong tournament, it also bounds the self-converse subclass from below.

**Theorem.** For every n>=3 with 6 not dividing n, the minimum H among self-converse strong order-n tournaments is exactly L(n), and a minimizer has a unique Hamiltonian cycle.

This is an existence and numerical minimum theorem, not a classification of all minimizers.

## 4. Multiples of six: a sharp family optimum, not a global classification

A separate reviewed Session197 result, `../session197-orbit-wildcard/EVEN_CONVERSE.md`, Section 1, proves H=1 mod 4 for every even-order self-converse tournament. Its frozen proof hash is `0fc0e7ca35a1d48626ffd1b49251b64dcc0564f66594bfd0e70307287c4c4cbb`; the independent report is `../session197-orbit-wildcard-graph-audit/EVEN_CONVERSE_REPORT.md`.

For n=6t, (3) equals `3*5^(2t-1)`, which is 3 mod 4. Thus every self-converse strong tournament has

    H >= 3*5^(2t-1)+2.                                (4)

At n=6 the single gap 5 gives H=17, attaining (4). Therefore the self-converse strong minimum at order six is exactly 17.

For t>=2 the palindromic list

    3 repeated (t-2) times, 2,2,3,2,2,
    3 repeated (t-2) times

has total gap sum 6t-1 and yields

    H = 405*25^(t-2).                                 (5)

For example order twelve has the explicit count 405. The global lower bound there is 377. **No equality between those two numbers is asserted.** In general the upper bound (5) is 27/25 times the unrounded Busch bound, hence is less than 27/25 times the lower bound (4).

Moreover (5) is the exact minimum within the palindromic-gap family defined in Sections 1–2. Here is an all-order proof of this narrower statement.

Write f(d)=1+2^(d-1). A side gap d>=5 can be replaced, together with its reflected partner, by the reflected pair of lists `3,d-3` and `d-3,3`. This preserves the order and palindrome but strictly decreases the product, since

    f(d)-f(3)f(d-3)=3*2^(d-4)-4 > 0.

A side gap 4 can be split as 2,2 without changing the product, because f(4)=f(2)^2=9. Repetition reduces all side gaps to 2 or 3.

The total gap sum 6t-1 is odd. Therefore the palindrome has one central odd gap, necessarily at least 3. A central odd gap d>=7 can be replaced by `2,d-4,2`, strictly improving the product:

    f(d)-f(2)^2 f(d-4)=7*2^(d-5)-8 > 0.

Its center is eventually 3 or 5. Let a count pairs of side gaps 2 and b count pairs of side gaps 3.

- Center 3 gives `2a+3b=3t-2`, so `a=2 mod 3`. For t>=2 the least feasible a is 2 and b=t-2. Increasing a by 3 and decreasing b by 2 multiplies the cost by `3^6/5^4=729/625>1`. The optimum is therefore `5*3^4*5^(2t-4)=405*25^(t-2)`.
- Center 5 gives `2a+3b=3t-3`, so `a=0 mod 3`. The optimum is a=0,b=t-1, with count `17*25^(t-1)=425*25^(t-2)`, strictly larger when t>=2.

At t=1 only the center-5 case is feasible, giving 17. The transformations show that every original palindromic gap list has product at least the reduced optimum, and the displayed lists attain it. This is a family optimum only; it does not close the global gap in (4)–(5).

## 5. Reproduction and evidence separation

`verify_self_converse_seams.py` is standalone standard-library Python and uses explicit exception checks, including under `-O`. It checks actual adjacency rows, root path-cover responses and Hamiltonian cycles for all gap compositions of sum at most eleven; reflection on palindromic cases; the displayed all-order formulas on a long finite range; and exhaustive palindromic inventories at small orders plus a separate dynamic-programming family minimum.

The verifier's finite checks are controls on the construction and arithmetic. The one-seam induction, reflection proof, published lower bound, and gap-reduction argument prove the all-order claims. The certificate records the finite ranges exactly. The complete `{7,21}` spectrum, all C=1/C=2 structural questions, and the global self-converse minimum at orders 6t>=12 remain open here.
