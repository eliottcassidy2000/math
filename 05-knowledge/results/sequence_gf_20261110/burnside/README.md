# A connected-diagram compiler for Burnside generating functions

**Status:** PROVED below as identities of formal series; FINITE-EXACT controls
in `cluster_compiler.out` and `cluster_order4.out`. This is a new research note,
not a canon promotion or a literature-priority claim.

**Result in one sentence.** Expand the existing odd-cycle Burnside enumerator
at `q=1+s`, take its logarithm, and a Taylor coefficient of order `r` uses only
connected graphs on at most `r+1` cycle vertices; a divisor/Möbius transform
turns every graph into explicit `lcm`-coupled Lambert sums of logarithms and
rational functions. This gives exact generating functions and convolution
recurrences for new marked-orbit sequences, in polynomial time at each fixed
Taylor order. It is not a faster full A000568 counter.

## 1. Inheritance and the changed question

The source is the existing Davis/Burnside formula, not a newly discovered
formula for unlabelled tournaments:

- [THM-514: unit tail and core kernel](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
- [THM-516: gcd/Euler-phi quadratic form](../../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md).
- [THM-4845: gcd Gaussian field and the non-D-finite boundary](../../../../01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md).
- [THM-5097: integral future-response recurrence](../../../../01-canon/theorems/THM-5097-factorial-gauge-integral-cycle-response-recurrence.md) already changes the rational response states to integer ones. This note does not repeat that change of basis.
- [THM-5098: precision-dependent cycle cutoff](../../../../01-canon/theorems/THM-5098-precision-cycle-cutoff-and-normalized-tournament-congruences.md) already bounds surviving cycles at fixed 2-adic precision. The bound here is instead on **connected cycle interactions at fixed marking order**.

The closest proved mechanism is the gcd divisor factorization. The hostile
example is the all-fixed-point cycle type, whose edge-orbit count is
`binom(n,2)` and which prevents any fixed Taylor jet near `q=1` from approximating
the full `q=2` count asymptotically. The corrected near miss is mistaking the
Gaussian expectation for a closed evaluation or a faster algorithm. The
underused sidecar is the full edge-orbit polynomial of THM-516, rather than
just its value at 2.

Concept board: odd-cycle gas; positive marked interactions; primitive-subset
Möbius inversion; divisor-to-lcm graph transform; exact coefficient recurrence;
Hamiltonian-path support, which is **not retained** by this scalar census.

## 2. The marked enumerator and the soluble base point

For an odd partition `lambda=(p_1,...,p_k)` put

    h(p) = (p-1)/2,
    e(lambda) = sum_i h(p_i) + sum_{i<j} gcd(p_i,p_j),
    z(lambda) = product_p p^m_p m_p!.

Define the ordinary size GF

    A(x,q,u) = sum_lambda q^e(lambda) u^ell(lambda) x^|lambda| / z(lambda).

The variable `u` marks permutation cycles, not tournament SCCs or vertices.
At `q=2,u=1` this is the A000568 GF. At `q=1`, interactions disappear:

    A(x,1,u) = exp(u artanh(x)) = ((1+x)/(1-x))^(u/2).

In particular

    [x^(2m)] A(x,1,1) = [x^(2m+1)] A(x,1,1)
                      = binom(2m,m)/4^m.

For the fixed-point-free odd core, restrict all `p_i` to odd `p_i>=3`:

    A_core(x,1,u) = exp(u(artanh(x)-x)).

At `q=2`, `[x^m u^t] A_core` is exactly the inherited kernel `B[m,t]`.

Set

    L_r(x,u) = [s^r] log A(x,1+s,u),
    B_r(x,u) = [s^r] A(x,1+s,u).

Then

    [x^n] B_r(x,u) = sum_(|lambda|=n) binom(e(lambda),r) u^ell(lambda)/z(lambda).

Thus `n! [x^n] B_r` is an integer polynomial in `u`: it counts an odd-order
permutation together with a choice of `r` distinct orbits of its action on
unordered vertex pairs. These are **marked edge-orbits**, not marked edges of
a tournament and not Hamiltonian paths.

## 3. Primitive-subset kernels

For `b>=1` define

    Psi_b(d) = sum_(c|d) mu(d/c) binom(c,b).

Then

    binom(g,b) = sum_(d|g) Psi_b(d),
    Psi_b(d) >= 0.

**Combinatorial proof.** `Psi_b(d)` counts `b`-element subsets `S` of
`{1,...,d}` with `gcd(d, all elements of S)=1`. Partition all `b`-subsets of
`{1,...,g}` by that gcd and divide their elements by it. This gives the first
identity, and Möbius inversion gives the formula and positivity. In
particular `Psi_b(d)=0` for `d<b`.

Only `b>=1` is used on graph edges. If one wants `b=0`, use
`gcd(d,empty)=d`, so `Psi_0(d)=1[d=1]`; no empty-edge factor is inserted.

The first cases connect directly to the inherited gcd field:

    Psi_1(d) = phi(d),
    Psi_2(d) = (J_2(d)-phi(d))/2,
    Psi_3(d) = (J_3(d)-3 J_2(d)+2 phi(d))/6,

where `J_j(d)=sum_(c|d) mu(d/c)c^j` is the Jordan totient. Thus the marking
compiler uses a **subset version of the Jordan totients**, not only Euler phi.

For odd `d` and `a>=0`, write `D=x d/dx` and set

    V[a,d](x) = sum_(p odd, d|p) binom((p-1)/2,a) x^p/p
              = binom((D-1)/2,a) artanh(x^d)/d.                 (1)

The binomial in the last expression is a polynomial differential operator.
For every fixed `a`, it is a finite linear combination of rational functions
of `x^d` and `artanh(x^d)`. For the core version, remove the `p=1` term:
only `V[0,1]` changes, to `artanh(x)-x`.

## 4. Connected-diagram theorem

Let `G` be a connected simple graph on labelled vertices `{1,...,k}`. Give
its edges positive integers `b_e>=1`, its vertices integers `a_i>=0`, and
require

    sum_e b_e + sum_i a_i = r.

For every edge choose an odd divisor label `d_e>=1`. At a vertex set

    D_i = lcm(d_e : e incident to i),

with the empty lcm equal to 1. Then, for every `r>=1`,

    L_r(x,u)
      = sum_(k=1)^(r+1) u^k/k!
          sum_(G connected on [k])
          sum_(a,b: total r, b_e>=1)
          sum_(d_e odd>=1)
            (product_e Psi_(b_e)(d_e))
            (product_i V[a_i,D_i](x)).                       (2)

For `k=1`, the sole graph has no edges, and its contribution is `u V[r,1]`.
Also `L_0(x,u)=u artanh(x)`.

All sums are interpreted coefficientwise in `x,s,u`. At fixed size `n`,
only finitely many labels can contribute. Formula (2) proves:

1. `L_r` is a polynomial in `u` of degree at most `r+1`.
2. All coefficients of `L_r(x,u)` are nonnegative rational numbers.
3. The whole coefficient of `s^r` is expressible using finitely many graph
   shapes and explicit divisor/lcm sums. There is no partition enumeration in
   this compiler.

**Proof.** Rewrite the partition sum as a labelled-cycle gas:

    A(x,q,u) = sum_(k>=0) 1/k! sum_(p_1,...,p_k odd)
        product_i (u x^p_i q^h(p_i)/p_i)
        product_(i<j) q^gcd(p_i,p_j).

The denominator `k!` turns repeated lengths into `m_p!`, so this is exactly
the original sum. Expand each pair factor as

    q^gcd(p_i,p_j) = 1 + (q^gcd(p_i,p_j)-1).

Choosing the second term on a set of pairs produces a simple graph on the
cycle labels. The elementary labelled exponential formula says that the
logarithm keeps its connected graphs. This is the classical Mayer/connected
cluster expansion, applied here to the inherited cycle gas.

Now put `q=1+s`. A present edge has expansion

    (1+s)^g-1 = sum_(b>=1) binom(g,b) s^b,

and a vertex contributes

    (1+s)^h = sum_(a>=0) binom(h,a) s^a.

A connected graph on `k` vertices has at least `k-1` edges, so a coefficient
of order `r` requires `k<=r+1` and at most `r` edges. Apply the primitive-subset
identity on each edge. Its divisor must divide both endpoint lengths. Hence
the length at vertex `i` must be divisible by the lcm of its incident divisors;
the length sums then separate and become (1). This proves (2).

**Complexity, precisely.** To obtain the prefix through `x^N` at fixed `r`,
there are only `C_r` graph/marking choices, where `C_r` depends on `r` but not
on `N`. Each has at most `r` edges. Restricting each edge divisor to `d<=N`
gives `O(N^r)` tuples, and each product of at most `r+1` truncated vertex
series can be formed in `O(r N^2)` rational operations. Consequently a naive
implementation costs `O(C'_r N^(r+2))` rational arithmetic operations; it can
prune any tuple with `sum_i D_i>N`, and cache equal lcm profiles. The diagram-stage
rational operands have polynomial bit length at fixed `r`: their denominators
divide a fixed power of `lcm(1,...,N)` times an `r`-dependent constant. The
subsequent coefficient recurrence has at worst factorial-size denominators,
also of polynomial bit length. Thus the construction also has polynomial
bit complexity. No uniform polynomial bound in **both**
`r,N` is asserted. The connected-graph constant grows rapidly with `r`.

This is a fixed-marking-order result, separate from THM-5098's fixed-adic-
precision result and THM-5097's full-count integer recurrence.

## 5. First two explicit generating functions

All divisor indices in this section are positive odd integers. The first
layer is

    L_1(x,u) = u V[1,1]
                 + u^2/2 sum_d phi(d) V[0,d]^2.             (3)

In elementary functions,

    L_1(x,1) = (1/2)(x/(1-x^2)-artanh(x))
                  + (1/2) sum_d phi(d)/d^2 artanh(x^d)^2.

The second layer is

    L_2(x,u) = u V[2,1]
       + u^2 sum_d phi(d) V[1,d] V[0,d]
       + u^2/2 sum_d Psi_2(d) V[0,d]^2
       + u^3/2 sum_(d,e) phi(d) phi(e)
                    V[0,lcm(d,e)] V[0,d] V[0,e].            (4)

The four terms are, respectively: two internal orbit marks on one cycle;
one internal and one cross-cycle mark; two marks on a single interaction;
and a path of two interactions on three cycles. This graph interpretation
explains the factors, instead of guessing them from initial coefficients.

A directly usable coefficient formula follows from (3). Let
`H_odd(k)=sum_(j=1)^k 1/(2j-1)`. Then

    [x^(2m+1)] L_1(x,1) = m/(2m+1),
    [x^(2m)] L_1(x,1)
       = (1/(2m)) sum_(d|m, d odd) phi(d)/d H_odd(m/d),       (5)

for `m>=0` in the first line, `m>=1` in the second. To prove the second line,
use `[x^(2k)] artanh(x)^2=H_odd(k)/k`, obtained by multiplying its two odd
power series and using `1/(a(2k-a))=(1/a+1/(2k-a))/(2k)`.

The first moment GF is now explicit:

    B_1(x,1) = sqrt((1+x)/(1-x)) L_1(x,1).

The second marking GF is likewise explicit:

    B_2(x,1) = sqrt((1+x)/(1-x)) (L_2(x,1)+L_1(x,1)^2/2).

These are ordinary GFs with rational coefficients. Multiplication of the
`x^n` coefficient by `n!` gives the corresponding integer sequence.

## 6. A triangular recurrence, two moment conventions, and bounds

Write `b[n,r]=[x^n] B_r(x,u)` and `l[m,j]=[x^m] L_j(x,u)`. Then

    b[0,0]=1,   b[0,r]=0 for r>0,
    n b[n,r] = sum_(m=1)^n sum_(j=0)^r
                   m l[m,j] b[n-m,r-j].                    (6)

This is simply `x partial_x A=A x partial_x log A`. It computes every marked
sequence through the chosen Taylor order once (2) supplies the kernels.
It is a convolution recurrence with arithmetic kernels, not a fixed-width
polynomial-coefficient recurrence for A000568.

For ordinary moments of `e`, substitute `q=exp(t)`. If

    log A(x,exp(t),u) = L_0 + sum_(m>=1) kappa_m(x,u)t^m/m!,

then the exact Stirling conversion is

    kappa_m = sum_(r=1)^m r! S(m,r) L_r.

The same conversion takes `binom(e,r)`-moment sequences to `e^m`-moment
sequences. Thus the Jordan-totient Gaussian/Wick viewpoint and the
primitive-subset positive-marking viewpoint are two bases of the same
compiler. They should not be counted as independent discoveries.

There is an exact full-count identity

    A(x,2,u) = sum_(r>=0) B_r(x,u),                           (7)

because `s=1`. For each `x^n`, the sum terminates at `r=binom(n,2)`. At `u=1`,
truncating it gives a certified coefficientwise lower bound on A000568.

**Necessary failure boundary.** No fixed number of Taylor layers approximates
the full count at large `n`. Indeed

    0 <= b[n,r] <= binom(binom(n,2),r) [x^n] A(x,1,1)
                  = O_r(n^(2r-1/2)),

whereas the identity permutation alone gives

    A000568(n) >= 2^binom(n,2)/n!.

Thus the ratio of every fixed Taylor truncation to A000568 tends to zero.
This is a precise boundary, not just a warning about possible numerical
accuracy. The marked sequences have useful explicit GFs even though their
finite collection cannot carry the full tournament count.

**Formal versus analytic.** All identities above are formal in size and
marking variables, so they do not require a convergent A000568 OGF. At each
fixed Taylor order, the polynomial coefficient bound gives radius at least
1. It does **not** follow that the full positive-temperature series is
analytic near `q=1` for a fixed real `x>0`; the unbounded all-fixed-point
subfamily already prevents that conclusion. No exchange of these two
summations is used.

## 7. Reproduction and finite examples

Run from the repository root:

    python3 05-knowledge/results/sequence_gf_20261110/burnside/cluster_compiler.py --n 18 --order 3
    python3 05-knowledge/results/sequence_gf_20261110/burnside/cluster_compiler.py --n 24 --order 4
    PYTHONDONTWRITEBYTECODE=1 python3 05-knowledge/results/sequence_gf_20261110/burnside/benchmark_jets.py

The main engine uses labelled connected graphs, primitive-subset weights,
lcm profiles, and exact Fraction series multiplication. A separate control
uses odd partitions directly. They agree for every reported size and layer.
Controls also cover (3), (4), (5), coefficient positivity, factorial
integrality, the `q=1` base, the core restriction, marked-cycle specializations
`u=2,3`, and complete Taylor reconstruction of `q=2` through `n=4`.

For `n=0..8`, the first marked integer sequence is

    n! b[n,1] = 0, 0, 1, 5, 22, 138, 927, 7431, 57564.

The second is

    n! b[n,2] = 0, 0, 0, 3, 23, 189, 1777, 17140, 152538.

The third is

    n! b[n,3] = 0, 0, 0, 1, 20, 200, 2399, 28440, 309452.

The first three logarithmic marking layers begin

    L_1 = x^2/2+x^3/3+x^4/3+2x^5/5+11x^6/30+...,
    L_2 = x^3/2+x^4/3+7x^5/10+9x^6/10+...,
    L_3 = x^3/6+2x^4/3+2x^5/3+74x^6/45+....

The explicit hostile control is order 3 at `n=4`:

    sum_(r=0)^3 b[4,r] = 37/12 < 4 = A000568(4).

Timing is implementation evidence only. The initial full jet through
`N=18,r=3` took about `0.05s` on this worker; a separate compiler-only run
through `N=80,r=3` took about `4.7s`, and through `N=200,r=2` about `1.2s`.
These are not comparisons with the full-count recurrence, which computes a
different quantity, and are not fastest-known claims.

## 8. What transfers, and what does not

**Source -> target:** the odd-cycle weighted permutation sum -> connected
cycle-interaction graphs -> edge-divisor labels -> vertex lcm filters.

**Preserved:** every coefficient of the `q=1` Taylor expansion, including cycle
count; the core restriction; factorial integrality of marked permutation
counts; exact reconstruction with the complete finite jet for a fixed size.

**Destroyed if one retains only a fixed jet:** large edge-orbit multiplicity;
full tournament isomorphism class information; Hamiltonian-path counts and
their support. The missing sidecar for Chan–Pak is still a faithful
Hamiltonian-path-value refinement. Neither the census nor its positive
marking layers prove that an odd value is attained.

**Bridge to other enumerators:** the connected-graph argument works for any
allowed cycle set and any nonnegative integral one-cycle/pair interaction.
For ordinary graphs use all positive lengths, `h(p)=floor(p/2)`, and the same
gcd interaction. Replace the vertex kernel by
`sum_(p>=1,d|p) binom(floor(p/2),a)x^p/p`; its base is
`-log(1-x^d)/d`, with the parity correction retained in higher vertex marks.
For directed graphs with two states per ordered pair, use
`h(p)=p-1` and pair interaction `2 gcd(p,q)`; replace `Psi_b` by the Möbius
transform of `binom(2d,b)`. These are proved substitutions of the same graph
identity, not an assertion that their vertex kernels are the odd one (1).

**Next useful experiments:** retain the coefficient of `u^k` instead of
specializing it; compile rooted/marked-path versions when a genuine
Hamiltonian-path observable is available; compare positive ordinary marking
layers with the signed connected-overlap layers of the heat lane. A negative
weight in the latter is an obstruction to identifying those two carriers,
not a contradiction to positivity here.

Novelty claim is local and limited: this note supplies the explicit
primitive-subset/divisor/lcm marking compiler, its first two closed kernels,
and the fixed-order complexity theorem. The exponential formula, Mayer
expansion, Möbius inversion, Jordan totients, and Davis/Burnside input are
classical; the Gaussian and response recurrences are inherited. No broad
literature search or claim of first discovery outside the repository is made.
