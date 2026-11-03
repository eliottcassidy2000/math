# Decimal carries, dyadic sections and exact source memory

**PROVED elementary identities and operation-relative minimality.** Local
labels D1--D5. Computations are controls, not the source of the all-size
claims. The prior [nine-cycle theorem](../session143_magic_square_cyclic_bridges_20261103/magic_cycle/THEOREM.md)
already identifies translation by four modulo nine and its ternary carry.

## D1 The displayed cycle has a full integer carry lift

Let s(n) be the decimal digit sum and f(n)=s(4n). The first ten positive
values are `4,8,3,7,2,6,10,5,9,4`. Reducing modulo 9, and choosing
representatives 2,...,10, gives exactly the user's nine-cycle. The reason
is 10=1 mod 9, so s(4n)=4n mod 9.

The full single-pass sum retains an additional integer coordinate. For
example 48 has sum 12 in the same phase as 3. For every n>=0, the exact
sequential recurrence is

`f(n+1)=f(n)+4-9 c(n)`,
`c(n)=sum_(j>=1) [4n mod 10^j >= 10^j-4]`.

Only finitely many terms contribute. Indeed
`s(x)=x-9 sum_(j>=1) floor(x/10^j)`; subtract at x+4 and x. Thus the
phase cycle and its actual decimal carry history coexist. At n=2499,
the step 9996->10000 has four carries and changes the sum by -32.
Retaining only the nine-cycle would lose that jump size.

## D2 Multiplication by four has a nearest-neighbour digit potential

Write n in decimal as digits d_1...d_L, including leading zeros if desired.
Then, for every n>=0,

`s(4n)=sum_i s(4d_i)-9 sum_(i=1..L-1) [d_i in {2,7}] [d_(i+1)>=5]`. (1)

Proof: put h(q)=[q=2 mod 5]. For c=0,1,2,3,

`s(4q+c)=s(4q)+c-9 h(q)[c>=2]`.

An addition of 1 cannot carry from the even units digit of 4q. Adding 2
or 3 carries exactly when that digit is 8. Its tens digit is then even,
because 4q is divisible by 4, so there is no second carry. For n=10q+d,
split 4d into its units digit and c=floor(4d/10). This yields

`f(10q+d)=f(q)+s(4d)-9 h(q)[d>=5]`,
`h(10q+d)=[d in {2,7}]`,

and proves (1) inductively.

Consequently `(f(n),h(n))` is the **smallest behavioural quotient** for
the following precise menu: append any decimal digit, and observe s(4n).
Equal pairs have equal responses after every future suffix. Conversely,
the empty suffix distinguishes different f values, and suffix 5
distinguishes different h values when f agrees. This quotient has an
unbounded integer accumulator and a two-state controller.

Its homogeneous linear dimension over Q is exactly three. Coordinates
`(f,h,1)` give an upper bound. States n=0,1,2 span all three coordinates;
the empty, 1 and 5 suffix observations span the dual. Hence the reachable
and observable linear representation cannot have smaller dimension.

## D3 An exact generating function and moment formulas

Let `S(z)=1+z+...+z^9` and `P(z)=1+z^2+z^3+...+z^10`. For every L>=1,

`sum_(0<=n<10^L) z^(s(4n)) = P(z) S(z)^(L-1)`.             (2)

Thus the length generating function, including L=0, is

`1 + t P(z)/(1-t S(z))`.

Here is a proof directly from the source-preserving recurrence. Use
states h=0,1 and attach weight z to each unit of the increment in (1).
The transition matrix, with rows the old state, is

```text
M00=1+z^2+z^3+z^4+z^5+z^6+z^7+z^9
M01=z^8+z^10
M10=2+z^3+z^4+z^7+z^-7+z^-3+z^-4
M11=z+z^8
```

These are Laurent polynomials because an individual carry can decrease
the running sum. The complete word weights have nonnegative exponents.
Direct multiplication gives det M=0 and tr M=S(z), so M^2=S M.
The initial row is (1,0), the terminal column is (1,1), and the one-digit
sum is P. This proves (2) for every length.

For uniform N in 0,...,10^L-1 it follows that

`E s(4N)=9L/2+9/10`,
`Var s(4N)=33L/4+99/100`,
`sum_(n<10^L) s(4n)=9(5L+1)10^(L-1)`.

The first-digit weights are `{0,2,3,...,10}`, with sum 54 and square sum
384. The nine nonzero weights recover the prior cycle/magic-square
labels and their mean 6. Equation (2) is a distribution identity;
it does not assert independence of the actual digits along each source.
The all-digit weighted sum has a shorter recurrence than the individual
digit operations. That is the same distinction between a weighted
attachment recurrence and a carrier closed under every labelled word.

## D4 Every fixed dyadic multiplier has finite-range decimal interactions

Let m=2^a, a>=1. In appending d to a prefix q, put c=floor(md/10). Then

```text
s(m(10q+d)) = s(mq)+s(md)-9 K_a(q,d),
K_a(q,d) = sum_(j=1..a-1)
           [(mq mod 10^j)+(c mod 10^j) >= 10^j].
```

This follows from the ordinary addition carry formula for mq+c.
There are no carries at a boundary j>=a: m divides 10^j, mq mod 10^j
is a multiple of m at most 10^j-m, and 0<=c<m. For j<a the residue
mq mod 10^j depends only on q mod 5^j, hence only on its last j decimal
digits. Thus the entire digit sum is an exact local potential involving
at most a consecutive input digits. For a=1 it is digit-additive; for
a=2 it is exactly (1). This is a sufficient locality bound, not a
minimality theorem for every a.

## D5 Halving restores a different, unbounded clock

The old [intrinsic halving decoder](../decoder_halving_20260925.md) writes
every n>0 as `2^v u`, u odd. Exactly v legal halvings reach u. That is
the height coordinate lost by the append observer: all n=10^k, k>=1,
have `(f,h)=(4,0)`, but their halving clocks are k.

More strongly, allow decimal appends and the Boolean test that another
halving is legal. **Every two distinct positive integers are separated
by a future context.** Let r=v2(m-n). Choose k>=1 with 5^k>=2^(r+1),
let t be the residue of -5^k n modulo 2^(r+1), and append the same
k-digit suffix a=2^k t to n and m (leading zeros allowed). Then

`v2(10^k n+a)>=k+r+1`, `v2(10^k m+a)=k+r`.

The first follows by construction. The difference of the two sources
has valuation k+r, giving the second. Also 0<=a<10^k, so this is an
actual decimal suffix. The legality of the (k+r+1)st halving separates
them. No finite quotient of all positive integers can be faithful for
this larger menu. This is a constructive guard result, not a claim
about computational hardness or Collatz convergence.

There is also an exact even-realm conjugacy. For each a>=1, on positive
multiples x of 2^a define

```text
C_a(x)=x/2       if v2(x)>a,
C_a(x)=3x+2^a    if v2(x)=a.
```

Then `C_a(2^a n)=2^a C(n)` for the ordinary Collatz map C. Native steps
and source integers are retained. At a=1 the literal orbit is

`6 -> 20 -> 10 -> 32 -> 16 -> 8 -> 4 -> 2`,

followed by the cycle `2 -> 8 -> 4 -> 2`. Thus 6,10,4 are not merely
residue labels here. Removing all powers of two from an unscaled integer
reaches the odd section; this conjugate stops its halving at height a
to remain in its declared even section. The maps are different charts
of the same dynamics, not independent convergence evidence.

For odd n, if b=v2(3n+1) and u=(3n+1)/2^b, then
`u=2^(-b)=(-1)^b mod 3`. The ternary phase therefore remembers the
parity of the halving clock, while higher powers of 3 retain more of
the source. Decreasing height within one tower does not prove descent
across successive odd sections.
