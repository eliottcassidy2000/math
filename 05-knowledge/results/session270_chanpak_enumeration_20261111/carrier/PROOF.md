# One marked carrier, two exact generating functions

**Status: PROVED IN NOTE.** Burnside and path-cover bookkeeping are inherited
mechanisms. The common marked-core carrier and formulas below synthesize them;
no external priority claim. This is not the full unmarked A000568 census at
order k+m and does not prove Chan–Pak.

## The carrier and the quotient

Fix k>=1. A structure has k core vertices with an arbitrary tournament C,
and m further vertices with a specified transitive order. Every core-spine
arc is arbitrary. Isomorphisms preserve the designated core subset and fix
each position of the ordered spine. They may permute core vertices.
Let B(k,m) count these ordinary **marked** classes.

For an odd partition lambda of k, write t(lambda) for its number of parts,
z(lambda) for its permutation centralizer and e(lambda) for its pair-orbit
count. Then

    B(k,m)=sum_(lambda odd partition of k) 2^e(lambda) 2^(m t(lambda))/z(lambda).  (1)

**Proof.** Apply Burnside to S_k acting on the core and all cross arcs.
An even permutation cycle reverses a core pair and fixes no core tournament.
For an odd cycle type, 2^e core tournaments are fixed. For each ordered spine
position, the cross orientations must be constant on each core cycle, giving
2^t choices independently at each of the m positions. Divide by z and sum.

Thus B(k,0)=A000568(k), and for fixed k

    sum_(m>=0) B(k,m) z^m = sum_lambda (2^e(lambda)/z(lambda))/(1-2^t(lambda) z). (2)

Every t=k,k-2,...,1 or 2 occurs with positive coefficient. Therefore the
minimal constant-coefficient recurrence over Q has characteristic polynomial

    product_(1<=t<=k, t=k mod2) (E-2^t),                       (3)

of order ceil(k/2). Distinct exponential sequences are linearly independent,
so no root can cancel. This is an exact recurrence in **spine length m**,
not a forbidden finite recurrence in unrestricted core size k.

Examples:

    B(2,m)=4^m,
    B(3,m)=(2*2^m+4*8^m)/3,
    B(4,m)=(4*4^m+8*16^m)/3.

Once the k-dependent coefficients are precomputed, a query uses O(k) integer
powers or an order-ceil(k/2) recurrence. Output length still grows with km.
The expensive core inventory has not disappeared; evaluation at m=0 only
returns the A000568(k) input.

## Labelled path mass becomes a binomial moment

Let Z(k,m) sum directed Hamiltonian path counts over all labelled core
orientations and all cross orientations, with the spine order fixed.
For k>=1 and m>=0,

    Z(k,m)=k! 2^(binom(k,2)-2k) sum_(j=0)^(k+1) binom(k+1,j)(2^k j)^m,          (4)

where 0^0=1. The prefactor may be rational, but the whole expression is an
integer. Equivalently, in the uniform labelled ensemble,

    E H = k! 2^(1-k) E[X^m],     X~Binomial(k+1,1/2).          (5)

**Proof.** Let q_r(C) count ordered covers of C by r nonempty directed paths.
The `chanpak/` averaged-spine theorem gives its exponential generating function

    sum_m Z_C(m) z^m/m!
      = (e^(2^k z)+1)^2/4 * sum_(r=1)^k q_r(C)((e^(2^k z)-1)/4)^(r-1).

Independently, sum_C q_r(C)=k! binom(k-1,r-1)2^(binom(k,2)-k+r): choose the
ordered vertex permutation and r-1 cuts; the k-r path arcs are forced and
all other core arcs are free. Substitute and apply the binomial theorem.
The result is k!2^(binom(k,2)-2k)(e^(2^k z)+1)^(k+1), proving (4).
There are 2^(binom(k,2)+km) labelled structures, giving (5).

For k=2 and one fixed orientation of the labelled core, the H mass of the
4^m words is

    Z_C(m)=Z(2,m)/2=(12^m+3*8^m+3*4^m+0^m)/8.

For m>=1 the zero-base term vanishes. In contrast the marked class count is
4^m. Therefore the mean is (3^m+3*2^m+3)/8 for m>=1.
Large mass or large class count does not imply all target H values appear.

## Why this is a genuine bridge and where it stops

Source: arbitrary core plus ordered transitive spine. Two observers:
(1) marked class count, obtained by Burnside's cycle-number spectrum;
(2) labelled Hamiltonian path mass, obtained by the ordered-cover polynomial.
Both produce finite exponential sums in the same extension parameter m.
The first observer has bases 2^t; the second has bases 2^k j. They are
not the same recurrence, even though their underlying graphs coincide.

Forgotten information is different: (1) erases each class's path count;
(2) erases graph isomorphism and the distribution of path counts. Summing
ordinary core classes uniformly instead of labelled cores would be wrong:
labelled mass requires weights k!/|Aut(C)|. Forgetting the designated core
also identifies many marked structures and does not preserve B(k,m).
The two-core witness compiler exploits this marked language for existence;
it never treats its words as distinct unmarked classes.

Controls in `check.py` enumerate the entire marked universe for small k,m,
canonicalize under actual core permutations, and independently count paths
by subset DP. The formulas are also tested at larger m without graph listing.
