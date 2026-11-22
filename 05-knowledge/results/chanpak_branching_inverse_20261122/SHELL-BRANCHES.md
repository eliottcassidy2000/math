# A two-port binary code with one ordinary Hamiltonian-count readout

**Status: PROVED.** This gives actual disjoint ordinary count images for
every fixed-width family of binary words, not merely decoded counts
or a vector of observations. It does not make the image an interval or
give a Chan–Pak descent for arbitrary targets.

## 1. Two branches, with a local inverse in the slot profile

Write `Z_U(r)` for the number of assignments to r labelled path
slots in a tournament U; slots may be empty, and each nonempty slot
has one directed Hamiltonian path on its assigned vertices. Put
`Z_U(0)=0` for nonempty U. By the odd-cycle collection formula (the
principal specialization of the Rédei–Berge function), `Z_U(q)`
as a polynomial for an order-n host has nonnegative integer
coefficients and exponents in `n,n−2,…`, ending at 1 for odd n
or 2 for positive even n.

Add two outer ports x and y with `x → U → y`. If `x → y`,
call the extension M; if `y → x`, call it S. Thus the two branches
differ at a single actual arc, and the port pair is the same in each.
These operations act on the slot polynomial `f=Z_U` by

```text
(M f)(q)=q² f(q);
(S f)(q)=q² f(q)+q[f(q+1)−f(q−1)].
```

Indeed, in different slots x and y do not constrain the old assignments.
In the same slot the ordinary x-before-y block survives in S iff its
old section is nonempty. The y-before-x orders in S split the old
section into the parts before y and after x, giving one extra labelled
slot. Counting the q common-slot choices proves the displayed
identity. The S extension is strong for nonempty U. Both operations
preserve self-converseness: extend an anti-involution of U by swapping
the two outer ports. These operators are inherited from the
[shell follow-up](../chanpak_matching_folds_branching_20261122/SHELL-INFORMATION.md).

On nonzero even slot polynomials from nonempty even-order hosts, M raises
the order at `q=0` by 2, while S resets it to exactly 2: the
coefficient at `q²` is `2f′(1)>0`. On odd profiles, S instead
resets the order to 1, with coefficient `2f(1)>0`, since
`f(q+1)−f(q−1)=f(1+q)+f(1−q)`. Both operators are injective
on polynomials: M is a degree shift, and the highest degree of
`(Sf)(q)` is two more than the highest degree of f, with the same
leading coefficient, since `q[f(q+1)−f(q−1)]` has degree at most
deg f. Consequently the last branch is read at `q=0`, and
the known branch can be stripped by exact triangular elimination.
Starting with the transitive two-vertex seed, every binary word of
length m therefore gives a distinct profile at order `2m+2`, with
an explicit inverse. The rest of this note converts this local
polynomial inverse into a single genuine count readout.

## 2. A single strong observation separates all fixed-order profiles

**Theorem.** For each n≥3 there is a common outer-port wrapper
W adding `O(n³ log n)` vertices such that from the single integer
`H(W(U))` one can recover the full slot polynomial `Z_U` of
**any** order-n host U. All output tournaments are strong. If U is
self-converse, so is `W(U)`. For n=1,2 the profile is fixed.

Let `R=ceil(n/2)`, `k=R−1`, and let t be a positive
integer specified below. Define `W=(S∘M^t)^k`, read as iterated
physical outer-pair extensions. All these operations retain and
reuse the outer-pair interface; an active S step is performed last.
The order bound will follow from `k=O(n)` and `t=O(n² log n)`.

View an evaluation of S at positive integer r as a signed sum of
three reads of its input:

```text
r → r   with factor r²;
r → r+1 with factor r;
r → r−1 with factor −r.
```

Since M has factor `r²`, expanding the k stages from the final
read `r_0=1` gives

```text
H(W(U)) = sum_(j=1)^R C_j Z_U(j),
C_j = sum_(gamma ending at j) alpha_gamma A_gamma^t,
A_gamma= product_(i=1)^k r_i².
```

Here `γ=(r_0,…,r_k)` is a positive integer walk with
steps in `{−1,0,1}`, and `α_γ` is the product of
the signed transition factors above. Walks touching zero vanish,
because their M multiplier at zero is zero (t>0). In particular
`α_γ` is a nonzero integer on each retained walk.
All coordinates are at most R, and the total absolute coefficient
weight `sum_γ |α_γ|` is at most `C=3^k R^(2k)`.

**Dominating walk.** For any `1≤j≤R`, the walk

```text
r_i*=min(1+i, j+k−i)  (0≤i≤k)
```

is the unique walk which is coordinatewise maximal among all walks
ending at a value at most j. The initial-point and endpoint
step bounds give the two displayed upper bounds, and their minimum
is a positive walk with allowed steps and ending exactly at j.
If `B_j=A_(γ*)`, then every other such walk lowers at least
one coordinate and raises none, so its multiplicative ratio is
at most `(1−1/R)²≤exp(−2/R)`.

Every `Z_U(j)` lies in `[0,E]` for `E=n! R^n`: choose arbitrary
vertex-to-slot assignments and then count no more than n! possible
orders across the slots. Take for instance

```text
t = R ceil(log_2(8 C E)).
```

Then `CE exp(−2t/R)≤1/8`. The term `B_j^t`
with nonzero integral factor from the dominating walk cannot be
erased by all other walks ending at j:
`|C_j|≥7B_j^t/8`. Moreover, for each j, all contributions from
lower endpoints have total absolute size at most `B_j^t/8`.
Crucially this comparison is for a fixed physical readout and a
uniform bound for *all* possible host profiles, not a profile-specific
choice after the difference is known.

**Decoding.** Begin at `j=R`. Round (H(W(U))/C_R) to the
nearest integer; the result is exactly `Z_U(R)` because lower
contributions are less than `|C_R|/2`, including when `C_R<0`.
Subtract `C_R Z_U(R)`, repeat at `j=R−1`, and continue down.
The constants `C_j` are computable without enumerating walks,
by a k-step signed tridiagonal transfer from 1; its transition
into state r is multiplied by `r^(2t)`.
Finally, values at `1,…,R` determine `Z_U`: divide by q
for odd order, or by `q²` for even order, and interpolate a
polynomial in `q²` of degree at most `R−1`. All arithmetic
is exact on integers or rationals of polynomial bit length in n.

The number of added vertices is `2k(t+1)=O(n³ log n)`, as claimed.
The last S gives strongness independently of the earlier neutral
outer pairs. ∎

## 3. Genuine branching and the local-to-height lesson

Starting with a transitive two-vertex seed, add a free sequence of
m pairs, each either M or S. Then apply the common wrapper for
`n=2m+2`. The `2^m` ordinary Hamiltonian counts are now
pairwise distinct and exactly decodable, realized by strongly
connected self-converse tournaments at one common order
`O(m³ log m)`. In particular the first digit, and each subsequent
digit conditioned on preceding digits, give disjoint images at the
same order. This is not the unrestricted odd-integer spectrum:
the padding is expensive and the images need not be intervals.
For a prescribed width, the decoder is a recognizer: on an arbitrary
integer, recover candidate evaluations, interpolate, and strip exactly
m branches. Accept only if the permitted seed is reached and the
re-encoded count agrees with the input. This inverts targets inside the
constructed images, not outside them.

One can cover both odd residue sheets with equally many codewords, at
one common order. For odd starting order take either the transitive
triple (H=1) or the directed triangle (H=3), then all m branch words
and the same wrapper for n=2m+3. The inverse recovers the word and the
seed. Both M and S preserve H modulo four for any host of order at
least two. M preserves H; S adds Z_U(2), divisible by four because
its empty odd-cycle collection contributes 2^n and any nonempty
collection contributes at least 2 times a positive power of 2.
Thus each seed generates 2^m distinct counts on its own residue
sheet. This is disjoint mod-four *branching*, not completeness on
either sheet.

There is also a literal separation of the fixed and paired channels.
On an odd-order host with fixed center c, either outer extension
preserves the reversal-fixed Hamiltonian-path count F. The new pair
cannot be internal in the first half (an incoming old arc would
force y, an outgoing arc would force x); it cannot be last there
(the forced y receives its arc from c). Thus it must be first,
contributing x followed by an old fixed half-path, a bijection.
The choice of its internal arc has no effect. Consequently all
2^m words on each of the two three-vertex seeds above have **F=1**,
including throughout the final wrapper, and their separations in
ordinary H are carried entirely by the number R=(H−F)/2 of two-element
path orbits. In fact this works with **any** prescribed positive odd
F: Session29 realizes F at odd order, and THM-5136(Q)(c) allows either
seed value of H modulo four by choosing the internal mate arcs,
which no fixed half-path can use. For each resulting seed the same
construction gives binary branching of reflected-pair counts with
F held fixed. The final active shell retains strongness.

For a coherent repeatable check, the companion script tests literal
subset dynamic programs on small actual tournaments and separately
tests the exact padded decoder through eight free port pairs. The
literal checks use small pads and do not assert that the loose general
padding constant is necessary.

**Analogy, not proof dependency.** [Pasten–Vidaux, §4](https://arxiv.org/pdf/2611.24748)
uses the formal-group leading term to arrange simultaneous local
approximation, then compares valuation information against global
height (§§2,5–7). Our local read is q-adic order of a slot
polynomial; the exact global read is recovered only after the
walk-weight gap beats the expressly bounded integer coefficient
height. Here the gap is real and the argument unconditional.
The [RH/GRH announcement remark](https://alpo.ge/remark.pdf)
emphasizes iterating lifts to gain multiplicity; this is an
analogy for repeating the same observation interface rather than
assuming that a single unweighted shell separates all profiles.
Nothing about the announced proof of RH/GRH is needed or asserted.

Reproduce the new finite checks with `python3 check_shell_branches.py`.
