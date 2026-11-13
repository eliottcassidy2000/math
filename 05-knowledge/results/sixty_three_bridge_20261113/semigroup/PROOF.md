# Activation ideals: what 63 repairs, and what localization forgets

**Status: PROVED elementary algebraic statements and source-bound tournament
applications.** The prime/seven-prime reduction, the arithmetic atoms, and
the smooth-sector constructions are **INHERITED**, not new results. The
countermodels below are abstract multiplicative monoids, **not** asserted
alternative tournament spectra. Chan–Pak completeness remains **OPEN**.

## 1. The inherited mechanism and its precise boundary

Let \(\mathcal O\) be the positive odd integers, and write

\[
\mathcal H=\{H(T):T\text{ a finite tournament}\},\qquad
M=\mathcal O\setminus\{7,21\}.
\]

The actual spectrum \(\mathcal H\) is a multiplicative monoid: the ordinal
join has a literal bijection
\(\operatorname{Ham}(T\triangleright U)
\cong\operatorname{Ham}(T)\times\operatorname{Ham}(U)\).
It contains \(1,3,5,35,49,63,343\), and is contained in \(M\).

The source dependencies are:

* [THM-4094, §§5–7, matching deficit and two-prime-lane completeness](../../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md): exact multiplication, the two infinite lanes, and the finite carry sources \(49,63,343\).
* [Session193, prime-atom induction, §2](../../session193_triangulated_proof_20261105/prime_induction/PROOF.md): the exact arithmetic atoms of \(M\); a least missing allowed value must be such an atom.
* [Session241, provenance, §3](../../session241_prime_relations_20261110/provenance/README.md): the older \(\{3,7\}\)-smooth monoid and its factorization laws, credited there to Session160 and the ordered-response normal forms.
* [Session219, Ford lane](../../session219_fermat_totient_chanpak_20261107/ford/PROOF.md) and [Session266, smooth lane](../../session266_fiber_radix_20261111/smooth/PROOF.md): sharper component-count and smooth-factorization information is already available.
* [Session299, §§2–3](../../session299_three_and_two_holes_20261113/README.md): ordinal multiplication by three, local triangle substitution, and arithmetic divisibility are different operations. Its reduction to 3-free targets is inherited here.
* [THM-2500, finite-hole relative atoms, §§3–4](../../../../01-canon/theorems/THM-2500-finite-hole-relative-atoms-prime-rays-and-cover-bounds.md): the general divisor-cage/prime-ray theorem and sharp cubic bound on the remaining weak atoms already exist. The two-lane geometry is a specialization of that older result, not a new generic finite-hole theorem here.
* [Session275, normalized cofactor descent](../../session275_exact_selection_20261111/descent/PROOF.md): a physical proper-graft cofactor can be 21, although its normalized gadget packet need not be physically realizable. Section 6 gives a specific consequence at total 63.

The arithmetic atoms of \(M\) are

\[
\mathcal A=
\{p:p\text{ odd prime},\ p\ne7\}
\cup\{7p:p\text{ prime},\ p\ge5,\ p\ne7\}
\cup\{49,63,343\}.                                      \tag{1.1}
\]

They are indecomposable **values in \(M\)**. They are not the same thing as
strong tournament classes. Every realized atom requires one nontrivial
strong component; the converse fails, for example at a strong tournament
of count nine.

In particular \(63\) is an atom of \(M\), although it is composite as an
ordinary integer. A nontrivial proper factorization into odd integers uses
7 or 21, so cannot be an ordinal factorization into actual tournament
values. Thus a tournament of count 63 has no triangle SCC to remove, even
though \(v_3(63)=2\). The valuation records a numerical factor, not a
detachable physical component.

The concept board is: ordinal products / valuation addresses / activation
ideals / conductor multipliers / localization / physical division. The
canonical hostile is \(7,21\notin\mathcal H\), \(63\in\mathcal H\).
The unused sidecar developed below is the first **actual source** at each
minimal activation address.

This is a different process from the prior
[strict synchronous product-saturation classification](../../multiplicand_saturation_chanpak_20261111/classification/PROOF.md).
That process forbids equal-valued parents and counts rounds. Ordinal joins
allow separate copies of the same tournament, and this note concerns
existence without a round budget. Its rank-one saturation criterion must
not be imported into the present monoid.

## 2. A finite activation certificate for any fixed multiplier bank

Let \(S\subseteq\mathcal O\) be any multiplicative monoid. Choose attained
multipliers \(q_1,\ldots,q_d\in S\), all greater than one, and an arbitrary
positive odd target \(b\), which need not belong to \(S\). Define

\[
I_b=\left\{e\in\mathbb N^d:
          b\prod_{i=1}^d q_i^{e_i}\in S\right\}.          \tag{2.1}
\]

**Theorem 2.1 (activation ideal and finite-box certificate).**

1. \(I_b\) is upward closed in the coordinatewise order. It has finitely
   many minimal addresses, possibly none, and is their upward closure.
2. The complement \(\mathbb N^d\setminus I_b\) is finite if and only if
   each coordinate axis has an activation: there are \(r_i\ge0\) with
   \(r_i e_i\in I_b\), where \(e_i\) is the coordinate unit vector.
3. Given such activations, every unresolved address is in the explicit box
   \(0\le e_i<r_i\). Positive witnesses for the axis activations and exact
   membership certificates inside that box determine all of \(I_b\).
4. If \(S=\mathcal H\), retain an actual tournament for every selected
   activation and for every \(q_i\). Every positive address above a selected
   activation then has an explicit ordinal-join witness. The construction
   begins from that attained activation, never from an absent \(b\).

**Proof.** Upward closure is multiplication by the known \(q_i\).
For finiteness of minimal addresses, every infinite sequence in
\(\mathbb N^d\) has a componentwise nondecreasing subsequence: in one
coordinate either some value occurs infinitely often, or an increasing
subsequence can be chosen; repeat for the finitely many coordinates.
Thus an infinite antichain is impossible. Every address in \(I_b\) lies
above a minimal one, by taking a minimal element of the finite box below
it. This proves item 1, the elementary form of the monomial-ideal finiteness
principle.

If the complement is finite, each infinite axis eventually meets \(I_b\).
Conversely, an address outside the displayed box has \(e_i\ge r_i\) for
some \(i\), and is above the corresponding axis activation. If some
\(r_i=0\), then \(0\in I_b\), the whole space is active, and the box is
empty. Finally, if \(a\le e\) is a certified active address, join its
tournament with \(e_i-a_i\) copies of the \(q_i\) witness for each \(i\).
The count is exactly \(b\prod q_i^{e_i}\). ∎

The address map need not be injective: the \(q_i\) may share primes or
satisfy multiplicative relations. This causes no problem for (2.1), but
different minimal addresses may encode the same integer. The theorem
neither gives an a priori bound on unknown minimal addresses nor makes
global nonattainment easy to certify. It is a finite certificate **once
the axis witnesses and the remaining finite-box decisions are supplied**.

Equivalently identify an upward set with its monomial ideal in
\(\mathbb Z[x_1,\ldots,x_d]\). Multiplying the base target obeys the exact
shift law

\[
I_{bq_i}=I_b:x_i
       =\{e:e+e_i\in I_b\}.                              \tag{2.2}
\]

This is an identity of target predicates. It is not a graph operation
which divides the count of a supplied witness.

## 3. The exact 7 → 21 → 63 activation boundary

Take the actual tournament spectrum, base \(b=7\), and multiplier bank

\[
(q_1,q_2,q_3)=(3,5,49).
\]

The inherited sources \(63,35,343\) certify respectively the addresses
\((2,0,0),(0,1,0),(0,0,1)\). Their axis box consists only of
\((0,0,0)\) and \((1,0,0)\), which represent the known forbidden totals
7 and 21. Therefore, at every order and all addresses,

\[
\boxed{I_7=(x^2,y,z).}                                   \tag{3.1}
\]

This is an exact theorem about actual existence, not an extrapolation from
small exponents. Its positive recipe is explicit:

\[
7\,3^a5^b49^c=
\begin{cases}
35\,3^a5^{b-1}49^c,&b\ge1,\\
343\,3^a49^{c-1},&b=0,\ c\ge1,\\
63\,3^{a-2},&b=c=0,\ a\ge2.
\end{cases}                                             \tag{3.2}
\]

All factors on the right are actual attained values. These formulas
repackage the inherited smooth construction; the added interface is the
finite activation certificate and the identity of the first physical
source used by each address.

The same multiplier bank and (2.2) give

\[
I_{21}=(x,y,z),\qquad I_{63}=(1).                         \tag{3.3}
\]

In the one-coordinate bank \(q=3\), the exact first activation of 7 is at
exponent two. For an arbitrary multiplicative \(S\) containing 3, define
on 3-free odd integers

\[
\tau_S(m)=\min\{k\ge0:3^km\in S\}\in\mathbb N\cup\{\infty\}.
\]

Then a fibre is either empty or an integer tail, and

\[
\tau_S(mn)\le\tau_S(m)+\tau_S(n).                        \tag{3.4}
\]

Conversely any such subadditive function with \(\tau(1)=0\) defines a
multiplicative monoid containing 3 by taking the corresponding tails.
This follows immediately by adding the two 3-adic exponents. For the
actual spectrum, \(\tau_{\mathcal H}(7)=2\). Completeness would additionally
require \(\tau_{\mathcal H}(m)=0\) for every other 3-free odd \(m\); that
requirement is still open. The activation of one fibre does not establish
those zero-delay claims.

## 4. Exact conductor and the loss under localization

For \(S\subseteq\mathcal O\), define its conductor into all odd integers by

\[
\mathfrak c(S)=\{a\in\mathcal O:a\mathcal O\subseteq S\}. \tag{4.1}
\]

It is an ideal for multiplication by arbitrary odd integers, and is a
subset of \(S\). For the candidate target monoid,

\[
\boxed{\mathfrak c(M)=\mathcal O\setminus\{1,3,7,21\}
                         =M\setminus\{1,3\}.}            \tag{4.2}
\]

Indeed a multiplier fails precisely if it divides 7 or 21; their odd
divisors are exactly \(1,3,7,21\). This computation concerns \(M\), not
the still-unknown full \(\mathcal H\).

For a monoid containing \(a>1\), write
\(S[a^{-1}]=\{s/a^k:s\in S,k\ge0\}\). Localization permits formal division
and can erase the activation delay. For example

\[
7=63/3^2,\qquad21=63/3
\]

belong to \(\mathcal H[3^{-1}]\) even though neither belongs to
\(\mathcal H\). Formula (2.2) records the same loss as repeated colon
shifts. A localized representation retains a numerator and denominator;
recovering a tournament with the unscaled target still requires a physical
operation that removes the denominator factor. A rational identity alone
does not provide one.

## 5. Universal hostile: infinitely many primitive holes survive all these tests

The following exact family shows how much stronger the missing physical
step must be. Let

\[
\mathcal A_*=
\{p:p\ge5\text{ prime},\ p\ne7\}
\cup\{7p:p\ge5\text{ prime},\ p\ne7\}.
\]

Choose **any** subset \(E\subseteq\mathcal A_*\), finite or infinite, and
put

\[
S_E=M\setminus E.                                       \tag{5.1}
\]

**Theorem 5.1 (deleted-atom countermodels).**

1. \(S_E\) is multiplicatively closed and contains \(1,3,49,63,343\).
2. Every positive odd multiple of 3 except 21 belongs to \(S_E\). In
   particular \(9\mathcal O\subseteq S_E\); also
   \(49\mathcal O\subseteq S_E\).
3. For **every** positive odd integer \(n\), both \(n^2\) and \(n^3\)
   belong to \(S_E\). Thus its power-saturation
   \(\{n\in\mathcal O:n^k\in S_E\text{ for some }k\ge1\}\)
   is all of \(\mathcal O\), and closing it under the square/cube rule
   \(n^2,n^3\Rightarrow n\) fills every hole in one step.
4. For every \(a\in S_E\), \(a>1\),
   \[
   S_E[a^{-1}]=\mathcal O[a^{-1}].                        \tag{5.2}
   \]
   More concretely, \(na^2\in S_E\) for every positive odd \(n\).
5. On 3-free odd inputs, the activation function is exactly
   \[
   \tau_{S_E}(m)=
   \begin{cases}2,&m=7,\\1,&m\in E,\\0,&\text{otherwise}.
   \end{cases}                                          \tag{5.3}
   \]

**Proof.** A product of two nonunits of \(M\) cannot be an arithmetic atom
of \(M\). Deleting any collection of its atoms therefore preserves
closure; identity multiplication causes no problem. None of the retained
finite sources is in \(\mathcal A_*\), proving item 1.

Every deleted value is 3-free, proving the first part of item 2. No deleted
value is divisible by 49: an ordinary prime is not, and \(7p\), with
\(p\ne7\), has seven-adic valuation one. The two original holes are not
divisible by 49 either. This proves the rest.

If \(n\in M\) is greater than one, its square and cube are products of
nonunits of \(M\), so are not deleted atoms. The exceptional input 7 has
square 49 and cube 343, both retained. For the input 21, its square and cube
belong to \(M\), have a factor 3, and are not in \(E\). The unit is
immediate. This proves item 3.

For item 4, if \(na\in M\), then \(na^2=a\cdot na\) is a product of two
nonunits of \(M\) and cannot be deleted. The case \(na=7\) is impossible
for \(a\in S_E\), \(a>1\). The only possible case \(na=21\) is
\(a=3,n=7\); then \(na^2=63\), which is retained. Hence
\(n/a^k=(na^2)/a^{k+2}\) lies in the left side of (5.2). The reverse
inclusion is tautological. Finally each deleted atom becomes a nondeleted
multiple of 3 at its first step, whereas 7 first enters at 63. This proves
item 5. ∎

The conductor retains more information than its mere nonemptiness. Put
\(P_E=\{p:7p\in E\}\). The exact formula is

\[
\mathfrak c(S_E)=\mathcal O\setminus
                 \bigl(\{1,3,7,21\}\cup E\cup P_E\bigr). \tag{5.4}
\]

This follows by listing the odd divisors of the missing values: a deleted
prime \(p\) contributes \(1,p\), and a deleted \(7p\) contributes
\(1,7,p,7p\). Thus full conductor equality with (4.2) would detect the
added holes; we do **not** claim otherwise. The hostile shows that even
the universal conductor multipliers 9 and 49, all square/cube attainments,
and all the displayed localizations fail to determine those holes.

For any prescribed finite correct prefix \(M\cap[1,R]\), enlarge \(R\)
to at least 7 if necessary, and choose

\[
E_R=\{p,7p:p>R\text{ an odd prime},\ p\ne7\}.            \tag{5.5}
\]

Then \(S_{E_R}\) agrees with that prefix and all its multiplicative
consequences, passes every property in Theorem 5.1, yet omits infinitely
many members of **both** primitive lanes. It also agrees with the entire
correct 3-divisible slice, not merely a finite portion of it. This is a
strictly scoped impossibility result for inference from those algebraic
properties. Actual tournament insertions and marked interfaces supply
extra structure not imposed on \(S_E\).

## 6. A physical normalization boundary specifically at 63

For a nonempty tournament gadget \(B\), put

\[
Q(B)=\sum_{X\subseteq V(B)}H(B[X])H(B[V(B)\setminus X]),
\qquad H(\varnothing)=1.
\]

For a rooted source \((A,a)\) with \(a\) not a sink, its **proper graft**
keeps the source and gadget internal arcs and adds
\(A-a\to B\to a\). Write \(N\) for its path count and
\(g=\gcd(H(B),Q(B))\). Session275's proved cofactor classification says

\[
N/g=21
\quad\Longleftrightarrow\quad
Q(B)=20H(B)\text{ and }A\text{ is transitive with }a
\text{ penultimate}.                                   \tag{6.1}
\]

**Corollary (PROVED, this proper-graft language only).** No proper graft
of count 63 has normalized cofactor 21.

Indeed (6.1) would give \(g=H(B)=63/21=3\) and \(Q(B)=60\).
Every tournament of count three is a single triangle SCC with \(t\ge0\)
singleton SCCs before and after it. Ordered two-slot covers multiply over
ordinal components; the triangle has \(Q=12\) and a singleton has \(Q=2\).
Consequently every such gadget has

\[
H(B)=3,\qquad Q(B)=12\cdot2^t,                           \tag{6.2}
\]

which cannot equal 60. ∎

This is stronger than merely noting that 21 is a forbidden ordinary count,
because a normalized cofactor of 21 is **not** universally forbidden.
Session275 gives an actual gadget with
\((H,Q)=(29,580)\) and a proper graft of count
\(609=21\cdot29\). Its normalized packet \((1,20)\) is nonphysical:
count one forces a transitive tournament, whose \(Q\) is a power of two.
Thus even a positive physical graft can have a forbidden normalized count.
The extra constraint at 63 is the small integral content \(g=3\), which
forces the incompatible packet (6.2). This corollary excludes neither
other constructions of 63 nor arbitrary occurrences of a gcd-three
observer; it concerns exactly the normalized cofactor in (6.1).

## 7. Operational consequence and exact remaining question

The positive finite-bank procedure is:

1. Specify attained multipliers and retain literal witnesses for them.
2. Construct actual sources at axis activation addresses.
3. Resolve the remaining finite box, keeping global nonattainment proofs
   distinct from bounded searches.
4. Compile every certified address by ordinal joins from one of its
   attained sources.

For \(7\) and \((3,5,49)\), this procedure is complete and explains why 63
is the first ternary activation. It is a compact witness compiler for an
infinite smooth family already known to exist.

What does **not** follow is an operation cancelling 3 from the first
attained member of an arbitrary fibre. The countermodels (5.1) pass even
the universal assertions that every permitted 3-divisible target is
attained and every odd square and cube is attained, while they can still
miss infinitely many 3-free primitive targets.

A useful new constructive question is therefore explicitly source-bound:
given a marked witness of count \(3p\) or \(21p\) in a specified family,
can one produce a witness of count \(p\) or \(7p\), with a checkable
interface guard, for every unseeded prime \(p\)? Such a cancellation
cannot follow from divisibility, localization, or SCC removal alone.
The strong source and its marked response data are the missing sidecar.
This question is **OPEN**, and this packet supplies no uniform cancelling
graph operation.

## 8. Controls

Run

```bash
python3 -B 05-knowledge/results/sixty_three_bridge_20261113/semigroup/check.py
```

The companion [check.py](check.py) and [receipt.json](receipt.json) use only
exact integer arithmetic. They compare bounded arithmetic atom tests,
deleted-atom membership, conductor divisors, square/cube membership,
localization numerators, ternary activation delays, and the explicit
finite-bank source recipes. They do not enumerate tournaments and do not
present any abstract countermodel as an actual spectrum. The all-parameter
claims follow from the proofs above, not the bounded controls.
