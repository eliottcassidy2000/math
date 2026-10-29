# Multiplier-resolved recurrence and a sharp port threshold

**PROVED scoped below; FINITE-EXACT checks in charged.py and jets.py.
Same-agent proof audit and alternate computations, not an independent review.**
Packet-local labels F1--F4 are not reserved global theorem IDs.

## Objects and inherited input

Let p>3 be prime, p=3 mod 4, h=(p-1)/2, and H the square multipliers in
F_p^*. On the p-ary rooted tree of depth k let A_k be the iterated affine
wreath action with root group {x -> ax+b : a in H}. Put A_0=C_0=1 acting
on one point. Let t negate every coordinate and C_k=C_(A_k)(t).
The underlying tournament is the lexicographic power of the Paley atom;
vertices are words, and their first differing coordinates determine the arc.
This is an actual tournament action, without a tie-breaking gauge.

[PORT_LAWS P4--P5](../complement_port_crt_20261028/PORT_LAWS.md) proves
that the square-root marking loop subgroup is

    Hol_k = ker(chi:C_k -> H^k),

where chi_j is the product of local multipliers at depth j. It also proves
surjectivity, and identifies the full tournament automorphism action for
the checked prime bases 7 and 11. For other primes the statements here concern
the specified affine wreath action; no larger full-Aut identification is used.
The p=3 case has H=1 and Hol=C, as inherited.

The uniform source is a group element; ports are ordered vertices with
repetitions permitted. Burnside identifies E fix(g)^r with r-port orbit count.
An arity r below refers to this observable, not to arbitrary legal contexts.

## F1. The full charged cycle inventory

Use group-algebra symbols e_u for u in (Z/h)^k and cycle variables x_i.
Let A_k and C_k in this section denote normalized inventories
E[e_chi(g) product_i x_i^(c_i(g))]. Define

    Psi_d: x_i -> x_(di),
    Theta_d: e_u -> e_(2u), x_i -> x_(di)^2.

Let e_a prepend the exponent of a in a chosen generator of H, and let
d(a) be its multiplicative order. With A_0=C_0=x_1,

    A_k = [ e_1 (A_(k-1)^p + (p-1) Psi_p A_(k-1))
             + p sum_(a in H, a!=1)
                 e_a A_(k-1) (Psi_d(a) A_(k-1))^((p-1)/d(a)) ] / (ph),

    C_k = (1/h) sum_(a in H) e_a C_(k-1)
                          (Theta_d(a) A_(k-1))^(h/d(a)).

Here e_1 prepends the identity charge, not exponent one. Products add
suffix charges modulo h. The ordinary Hol cycle index is

    Z_Hol_k = h^k [e_0] C_k.                              (F1)

Proof: along a root permutation cycle of length d the product of independent
child actions is uniform in A_(k-1), and a child cycle of length i becomes
one of length di. The zero child in C is independent. The nonzero root
children occur in h mirrored pairs; corresponding product actions are
conjugate by t, hence have identical cycles and multiplier characters.
They therefore duplicate cycle multiplicities and double charges. Since h
is odd this doubling is invertible. Surjectivity of chi gives the factor
h^k in (F1). This proves the recurrences including all multiplicities.

`charged.py` uses unnormalized integer inventories: a root cycle product
has |A_(k-1)|^(d-1) preimages. It independently enumerates all 83,349 literal
permutations of C_2 at p=7 and matches every charged cycle cell.

## F2. Only the deepest nontrivial Fourier character matters for fixed ports

For lambda in the dual of H^k define

    F_(k,lambda)(z) = E_(A_k) lambda(chi(g)) exp(z fix(g)),
    G_(k,lambda)(z) = E_(C_k) lambda(chi(g)) exp(z fix(g)).

These are formal exponential generating series. Write F_k,G_k for trivial
lambda. For nontrivial lambda let j be the deepest nontrivial coordinate,
numbered from root depth zero. The series depend only on k,j, not on the
other character coordinates or the particular nontrivial character there.
Write them f_(k,j),g_(k,j). There are h^j(h-1) such characters.

The ordinary recurrences, with F_0=G_0=exp(z), are

    F_k = [F_(k-1)^p + p(h-1)F_(k-1) + p-1]/(ph),
    G_k = G_(k-1) [F_(k-1)(2z)^h+h-1]/h.              (F2a)

The nontrivial recurrences are

    f_(k,0) = [F_(k-1)^p-pF_(k-1)+p-1]/(ph),
    g_(k,0) = G_(k-1)[F_(k-1)(2z)^h-1]/h,

    f_(k,j) = f_(k-1,j-1)^p/(ph),                     j>=1,
    g_(k,j) = g_(k-1,j-1) f_(k-1,j-1)(2z)^h/h.        j>=1.  (F2b)

Consequently the all-order Hol port EGF is

    B_k(z) = G_k(z) + sum_(j=0)^(k-1) h^j(h-1)g_(k,j)(z).   (F2c)

Proof: a nontrivial suffix character averages to zero on a uniform child
group. Every nonidentity root affine permutation has a nonfixed cycle,
so such a suffix kills its contribution to F. Only the identity root
remains. For a trivial suffix, summing a nontrivial root character over
H minus its identity gives -1, whereas the trivial sum is h-1. This gives
F2a and the first f recurrence. The same reasoning on mirrored cycles
gives the G recurrences. Doubling a character preserves its deepest
nontrivial coordinate because h is odd; induction proves the stated
collapse. Finally sum all characters to project onto chi=1.

The collapse is for fixed-point observations. F1 retains more cycle data;
neither formula identifies individual transporters. The inherited parallel
and mirrored C3 actions have the same full cycle index and different joint
transporter counts 1 and 0; charged.py actively retests them.

## F3. Exact first visibility at arities 1,p,p^2,...

Let K_(k,l) be the subgroup of C_k on which the first l multiplier
coordinates vanish, 0<=l<=k. Thus K_(k,0)=C_k and K_(k,k)=Hol_k.
Its port EGF is F2c truncated to j<l. For j<k, imposing the next coordinate
changes no r-port count when r<p^j, and changes every r-port count when
r>=p^j. In particular its first nonzero difference is

    orbit_(K_(k,j+1))(p^j) - orbit_(K_(k,j))(p^j)
      = (h-1) 2^(p^j) (p^j)! / (ph)^E_j,

    E_j = ( ((p^j-1)/(p-1)) - j )/2.                  (F3)

This number is independent of k as long as k>j.

Proof: F_k(0)=G_k(0)=1 and F'_k(0)=1. Put u=F_(k-1)-1.
The numerator for f_(k,0) is sum_(i=2)^p binom(p,i)u^i;
it has nonnegative coefficients and leading term binom(p,2)z^2.
Therefore f_(k,0)=z^2+O(z^3), and g_(k,0)=2z+O(z^2).
Induction in F2b gives

    ord_z f_(k,j)=2p^j,
    [z^(2p^j)] f_(k,j)=(ph)^(-((p^j-1)/(p-1))),
    ord_z g_(k,j)=p^j,
    [z^(p^j)] g_(k,j)=2^(p^j)/(h^j (ph)^E_j).

All coefficients at and above the leading degree are positive. Multiplying
by h^j(h-1) and by r! proves F3 and the strictness claim. This positivity
is derived after character cancellation; it is not an assumption that
different Fourier contributions cannot cancel.

For p=7 the first gaps are 4 at one port and 1,290,240 at seven ports.
For p=11 the second gap is 326,998,425,600 at eleven ports. The full rational
jets verify the next threshold 49 for p=7 at both depths three and four.
For any fixed arity cutoff R only j<=floor(log_p R) contribute. This is an
exact way to omit deeper charge information in this counting task.

## F4. Low-order consequences and attachment use

For k>=1 the one- and two-port Hol counts are

    b_(k,1)=p+2k-2,
    b_(k,2)=8k^2+(10p-30)k+p^2-10p+22.                (F4)

Only the j=0 character block can contribute through arity two, giving F4
by differentiating F2. For p=7 these start (7,49),(9,113),(11,193).
The occurrence of 113 has no established map to the Borwein denominator
113; RECOVERY records the distinct reasons for those values.

For a fixed reversal marking, q choices of an attachment profile per
vertex give Z_C(q,q,...) or Z_Hol(q,q,...) orbits, according to whether all
centralizer maps or only loop-generated corrections are allowed. At p=7,
k=2, binary profiles give respectively 8,951,645,952 and 61,972,938,752.
These are marked-profile counts, not A000568 values or a reconstruction of
all general strong atoms. Actual native complements still need the ordered
path data and legality guards from the inherited composition carrier.
