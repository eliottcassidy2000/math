# A complete equality graph can recover the discarded clock

Let V be a separately verified integer potential for the height-adaptive
E,R graph, with

    V(t)>=V(s)+8*c(s,a,t)-5*a,
    V(entry)>=8*S_E(1^k a)-5*(k+a).

Assume all retained potentials are finite, and fix the terminal state
(z,m)=(0,1). Write V(terminal)=-J. A non-one actual word with endpoint
identity modulo5^E and H_E=J/8 must trace a path on which every displayed
inequality is equality. This follows by telescoping the nonnegative
integer slacks; it is not an inference from one extremizing path.

## All entry phases, including omitted long prefixes

At its first non-one letter a, write the preceding ones length as
k=b*T_ell+s, with0<=s<T_ell and ell=max(R,h(z_entry)). The entry state
is unchanged on replacing k by s. The floor-sum entry argument gives

    S_E(1^k)-S_E(1^s)<=lambda_E*b*T_ell.

Consequently its integer excess decreases by at least

    (5-8*lambda_E)*b*T_ell=b*T_ell/5^(E-1)>0

when b>0. Since the excess difference is an integer, it decreases by at
least one. Such an entry cannot lie on an equality path. Thus only the
least phase representatives used by the certificate need be enumerated.

If Vmin and Vmax are the certified range, every letter a>=65 has
transition excess at most Vmax+8E-325, and every omitted entry has excess
at most8E-325. When both are strictly below Vmin, no such letter is tight.
The strict version is required here; non-strict tail validation alone
would not justify removing a potentially tight omitted edge.

## Exact reverse edges

At target(z,m), for each1<=a<=64, an incoming state necessarily has

    z0=(2^a*(z-1)+2)/3 mod5^E.

Let M0=5^max(R,h(z0)), M=5^max(R,h(z)). Its possible retained multipliers
are exactly all unit m0 modulo M0 satisfying

    m0=m*2^a/3 mod min(M0,M).

Every such state is tested against the equality

    V(z0,m0)+8*min(h(z),v5(m-1))-5*a=V(z,m).

This includes every compatible lift, rather than picking one inverse.
Starting at the terminal and repeatedly taking every equality predecessor
therefore gives the complete backward equality graph. Each canonical entry
is tested independently against its exact floor-sum charge.

If the resulting finite graph is acyclic, topological dynamic programming
computes the maximum original clock over every equality path. Each entry
starts with clock k+a, and each edge adds its actual letter a. Acyclicity
is checked, not inferred from a visually short witness. If it fails or a
resource guard stops enumeration, no complete maximum is reported.

## Semantic exclusion

Suppose the complete graph has maximum clock K. Every actual equality word
then has clock A<=K, length1<=L<A, and multiplier3^L/2^A=1 modulo5^E.
A literal exclusion of all these finitely many pairs excludes actual
equality. The integer lattice then strengthens the actual bound to

    H_E >= (J+1)/8.

This conclusion is about actual words, although the input potential and
tight graph admit artificial lifted paths. It uses the restored original
clock precisely where the abstraction had forgotten the multiplier's
higher digits. It does not assert that every surviving multiplier pair is
realizable, nor that a single abstract extremizer classifies equality.

The E7 pilot has84 states,83 edges,five tight canonical entries and
clocks26..42. The completed E10 refinement through two slack layers appears
in `NEAR_TIGHT_CLOCK_THEOREM.md` and closes E14 as proved in `E14_CLOSURE.md`.
