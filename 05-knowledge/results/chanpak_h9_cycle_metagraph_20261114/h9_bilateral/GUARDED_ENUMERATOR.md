# A guarded exact labelled and unlabelled lift from paths to cycles

A **core** D is a nonempty tournament whose first and last
strongly connected components are nonsingleton. It need not
be strong; a single nontrivial strong component is allowed.
A **marked core** has an apex outcontact set U, considered up
to the automorphism group of D. Take the ordinal puncture
TT_l → D → TT_m with both l,m≥2. Add one apex. Require the
first two tape vertices to be beaten by the apex and the last
two tape vertices to beat the apex; all other tape contacts
are independent binary choices.

1. The apex is the unique vertex whose deletion leaves the
   tournament non-strong. After deleting any other vertex
   the remaining puncture still has a Hamiltonian path (by
   the usual tournament insertion argument), whose extreme
   surviving singleton vertices have their closing contacts
   by the two guards. Adding the apex closes it to a spanning
   cycle, witnessing strong connectivity.
2. The apex deletion and its ordered SCC factorization recover
   both tapes and the entire marked core. Core marks are
   identified *only* by core automorphisms. Conversely every
   core automorphism preserving U extends by fixing the
   apex and the two tapes. Hence the resulting automorphism
   group is exactly Aut(D,U).
3. The number of Hamiltonian cycles, counted up to rotation,
   is **H(D)**. Deleting the apex bijects such cycles with
   puncture Hamiltonian paths: every such path starts/ends
   at the extreme singleton vertices, and therefore closes.
   The two ordinal transitive tapes do not change H(puncture).

Write m(D) for the number of Aut(D)-orbits of all subsets U.
Use z for order and v to record the cycle count. For this
recognizable chart (not the class of all tournaments), the
exact ordinary isomorphism-class series is
$$
 \frac{z^5}{(1-2z)^2}
 sum_{[D]	ext{ trimmed core}}m(D)z^{|D|}v^{H(D)}.
$$
Indeed each tape has generating series z²/(1−2z), and the apex
contributes one vertex. The labelled exponential generating
series has the same expression with m(D) replaced by
2^|D|/|Aut(D)|: sum 1/|Aut(D,U)| over U-orbits and use
orbit-stabilizer. No assertion that a general core is rigid
is needed.

This can also be written in terms of just strong core classes.
Use a formal variable [h] with [h][k]=[hk]; at the end replace
[h] by v^h *as a linear operation*, not before multiplying.
Let A=Σ_C m(C)z^|C|[H(C)], summing over nonsingleton strong
tournaments C. Then the core series equals
$$
        A+Aigl(1-2z-Aigr)^{-1}A.
$$
An ordinal SCC word of length one contributes A; if longer,
its two ends contribute A and its internal word may use
either a nontrivial strongly connected marked block or a
singleton with one of two contact bits. Automorphisms
cannot permute different SCC positions since they form a
strict linear order. Path counts multiply because every
Hamiltonian path runs through the components in order.
The same word construction works for the exponential series
by using 2^|C|/|Aut(C)| as the atom weight.

Converse acts by reversing the tapes, swapping their roles,
and complementing contacts under a core anti-isomorphism.
A self-converse chart therefore needs equal widths and its
two addresses equal under the ordinary-left/reverse-
complement-right convention. It also needs a core
anti-isomorphism exchanging U with its complement,
impossible when |D| is odd. These constraints apply
inside this canonically recognizable chart.

## Exact converse and cycle-marked enumeration in the chart

For a self-converse trimmed core D, choose any anti-isomorphism a.
On its subset-orbit classes the involution sends U to
a(D\U). It is independent of the chosen a, since the ratio of
two anti-isomorphisms is an automorphism. Write m^+(D) for its
number of fixed orbits. The ordinary self-converse class series
of the chart is exactly
$$
  \frac{z^5}{1-2z^2}
       \sum_{[D]\ {m self\mbox{-}converse}}
                     m^+(D)z^{|D|}v^{H(D)} .
$$
There are 2^(l−2) left addresses of width l obeying the guard.
The right address is then forced to agree, and the core test is
exactly the fixed-orbit condition. Odd cores have m^+(D)=0.
Converse maps the entire chart to itself, so the merged class
series is half the sum of the ordinary and self-converse series.

Every automorphism fixes the unique apex. The induced action on
directed Hamilton cycles is therefore free: a rotation of a
cycle fixing the apex is the identity on all vertices. Thus a
marked graph contributes H(D)/|Aut(D,U)| unrooted cycle-marked
classes, not H(D) unless rigid. Starting from the EGF by labelled
graphs above, applying v∂/∂v produces the ordinary unrooted
cycle-marked series (the coefficient weights coincide by
orbit–stabilizer). Then z∂/∂z gives the ordinary rooted
closed-cut series, since the cycle stabilizer is trivial.
This makes the passage from a core path count to metagraph
presentations exact within the chart.

For any attained path count h>1, trim the initial and
terminal singleton SCCs of a witness; the surviving core
has the same H. Thus the construction converts that path
count into the same cycle count on a canonically rooted
family at every sufficiently large order, with exponentially
many labelled or unlabelled variants as above. For h=1
one can close an ordinary transitive spine by one apex to
get a one-cycle witness; the empty core falls outside this
unique-split enumerator. This transfer does not prove new
values occur as path counts, and a cycle count can be
realizable outside the image of the chart.
