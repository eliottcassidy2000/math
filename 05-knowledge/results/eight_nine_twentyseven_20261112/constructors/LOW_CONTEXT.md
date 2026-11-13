# Exact low uniform-context spectrum and a primitive source with a hole

**PROVED corollary of inherited active/inactive module bounds.** The new
application distinguishes the actual positive context spectrum from its
generated integer ideal. It does not establish a new hole in the unrestricted
tournament spectrum or resolve Chan–Pak.

## 1. Recovered theorem and the exact initial segment

Let M be any nonempty tournament, let h=H(M), and let F_M(2) count ordered
spanning two-slot path covers allowing either slot empty. Write

    L(M)=h+F_M(2).

Let S be the Hamilton-path count spectrum of all nonempty tournaments, and
let U(M) be the spectrum of all uniform contexts containing the fixed module
M, including the empty exterior. Then

    U(M) intersect [1,L(M))
       = (h S) intersect [1,L(M)).                      (1)

Moreover L(M) is attained by the cyclic triangle wrapper
`M→a→b→M`, and is the **minimum** count among contexts in which the
collapsed module vertex lies in a nontrivial strong component.

This is a direct spectral consequence of
[Session142, path-cover amplification, section3](../../session142_tournament_odd_spectrum_20261102/path_cover_amplification/PROOF.md),
which already proves both cases and the sharp bound. The underlying
maximal-run criterion and floor are in
[Session136, gluing, sections1 and3](../../session136_full_tournament_spectrum_20261102/gluing/PROOF.md).
The following short proof fixes the exact quantifiers and conventions; no
new priority is claimed for that mechanism.

Collapse M to one quotient vertex c. There are two exhaustive cases.

* If c is a singleton strong component, the ordered strong-component
  decomposition has every earlier component beating M and every later
  component losing to M. Every Hamilton path visits these components in
  order, so `H(E[M])=h H(E)`. Thus the multiplier is an actual tournament
  count. When E is empty the multiplier1 is also realized by a nonempty
  transitive tournament.
* If c lies in a nontrivial strong component, it belongs to a directed
  triangle. Indeed, along a directed cycle through c, there is a transition
  from an out-neighbour of c to an in-neighbour; those two vertices and c
  form the triangle. Its lift is an induced `M→a→b→M`. The wrapper has
  count `3h+C_2(M)=h+F_M(2)`: three one-run arrangements contribute3h,
  and the only two-run arrangement is `M_1,a,b,M_2`, counted by C_2(M).
  Hamilton counts are monotone under induced tournament extension, so
  every such context has count at least L(M).

For the monotonicity just used, deterministically insert a new vertex in
a legal position of each old Hamilton word. A legal position exists by
the usual first/last or consecutive transition argument. Deleting the new
vertex recovers the old word, so the map is injective; iterate over the
added vertices. This gives the lower bound without a strict-inclusion
assumption.

These cases prove one inclusion in (1). Conversely, for every actual
tournament V, the ordinal join `M→V` is a uniform context with H equal
to hH(V), proving the other. The triangle wrapper itself attains the
threshold. No classification of equality cases in larger contexts is
needed or claimed.

## 2. The two eight-vertex examples

For the original asymmetric source C in [PROOF.md](PROOF.md),

    h=513,    F_C(2)=7128,    L(C)=7641.

Since7 is not an unrestricted tournament count and
`7·513=3591<7641`,

    3591 is not in U(C), although 27 divides3591.        (2)

Thus even a positive odd member of its exact ideal need not be realized
by a uniform context. This is an actual all-context absence theorem,
not a missing value in a finite menu.

The primitive escape C4 has

    h=517,    F_C4(2)=7144,    L(C4)=7661,
    uniform-context ideal generator1.

Nevertheless

    3619=7·517 is not in U(C4).                         (3)

Primitivity is therefore insufficient for actual positive context
coverage, even for the concrete primitive source obtained by one of the
two arc switches. The proof applies to arbitrary exterior size and arcs.

Using the inherited attained spectrum below15, equation(1) gives the
entire initial segments explicitly:

| module | all uniform-context counts strictly below its active threshold |
|---|---|
| C | 513,1539,2565,4617,5643,6669 |
| C4 | 517,1551,2585,4653,5687,6721 |

These are respectively h times `1,3,5,9,11,13`. The next odd multiplier15
already exceeds each displayed threshold. The positive multipliers have
small explicit sources: singleton; directed triangle; `(2,4,9,3)` for5;
the ordinal join of two directed triangles for9; `(24,1,3,6,14)` for11;
and `(18,8,3,5,14)` for13. Ordinal join with the fixed module realizes
every table entry. The all-order omission of7, and of21 used below, is
inherited from **THM-1370**,
[h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
Only its proved omissions and small positive spectrum are used, not the
open claim that these are the only omissions.

## 3. Size-dependent inherited holes

For M of order n, cut each Hamilton word in any of its n+1 positions,
including the empty prefix and empty suffix. Concatenation recovers the
word and the cut, so these `(n+1)h` objects inject into its ordered
two-slot covers. Therefore

    F_M(2)≥(n+1)h,       L(M)≥(n+2)h.                  (4)

This is precisely the inherited path-cut incidence bound in Session142;
its even-order parity sharpening is unnecessary here.

More generally, if an integer q is absent from S and `qh<L(M)`, then
`qh` is absent from U(M). In particular the following order-only conditions
are sufficient:

    n≥6  => 7H(M) is not a uniform-context count of M;
    n≥20 => 21H(M) is not a uniform-context count of M. (5)

The exact criterion `qh<L(M)` can be stronger for a given source. The
orders in (5) are sufficient bounds, not asserted minimal. None of these
claims says7H(M) or21H(M) is absent from the unrestricted spectrum:
a different source or a nonuniform construction can realize such values.

## 4. Consequence for the current constructors

The context ideal controls signed integer combinations of responses.
Equation(1) supplies a separate, exact positive threshold. It explains
why the primitive C4 exit and the modulo8 selector are usable resources
without themselves constituting an arbitrary-target proof.

The root-child return in [PROOF.md](PROOF.md), section4, crosses the
required boundary explicitly: its first child sees the selected vertex
of C differently from the other seven vertices. Hence the original C
ceases to be a uniform module. That physical exception is outside the
scope of (1), rather than a violation of it. Future exact-target synthesis
must similarly retain the selected incidence or permit actual source
replacement; another uniform wrapper cannot remove the low-spectrum
holes just proved.
