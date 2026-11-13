# An exact selector for an unbounded repair language

**PROVED reduction and termination; FINITE-EXACT complete head.** This is
an exact inverse for the guarded23-core prefix language. Membership in
that language is a sufficient condition for a tournament target, not a
necessary condition for the full Hamilton spectrum.

Define

    A23={23H(L)+60P_L(2): L nonempty, gcd(H(L),30P_L(2))=1}.

Each member comes with the [same-target preparation and indefinite chain
entry](REPAIR23.md). A selector must return L and the literal core/port
program; a pair of feasible integers(h,f) is not enough.

## 1. The exact reduction and its finite height

For any tournament L of order n, every ordered bipartition has at least
one path in each nonempty part. Also, splitting every Hamilton path at
each of its n+1 cuts gives distinct ordered path covers. Consequently

    f=P_L(2)≥2^n,       f≥(n+1)H(L).                  (1)

For target N=23h+60f, these give

    h≤N/143,
    n≤floor(log_2((N-23)/60)).                       (2)

N is odd and prime to15 under the guard. Since f is even,

    h≡47N (mod120),       f=(N-23h)/60,              (3)

because23·47≡1 mod120. These are inexpensive necessary filters, not
sufficient realizability conditions. In particular the least target is143,
attained by the singleton prefix.

Equations(2) yield a terminating exact algorithm for **any** input N:
enumerate all tournaments through the finite n bound, compute h,f, and
test equality and the guard. Return the first actual L, or a complete
rejection from A23 when the finite search exhausts. The bound is logarithmic
in N, but unrestricted tournament enumeration is not polynomial in N;
no fast arbitrary-N implementation or global existence theorem is claimed.

## 2. SCC multiplication compresses the search without losing responses

Write L=S1→...→Sj for its actual SCC order. Then

    (H(L),P_L(2))=product_i(H(Si),P_Si(2))             (4)

coordinatewise. Restricting the two paths to each component and
concatenating proves this, including empty slots. A singleton contributes
(1,2). A finite target bound permits only finitely many strong generators,
and each has a known order bound. Products failing gcd(h,30f)=1 can be
discarded permanently, since multiplication cannot remove a common factor.

For this observer, different SCC orders and singleton-gap placements can
have the same pair. One actual ordered representative suffices to build
the final target because the root is in the separate terminal23-core.
This quotient does not classify L itself or preserve arbitrary roots inside L.
Retaining one full prefix graph and its component word is the lifting data.

## 3. A complete implementation below26521

The bank contains every ordinary strong tournament through order7. Suppose
an admissible prefix contains a strong component S of order at least8.
Busch's inherited minimum gives H(S)≥45; the guard excludes multiples
of3 or5, so H(S)≥47. By(1), P_S(2)≥9·47=423, and this count is even,
so P_S(2)≥424. Other SCC factors cannot lower either coordinate. Therefore

    N≥23·47+60·424=26521.                            (5)

Thus below26521 the small strong bank and arbitrary products of its pairs
are **complete at all prefix orders**. The implementation enumerates the
bounded product monoid, retains actual component words, and independently
recounts every prefix and output graph. It includes the singleton generator
but excludes the empty prefix. No near-minimum conjecture is used in(5).

Run `python3 -B 05-knowledge/results/surplus_selector_20261113/prefix_selector.py`.
The receipt lists every accepted target and every missing odd target prime
to15 in this interval. Its negative answers are only statements about A23.
For example881 has no admissible h in(2)–(3), while881 is already known
to be a tournament count. This is a hostile control against mistaking the
repair language for the full spectrum.

## 4. What the root grammar still needs

The [atomic grammar](../ternary_core_grammar_20261113/GRAMMAR.md) gives an
inverse parse once an actual strong graph is supplied. Its exact surplus
identity needs a smaller parent and an actual cross matrix. Here (1)–(5)
give target-dependent finite search and reconstructible witnesses for one
unbounded sublanguage, with a proved endpoint continuation after success.
What remains open is a theorem that the union of such physically closed
languages contains every permitted odd target. Merely solving(3), finding
a primitive formal chart point, or covering all residues does not prove it.
