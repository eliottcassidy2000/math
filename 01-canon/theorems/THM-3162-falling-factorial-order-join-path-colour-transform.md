---
id: THM-3162
title: "Falling-factorial order-join path-colour transform"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The complete
  spanning path-cover profile of a tournament becomes multiplicative under
  order-join after the falling-factorial transform
  Q_T(t)=sum_d pc_T(d)(t)_d.  Equivalently, order-join is the bipartite
  matching convolution on path counts.  Every fixed-depth profile along an
  arbitrary SCC chain is an exact finite difference of products Q_i(j), and
  repeated joins give finite exponential sums and constant-order recurrences.
  The transform forgets SCC order and does not extend to cyclic substitution
  or make arbitrary tournament path counting easy.
audit: >
  The canonical subset Hamilton-path/path-cover engine checks all 1,099
  labelled tournaments through order five, all 5,625 ordered pairs through
  order four, SCC factorization, 210 Laguerre kernels, and 4,275 repeated-join
  coefficients.  An independent permutation-cut engine checks 761 direct
  joins through selected order eight, 528 ordinary-power falling-factorial
  identities, fixed-depth exponential forms, and both hostile boundaries.
  Normal and optimized transcripts byte-match.
source: root/frontier-synthesis/tournament-sequence-transform/2026-08-02
depends_on:
  - THM-1862-order-join-reduction-principle
  - THM-3121-path-cover-walk-content-substitution-kernel
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
related:
  - THM-2412-delta-exponential-and-central-newton-layer-split
script: 04-computation/tournament_order_join_falling_factorial_transform_thm3162.py
output: 05-knowledge/results/tournament_order_join_falling_factorial_transform_thm3162.out
script_sha256: ceb1e3901e12940519ec15f45c325386b38b17db853a6f0716e23404ca9ac7fa
output_sha256: 00f59a138654231ff257b1e9f2e0a7d84010b1b55f27309930a81e090f905aa4
independent_script: 04-computation/tournament_order_join_falling_factorial_independent_audit_thm3162.py
independent_output: 05-knowledge/results/tournament_order_join_falling_factorial_independent_audit_thm3162.out
independent_script_sha256: b8b2b7ff8b61470e844e2e0130b842b64cde9e16aa4425ba9600b2a5d32dc564
independent_output_sha256: 560a57a1a42651d9936f4c0e143406a51b79fc442fd217f0df757541aeaa1c94
hash_basis: LF-normalized bytes
---

# THM-3162 -- falling-factorial order-join path-colour transform

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Multiplicative path-colour transform

For a nonempty tournament (T), let

\[
 p_T(d)=\#\{\text{spanning covers of (T) by (d) unordered directed
 paths}\}.                                                   \tag{1}
\]

Thus (p_T(1)=H(T)), the Redei Hamiltonian-path count.  Put (p_T(0)=0)
and define the **injective path-colour polynomial**

\[
 Q_T(t)=\sum_{d=1}^{|T|}p_T(d)(t)_{\underline d},\qquad
 (t)_{\underline d}=t(t-1)\cdots(t-d+1).                    \tag{2}
\]

For disjoint tournaments (A,B), write (A\triangleright B) for the
order-join in which every vertex of (A) beats every vertex of (B).  Then

\[
 \boxed{Q_{A\triangleright B}(t)=Q_A(t)Q_B(t).}              \tag{3}
\]

More explicitly, the untransformed profile satisfies

\[
 p_{A\triangleright B}(d)=
 \sum_{\substack{a,b\ge1,\ 0\le k\le\min(a,b)\\a+b-k=d}}
 p_A(a)p_B(b)\binom ak\binom bk k!.                         \tag{4}
\]

Equations (3)--(4) are all-order operation-response laws.  THM-1862's
Hamiltonian multiplicativity is their depth-one slice, not a separate
phenomenon.

## 2. Coloured-cover bijection

For every nonnegative integer (m), (Q_T(m)) counts spanning directed-path
covers whose components receive distinct colours from an (m)-colour set.
Indeed, a (d)-component cover has exactly ((m)_{\underline d}) injective
colourings.

Choose independently an injectively coloured cover of (A) and one of
(B).  If a colour occurs on both sides, concatenate the unique (A)-path
of that colour with the unique (B)-path of that colour.  The joining arc is
always directed from (A) to (B).  Colours occurring on only one side are
left unchanged.

Conversely, a directed path in (A\triangleright B) crosses the block cut at
most once: after it enters (B), it cannot return to (A).  Splitting every
coloured component at that cut recovers the two independently coloured
covers.  This is a bijection, so (3) holds for every integer (m\ge0), hence
as a polynomial identity.

If the two uncoloured covers have (a,b) paths and exactly (k) colours are
shared, the shared paths can be selected and paired in

\[
 \binom ak\binom bk k!                                      \tag{5}
\]

ways, and the joined cover has (a+b-k) components.  Forgetting the colours
in the same bijection proves (4).

## 3. Falling-factorial and Laguerre conjugacy

On formal profile monomials define

\[
 x^a\star x^b=
 \sum_{k=0}^{\min(a,b)}\binom ak\binom bk k!x^{a+b-k}.       \tag{6}
\]

The elementary product identity

\[
 (t)_{\underline a}(t)_{\underline b}
 =\sum_{k=0}^{\min(a,b)}
 \binom ak\binom bk k!(t)_{\underline{a+b-k}}               \tag{7}
\]

says precisely that the map (x^d\mapsto(t)_{\underline d}) conjugates
(star) to ordinary multiplication.  One can read (7) by counting two
injections according to the size (k) of their common image, or prove it by
coefficient extraction.

There is also a literal one-polynomial closed form.  If (a\le b), then

\[
 x^a\star x^b=a!x^bL_a^{(b-a)}(-x),                         \tag{8}
\]

because

\[
 L_a^{(b-a)}(-x)=
 \sum_{j=0}^a\binom b{a-j}\frac{x^j}{j!}.                  \tag{9}
\]

Thus the matching response is simultaneously a rook polynomial, a Laguerre
polynomial, and a falling-factorial product.  These are three coordinates on
one exact kernel, not analogies between unrelated counts.

## 4. SCC products and exact fixed-depth sequences

Let (S_1,\ldots,S_r) be the strong components of (T) in condensation
order.  A tournament condensation is transitive, so

\[
 T=S_1\triangleright\cdots\triangleright S_r,qquad
 Q_T(t)=\prod_{i=1}^r Q_{S_i}(t).                            \tag{10}
\]

The Newton inversion formula gives

\[
 p_T(d)=\frac{\Delta^dQ_T(0)}{d!}
 =\frac1{d!}\sum_{j=0}^d(-1)^{d-j}\binom dj
   \prod_{i=1}^rQ_{S_i}(j).                                 \tag{11}
\]

For integer (0\le j\le d), terms of depth greater than (j) vanish in
(2).  Consequently the depth-(d) answer in (11) consumes only

\[
 \{p_{S_i}(a):1\le a\le d\},                               \tag{12}
\]

not the full profiles of arbitrarily large components.  After these short
jets are supplied, (11) costs (O(rd^2)) elementary arithmetic, or
(O(rd)) with the (Q_{S_i}(j)) table precomputed.  Transforming and
multiplying complete profiles by signed Stirling tables costs (O(N^2))
arithmetic for total order (N), replacing the direct cubic matching
convolution.

For the repeated join (T_r=A^{\triangleright r}), (11) becomes the finite
exponential sum

\[
 \boxed{
 p_{T_r}(d)=\frac1{d!}\sum_{j=0}^d(-1)^{d-j}\binom dj
 Q_A(j)^r.}                                                \tag{13}
\]

For fixed (d), the sequence in (r\) therefore satisfies a constant-order
linear recurrence whose distinct characteristic roots lie in

\[
 \{Q_A(0),Q_A(1),\ldots,Q_A(d)\}.                           \tag{14}
\]

Its order is at most (d+1), regardless of the number of vertices in
(T_r).  At (d=1), (13) reduces to

\[
 H(A^{\triangleright r})=H(A)^r.                            \tag{15}
\]

This is the promised closed-form sequence law: it is derived from the
operation, not fitted to initial values.

## 5. Exact information boundary

The transform deliberately forgets condensation order because ordinary
multiplication is commutative.  This loss occurs already on four vertices.
Let (C_3) be the directed triangle.  The tournaments

\[
 K_1\triangleright C_3,qquad C_3\triangleright K_1          \tag{16}
\]

have different score sequences

\[
 (1,1,1,3)\ne(0,2,2,2),                                   \tag{17}
\]

so they are nonisomorphic, but both have path-cover profile

\[
 (p(1),p(2),p(3),p(4))=(3,9,6,1).                          \tag{18}
\]

Equation (10) retains the strong factors and every path-cover count, but not
their linear order.  The ordered SCC list is the required sidecar for any
consumer that needs the original tournament.

Nor does (3) extend to a cyclic quotient.  The transitive three-vertex
tournament is (K_1^{\triangleright3}) and has

\[
 Q_{T_3}(t)=t^3,                                           \tag{19}
\]

whereas

\[
 Q_{C_3}(t)=t^3+2t.                                        \tag{20}
\]

The extra (2t) is the cyclic quotient-walk contribution isolated by
THM-3121/3134.  Cyclic substitution needs the complete run sidecar; it is not
an order-join product in disguise.

Finally, the transform does not produce an arbitrary input profile: its first
falling-factorial coefficient is already the Hamiltonian-path count.
Equations (3), (11), and (13) are efficient response and sequence laws once
the factor jets are supplied; they assert no arbitrary-input complexity
collapse.

## 6. Connection contract and exact replay

The connection is typed as follows.

~~~text
source:     spanning path-cover profiles of order-join factors
target:     the joined profile and its fixed-depth sequences
map:        colour equal path components, then concatenate across the cut
preserved:  every unordered spanning path-cover count
lost:       order of the strong components
sidecar:    ordered SCC list when a later consumer needs the tournament
test:       falling-factorial product and finite-difference inversion
~~~

The falling-factorial coordinate is the same discrete basis appearing in
THM-2412, but no factorial or Gaussian predicate transfers here.  What
transfers is the exact lowering/product calculus.

The canonical companion uses a subset Hamilton-path DP followed by a
least-vertex set-partition DP.  The independent companion instead enumerates
permutations by backward-adjacency cuts.  Together they check all labelled
tournaments through order five, thousands of direct joins, SCC products,
Laguerre/falling-factorial identities, repeated closed forms, and the two
hostile boundaries.  Run

~~~text
python3 04-computation/tournament_order_join_falling_factorial_transform_thm3162.py
python3 -O 04-computation/tournament_order_join_falling_factorial_transform_thm3162.py
python3 04-computation/tournament_order_join_falling_factorial_independent_audit_thm3162.py
~~~

and compare with the declared transcripts.

**End of proof.**
