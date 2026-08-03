---
id: THM-3181
title: "Tournament half-grid reciprocity and repeated-join recurrence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The tournament
  path-colour polynomial has pure order parity: Q_T(-t)=(-1)^n Q_T(t), so
  the negative-colour coordinate of THM-3166 is exactly the ordinary positive
  path-colour coordinate.  Its nonnegative-integer generating function is
  x B_T(x)/(1-x)^(n+1), giving an explicit binomial inverse and recovery of
  the complete backward-adjacency distribution and path-cover profile from
  only Q_T(1),...,Q_T(ceil(n/2)).  The older tournament Worpitzky sequence is
  Delta Q_T.  Along SCC chains every fixed extreme coefficient uses only a
  finite product jet, and repeated joins give an exact polynomial-exponential
  formula with a minimal constant-coefficient recurrence of triangular order.
  A truncated endpoint-deletion DP computes fixed jets without computing the
  full distribution.  SCC order and cyclic-substitution run data remain lost;
  no arbitrary-tournament complexity collapse is claimed.
audit: >
  The direct-permutation companion checks all 1,099 labelled tournaments
  through order five, 8,792 instances each of parity, positive binomial
  reciprocity, the rational generating function, and Worpitzky integration,
  plus 1,099 half-grid inversions, SCC hostiles, 40 repeated-join tails through
  total order eight, the C3 order-six annihilator, and the transitive Eulerian
  boundary.  An independent engine uses endpoint subset DP and a separate
  Held--Karp/canonical-set-partition path-cover computation on 283 fixed cases
  through order eight and 803 truncated jets, then checks an SCC chain and the
  nine-vertex cyclic-substitution hostile 3159 versus 33.  Normal, optimized,
  and stored outputs match byte for byte.
source: root/frontier-synthesis/tournament-half-grid/2026-08-02
depends_on:
  - THM-3166-falling-factorial-order-join-path-colour-transform
  - THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform
  - THM-087-worpitzky-tournament
related:
  - THM-3121-path-cover-walk-content-substitution-kernel
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
script: 04-computation/tournament_half_grid_reciprocity_repeated_join_thm3181.py
output: 05-knowledge/results/tournament_half_grid_reciprocity_repeated_join_thm3181.out
script_sha256: fb66d1a3d6ff32060a5cfb901b0bd18fa1c31d726a8bd14af6b8eb1fcd3fe87d
output_sha256: 2393c07e7a2767a16274db3fe6899f50fd7db3134b9eb11a5c0ee0c127859a41
independent_script: 04-computation/tournament_half_grid_reciprocity_independent_audit_thm3181.py
independent_output: 05-knowledge/results/tournament_half_grid_reciprocity_independent_audit_thm3181.out
independent_script_sha256: 4826061239748b2d765d4b0060616733fcd34b186c2bdc902a49f97a7c173f72
independent_output_sha256: 9a5ab8d3203acd78b19b1813c7b9f7861e6607a30ccb9ca5502d7e2c150fed3f
hash_basis: LF-normalized bytes
---

# THM-3181 -- tournament half-grid reciprocity and repeated-join recurrence

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3166 makes the falling-factorial path-colour polynomial multiplicative
under order join and interprets its negative integer values positively.  The
permutation-reversal sidecar in THM-3134 strengthens that result: the
negative coordinate is the positive coordinate itself, and only half of the
integer grid is needed to reconstruct the complete profile.  At fixed extreme
depth, repeated order joins then become exact constant-recurrence sequences.

## 1. Objects and two binomial forms

Let `T` be a tournament on `n>=1` vertices.  Write

\[
 p_T(d)=\#\{\text{spanning covers of }T\text{ by }d\text{ unordered
 directed paths}\},
\]

and define the THM-3166 path-colour polynomial

\[
 Q_T(t)=\sum_{d=1}^n p_T(d)(t)_{\underline d}.             \tag{1}
\]

For a vertex permutation `pi`, let `b_T(pi)` be the number of consecutive
pairs oriented backward along `pi`, and put

\[
 A_T(b)=\#\{\pi:b_T(\pi)=b\},\qquad
 B_T(x)=\sum_{b=0}^{n-1}A_T(b)x^b.                         \tag{2}
\]

Reversing a permutation exchanges its forward and backward consecutive
adjacencies.  Therefore

\[
 A_T(b)=A_T(n-1-b),\qquad x^{n-1}B_T(x^{-1})=B_T(x).       \tag{3}
\]

The path-colour polynomial has the two exact binomial forms

\[
 \boxed{Q_T(t)=\sum_{b=0}^{n-1}A_T(b)
                    \binom{t+n-1-b}{n}
             =\sum_{b=0}^{n-1}A_T(b)\binom{t+b}{n}.}       \tag{4}
\]

These are polynomial identities in `t`, with generalized binomial
coefficients.

### Proof

Put `F_T(d)=d!p_T(d)`.  Concatenating the ordered paths of a cover gives a
permutation with cuts.  For a permutation having `b` backward and
`g=n-1-b` forward adjacencies, every backward position must be cut and an
arbitrary subset of the forward positions may be cut.  Its contribution to
`Q_T(t)=sum_d F_T(d) binom(t,d)` is

\[
 \sum_{s=0}^g\binom gs\binom{t}{1+b+s}
   =\binom{t+g}{n},                                      \tag{5}
\]

by Vandermonde.  Summing proves the first form of (4), and (3) proves the
second.  This also gives a self-contained proof that the cut construction is
valid for polynomial, not only nonnegative integer, colours.

## 2. Pure parity and positive-negative identification

The generalized-binomial reflection identity gives

\[
 \binom{-t+b}{n}=(-1)^n\binom{t+n-1-b}{n}.
\]

Applying it to the second form of (4) yields

\[
 \boxed{Q_T(-t)=(-1)^nQ_T(t).}                            \tag{6}
\]

Thus `Q_T` contains only ordinary powers of the same parity as `n`.  More
specifically, the positive reciprocity coordinate from THM-3166 satisfies

\[
 R_T(m):=(-1)^nQ_T(-m)=Q_T(m)
     =\sum_{b=0}^{n-1}A_T(b)\binom{m+b}{n}                \tag{7}
\]

for every integer `m`.  Negative colours do not supply a second independent
coordinate; they reveal a positive interpretation of the existing one.

## 3. Rational integer-value series and half-grid inverse

For `0<=b<n`,

\[
 \sum_{m\ge0}\binom{m+b}{n}x^m
   =\frac{x^{n-b}}{(1-x)^{n+1}}.                          \tag{8}
\]

Using (3) and (7) gives the rational generating function

\[
 \boxed{\sum_{m\ge0}Q_T(m)x^m
       =\frac{xB_T(x)}{(1-x)^{n+1}}.}                     \tag{9}
\]

Extracting the coefficient of `x^k` after multiplying by
`(1-x)^(n+1)` gives, for `1<=k<=n`,

\[
 \boxed{A_T(k-1)=A_T(n-k)=
   \sum_{m=0}^k(-1)^{k-m}\binom{n+1}{k-m}Q_T(m).}         \tag{10}
\]

Since `Q_T(0)=0` and `A_T` is palindromic, the values

\[
 Q_T(1),Q_T(2),\ldots,Q_T(\lceil n/2\rceil)              \tag{11}
\]

recover the complete backward-adjacency distribution.  They then recover the
complete path-cover profile through the THM-3134 cut dictionary

\[
 d!p_T(d)=\sum_{b=0}^{d-1}
   \binom{n-1-b}{d-1-b}A_T(b).                            \tag{12}
\]

Thus the `n` negative evaluations in THM-3166 are sufficient but redundant:
`ceil(n/2)` ordinary positive evaluations suffice.  No sharp minimality claim
among the set of realized tournament polynomials is made.

## 4. THM-087 Worpitzky is the discrete derivative

Let `a_m(T)` be the degree-`n-1` Worpitzky sequence of THM-087.  Its
forward-adjacency polynomial equals `B_T` coefficientwise by (3), and THM-087
states

\[
 \sum_{m\ge0}a_m(T)x^m=\frac{B_T(x)}{(1-x)^n}.             \tag{13}
\]

Comparing (9) and (13) proves

\[
 \boxed{Q_T(0)=0,\qquad Q_T(m+1)-Q_T(m)=a_m(T).}          \tag{14}
\]

The path-colour polynomial is therefore the canonical discrete
antiderivative of the older tournament Worpitzky polynomial.  At the
transitive boundary, `a_m=(m+1)^n-m^n` and `Q_T(m)=m^n`.

## 5. SCC-chain finite jets

Let the strong components of `T`, in source-to-sink order, be
`S_1,...,S_r`; then `T=S_1▷...▷S_r`.  THM-3166 gives

\[
 Q_T(m)=\prod_{i=1}^rQ_{S_i}(m).                          \tag{15}
\]

Substitution in (10) yields, for `1<=k<=ceil(n/2)`,

\[
 \boxed{A_T(k-1)=A_T(n-k)=
 \sum_{m=1}^k(-1)^{k-m}\binom{n+1}{k-m}
       \prod_{i=1}^rQ_{S_i}(m).}                          \tag{16}
\]

Given the component `k`-jets, one selected extreme coefficient takes
`O(rk)` arithmetic after binomial preparation; every coefficient through
depth `k` takes `O(rk+k^2)`.  The product in (15)--(16) is commutative and
therefore still destroys SCC order.  An ordered SCC list is a mandatory
structural sidecar.

## 6. Repeated joins: closed form and minimal recurrence

Let `A` be a tournament of order `a`, put

\[
 T_r=A^{\triangleright r},\qquad \lambda_m=Q_A(m),         \tag{17}
\]

and fix `k>=1`.  For every `r` with `ar>=k`, define the symmetric extreme
coefficient

\[
 s_r^{(k)}=A_{T_r}(k-1)=A_{T_r}(ar-k).                    \tag{18}
\]

Equation (16) becomes the exact polynomial-exponential formula

\[
 \boxed{s_r^{(k)}=
 \sum_{m=1}^k(-1)^{k-m}\binom{ar+1}{k-m}\lambda_m^r.}    \tag{19}
\]

The values `lambda_1,...,lambda_k` are positive and strictly increasing,
because (14) gives `lambda_(m+1)-lambda_m=a_m(A)>0`.  In (19), the
coefficient of `lambda_m^r` is a polynomial in `r` of exact degree `k-m`,
with nonzero leading coefficient.  Linear independence of distinct
polynomial-exponential sequences now proves that the reduced ordinary
generating-function denominator and the minimal constant-coefficient tail
recurrence polynomial are respectively

\[
 \boxed{\prod_{m=1}^k(1-\lambda_mz)^{k-m+1},\qquad
        \prod_{m=1}^k(E-\lambda_m)^{k-m+1},}              \tag{20}
\]

where `(Es)_r=s_(r+1)`.  Their order is

\[
 \sum_{m=1}^k(k-m+1)=\frac{k(k+1)}2,                     \tag{21}
\]

independent of the join level and the seed order.  The first three instances
are

\[
\begin{aligned}
 A_{T_r}(0)&=\lambda_1^r,\\
 A_{T_r}(1)&=\lambda_2^r-(ar+1)\lambda_1^r,\\
 A_{T_r}(2)&=\lambda_3^r-(ar+1)\lambda_2^r
             +\binom{ar+1}{2}\lambda_1^r.
\end{aligned}                                             \tag{22}
\]

For `A=K_1`, `lambda_m=m`, and (19) is exactly the classical closed form for
the Eulerian number with `k-1` descents.  Equation (19) is therefore a
seed-tournament Eulerian generalization, not a fitted recurrence.

## 7. Endpoint deletion and truncated computation

For a nonempty vertex subset `S` and `v in S`, let `D_(S,v)(z)` count
permutations of `S` ending at `v`, by backward consecutive adjacencies.  Then

\[
 D_{\{v\},v}(z)=1,\qquad
 D_{S,v}(z)=\sum_{u\in S\setminus\{v\}}D_{S\setminus\{v\},u}(z)
                   z^{\mathbf1[v\to_Tu]},                 \tag{23}
\]

and

\[
 B_T(z)=\sum_{v\in V(T)}D_{V(T),v}(z).                    \tag{24}
\]

Indeed, appending `v` after `u` adds one backward adjacency exactly when
`v ->_T u`.  Truncating every polynomial modulo `z^k` computes
`A_T(0),...,A_T(k-1)` in

\[
 O(kn^2 2^n)\text{ arithmetic},\qquad O(kn2^n)\text{ storage}. \tag{25}
\]

The first form of (4), evaluated at the positive integer `k`, truncates to

\[
 \boxed{Q_T(k)=\sum_{b=0}^{k-1}A_T(b)
                  \binom{n+k-b-1}{n}.}                    \tag{26}
\]

Thus a fixed-colour evaluation needs no full backward distribution.  At
`k=1`, however, this computation already contains
`Q_T(1)=A_T(0)=H(T)`.  The exponential endpoint DP is a parameterized exact
reduction, not a polynomial-time algorithm for arbitrary tournaments.

## 8. Hostile boundaries and destroyed information

### SCC order is not recovered

The nonisomorphic tournaments `K_1▷C_3` and `C_3▷K_1` have score sequences

\[
 (3,1,1,1)\qquad\text{and}\qquad(2,2,2,0),                \tag{27}
\]

but their path-colour polynomials and complete backward-adjacency
distributions agree.  Equations (15)--(16) cannot recover the order of strong
components.

### Cyclic substitution has no scalar plethystic law

For the directed triangle, `Q_C3(t)=t^3+2t`.  Independent direct computation
gives

\[
 H(C_3[C_3,C_3,C_3])=3159,                               \tag{28}
\]

whereas the naive scalar plethysm gives

\[
 Q_{C_3}(Q_{C_3}(1))=Q_{C_3}(3)=33.                      \tag{29}
\]

Thus (15) is an order-join law, not a lexicographic-substitution law.
THM-3121's quotient-walk content and the complete factor path-cover profiles,
equivalently THM-3134's endpoint jets, remain necessary for cyclic
substitution.  The present theorem does not compute arbitrary base profiles,
turn scalar Hamiltonian data into cyclic run data, or collapse `#P`-hard
path counting.

## 9. Exact verification

Run

```bash
python3 04-computation/tournament_half_grid_reciprocity_repeated_join_thm3181.py
python3 -O 04-computation/tournament_half_grid_reciprocity_repeated_join_thm3181.py
python3 04-computation/tournament_half_grid_reciprocity_independent_audit_thm3181.py
python3 -O 04-computation/tournament_half_grid_reciprocity_independent_audit_thm3181.py
```

The primary direct-permutation engine checks:

- all `1,099` labelled tournaments of orders `1..5`;
- `8,792` instances each of parity, the positive binomial formula, the
  rational generating function, and Worpitzky integration;
- `1,099` complete half-grid inversions;
- SCC product and SCC-order-loss controls;
- `40` repeated-join extreme coefficients through total order eight;
- the `C_3`, `k=3` order-six annihilator with roots `(3,12,33)` and
  multiplicities `(3,2,1)`; and
- the transitive Eulerian boundary through order ten.

The independent engine uses no permutation enumeration.  It combines the
endpoint-deletion polynomial DP (23) with an independent Held--Karp induced
path engine and canonical least-vertex set-partition recurrence for path
covers.  It checks `283` fixed tournaments of orders `1..8`, `803` truncated
jets, an SCC chain, and the order-nine cyclic hostile (28)--(29).  Normal and
optimized transcripts byte-match the stored outputs, and the four artifact
hashes match the YAML header.

**End of proof.**
