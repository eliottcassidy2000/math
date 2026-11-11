# A positive Gram kernel inside literal cyclic tournament substitution

**PROVED IN THIS NOTE.** New factorization, rank, equality and observer
interpretation. The underlying quotient-walk formula is inherited from
THM-3121-path-cover-walk-content-substitution-kernel, rather than claimed new.
No independent canon audit or literature-priority assertion.

## 1. Source objects and the recovered kernel

For a nonempty finite tournament S, define

    f_S(r) = r! times the number of spanning covers by r unordered
             nonempty directed paths,  r >= 1,
    f_S(0) = 0.

Extend by zero past the order of S. Thus f_S(1)=H(S). Let C3[S,T,U] be
literal tournament substitution into a directed triangle S→T→U→S.
THM-3121 gives H as the sum of three diagonal terms and six nearest-shell
terms in the three ordered path-cover profiles. This is the correct
replacement for multiplying bare H values: a global path can revisit blocks.

For fixed U, put u_r=f_U(r), s_r=f_S(r), t_r=f_T(r). Regrouping those nine
terms gives the following exact identity:

    K_U(S,T) := H(C3[S,T,U])
      = Σ_(r>=1) u_r [
          (s_(r−1)+s_r)(t_(r−1)+t_r)
          +(s_r+s_(r+1))(t_r+t_(r+1)) + s_r t_r ].           (1)

Each coefficient is positive. In particular this is not a signed spectral
analogy: every K_U entry is the count of paths in an explicit tournament.

### Proof

The inherited diagonal coefficient is 3u_r s_r t_r. The terms with one
profile index raised by one give
u_r(s_(r+1)t_r+s_r t_(r+1))+u_(r+1)s_r t_r. The terms with two raised
indices give u_r s_(r+1)t_(r+1)+u_(r+1)(s_(r+1)t_r+s_r t_(r+1)). Expanding
(1), and shifting the first pair's index down once, yields precisely these
terms. Every sum is finite. This proves (1). QED.

## 2. Rank and exact equality

Let c=|U| and v_S=(s_1,...,s_(c+1)). Equation (1) is v_S^T W_U v_T, with

    (W_U)_(r,r) = u_(r−1)+3u_r+u_(r+1),
    (W_U)_(r,r+1) = u_r+u_(r+1),

for 1<=r<=c+1; other entries vanish and u_0=u_(c+1)=u_(c+2)=0.
Every u_r is positive for 1<=r<=c: split one Hamiltonian path of U into r
nonempty consecutive pieces. The classical existence of a Hamiltonian path
in a tournament is the only positivity input.

For any real vector z=(z_1,...,z_(c+1)), with z_0=0,

    z^T W_U z = Σ_(r=1)^c u_r [
       (z_(r−1)+z_r)^2 + (z_r+z_(r+1))^2 + z_r^2 ].       (2)

It vanishes only at zero: the terms z_r² force z_1,...,z_c=0, and the last
pair square forces z_(c+1)=0. Therefore W_U is positive definite.

Consequently, for every finite collection of tournaments S_i, the matrix
[K_U(S_i,S_j)] is positive semidefinite and

    K_U(S,T)^2 <= K_U(S,S) K_U(T,T).                       (3)

Equality holds exactly when the two truncated ordered path-cover vectors
are proportional. Distinct tournaments can have proportional profiles, so
this is not an isomorphism test. Because the first coordinates are positive,
the proportionality factor in the equality case is positive.

The kernel has rank exactly c+1 on the whole collection of tournaments.
Indeed, for the transitive tournament TT_a,

    f_(TT_a)(r) = r! S(a,r),                              (4)

where S(a,r) is a Stirling number of the second kind: each set partition
has a unique increasing directed path inside each block. The matrix with
rows a=1,...,c+1 and columns r=1,...,c+1 is triangular with nonzero diagonal
r!, hence invertible. Its W_U Gram matrix is therefore nonsingular. Every
(c+2)-row Gram determinant is zero, and the displayed (c+1)-row one is
positive. The dimension is sharp, not only an upper bound.

## 3. A finite separating bank of actual graph contexts

For fixed U, the c+1 counts

    H(C3[S,TT_a,U]),  a=1,...,c+1,                        (5)

determine the complete truncated profile (f_S(1),...,f_S(c+1)). This is a
linear inversion with invertible matrix (f_(TT_a)(r)) W_U. Thus two source
tournaments have the same response in every context C3[S,T,U] if and only
if they agree at this finite bank. This is a genuine contextual equivalence
certificate, parallel to the nine-response bank for the five-core apex
formula. It does not determine profile coordinates above c+1 or responses
with a larger observer U.

The map is S -> its truncated profile -> all observer-U cyclic responses.
It preserves exactly those responses, loses larger jets and vertex-level
structure, and restores the response information by the c+1 actual probes.
The cheapest decisive test is an exact recount of those probes, followed by
the rational linear inverse. No enumeration of every possible T is needed.

## 4. The one-vertex observer exposes the correct norm

When U=K1, write h=H(S), p=f_S(2), k=H(T), q=f_T(2). Then

    K_(K1)(S,T) = 2hk+(h+p)(k+q),                         (6)
    K_(K1)(S,S) = 2h²+(h+p)².                            (7)

The Gram determinant has the exact square-times-two form

    K(S,S)K(T,T)−K(S,T)^2 = 2(hq−kp)^2.                  (8)

This is a norm geometry over Q(sqrt(−2)), not automatically the ordinary
Pythagorean cone. In particular, if the two diagonal counts equal N, the
identity K(S,T)^2+d²=N² can have integer d only when hq=kp; there is no
nondegenerate primitive Pythagorean triple obtained from this two-channel
Gram area. The obstruction is arithmetic, not an absence of positivity.

This comparison is useful precisely because the older Gaussian/Farey
construction uses the ordinary norm x²+y². A norm-preserving recipe must
carry its quadratic form, not only two integer registers. Changing that
form silently would be the first invalid implication in a claimed transfer.
A possible new route is to use a larger fixed observer and study its
integral quadratic form and its *realized* path-cover vectors together.
That realization problem remains OPEN.

## 5. Controlled cycle packing and all-size recurrence

For S=TT_a, T=TT_b, U=TT_c, each directed cycle uses at least one vertex of
each block. Hence the maximum number of vertex-disjoint directed cycles is
at most min(a,b,c), and equals that minimum by using disjoint directed
triangles. Fixing c is an exact cycle-packing cap.

Equation (1) and (4) give every count with that cap. Since

    r! S(a,r) = Σ_(j=0)^r (−1)^(r−j) binom(r,j) j^a,

for positive a the sequence K_U(TT_a,T) obeys the constant-coefficient
recurrence

    [product_(j=1)^(c+1) (E−j)] K_U(TT_a,T) = 0,          (9)

where E shifts a to a+1. The order c+1 is minimal as a *common* recurrence
for all T: the c+1 transitive probes recover the full span of j^a for
j=1,...,c+1. An individual T can have a shorter recurrence. This fixed
observer theorem does not contradict THM-3202's moving-jet obstruction:
when all three block sizes grow, the required c and recurrence order grow.

For c=1, (6) becomes the inherited safe-inflation family

    H(C3[TT_a,TT_b,K1]) = 2+(2^a−1)(2^b−1).

For c=2, put P_a=2^a−1, R_a=3^a−2·2^a+1, V_a=2^a−2. Then

    H(C3[TT_a,TT_b,TT2])
       = 2+3P_a P_b+2R_a R_b+2V_a V_b.                  (10)

This is the controlled-packing-two continuation requested by Session251.
It gives H(C3[TT2,TT2,TT2])=45. Applying the intercyclic formula to that
same graph illegally gives29, because two disjoint cycles now exist.
Formula (10) includes them through the full profile kernel.

## 6. Chan–Pak contribution and scope

The new ingredient is a positive, sharp-dimensional geometry of *actual*
cyclic graph counts, plus finite contextual probes and fixed-observer
recurrences. It offers a bank of extra source graphs beyond the five-core
one-apex language, with constructive output and an exact cycle-packing cap.
One can search norm fibres, equal-response collisions and observer changes
without claiming that scalar H determines substitution.

It does not prove all odd values except7,21, a universal three-branch source
operation, a bounded observer suffices, or that arbitrary vectors in the
positive Gram space are realizable by tournament profiles. Chan–Pak and
those realization questions remain OPEN.

## Sources and reproducibility

- `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`:
  inherited exact nine-shell kernel and path-cover interpretation.
- `01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`:
  inherited path-cover/newton and order-join transforms.
- `01-canon/theorems/THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md`:
  inherited moving-observer boundary.
- Session251 `cycles/README.md`: safe inflation and the proposed packing-two
  departure; inherited value45 hostile.

The companion `check.py` tests only the new factorization and consequences,
using a small independently implemented subset-DP recount for actual graphs.
Its finite ranges are controls for the proofs, not an all-target census.
