# Sharp weighted stability at orders five and six

Status: PROVED, independently audited. The order-five (1,13/16) corner
has an exact all-orthant Bernstein certificate; the other principal claims
below have direct algebraic nonnegative decompositions.

## Definitions and low-order cycle expansion

A is a weighted tournament, H its weighted Hamilton-path sum,
u=det(I+A), v=(-1)^n det(A), V=sum_edges q_e with q_ij=A_ij A_ji,
Delta=H-u-|v|, and G_+=H-u+v, G_-=H-u-v.
Let C_j be the total weight of directed j-cycles, each cycle counted once
up to cyclic rotation, and put tau=C_3. Let M_ab denote total weight of
unordered disjoint cycle collections of lengths a,b (q_e is a 2-cycle
weight); similarly M_222 is the sum over three disjoint 2-cycles.
M_22 counts each unordered edge matching once, and M_33 counts each
unordered pair of directed triangles once.

For a triple R={i,j,k}, let K=A-A^T and

    s_R=K_ij+K_jk+K_ki,
    Q_R=sum_{e subset R} q_e.

The choice of cyclic orientation changes s_R at most by sign, and

    Q_R=2 tau_R+1/4-s_R^2/4.                         (T)

This follows by expanding
4 tau_R=1+K_ij K_jk+K_jk K_ki+K_ki K_ij and
4Q_R=3-K_ij^2-K_jk^2-K_ki^2.

## Order five: a sharp nonnegative decomposition

Determinant cycle expansion and the weighted odd-cycle formula give

    H=1+2tau+2C5,
    u=1-V+tau-C4+C5+M22-M23,
    v=M23-C5.

Thus

    G_-=V+tau+C4+2C5-M22,
    G_+=V+tau+C4+2M23-M22.

For each edge e its complement R is a triple. Since
sum_e q_e Q_(e^c)=2M22 and sum_e q_e tau_(e^c)=M23, identity (T) gives

    V/8-M22=-M23+(1/8)sum_e q_e s_(e^c)^2.

Consequently

    G_--3tau/4-7V/8
      =C4+2C5+sum_R tau_R(1/4-q_(R^c))
           +(1/8)sum_e q_e s_(e^c)^2 >=0,

    G_+-3tau/4-7V/8
      =C4+tau/4+M23+(1/8)sum_e q_e s_(e^c)^2 >=0.

Therefore the all-weighted inequality

    Delta >= (3/4)tau+(7/8)V                        (5A)

is analytically proved. In particular the optimal unreserved stability
constant at n=5 is 7/8, attained by two fair pairs and one singleton.
Because tau has positive coefficient in (5A), equality with V>0 in the
unreserved inequality forces tau=0, hence H=1, so H1_STABILITY.md gives
the complete equality classification.

An exact degree-two Bernstein certificate additionally proves

    Delta >= tau+(13/16)V.                         (5B)

The constant 13/16 is sharp: an ordinal sum of a fully fair 3-block and a
fair 2-block has (H,u,v,V,tau)=(3/2,3/8,1/16,1,1/4), hence
Delta=17/16. It also refutes the false proposed strengthening
Delta>=tau+7V/8 by the exact deficit -1/16.

The two valid corners (3/4,7/8) and (1,13/16), plus tau,V>=0, prove the
complete coefficient region for inequalities Delta>=a tau+b V:

    a <= 1,   b <= 7/8,   a/4+b <= 17/16.            (5C)

No lower bounds on a,b are needed. Necessity is witnessed respectively by
C3 ordinal-sum two singletons, two fair pairs plus a singleton, and the
fair 3-block plus fair 2-block. Sufficiency follows by convex interpolation
between (5A),(5B) and downward closure. For a<3/4 use (5A); for a>1 no
coefficient b works because the first witness has V=0.

The (5A) equality set is not discrete: for an ordinal sum of a fair pair
and a weighted triple, its G_- margin equals s_R^2/32. Thus every weighted
triple satisfying s_R=0 gives equality. These interior zeros explain why
an unsubdivided tensor Bernstein certificate for (5A) has negative entries
although (5A) is true; the analytic sum-of-squares is the correct repair.

## Order six: sharp triangle-plus-variance rectangle

The exact cycle expansions are

    H=1+2tau+2C5+4M33,
    u=1-V+tau-C4+C5-C6+M22-M23+M24+M33-M222,
    v=-C6+M24+M33-M222.

Hence

    G_+=V+tau+C4+C5-M22+M23+4M33,
    G_-=V+tau+C4+C5+2C6-M22+M23-2M24+2M33+2M222.

For an edge complement R of size four, sum (T) over its four triples.
Each edge occurs twice, so

    V_R=tau_R+1/2-(1/8)sum_{T subset R, |T|=3} s_T^2.

Let S=sum_e q_e sum_{T subset e^c,|T|=3} s_T^2. Summing the preceding
identity against q_e gives

    M22=M23/2+V/4-S/16.

Therefore

    G_+-tau-3V/4=C4+C5+M23/2+4M33+S/16 >=0,

    G_--tau-3V/4
      =C4+C5+2C6+M23/2-2M24+2M33+2M222+S/16 >=0.

For the last inequality use M24<=C4/4: every directed 4-cycle has a
unique complementary edge, whose q is at most 1/4. Thus C4-2M24>=C4/2.
All other displayed terms are nonnegative. We have the analytic theorem

    Delta >= tau+(3/4)V,                           (6A)

with sharp variance coefficient (three fair pairs attain it). Together
with a deterministic triangle and three singletons, this also proves the
complete coefficient region Delta>=a tau+bV at n=6 is the rectangle

    a<=1, b<=3/4.

Equality with V>0 for the unreserved Delta>=3V/4 again forces tau=0,
so the H=1 classification makes three fair pairs the unique block type.

## Exact certificate mechanism and reproducibility

For n=5,6 partition the weighted cube into binary orientation orthants
A_ij=(1+s_ij z_ij)/2 for i<j, with z_ij in [0,1] and s_ij in {+1,-1}.
Relabeling reduces these to all 12 and 56 tournament orbits respectively;
the enumerator explicitly partitions the 2^10 and 2^15 binary orientations.
Each signed target is degree at most two in every z variable.

`bernstein.py` expands H and both determinants by permutations, with all
coefficients scaled by 2^n. Its tensor transform maps power coefficients
(c0,c1,c2) on one axis to doubled degree-two Bernstein coefficients

    (2c0, 2c0+s c1, 2c0+2s c1+2c2).

Thus all coefficients remain integral. A polynomial with nonnegative
Bernstein coefficients is nonnegative everywhere on the closed orthant.
Every arithmetic path is bounded by 2^m times the initial coefficient l1
bound, safely within int64 at n=5,6; independent audit also reconstructs
coefficients by local arc factors and subset-DP rather than this transform.

    python3 -B bernstein.py 5 7/8 > bernstein5.jsonl
    python3 -B bernstein.py 6 3/4 > bernstein6.jsonl
    python3 -B reserved_cycles.py 5 13/16 > triangle5_13over16.jsonl
    python3 -B reserved_cycles.py 6 3/4 > triangle6.jsonl

Arrays are deliberately regenerated, not stored. The coefficient checks
cover the entire real cube through Bernstein positivity, NOT only a finite
set of weighted matrices. The order-five (5B) certificate is load-bearing;
other certificates independently check the separate analytic proofs.
