# Outer six-root K3 Picard number is at least nineteen

**Status: PROVED over C**, using the inherited K3/quotient-transfer theorem
and the standard Kodaira and Shioda–Tate inputs pinned below. All displayed
algebra is reproduced exactly by `check_lower_bound.py`; the checker is a
certificate for that algebra, not a replacement for the cited geometric
inputs. This lower-bound snapshot is independent of the separate finite-field
upper-bound investigation. It excludes outer rank 18 but alone does not
choose between ranks 19 and 20.

## 0. Inheritance, objects, and scope

Let S be the minimal resolution of the six-root surface obtained by deleting
E−x−z from the hourglass seven forms. The inherited deck packet identifies
the smooth model Z of its even-sign quotient by H=(C2)^4 with the resolution
of the six-line double-plane

    W² = E(E+x)(E−x)(E+z)(E−z)(E+x+z).

Both S and Z are K3. The quotient-transfer theorem in
`../deck-k3/HODGE_DECOMPOSITION.md`, Sections 1 and 4, identifies their rational
transcendental Hodge structures. Thus rho(S)=rho(Z). This is NOT an integral
lattice isometry, an isomorphism of surfaces, or a transfer of Kummer status.
The degree-sixteen pullback scales the pairing by sixteen.

Closest proved mechanism: the actual symplectic quotient and rational
transcendental transfer, rather than similarity of equations. Canonical
hostile: a lower-bound divisor count is not an equality at a special member.
Corrected near miss: checking only finite places of an elliptic equation
misses a rank-six fiber lattice at infinity. Least-used sidecar: a concrete
section, whose specialization supplies a nontorsion certificate.

## 1. Birational elliptic model and an actual section

On E=1 write t=x, d=1−t². Then

    w² = d(1−z²)(1+t+z).

The change

    X=−dz,   Y=dw,
    z=−X/d, w=Y/d                                      (1)

is invertible over Q(t) and gives

    Y²=(X−d)(X+d)(X−d(t+1)).                           (2)

Indeed the right side after X=−dz is
`d³(1−z²)(1+t+z)=d²w²`. The map is birational over the open set d!=0;
there is no claim that (1) extends as a coordinate isomorphism over t=±1.
The smooth minimal elliptic surface associated to (2), including its zero
section, is the same minimal K3 model Z, since smooth minimal K3 models
of the same surface function field are isomorphic.

The line z=−t has w=d and therefore gives the rational section

    P=(dt,d²).                                         (3)

Direct substitution verifies it: the three factors at X=dt multiply to d⁴.
The three points (d,0), (−d,0), (d(t+1),0) are rational two-torsion;
P is not to be confused with any of them.

## 2. Complete singular-fiber and minimality calculation

In `Y²=X³+a2 X²+a4 X+a6` the coefficients are

    a2=−d(t+1),   a4=−d²,   a6=d³(t+1).

The standard Weierstrass invariants, with Delta equal to sixteen times the
polynomial cubic discriminant, are

    c4=16 d²(t²+2t+4),
    c6=−64(t−2)(t−1)³(t+1)⁴(t+4),
    Delta=64 d⁶ t²(t+2)².                              (4)

At infinity use s=1/t and the integral coordinates

    X_inf=s⁴X,  Y_inf=s⁶Y.

The coefficients there are

    a2_inf=−s(s²−1)(s+1),
    a4_inf=−s⁴(s²−1)²,
    a6_inf=s⁵(s²−1)³(s+1).

In particular

    c4_inf=16s²(s²−1)²(1+2s+4s²),
    Delta_inf=64s⁸(s²−1)⁶(1+2s)².                     (5)

Here is the entire valuation table in characteristic zero:

| Place | v(c4) | v(c6) | v(Delta) | Kodaira type | m_v−1 |
|---|---:|---:|---:|---|---:|
| t=0 | 0 | 0 | 2 | I2 | 1 |
| t=−2 | 0 | 0 | 2 | I2 | 1 |
| t=1 | 2 | 3 | 6 | I0* | 4 |
| t=−1 | 2 | 4 | 6 | I0* | 4 |
| t=infinity | 2 | 3 | 8 | I2* | 6 |

Each displayed model is integral at its place and has v(Delta)<12, so no
integral admissible minimization can lower its discriminant valuation:
any decrease is a positive multiple of twelve. Thus the table is for
minimal models, not a nonminimal equation. Translation to short form has
coefficients −c4/48 and −c6/864, so the valuation criteria in the cited
Kodaira table apply. There are no additional singular places because (4)
has no other finite zero and (5) treats infinity. The Euler numbers sum to
`2+2+6+6+8=24`, independently consistent with K3.

The trivial lattice has root part A1+A1+D4+D4+D6, hence rank

    2 + 1 + 1 + 4 + 4 + 6 = 18.                       (6)

The leading 2 comes from the zero section and a fiber. Fiber multiplicity
is one because a section exists. By Shioda–Tate,

    rho(Z)=18+rank E(C(t)).                             (7)

## 3. Nontorsion, with an elementary denominator proof

The fiber at t=2 is smooth: d=−3 and Delta(2)!=0. The specialization is

    Y²=(X+3)(X−3)(X+9),  P_2=(−6,9).

The *integral translation* u=X+3 turns this into the short Weierstrass
model

    E_6: Y²=u³−36u,     Q=(−3,9).                     (8)

Its tangent slope at Q is −1/2, and

    2Q=(25/4,−35/8).                                  (9)

Standard Nagell–Lutz integrality already proves Q has infinite order: a
rational torsion point, including 2Q, would have integral coordinates on
(8). There is no generalized-Weierstrass or minimal-model loophole here.
However the following argument proves precisely what is needed without
importing that theorem.

For a finite point R on (8) with Y(R)!=0 the tangent formula simplifies to

    u(2R)=(u(R)²+36)²/[4u(R)(u(R)²−36)].               (10)

If v_2(u(R))=−2k with k>=1, then u(R)² has valuation −4k, while
v_2(36)=2. Therefore the numerator in (10) has valuation −8k and its
denominator has valuation 2−6k. Thus

    v_2(u(2R))=−2k−2.                                 (11)

The denominator is nonzero, since a quantity of negative valuation cannot
be 0 or ±6. Starting from (9), induction gives
`v_2(u(2^n Q))=−2n` for every n>=1. These are infinitely many distinct
points, proving Q nontorsion. The finite exact checker evaluates several
initial doubles as a positive control; the proof of the infinite assertion
is the valuation induction, not a bounded numerical test.

Specialization at a smooth fiber respects addition. If the section P were
torsion, its t=2 specialization Q would be torsion. Hence P is nontorsion
over Q(t), and remains nontorsion after extending constants to C.
Equations (7) and (11) prove

    rho(Z)>=19,  rho(S)>=19.                           (12)

No theorem asserting injectivity of the entire Mordell–Weil specialization
map is needed; just the fact that a group homomorphism sends torsion to
torsion. Nontorsion points can in general specialize to torsion, so the
reverse implication would be false.

## 4. Consequence for the seven-root hourglass

The inherited Hodge decomposition and central deletion calculation give

    rho(Y_7)=6 rho(S)−56.

Together with the K3 upper bound rho(S)<=20, (12) sharpens the previous
three-value set `{52,58,64}` to

    rho(S) in {19,20},   rho(Y_7) in {58,64}.           (13)

In particular 52 is excluded. Any exact choice between these values must
come from an additional upper bound, not from this section or genericity.

## 5. The exact Session57 congruent-number bridge

Equation (8) is precisely the congruent-number curve of area 6. The
classical rational map sends (u,Y), Y!=0, to signed sides

    a=(u²−36)/Y,  b=12u/Y,  h=(u²+36)/Y.

The curve equation proves directly `a²+b²=h²` and `ab/2=6`. At Q=(−3,9)
it gives (−3,−4,5), hence the positive 3–4–5 triangle. More specifically,
adding the rational two-torsion point (0,0) sends Q to (12,36), exactly the
point used in Session57's `congruent-curve/PROOF.md`. Both have the same
double (9). Its square-center image is

    c=4u(2Q)=25,  d_gap=24,
    (c−d_gap,c,c+d_gap)=(1,25,49).

This is an explicit curve map and a torsion-translation identification,
not a numerical analogy. It does NOT produce three distinct centers in
arithmetic progression, a magic square of squares, or a global result about
all rational points. Also, the special fiber's j=1728 does NOT imply CM of
the whole K3 surface or identify it with the central Kummer K3 (whose CM
calculation concerns its associated abelian surface, not an isolated fiber).

## 6. Primary pins, dependency and loss ledger

1. Matthias Schuett and Tetsuji Shioda, *Elliptic Surfaces*,
   https://arxiv.org/pdf/0907.0298 : Section 4.5, Table 1, printed p.14,
   gives the characteristic-not-2-or-3 valuation classification; Proposition
   6.6 and Corollary 6.13, printed pp.24 and 29, give trivial-lattice rank
   and Shioda–Tate. Imported only these standard geometric inputs.
2. Alice Garbagnati and Alessandra Sarti, *Kummer surfaces and K3 surfaces
   with (Z/2Z)^4 symplectic action*, Section 10.2, printed pp.54–55,
   https://sarti.pages.math.cnrs.fr/home/GarbagnatiSarti_RMJ.pdf , is the
   primary pin for the diagonal-six-root/six-line-quotient model. The actual
   rational transcendental transfer proof used here is already in the
   inherited `HODGE_DECOMPOSITION.md`, Sections 1 and 4.
3. Session57 internal source: `05-knowledge/results/
   session57_square_sum_magic_paths_20261025/congruent-curve/PROOF.md`,
   Sections 2–3. The elementary map and the exact point identification are
   re-proved above; no new external rank or generator statement is imported.

Connection contract: source S -> quotient model Z -> elliptic fibration E.
Maps: degree-sixteen quotient correspondence, then birational map (1).
Preserved: transcendental rank under the first map, surface function field
under the second. Lost: integral lattice index, arithmetic twists/descent
information in the Hodge transfer, and the exceptional-fiber coordinate
charts in the affine birational expression. Sidecars: inherited transfer
proof, full minimal-fiber table, and the section/smooth-specialization data.
Cheapest hostile tests: drop infinity (rank becomes 12, incorrectly), use
X=+dz (identity fails), or replace d² by d in P (membership fails).

Reproduction: `python3 check_lower_bound.py` and
`python3 -O check_lower_bound.py`. No Git operations are part of this packet.
