# Logarithm-free balanced lifts from orbit Gauss sums and legal menus

**PROVED scoped, using CITED Lovett--Meka partial colouring; same-agent
review.** This strengthens LIFTS.md's existential square-root-log bound.
It applies to the specified ramified tower with fixed p. It is not a
general growing-orbit discrepancy theorem or a complement-closed grammar.

Use all notation and hypotheses of [LIFTS.md](LIFTS.md). Put m=p^k,
N=m², and let a parent be fixed at level k-1, k>=2. Define

    D_p=36p ceil(log_2 p) sqrt(2 log(16p)).                    (1)

**Theorem L4.** Every fixed parent has a balanced, hence prime, lift with

    ||K_k|| <= max{p ||K_(k-1)||, D_p sqrt(N_k)}.             (2)

Consequently every base state extends to a compatible prime tower with

    ||K_k|| <= D_p sqrt(N_k), k>=2.                          (3)

The constants are deliberately loose and depend on p. Parent signs,
child ports and their majority relation remain exact. The spectral scale
sqrt(N) is optimal up to constants: skew adjacency has squared Frobenius
norm N(N-1), annihilates constants, and thus has norm at least sqrt(N).
We do not implement Lovett--Meka here or claim an improved bit complexity.

## 1. The pointwise estimate missing from the energy bound

For a new character xi, v=v_delta(xi) is0 or1. For an orbit O containing
x of valuation ell, write S=sum_(y in O) chi_xi(y). If r=v+ell<2k, then

    r even: S=0;
    r odd:  |S|²=p^(k+(r-1)/2+1-2 ceil(ell/2)).              (4)

If r>=2k, chi_xi is identically one on O and its sine coefficient is zero.
In particular the real matrix entries in LIFTS.md satisfy

    |C_(xi,O)| <= 2 p^v p^((k-1)/2)
               <= 2 sqrt(p m) =: L.                        (5)

The full valuation pattern (4) is stronger than the row variance. A row
may have many nonzero entries, so its l1 norm does not give (3).
It also gives an exact two-channel decomposition: unit new characters
see only odd-valuation orbit signs, while valuation-one new characters
see only even-valuation signs. The vanishing-product case contributes
zero to both sine channels. This parity statement concerns arithmetic
Fourier responses, not parity independence in a Collatz orbit.

### Proof of the orbit-sum formula

Set D=a²+4=p d, where d is a p-adic unit. The Cayley parametrization

    h(t)=(1+t delta)/(1-t delta), t in Z/p^k Z,               (6)

gives exactly H=<g>. Its p^k values are distinct, its image is a group
(the parameter law is (s+t)/(1+D st)), and it contains g: the parameter
(g-1)/(delta(g+1)) is fixed by quadratic conjugation and is integral.
Since g has order p^k, equality follows.

Let lambda extract the coefficient of1 in the (1,mu) basis. Then
chi_xi(y)=exp(2pi i lambda(xi y)/p^k). In the delta basis,
lambda(A+B delta)=A-aB. If z=A+B delta, direct multiplication gives

    f_z(t)=lambda(z h(t))
      =[(A-aB)(1+D t²)+2t(BD-aA)]/(1-Dt²).                 (7)

For r=2s<2k, write xi x=p^s z with z a unit. Then A is a unit, and
f_z modulo p is affine linear with nonzero slope -2aA. Complete additive
sums modulo every positive p-power vanish: sum over the last digit,
using that the derivative remains nonzero modulo p. Therefore the full
H sum vanishes, and so does S.

For r=2s+1<2k, write xi x=p^s z with z=delta times a unit. Then p divides
A and B is a unit. Formula (7) becomes

    f_z(t)=f_z(0)+p F(t),
    F(t)=[2d(A-aB)t²+2(Bd-aA/p)t]/(1-Dt²).                 (8)

Modulo p, F is a quadratic with unit leading coefficient. Its derivative
has exactly one simple zero modulo p, lifting uniquely to a p-adic
critical point. The complete exponential sum of F modulo p^j has
absolute value p^(j/2). A short stationary-phase proof is sufficient:
for j>=2, last-digit cancellation removes all noncritical residue
classes; shift to the critical point and put t=c+p u. Dividing the
nonconstant part by p² leaves another function congruent to a quadratic
with unit leading coefficient. The sum becomes a unit phase times p
times the analogous sum at j-2. The bases j=0 and1 have magnitudes1
and sqrt(p); the latter is the elementary complete quadratic Gauss sum,
whose squared magnitude follows by expanding and summing the difference
variable. This proves the assertion without assuming a heuristic random
phase or invoking an unproved cancellation rate.

There are p^(s+1) repeated parameters when (8) is summed over t modulo
p^k. Thus the full H sum has absolute value

    p^(s+1) p^((k-s-1)/2).

The stabilizer of x has size p^ceil(ell/2). Divide by that stabilizer
and square to obtain (4). If v=0, r odd forces ell=2s+1 and the largest
magnitude occurs at ell=1. If v=1, it forces ell=2s and the largest
occurs at ell=0. These two cases give (5). QED.

## 2. Rounding whole legal choices, instead of individual signs

We use the following inherited consequence of Lovett--Meka,
[Theorem4](https://arxiv.org/html/1203.5747#S2): for an M-by-G real matrix
B, M>=G, with entries bounded by L', and any x in [-1,1]^G, some
eta in {-1,1}^G satisfies

    ||B(eta-x)||_infinity <=18L' sqrt(G log(16M/G)).           (9)

The [previous packet, Theorem A](../pinned_discrepancy_reset_tails_20261029/README.md)
derives this exact form from the arbitrary-start partial-colouring lemma
by repeated halving. Its assumptions are checked here for the NEW menu
matrix; no bounded-orbit hypothesis is silently reused.

For each parent group, label its p children cyclically. Take p legal
balanced patterns: the p cyclic shifts of a fixed block of (p-1)/2
minority signs. Each has signed sum equal to the parent sign; their
uniform mean on each child is parent_sign/p. For each kernel bit take
the two legal choices plus and minus with equal weights. These are
independent menus, not extra graph edges. Their total number is

    G=p^(k-1)+p-2 <=2m/p.

There are M=(p-1)p^(k-1) new absolute character rows, and G<=M.
The sum of the menus' mean Fourier contributions is zero, as proved
in LIFTS.md: the expected nonkernel symbol is constant on each p²-fibre.
Every menu value in a row has absolute value at most pL by (5).

Split each unfinished menu into two parts of as equal size as possible.
Let its two conditional mean vectors be u_i and v_i, and let alpha_i
be the proportion in the first part. Choosing a part replaces
alpha_i u_i+(1-alpha_i)v_i by u_i or v_i. Define

    B_i=(u_i-v_i)/2,  x_i=2alpha_i-1.

Every entry of B_i is bounded by pL. Apply (9) jointly across all menus
to choose their parts, paying at most

    18pL sqrt(G log(16M/G))
      <=36p m sqrt(2 log(16p))                              (10)

in every new row. Finished menus have zero columns and can be padded
to keep G fixed. Repeat on the selected submenus. At most ceil(log_2 p)
rounds leave a single legal pattern in each menu. Telescoping from the
zero total mean proves max_new |lambda|<=D_p m.

All old modes are unchanged by choices within these balanced menus:
they remain p times the parent's modes. This proves (2). The base norm
is at most p²-1, so its norm divided by sqrt(N_1) is less than p<D_p.
Induction proves (3). No source-port condition is weakened to an average;
every final menu leaf satisfies the exact parent constraint. QED.

## 3. Scope, checks and what remains

`check_lifts.py` computes S times its conjugate as an integer cyclotomic
polynomial and checks (4) for every new-row/orbit pair at golden levels
k=1,2,3 and at (a,p,k)=(3,13,1). This is separate from its row-energy
calculation and interval spectral search. The finite candidates themselves
are not produced by the partial-colouring theorem, and their small
constants do not establish a sharper uniform constant.

The stronger theorem answers the growing-orbit logarithm question ONLY
for these ramified compatible lifts. Arbitrary pinned signs may carry
an unavoidable parent spectral debt, explicitly retained in (2).
Uniformity as p varies, an efficient exact implementation, sharp constants,
and closure under native path complements remain open. In particular,
this theorem does not transfer a random or discrepancy estimate to a
deterministic Collatz orbit. The common research move is to preserve
legal finite choices while rounding their response, not to identify the
two mathematical dynamics.
