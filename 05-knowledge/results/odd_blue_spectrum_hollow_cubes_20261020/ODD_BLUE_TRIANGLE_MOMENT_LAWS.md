# Triangle polynomial and moment laws for hollow-frame cubes

Status: **PROVED addendum**, for the auxiliary class cubes of
`ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md`. Exact checks cover every packet
at `n=3,5,7,9`. This does not modify the original proof packet or its scope.

## 1. Exact packet polynomial

Fix an odd order `n=2k+1>=3`, a hollow reflection frame `D`, its odd
frame group `H`, and its pair-orbits `O_1,...,O_r`. In the notation of the
main theorem, complete the matching with bits `a_i=[x_i -> y_i]`, where
`x_i -> f -> y_i`. Let `d_D^+(x_i)` exclude the deleted matching arc and put

    delta_i = 2k-1-2d_D^+(x_i).

The proved integer matching-cube law is

    c3(T_a) = c3(T_0) + sum_i delta_i a_i.

Each `delta_i` is odd and is constant on each `H`-orbit. At a chosen
completion `a`, toggling the whole orbit `O_j` changes `c3` by

    Delta_j = delta_(O_j) (|O_j|-2w_j),
    w_j = sum_(i in O_j) a_i.                               (1)

Both factors in (1) are odd. Thus `Delta_j` is nonzero and odd.
Disjoint orbit supports make the `r` changes independent: for any subset
`J` of coordinates, the corresponding cube corner has triangle count
`c3(T_a)+sum_(j in J) Delta_j`.

Let `C` be this free auxiliary `Q_r` packet, and define

    c_min = c3(T_a) + sum_j min(0,Delta_j),
    d_j = |Delta_j|.

These definitions choose the corner minimizing the triangle count, rather
than imposing a sign convention on the orbit labels. The exact generating
polynomial, with multiplicity for classes sharing a triangle count, is

    F_C(q) := sum_(T in C) q^c3(T)
            = q^c_min product_(j=1,...,r) (1+q^d_j).        (2)

Every `d_j` is a positive odd integer, and `c_min>=0` because it is the
actual triangle count of a packet member. The free action proved in the
main theorem ensures that (2) counts `2^r` distinct tournament classes,
not merely `2^r` labelled descriptions.

The multiset of positive steps in (2) does not depend on the initial
corner: changing the corner only negates the corresponding signed steps.
It is a packet invariant, up to permuting coordinates.

## 2. Exact vanishing order and signed moments

Each factor `1+q^d_j` has a simple zero at `q=-1`, since `d_j` is odd.
Consequently `F_C(q)` is divisible by `(1+q)^r`, with **exact** multiplicity
`r` at `q=-1`.

For every integer `0<=s<r`,

    sum_(T in C) (-1)^c3(T) c3(T)^s = 0.                   (3)

Equivalently, the even- and odd-triangle shores have equal triangle
moments of every degree below `r`, including equal counts at degree zero.
This also holds with any polynomial in `c3` of degree below `r`.

One proof applies `(q d/dq)^s` to (2) and evaluates at `q=-1`. For an
explicit sharpness formula, substitute `q=-exp(t)` in (2):

    F_C(-exp(t))
      = (-1)^c_min exp(c_min t) product_j (1-exp(d_j t)).

The first nonzero Taylor coefficient has degree `r`. Hence

    sum_(T in C) (-1)^c3(T) c3(T)^r
      = (-1)^(c_min+r) r! product_j d_j,                  (4)

which is nonzero. The moment-vanishing range is therefore sharp on each
individual packet; additional cancellation can occur when packets are
summed and is not assumed.

## 3. Fixed-blue-degree consequences

Every packet has constant blue degree `B`. For an attained blue value `b`,
let `r_b` be the minimum dimension of a packet with `B=b`. Then the
blue-filtered polynomial

    F_(n,b)(q) = sum_(T in SC_n, B(T)=b) q^c3(T)

is divisible by `(1+q)^r_b`. Its signed triangle moments vanish in every
degree below `r_b`. This is a lower bound on its root multiplicity, not
an assertion of exact multiplicity after different packets are summed.

At `n≡1 (mod4)`, `n>=5`, the main theorem gives positive even cube
dimension, hence `r_b>=2`. Thus for **each fixed blue degree** `b`, both

    #{T : B(T)=b, c3(T) even}
       = #{T : B(T)=b, c3(T) odd},

and

    sum_(B(T)=b, c3(T) even) c3(T)
       = sum_(B(T)=b, c3(T) odd) c3(T).                    (5)

In particular the two uniform-class means of `c3`, conditioned on `B=b`,
are equal whenever that blue value occurs. This is a degree-conditioned
first-moment law, not equality of the full integer triangle distributions.

## 4. Small sharp boundaries

At `n=3`, the unique packet has

    F_C(q)=1+q.

Its shore counts agree, but the shore triangle totals are `0` and `1`.
The first signed moment is `-1`, so an unrestricted all-odd-order
first-moment conclusion would be false.

At `n=5`, the two exact packets have

    B=1: F_C(q)=(1+q)(1+q^3)=1+q+q^3+q^4;
    B=3: F_C(q)=q^3(1+q)^2=q^3+2q^4+q^5.

For `B=1`, each shore has two classes and triangle total `4`; for `B=3`,
each shore has two classes and triangle total `8`. But at `B=1` the even
shore's second moment is `16`, versus `10` on the odd shore. Its signed
second moment is `6`, exactly (4), so the `n≡1 mod4` first-moment law
does not automatically extend to second moments.

## 5. Exact verification and preservation boundary

`verify_odd_blue_moments.py` reads the frozen `pairing_n*.json` files and
imports no producer code. It independently reconstructs adjacency masks,
triangle counts and matching slopes from the canonical frame mask. For
every packet it checks:

- the whole polynomial identity (2), coefficient by coefficient;
- positive odd steps and the actual minimum corner;
- exact division by `(1+q)^r`, and failure of an additional division;
- every vanishing signed moment below `r` and the nonzero value (4);
- each fixed-`B` shore equality implied by its minimum packet dimension;
- the `n=3` first-moment and `n=5,B=1` second-moment hostiles.

All checks pass for:

| n | SC classes | auxiliary cube packets |
|---|---:|---:|
| 3 | 2 | 1 |
| 5 | 8 | 2 |
| 7 | 88 | 14 |
| 9 | 2752 | 184 |

Reproduce alongside the original data files:

    python3 verify_odd_blue_moments.py

Or use `--data-dir PATH` to select their directory. The output is
`odd_blue_moment_verification.json`; it includes input hashes, every
packet's minimum and odd step multiset, its exact signed moments, and all
fixed-blue-degree shore counts and triangle totals.

These are consequences of the previously audited auxiliary cubes.
They do not assert that cube edges are literal antipodal blue edges,
that integer `c3` is preserved by a toggle, or that the actual metagraph
is determined by these polynomial constraints. No even-order extension,
new blue-degree bound, or new exploratory classification is claimed.
