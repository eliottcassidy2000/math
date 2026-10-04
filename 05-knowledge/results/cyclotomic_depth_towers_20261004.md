# Cubic depth towers: missing binary phases and paired golden ideals

2026-10-04. **PROVED** for the all-height identities and obstructions below;
**FINITE-EXACT** for the six field-order computations, four factored binary
layers and eight golden layers. **OPEN** for all-height primitivity of the
normalized subfield element; no such assertion is needed. These are arithmetic
constructions, with no asserted identification with integer Collatz dynamics.
No historical novelty claim is made.

## Inheritance, target and hostile controls

The incoming [sixth-clock note](sixth_clock_branches_20261004.md), sections1,
3 and5, already proves the exceptional prime support of63, the norm19 pair,
and the order63 of `1+t` in `F2[t]/(t^6+t^3+1)`. Its binary prime-power lifts
keep the degree-six number field fixed. This note changes a different scale:
it takes successive *cubic field extensions*, of dimensions2,6,18,54,... .
Those two towers must not be conflated.

Closest proved mechanism: Frobenius sends a primitive ninth root to its
inverse halfway around its six conjugates. Corrected near miss: a generator
of a field need not generate its multiplicative group. Canonical hostile:
`1+t` at dimension18 misses19. Least-used sidecar: the subfield component
remaining after a cyclotomic phase is removed. The golden comparison inherits
[the prime-ideal versus norm distinction](golden_primes_carries_route_compiler_20261003.md)
and [the finite golden phase decoder](difference_families_20261004.md).

Working board: primitive prime *depth*; exact residue order; relative norm;
half-degree subfield; ideal branch; and retained Collatz address. The first
five have explicit maps here; the last is compared only through the typed
information-loss table, not through an invented dynamical equivalence.

## 1. A norm-compatible cubic tower that stops being primitive

For every integer `a>=1`, put

\[
 n=3^a,\quad m=3^{a-1},\quad q=2^m,\qquad
 K_a=\mathbb F_2[\zeta]/(\zeta^{2m}+\zeta^m+1).
\]

Here the quotient notation names the image of the indeterminate as `zeta`.
The polynomial is `Phi_(3^a)`. The elementary order formula

\[
 \operatorname{ord}_{3^a}(2)=2\,3^{a-1}
\]

follows from `v3(4^s-1)=1+v3(s)` and the evenness of a return exponent.
Consequently the primitive `3^a`-th roots form one Frobenius orbit: the
polynomial is irreducible and `K_a` is the field with `q^2` elements.
Choose embeddings by `zeta_(a-1)=zeta_a^3`; each step has degree3.

Set `b_a=1+zeta`. Because `2^m=-1 mod3^a`,

\[
 \zeta^q=\zeta^{-1},\qquad b_a^q=b_a/\zeta,
 \qquad \boxed{b_a^{q-1}=\zeta^{-1}}.                 \tag{1}
\]

For `h=(3^a-1)/2`, define the normalized element

\[
 c_a=b_a\zeta^h.
\]

Equation(1) and `2h=-1 mod3^a` give `c_a^q=c_a`, so `c_a` belongs to
the multiplicative group of the half-degree subfield `F_q`. The factor
`zeta^(-h)` has exact order `3^a`, coprime to `q-1`. Hence

\[
 \boxed{\operatorname{ord}(b_a)=3^a\operatorname{ord}(c_a)
          \mid 3^a(2^{3^{a-1}}-1)}.                 \tag{2}
\]

The index of this orbit in the full nonzero field is therefore a multiple of

\[
 D_a=\frac{2^{3^{a-1}}+1}{3^a}.                    \tag{3}
\]

`v3(2^m+1)=1+v3(m)=a` shows that this is an integer not divisible by3.
It is1 for `a=1,2` and greater than1 for every `a>=3`. Thus:

**PROVED.** `1+zeta_a` fails to be a multiplicative generator at every
level `a>=3`, although it generates the entire field over `F2` at every
level, because `zeta_a=b_a+1`.

**FINITE-EXACT.** The normalized element `c_a` is primitive in `F_q`
for `1<=a<=6`; at `q=2` this means order1. No all-height statement follows
from these six checks.

| a | Field dimension | Order of b_a | Full-group index |
|---|---:|---:|---:|
|1|2|3|1|
|2|6|63|1|
|3|18|13797|19|
|4|54|10871635887|1657009|
|5|162|587537948332709778907201293|9950006745799417075771|
|6|486|`729*(2^243-1)`|`(2^243+1)/729`|

At dimension18, the full order is `262143=3^3*7*19*73`; the chosen
element has order `13797=3^3*7*73`. It retains73 and misses19.

The final bounded probe reached a=6 with a complete exact certificate:

```text
2^243-1 = 7*73*487*2593*71119*262657*97685839
          *16753783618801*192971705688577*3712990163251158343.
```

All factors occur once. The script supplies explicit Pocklington witnesses
for the last three primes and exhaustive trial division for their supporting
primes; it uses no probable-prime assumption. For clarity, the certificate
principle needed here has a short proof. If `F|p-1`, `F^2>p`, and for each
prime q dividing F a witness satisfies `a^(p-1)=1 modp` and
`gcd(a^((p-1)/q)-1,p)=1`, then every prime divisor r of p has `F|r-1`.
It must exceed sqrt(p), so p is prime. Complete orders are checked by return
and every prime-divisor shortening. This certifies c_6 without factoring the
unused `2^243+1`; the large missing phase index is
`19389268200585836264288587113776883575610248525384021488302948711030121`.

The successful recursive predicate is a relative norm identity:

\[
 \boxed{N_{K_a/K_{a-1}}(b_a)=b_{a-1}}\quad(a\ge2). \tag{4}
\]

Indeed the three conjugates of `zeta_a` over the preceding field are
`zeta_a,omega*zeta_a,omega^2*zeta_a`, where `omega` has order3.
Their `1+` product is `1+zeta_a^3` in characteristic2. Norm compatibility
survives the loss of multiplicative primitivity. Its fibres are not labelled
without retaining a choice upstairs.
The normalized elements are norm-compatible as well: `N(c_a)=c_(a-1)`.
Indeed `N(zeta_a)=zeta_(a-1)` and the difference between their normalizing
exponents is `3^(a-1)`, the order of the lower root. Restricted to the
half-degree subfields this is their degree-three norm; its three Frobenius
conjugates are the same set, possibly in reversed order.

### A lossless repair: retain the missing phase register

The obstruction also supplies a canonical replacement representation.
Let `mu_e` denote the subgroup of e-th roots of unity in K_a. The numbers
`q-1,n,D_a` are pairwise coprime and their product is `q^2-1`. Therefore

\[
 K_a^*\simeq\mathbb F_q^*\times\mu_{3^a}\times\mu_{D_a},
 \qquad(r,v_n,v_D)\longmapsto r v_n v_D.             \tag{4a}
\]

This is a product of specified subgroups, not a choice of primitive element.
For any nonzero x, first compute

\[
 r=(x^{q+1})^{q/2},\qquad v=x/r.
\]

The square of r is the relative norm of x, so `v^(q+1)=1`. For `D_a>1`,
the two CRT projectors are

\[
 v_n=v^{D_a(D_a^{-1}\bmod n)},\qquad
 v_D=v^{n(n^{-1}\bmod D_a)}.
\]

For `D_a=1`, set `v_n=v,v_D=1`. These formulas reconstruct x exactly.
In the case of b_a the three components are

\[
 (r,v_n,v_D)=(c_a,\zeta^{-(3^a-1)/2},1).
\]

Thus the missing phase register is identically1 throughout the b_a orbit.
Allowing r to range over the entire subfield group and retaining both phase
registers restores all nonzero field states. Adding only the last register
to the b_a orbit reaches `<c_a> x mu_n x mu_D`; that is the full group only
when c_a is primitive. The middle register has a ternary address
of length a. Equation(6) below decomposes the order of the third register
into successive primitive-prime layers. In dimension6 its order is1;
in dimension18 it has order19. This is an explicit repair of the recursive
storage scheme at the first failed level.

## 2. The missing cofactor consists of later primitive-prime layers

For `j>=1`, define

\[
 C_j=\Phi_{2\cdot3^j}(2)
     =2^{2\cdot3^{j-1}}-2^{3^{j-1}}+1.
\]

Every `C_j` has exactly one factor3. Substituting
`x=2^(3^(j-1))=-1+3s` gives `x^2-x+1=3(1-3s+3s^2)`.
If a prime `p!=3` divides `C_j`, then `x` has order6 modulo p:
`x^3=-1`, while `x=-1` would force `p=3`. If `r=ord_p(2)`, then
`r/gcd(r,3^(j-1))=6`, forcing

\[
 \boxed{\operatorname{ord}_p(2)=2\cdot3^j}.        \tag{5}
\]

Thus every prime of `C_j/3` is new at exactly that index in `2^s-1`.
For `j=1` the quotient is1; for every `j>=2` it is greater than1.
This proves fresh prime existence on this particular tower without invoking
a global primitive-divisor theorem.

The factorization of `x^3+1` now gives the exact telescope

\[
 2^{3^{a-1}}+1=3\prod_{j=1}^{a-1}C_j,
 \qquad
 \boxed{D_a=\prod_{j=2}^{a-1}\frac{C_j}{3}}.        \tag{6}
\]

Empty products are1, including `a=1,2`. Consequently every successive
primitive-prime layer in(5) becomes an unavoidable missing factor in the
index of the later `b_a` orbit. The factor3 carries depth; the other primes
mark new return orders. This is an all-height mechanism for the distinction.

| j | C_j/3 | Exact order of 2 at each displayed prime |
|---|---:|---:|
|1|1|No new prime|
|2|19|18|
|3|87211|54|
|4|163*135433*272010961|162|

The owner's `63=3*19+6` is the first case `Phi_18(2)=57`; equation(6)
explains why the same19 is precisely the first phase factor missing from
the cubic extension of the 63-state construction. This does not mean the
larger field lacks elements of order19: it means this specified element
cannot visit those multiplicative phases.

## 3. A golden tower of pairwise coprime norms

Let `O=Z[phi]`, with `phi^2=phi+1`. For `k>=1`, put

\[
 d=3^{k-1},\quad u=\phi^d,\quad
 A_k=u^2+u+1,\quad B_k=u^2-u+1.
\]

These are `Phi_(3^k)(phi)` and `Phi_(2*3^k)(phi)`. Since d is odd,
`u'=-u^(-1)` and `L_d=u-u^(-1)` is the odd-index Lucas number. Direct
multiplication, retaining the signed quadratic norm, gives

\[
 N(A_k)=N(B_k)=L_d^2+3=:N_k,\quad
 A_kB_k=N_ku^2,\quad B_k=u^2\overline{A_k}.         \tag{7}
\]

The odd-index identity `L_(3d)=L_d^3+3L_d=L_d*N_k` yields

\[
 \boxed{N_1=4,\qquad N_{k+1}=N_k^3-3N_k^2+3}.     \tag{8}
\]

All terms are1 modulo3. For `l>k`, recurrence(8) gives
`N_l=3 modN_k`: the first step evaluates the polynomial at0, and all
subsequent steps at3, a fixed point. Therefore

\[
 \boxed{\gcd(N_k,N_l)=1\quad(k\ne l)}.            \tag{9}
\]

Every level supplies new rational prime factors. The levels need not be
prime:

\[
 4,\quad19,\quad5779,\quad
 192900153619=3079\cdot62650261,\quad\ldots
\]

The prime assertions in this display, including the two factors of the
fourth term, are certified by exhaustive trial division in the experiment.

## 4. The norm splits into two ideal branches at every level of odd norm

For `k>=2`, N_k is odd and4 modulo5, as follows by induction in(8).
Since `A_k-B_k=2u` and u is a unit, the ideals `(A_k)` and `(B_k)`
are coprime: a common prime ideal would lie over2, impossible for their
odd norms. Equation(7) therefore proves

\[
 (A_k)(B_k)=(N_k),\qquad
 O/(N_k)\simeq O/(A_k)\times O/(B_k).              \tag{10}
\]

This generalizes the incoming norm19 ideal equality to all k, including
composite N_k. It is stronger than recording only the product norm.

Write `A_k=a+b*phi`. Then `gcd(b,N_k)=1`. Otherwise a rational prime p
dividing both would divide a as well, by the norm equation. It would divide
both coefficients of `A_k`, hence also those of `B_k=u^2*conjugate(A_k)`,
contradicting coprimality. Define

\[
 r_k=-a b^{-1}\pmod{N_k},\qquad s_k=1-r_k.
\]

Both solve `X^2-X-1=0 modN_k`, and coefficient reduction gives the explicit
ring decoders

\[
 O/(A_k)\simeq\mathbb Z/N_k,\quad \phi\mapsto r_k;
 \qquad O/(B_k)\simeq\mathbb Z/N_k,\quad\phi\mapsto s_k.
                                                               \tag{11}
\]

The kernel equality uses inclusion plus equal index N_k. The roots differ
by a unit: `(r_k-s_k)^2=5 modN_k` and N_k is coprime to5. Equations(10–11)
are therefore an explicit two-branch CRT representation.

For each rational prime p dividing N_k, the reduction of r_k is a root of
`Phi_(3^k)`, in characteristic different from3. Its exact order is `3^k`.
Similarly s_k has order `2*3^k`. This proves both the splitting of p in O
and the congruence

\[
 p\equiv1\pmod{2\cdot3^k}.
\]

The same two exact orders hold modulo the possibly composite N_k. The
cyclotomic equations give divisibility of the orders, and reduction at any
prime divisor rules out a proper divisor. The range `k>=2` is essential:
N_1=4 belongs to the dyadic base, where the two ideals are not coprime.
In fact `A_1=2*phi^2` and `B_1=2`, so both ideals are `(2)`. The prime2
is inert in O, as `X^2-X-1` is irreducible modulo2; this coincidence of
branches is not ramification of2 in the golden quadratic field.

**Hostile to a stronger binary identification.** At k=2, both the binary
quotient `C_2/3` and the golden norm N_2 equal19. At k=3 they are87211
and5779. Moreover

\[
 \operatorname{ord}_{5779}(2)=5778,\qquad
 \operatorname{ord}_{5779}(r_3)=27,\qquad
 \operatorname{ord}_{5779}(s_3)=54.
\]

The first equality is certified by the factorization `5778=2*3^3*107`
and modular powers; in particular `2^54=2944 mod5779`. Thus the19 match
does not extend to an equality of binary and golden clocks. The survivor is
the shared cyclotomic construction with its base element explicitly named.

## 5. Typed transfers and research consequences

| Source -> target | Map / preserved predicate | Lost information / sidecar | Decisive control |
|---|---|---|---|
| K_a* -> preceding field | Relative norm; sends b_a to b_(a-1) | Kernel element; retain an upstairs choice | Equation(4) and all six finite orders |
| b_a -> half-field and phase | b_a=zeta^(-h)c_a; exact product order | No loss with both factors; c alone forgets phase | Equation(1), ord(zeta)=3^a |
| Full nonzero field -> b orbit | Generated subgroup | At least D_a cosets inaccessible | First failure is index19 in dimension18 |
| Golden ideal pair -> norm | Product gives (N_k) | Which branch has odd/even clock length | Roots r_k,s_k, orders3^k and2*3^k |
| Binary layer -> golden layer | Same cyclotomic indices, changed evaluation base | No clock-preserving map established | k=3:87211 versus5779 |
| Collatz ray -> residue address | Separate inverse-ray construction | Higher ternary digits and certified parent | [Proof-carrying addresses](inverse_ray_ternary_addresses_20261004.md) |

The promising recursive target is now precise: retain *relative-norm
compatibility*, and carry the phase choice separately. Requiring the same
`1+zeta` element to generate every larger multiplicative group is impossible.
For golden storage, use a norm and one of its two ideal branches; the norm
alone merges clocks of different lengths. Neither arithmetic carrier supplies
an integer's missing Collatz return certificate.

## Reproduction and scope

Run:

```text
python -X utf8 -B 04-computation/experiments/cyclotomic_depth_towers_20261004.py
python -O -X utf8 -B 04-computation/experiments/cyclotomic_depth_towers_20261004.py
```

Output: [cyclotomic_depth_towers_20261004.out](cyclotomic_depth_towers_20261004.out).
The standard-library script uses explicit runtime checks, so optimization
cannot delete verification. Its full-group-factorization range is `a=1..5`;
Rabin polynomial tests provide an independent irreducibility path, and full
factorizations are certified by products and trial-division primality.
The separate a=6 probe verifies irreducibility, the complete normalized
subfield order and b_6 order, and relative norms, using the displayed exact
Pocklington certificates for its three larger subfield factors.
Order certificates check return and every prime-divisor shortening.
The three-factor decoder is exhaustive in dimensions2 and6, and is checked
on the first16 nonzero polynomial states and b_a at dimensions18,54,162.
The binary factorizations cover `j=1..4`; golden identities, coprimality,
linear root decoders and branch orders cover `k=1..8`, with the dyadic
exception explicitly separated. All-height conclusions above are proved
algebraically rather than extrapolated from these finite ranges.
