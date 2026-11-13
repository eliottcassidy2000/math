# An exact inverse at the fixed prime cut, with all 2–3 residues

Keep the actual source of the [prime-derived Session 321 ray](../session321_selective_returns_and_fano_wiring_20261113/SELECTIVE_RAYS.md)
and its cut, rather than switching sources. Its output is
```
N=4595601 Q_V(2) − 1915 H(V).
```
Use an ordinal reservoir of b directed triangles, e copies of the
dyadic factor H=5, f copies of H=19, g copies of H=761, and transitive
padding. The last three factors are in the [witness bank](DYADIC_INVERSE.md).
The elementary factor profiles give
```
N = 3^b (4595601·2^m − 1915 D),
D = 5^e 19^f 761^g,
m >= 2b+5e+8f+14g.
```
Each unit of remaining cost is one singleton. In particular the 761
factor adds a 3-adic sign choice without introducing any odd factor
in Q2. Strongness belongs to the whole output by the fixed quotient.

## Absolute decoder

For exact positive N, b=v_3(N), since 3 divides 4595601 and the
intercept product is a 3-adic unit. Divide N by 3^b. The same strict
shell inequalities as before give the unique m, because D<=2^(m−1).
Check whether (4595601·2^m−N/3^b)/1915 is a positive integer and
factors into the three permitted prime powers. Then check the cost.
If it passes, the test constructs V; if not, this family cannot
produce N. This recovers both the ternary valuation and the dyadic
scale before trying to build the reservoir.

## Modular completeness of this same fixed cut

For every r,s>=0 it reaches every residue of an odd integer modulo
2^r 3^s. Here is an explicit construction; unlike the absolute
decoder, it can use arbitrarily much padding.

Take b=min(v_3(N),s), and require m>=r. When b<s, require also
m≡0 modulo phi(3^(s−b)), so 2^m≡1 there. The required conditions on D
are then
```
D ≡ −(N/3^b)/1915       (mod 2^r),
D ≡ (4595601−N/3^b)/1915 (mod 3^(s−b)).
```
Both prescribed values are units when the respective moduli are
nontrivial. Powers of 5 and 19 first realize the 2-adic residue,
as in the dyadic inverse note. If the sign modulo 3 must be changed,
multiply by 761 times a suitable **even** power of 5 to preserve
the 2-adic residue: 761≡1 mod 8, while even powers of 5 fill that
subgroup. This changes the sign because 761≡−1 mod 3 and the
compensation has sign +1. Finally powers of 5^L, with
L=2^max(r−2,1), fix the 2-adic residue and fill all principal units
at 3. Indeed v_3(5^L−1)=1 by v_3(25−1)=1 and repeated squaring.
This settles D using nonnegative exponents. If b=s there is no
3-adic unit condition. Now choose m in its prescribed phase large
enough to pay the cost. This proves the claim without a density
heuristic.

Through X the number of realized outputs here is O((log X)^5):
the coordinates b,e,f,g,m are all O(log X) by the size inequality
and the cost. Hence this single fixed cut is simultaneously locally
complete at 2 and 3 and absolutely sparse. The [general modular
aperture](RESIDUE_APERTURE.md) allows the modules to change more
freely and reaches arbitrary moduli; it does not erase the
absolute cost condition.

Two absolute examples illustrate target-first selection. At
N=9411754463 (prime), the decoder forces b=0, m=11 and D=19,
so it chooses the order-six factor and three singleton pads.
At N=4705885849=7·672269407, with prime cofactor, it forces
m=10 and D=5, so it chooses the order-four factor and five
singleton pads. These are exact choices in the family, not
a claim about minimal order or the novelty of these spectrum values.
