# Prime layers, shared exponent phases, and signed affine construction

**PROVED IN NOTE:** a recursive relation-lattice compiler at every modulus
coprime to6; its exact common-quotient gluing rule; a one-digit prime-power
lift; the remaining ternary contraction layer at every odd modulus; and
the signed extension of the inherited guarded-ancestor construction.
**FINITE-EXACT:** the declared controls in `receipt.json`.
**OPEN:** descent from every supplied positive source and exhaustion of
the integer cycles on either sign sheet. No literature-priority claim.

The positive construction replaces an entire multiplier table by three
integers that retain its synchronized 2/3 exponent relations. Prime powers
are built one digit at a time. Coprime factors are glued over an explicit
shared phase, rather than treated as independent choices. The affine carry
is then retained as a separate coordinate. This goes beyond the inherited
row-image, entropy, central-block and Smith-factor formulas.

## 1. Inheritance and the useful older mechanisms

Closest proved mechanism: [the October4 guarded affine lift](../../collatz_affine_guarded_lifts_20261004.md)
identifies the correct finite affine group and realizes its elements by
legal inverse paths from a supplied hub. Its canonical hostile is95:
the two local multiplier groups have product order72, but the synchronized
global multiplier group has order36. The corrected near miss is that CRT
on residue rings does not make a chosen pair of generators independent.
The least-used sidecar is the lattice of **one common pair of exponents**.

The live board is prime-power digits / common-quotient gluing / translated
carry / exact source guard / signed sheet / physical operation.

Several older, less-direct lanes suggest the same constructive repair:

* [September21 return-period seams](../../arithmetic_seams_20260921_primes.md)
  retain first-hit offset and prime-power clock, not only prime support.
  Its sharp fifteen-term coprimality window comes from actual modular
  targets; the number fifteen alone does no work.
* [The Monster shared-quotient note](../../cyclic_attachment_lattices_20261023/MONSTER.md)
  retains one common C2 value when gluing two surjective local actions.
  The map below is a literal arithmetic fibre product, not a claim that
  the Monster acts on Collatz trajectories.
* [Normal Frobenius-port transport](../../relation_tower_frobenius_ports_20261028/THEORY.md),
  T2–T3, recovers joint phases from an injective correlation carrier and
  glues along consistent shared ports. Here relation lattices provide
  the shared phases; no finite-field normality assumption is imported.
* [The October19 Carmichael/carry bridge](../../carmichael_collatz_carry_20261019/INHERITANCE_AND_BRIDGES.md)
  proves that scalar group characters erase the target-bearing affine
  translation. Our compiler therefore keeps that translation explicitly.
* [The prime-batch triangle note](../../prime_batch_triangle_20261020/README.md)
  distinguishes synchronized cyclic actions from independent factors.
  Its ordered-block envelope also retains the actual source and clock;
  the legal construction in section6 keeps the corresponding guards.

The [current center/entropy packet](../../modular_centers_collatz_pell_20261111/PROOFS.md)
and [golden content packet](../../central_modular_chanpak_20261110/information/PROOF.md)
already handle the central entries, singular rows, general CRT tables,
and Smith information loss. Those results are inputs, not new claims here.

## 2. A three-integer object for synchronized multipliers

Let m>=1 with gcd(m,6)=1, and define

    H_m = <2,3> <= (Z/mZ)^*,
    L_m = {(u,v) in Z² : 2^u 3^v =1 mod m}.         (1)

Negative powers mean modular inverses. The case m=1 is trivial. The
homomorphism from Z² onto H_m has kernel L_m, so

    H_m ≅ Z²/L_m.                                  (2)

The lattice has a unique basis of the form

    (h,0), (c,v),   h,v>0, 0<=c<h.                 (3)

Thus three integers `(h,c,v)` retain **all** exponent relations, not
just the order `hv=|H_m|`. An exponent pair (u,w) belongs to L_m iff

    v divides w,   h divides u−c(w/v).              (4)

The representative set `0<=u<h, 0<=w<v` gives every generated multiplier
exactly once. To reduce any exponent pair, first subtract `(c,v)` to
bring its second coordinate into that interval, then subtract `(h,0)`.

### Constructing the first prime layer

For a prime p>=5, let h=ord_p(2). Let v be the least positive integer
with `3^v in <2> mod p`, and choose c in [0,h) with `2^c 3^v=1`.
Then (3) is exactly L_p. Indeed any relation's second coordinate is a
multiple of v in the quotient `H_p/<2>`; subtracting that multiple of
(c,v) leaves a relation involving only2, hence a multiple of(h,0).
This is a finite prime-field calculation. It does not assume that2 or3
individually generates the whole unit group.

### One higher prime digit is one linear equation

Suppose L_(p^e) is already known, e>=1. On that lattice define

    λ_e(u,v) = (2^u 3^v−1)/p^e mod p.               (5)

The numerator is interpreted by evaluating modulo p^(e+1), then taking
its divisible difference from1. This makes (5) well-defined even with
negative exponents. Products `(1+p^e a)(1+p^e b)` show that λ_e is an
additive homomorphism to F_p. Therefore

    L_(p^(e+1)) = ker λ_e within L_(p^e).           (6)

In a basis b1,b2, compute the two digits λ1,λ2. If both vanish, keep the
basis. Otherwise interchange them if necessary so λ1!=0, put
`k=λ2/λ1 mod p`, and use

    p b1,   b2−k b1.                               (7)

These generate exactly the kernel. Normalizing them back to (3) completes
the recursive step. Each digit either preserves the group order or
multiplies it by p. The *direction* of the new constraint is retained.

An exact example shows why that direction matters:

    L_1093       = (364,104,1),
    λ_1         = (0,781) on its displayed basis,
    L_(1093²)   = (364,104,1093).                   (8)

The order of2 stays364 while the joint multiplier group grows by1093.
Watching only the base2 order would miss the new phase completely.
At5 the first four layers are `(4,1,1),(20,13,1),(100,93,1),(500,393,1)`.
These are literal examples of (6), not a universal ordinary-lift assumption.

There is also a compact order formula. Put

    t_p=min(v_p(2^(p−1)−1), v_p(3^(p−1)−1)).

Then

    |H_(p^e)| = |H_p| p^max(0,e−t_p).              (9)

To prove this, raising every element of H_(p^e) to p−1 extracts exactly
its principal-unit subgroup: the prime-to-p quotient has order dividing
p−1, and the power is invertible on the p-group kernel. The image is
generated by `2^(p−1),3^(p−1)`. Both are1 modulo p^t_p, and one has
exactly that depth. The binomial identity increases its valuation by
exactly one on every pth power, so it generates all p^max(0,e−t_p)
principal units at that depth. This proves (9). Formula (9) retains size;
the lattice additionally retains which exponent combinations realize it.

## 3. CRT is exact gluing over a computable common quotient

For coprime m,n, the complete recursion is

    L_(mn)=L_m intersect L_n.                      (10)

More precisely there is an exact sequence of abelian groups

    0 -> Z²/(L_m∩L_n)
      -> Z²/L_m × Z²/L_n
      -> Z²/(L_m+L_n) -> 0.                        (11)

The first map uses the same exponent pair in both factors. The final map
sends `(z_m+L_m,z_n+L_n)` to `z_m−z_n+(L_m+L_n)`. It is well-defined
and onto. Its kernel consists exactly of pairs whose representatives
can be altered by local relations to become equal, proving exactness.

Thus two local multiplier choices have a common global realization iff
their images in

    Q_(m,n)=Z²/(L_m+L_n)                            (12)

agree. The compatibility quotient's order is the exact index of the
actual global subgroup in the product of the local groups:

    |H_m||H_n| / |H_(mn)| = |Q_(m,n)|.             (13)

This is a constructive theorem: if `z_m−z_n=l_m+l_n`, then
`z=z_m−l_m=z_n+l_n` is an actual common exponent pair. Retaining the
decomposition gives a witness, not merely a cardinality test.

### Explicit lattice intersection without a large table

Let the two records be `(h1,c1,v1),(h2,c2,v2)`. Put

    l=lcm(v1,v2), d=gcd(h1,h2),
    δ=c1(l/v1)−c2(l/v2),
    v=l d/gcd(d,δ), h=lcm(h1,h2).                  (14)

Choose c modulo h solving

    c=c1(v/v1) mod h1,
    c=c2(v/v2) mod h2.                             (15)

The congruences are compatible by (14), and the new record is `(h,c,v)`.
To prove minimality of v, any common second coordinate is lt, and its
two first-coordinate congruences are compatible exactly when d divides
δt. Equations (14)-(15) give the least such t and all remaining relations.

The invariant factors of Q are equally small to compute. Form the four
columns `(h1,0),(c1,v1),(h2,0),(c2,v2)`. Let D1 be the gcd of all their
entries and D2 the gcd of all their2-by2 minors. Then

    Q ≅ C_(D1) × C_(D2/D1).                        (16)

This is the rank-two Smith calculation. Trivial C1 factors are omitted.
It is an information-bearing quotient, not an extra independent factor.

### A two-phase and a six-phase example

At95 the local records are `(4,1,1)` for5 and `(18,5,1)` for19.
Their common quotient is C2, the merged record is `(36,5,1)`, and the
actual group has36 elements rather than72. This recovers the older
quadratic-character obstruction by a general recursive procedure.

At1729=7·13·19, the local records are

    L_7=(3,2,2), L_13=(12,8,1), L_19=(18,5,1).

The first two glue freely to `(12,0,6)`, of index72. Gluing19 then
requires agreement in C6 and gives

    L_1729=(36,12,6), |H_1729|=216,
    |H_7||H_13||H_19|=1296.                        (17)

The six is a cyclic shared phase. It is not inferred from3!, a sixfold
picture, or the number of known cycles. Here its exact source, map and
compatibility predicate are all specified. Likewise the recovered
fifteen-prime/Ogg material motivates keeping levels and joint lifts;
it does not provide a fifteen-prime cutoff for these recursions.

## 4. Translation is the coordinate that makes this an affine compiler

For gcd(m,6)=1, the inherited finite group theorem gives

    G_m = (Z/mZ,+) semidirect H_m.                  (18)

Represent an affine map by `(u,v;b)`, meaning
`x -> 2^u 3^v x+b mod m`, and reduce (u,v) using L_m. Composition is

    (u,v;b) ◦ (u',v';b')
      = (u+u',v+v'; b+2^u3^v b') mod m.            (19)

The multiplier relation and translation are separate. Under CRT the
translation coordinates glue freely, while multiplier coordinates must
pass (12). Formula (19), prime-digit lifting, and the lattice intersection
form an exact recursive object closed under the actual finite affine
operations. They retain the order-dependent intercept which scalar
characters of G_m erase.

For a binary shortcut word w of length L with k ones, and sign ε=±1,
write its affine action as

    T_(ε,w)(x)=(3^k x+εB_w)/2^L,
    B_empty=0,
    B_(w0)=B_w, B_(w1)=3B_w+2^length(w).            (20)

Its state in (19) is `(-L,k; ε2^(−L)B_w)`. At a cycle clock
`m=|2^L−3^k|>1`, the slope is automatically the identity. The cycle
divisibility question lies entirely in the retained translation.
The lattice alone consequently cannot exclude cycles at their own
clocks. It makes their synchronized modular programs constructible
without pretending that an erased coordinate was irrelevant.

## 5. The exceptional prime3 completes the object at every odd modulus

Every odd modulus is uniquely `m=3^a q`, with gcd(q,6)=1. On the q part,
use the lattice-plus-translation compiler. On the3^a part, (20) has a
controlled contraction. Put d=3^min(k,a). Then its image is exactly

    {y mod m : 2^L y=εB_w mod d},                  (21)

and every attained value has exactly d preimages modulo m.

Indeed `3^k x=2^L y−εB_w mod m` is a linear congruence with gcd(3^k,m)=d.
Its image criterion and fibre size give (21); CRT gives the stated
separation. After k>=a odd branches, the3-primary output is independent
of the starting residue, while the q-primary map remains a permutation.
The actual surviving3-primary value still depends on the ordered carry.

This is a whole-word recursion, distinct from inspecting one product row.
It distinguishes reversible prime directions from a specific repeatedly
consumed prime direction, at every odd modulus. Both signs have exactly
the same fibre cardinalities; their translated images retain the sign.
No comparison of these finite fibre counts proves that one integer sign
sheet has only a specified number of cycles.

## 6. A positive signed consumer: legal ancestors with any finite action

The October4 construction in section1 is stated for the plus map. Its
proof extends to both signs, which gives this compiler an actual integer
construction consumer rather than only a modular observer.

**Signed guarded-lift theorem.** Fix ε in{−1,+1}, a positive odd hub u
coprime to3, a modulus m coprime to6, and any g∈G_m. There are infinitely
many positive odd ancestors n>u whose legal inverse shortcut route from
u to n has affine action g modulo m. They may be chosen coprime to3,
with every odd vertex before the final forward arrival at u greater than u.
Any supplied first-hit route from u to a specified root or cycle splices
to these ancestors. The theorem does not supply that route for u.

Here are the details altered or retained by the sign. Use inverse letters

    D(x)=2x,   E_ε(x)=(2x−ε)/3.

The finite affine group they generate is still (18), since their
commutator is translation by−ε/3, a generator of the additive residue
group. A finite positive inverse word with length L, e E-letters and
carry B has the exact form

    W(x)=(2^L x−εB)/3^e,
    D: B<-2B;  E: B<-2B+3^e, e<-e+1.               (22)

Its legal integer-input guard is `2^L x=εB mod3^e`. Reducing this guard
to each prefix proves intermediate integrality. Both signs preserve
positivity on their legal positive inputs; E_ε(x)>=x/3. When e>0 the
guard residue is a3-unit.

As in the inherited proof, choose an Euler circuit in the positive
D/E_ε Cayley digraph of G_m. Its E-ending prefixes visit every group
element, and one final extra E supplies a3-unit endpoint guard. Since2
generates the units modulo every3^e, a preliminary doubling count k
can make2^k u satisfy the whole-word guard. Add a joint multiple of
the ternary exponent period and ord_m(2) until k>=2L. All subsequent
states then exceed u, since2^k/3^L>1. Right multiplication by the
initial D^k merely permutes the group elements represented by the
marked prefixes, so one prefix supplies the prescribed g. Increasing k
by its joint period produces infinitely many ancestors with the same
finite action and first-hit property.

Thus the recursive lattice solver selects compatible finite actions,
and a sign-correct guard realizes them on actual integer paths. What it
constructs is an ancestor of the supplied hub, often very large. It does
not turn an arbitrary supplied unresolved integer into that ancestor or
establish a decreasing step from that integer.

## 7. Reproduction, exact scope, and the next constructive question

```bash
python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/prime_layers/compile.py
```

The program uses integer arithmetic only. It implements the prime-field
record, the linear digit lift, rank-two normalization, exact CRT lattice
intersection, and the shared quotient's invariant factors. Independent
direct multiplier generation checks every admissible odd modulus5..199.
There are19074 explicit signed-exponent membership comparisons on
`[-8,8]^2`, plus five named prime-power towers. All binary words of
lengths1..6, both signs, on every odd modulus3..45 give5544 independent
literal-map checks of (21). Source moduli, lattice steps, digit
functionals and quotient invariants are saved in `receipt.json`.

The source-target map is now explicit: prime-field relations -> one-digit
linear constraints -> compatible global exponent classes -> affine carry
state -> guarded ancestor of a supplied hub. Preserved predicates are
the exact finite action and legal path. The essential data not supplied
by these operations are the prescribed starting integer and a decreasing
height for it. A useful next theorem would couple the lattice-selected
action to a source-preserving descent or a bounded-cost reset. Another
count of prime rows, a center identity, or a list of known cycles does
not provide that missing operation.
