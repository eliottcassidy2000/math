# Fixed-count growth concentrates on long extension lifetimes

Date: 2026-11-11. **PROVED:** the exact H9 ordinary/labelled lifetime
fractions, the fixed-power concentration theorem, and its first asymptotic
correction. **INHERITED:** the complete H9 strong alphabet, the SCC word
enumerator, and the all-gap triangle filtration. No new tournament census
or universal target-coverage claim is made.

## Three proved inputs meet at an actual constructor

1. [PROOF.md](PROOF.md), §2, counts fixed-H ordinary classes and labelled
   mass by ordered non-singleton strong-component words. The word of r
   components and total size m contributes `binom(n−m+r,r)` ordinary classes,
   and this quantity divided by the product of its component automorphism
   orders to the labelled mass `L_n(H)/n!`.
2. Incoming checkpoint `b305e1d03f` proves the complete triangle-chain law in
   [TRIANGLE_ALL_K_PROOF.md](../../berggren_target_return_20261111/arithmetic/TRIANGLE_ALL_K_PROOF.md),
   with the constructor and graph-side interpretation in
   [TRIANGLE_FILTRATION.md](../../berggren_target_return_20261111/fibres/TRIANGLE_FILTRATION.md).
   For a triangles and singleton gaps, let ν be the maximum matching size
   of the path whose edge i exists iff internal gap i is zero. Then

       C_k(T)=3^max(a−k−min(k,ν),0),    lifetime L(T)=a−ν.       (1)

3. [KERNEL_SWITCH.md](../signed_arcs/KERNEL_SWITCH.md) recovers the two
   distinct strong H9 classes in the first25-prime graph's completion cube.
   The same root cut extends them to actual counts23 and25, respectively,
   with identical `(E,U,d)=(2,50,9)`. Each child count is coprime to3.

Here `C_k(T)` is the gcd of H(U) over **all** k-vertex induced extensions U
that retain T, and `L(T)=min{k:C_k(T)=1}`. These are the inherited P13
definitions, not the scalar minimum among a selected extension menu.
When H(T)=3^a, adjoining k sinks supplies the value3^a, so C_k divides3^a.
Consequently L(T) is exactly the first number of new vertices permitting
an extension count not divisible by3. Every extension with fewer new vertices
still has a factor3; at the lifetime there is an actual factor-removing child.

The live comparison is thus ordinary fibre growth / automorphism weight /
ordered gaps / extension gcd / literal one-birth switch. The hostile source
is two triangles with a positive singleton gap: H=9 alone cannot distinguish
it from either strong H9 core, but their first extension contents differ.

## 1. Almost every fixed-3^a source has exact lifetime a

Fix a positive integer a. At order n, let `P^ord_n` mean uniform ordinary
isomorphism classes with H=3^a. Let `P^lab_n` mean uniform tournaments on
the fixed label set[n] with that count. Converse classes are not identified.

**Concentration theorem.** For either convention,

    P_n[L(T)=a] →1,    and therefore P_n[L(T)≥a] →1.            (2)

The exact known subfamily witnessing both limits consists of

    TT_(g0) → C3 → TT_(g1) → ... → C3 → TT_(ga),
    g1,...,g_(a−1)≥1,   g0,ga≥0.

Its number of ordinary classes is

    S_a(n)=binom(n−3a+1,a),      n≥4a−1,                       (3)

and zero below that order. Its labelled mass is `S_a(n)/3^a`, or
`n! S_a(n)/3^a` actual labelled tournaments. Every member has lifetime a.

*Proof.* The non-singleton strong components of a count3^a tournament have
counts3^e with e≥1. Their number is at most a. Equality forces every exponent
to be1, and the unique strong count3 component is C3. Thus the only
maximal-length core word is a ordered triangles. The finite-word formula
from input1 gives

    N_n(3^a)=n^a/a!+O_a(n^(a−1)),
    L_n(3^a)/n!=n^a/(3^a a!)+O_a(n^(a−1)).                    (4)

For(3), reserve one singleton in each of the a−1 interior gaps. Distribute
the remaining n−4a+1 singletons among the a+1 ordered gaps. Stars and bars
gives `binom(n−3a+1,a)`. Unique SCC order prevents collisions. Each triangle
has automorphism order3, and no automorphism permutes SCC positions, so the
total automorphism order is3^a. Every interior gap is positive, hence ν=0
and(1) gives exact lifetime a. Finally(3) has the same leading coefficients
as(4) in both counting conventions. Its proportion tends to1; the two
events in(2) contain it and have probability at most1. QED.

This proves concentration for fixed a as n grows. It makes no uniform
assertion when a grows with n, and does not bound the lifetime of every
exceptional source by a.

## 2. A sharper law: 3a(a−1) versus19a(a−1)

For each fixed a≥2, the same merger gives

    P^ord_n[L(T)=a] =1−3a(a−1)/n+O_a(n^−2),
    P^lab_n[L(T)=a] =1−19a(a−1)/n+O_a(n^−2).                 (5)

Both formulas remain true with the event L(T)≥a. In particular
`P_n[L(T)>a]=O_a(n^−2)` in either convention. This last estimate is a
proportion statement, not a universal lifetime bound.

*Proof.* Write

    T_a(n)=binom(n−2a,a),       B_a(n)=binom(n−2a,a−1).

T_a counts all a-triangle words, including zero internal gaps. A word with
a−1 non-singleton components and product3^a must have exactly one H9
component and a−2 triangles. There are two strong H9 classes, both rigid,
and a−1 possible positions. Its total strong size is3a−1. All remaining
core words have at most a−2 components, hence total contribution O_a(n^(a−2)).
The word formula therefore gives

    N_n(3^a)=T_a(n)+2(a−1)B_a(n)+O_a(n^(a−2)),
    L_n(3^a)/n!=T_a(n)/3^a
             +2(a−1)B_a(n)/3^(a−2)+O_a(n^(a−2)).          (6)

Every a-triangle source outside(3) has a zero internal gap, so ν≥1 and its
lifetime is at most a−1 by(1). Every one-H9 word also has lifetime at most
a−1: extend that component by its explicit23 or25 child, and extend each
of the a−2 triangles to the strong order4 tournament of count5. Retain all
singleton gaps and all ordinal incidences. These a−1 births retain the
original source and have count23·5^(a−2) or25·5^(a−2), coprime to3.
Thus none of the degree-(a−1) exceptional families belongs to L≥a.

The difference between T_a(n) and S_a(n) has leading coefficient

    T_a(n)−S_a(n)
      =(a−1)n^(a−1)/(a−1)!+O_a(n^(a−2)).                  (7)

The H9 words add twice this leading coefficient in ordinary counting.
For labelled mass they carry weight3^−(a−2), which is **nine times** the
triangle-word weight3^−a. Hence the respective leading deficits are

    3(a−1)n^(a−1)/(a−1)!,
    19(a−1)n^(a−1)/(3^a(a−1)!).                          (8)

Divide by(4). Sources with shorter core words contribute only O(n^−2) to
either conditional probability; their individual lifetimes need not be
classified. This proves(5) for both events, and their difference gives the
claimed O(n^−2) upper tail. QED.

The constants3 and19 have an explicit origin: one lost-gap contribution
plus two H9 core choices; their labelled weights turn `1+2` into `1+2·9`.
They are not fitted numerical rates.

## 3. Exact H9 fractions at every order n≥7

The complete inherited H9 alphabet has just three core words:
the two rigid strong five-vertex cores and the word(C3,C3). Thus

    N_n(9)=2(n−4)+binom(n−4,2)=(n−4)(n−1)/2,
    L_n(9)/n!=2(n−4)+binom(n−4,2)/9
             =(n−4)(n+31)/18.                           (9)

For either strong core, its23 or25 child remains a valid one-birth extension
after arbitrary uniform transitive prefix and suffix padding. A sink birth
also retains count9. Their gcd is1, so all such padded sources have lifetime1.

For two triangles, the incoming law says lifetime1 when the interior gap
vanishes and lifetime2 when it is positive. Consequently **exactly**
`binom(n−5,2)` ordinary classes have lifetime2, each of automorphism order9.
All other H9 classes have lifetime1. Hence, for every n≥7,

    P^ord_n[L=2]=(n−5)(n−6)/((n−4)(n−1)),
    P^lab_n[L=2]=(n−5)(n−6)/((n−4)(n+31)).               (10)

The exact complementary counts are

    ordinary lifetime1: 3n−13,
    labelled lifetime1: n!(19n−77)/9.                   (11)

For orientation, at n=7 the fractions in(10) are1/9 and1/57; at n=8 they
are3/14 and1/26. Both tend to1, with deficits6/n and38/n respectively,
as predicted by(5). The labelled population initially favours the rigid
strong-core paddings, but its degree-two triangle population ultimately
dominates as well.

## 4. Meaning for constructor selection

An abundance of sources at a fixed count is compatible with increasing
concentration on sources whose **every** small induced extension retains
an unwanted divisor. At count3^a, a typical large-order representative
requires a new vertices before this divisibility barrier can disappear.
The rare rigid H9 cores are nevertheless useful: their distinct overlap
graphs expose the same one-arc23/25 actuator and generate the cheaper
one-birth factors used in the proof of(5).

This is a concrete growth-versus-usability law. The map is ordered SCC
word→gap vector→zero-gap matching→exact extension-gcd lifetime, with
automorphism weight retained for labelled counting. Scalar H loses the
gap vector; ordinary multiplicity loses the weight; the explicit child
witness supplies the constructive step that neither count alone implies.

The lifetime records a divisor-removal interface, not a prescribed target
count and not the separate 1-ready condition in[RIGIDITY.md](RIGIDITY.md).
The theorem does not assert that a gcd1 extension menu contains every desired
odd count. Chan–Pak completeness remains open.

## 5. Cheap exact controls

No graph census is required. The following self-contained algebra check was
run from the repository root. It compares the H9 formulas with the inherited
orders1–8 CSV, checks both rational identities for n=7,...,200, and enumerates
small gap vectors to verify the separated stars-and-bars numerator against
the matching lifetime. Explicit exceptions keep checks active under `-O`.

```python
import csv
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path

def require(p):
    if not p:
        raise RuntimeError('extension-lifetime control failed')

def compositions(s, k):
    if k == 1:
        yield (s,)
    else:
        for x in range(s + 1):
            for tail in compositions(s - x, k - 1):
                yield (x,) + tail

gap_cases = 0
for a in range(1, 5):
    for n in range(3*a, 4*a + 9):
        count = 0
        for gaps in compositions(n - 3*a, a + 1):
            run, nu = 1, 0
            for z in gaps[1:-1]:
                if z == 0:
                    run += 1
                else:
                    nu += run // 2
                    run = 1
            nu += run // 2
            count += a - nu == a
            gap_cases += 1
        expected = comb(n - 3*a + 1, a) if n >= 4*a - 1 else 0
        require(count == expected)

for n in range(7, 201):
    all_ord = 2*(n-4) + comb(n-4, 2)
    all_mass = F(2*(n-4)) + F(comb(n-4, 2), 9)
    sep = comb(n-5, 2)
    require(F(sep, all_ord) == F((n-5)*(n-6), (n-4)*(n-1)))
    require(F(sep, 9)/all_mass == F((n-5)*(n-6), (n-4)*(n+31)))
    require(all_ord-sep == 3*n-13)
    require(all_mass-F(sep, 9) == F(19*n-77, 9))

root = Path('05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity')
with (root/'spectrum_n1_n8.csv').open() as f:
    for row in csv.DictReader(f):
        if row['H'] == '9' and int(row['n']) >= 5:
            n = int(row['n'])
            require(int(row['ordinary_classes']) == (n-4)*(n-1)//2)
            require(int(row['labelled_tournaments']) == factorial(n)*F((n-4)*(n+31), 18))
print({'passed': True, 'gap_vectors': gap_cases, 'H9_orders': 194})
```

These controls verify the new counting assembly. The all-depth filtration
and completeness of the H9 alphabet are the named proved inputs, not claims
inferred from the small gap enumeration.
