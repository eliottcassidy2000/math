# Session 131: prime shells, Frobenius-odd sectors, and PRIMEGAME labels

## Results first

Let

- A = {2,3,5,7,11,13};
- Btail = {17,19,23,31,47};
- C = {29,41,59,71};
- D = {37,43,53,61,71}, exactly as supplied.

The first three sets are disjoint, and their established meaning is

    π(Fi22) = A,
    π(Baby Monster) = A ∪ Btail,
    π(Monster) = A ∪ Btail ∪ C.

Here π denotes prime divisors of the group order. Moreover A is precisely
{p : p² divides |Monster|}, while {2,3,5,7}, PRIMEGAME's four data labels,
is {p : p² divides |Baby Monster|}. These are exact consequences of the
[Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M),
[Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B), and
[Fi22](https://brauer.maths.qmul.ac.uk/Atlas/v3/group/Fi22) ATLAS orders.
They do not make PRIMEGAME's choice of labels intrinsic.

D is not disjoint from the first three sets: its intersection with their
union is {71}. If a complement among primes at most 71 was intended, that
complement is {37,43,53,61,67}. This is a conditional interpretation, not an
assumption that 71 was a typo. In fact the supplied D has a distinctive
geometric profile:

| p | 37 | 43 | 53 | 61 | 67 | 71 |
|---|---:|---:|---:|---:|---:|---:|
| g(X₀⁺(p)) | 1 | 1 | 1 | 1 | 2 | 0 |
| Frobenius-odd rank | 1 | 1 | 1 | 1 | 2 | 0 |
| Integral eigensum index | 2 | 2 | 2 | 2 | 4 | 1 |

The last two rows are explained in §3, not guessed from numerical labels.

## 1. What “absolute supersingular prime” means

Use the Ogg predicate: **every supersingular j-invariant in characteristic p
belongs to F_p**. This is equivalent to genus zero for
X₀⁺(p)=X₀(p)/⟨w_p⟩, and holds exactly for the fifteen Monster prime divisors
listed above. It is not the supersingular-reduction predicate for one fixed
elliptic curve over Q. The first six Ogg primes are indeed A.
[Duncan–Ono, introduction, equation (1.1)](https://arxiv.org/html/1411.5354v2).

Do not replace X₀⁺(p) by X₀(p): the latter has genus zero at prime levels
2,3,5,7,13, omitting 11. Our checker independently evaluates the classical
X₀ genus formula; at 11 it gives genus 1, while the Fricke quotient has
genus 0.

For p≥5, put s(p)=number of geometric supersingular j-values, f(p)=number
in F_p, and h(D)=proper class number of the imaginary quadratic order of
discriminant D. The exact formulas used are

    s(p) = floor(p/12) + ε(p),
    ε(p) = 0,1,1,2 for p ≡ 1,5,7,11 (mod 12),
    f(p) = (h(-4p) + [p≡3 mod4] h(-p))/2,
    g(X₀⁺(p)) = (s(p)-f(p))/2.

Equivalently, f(p) is h(-4p)/2, h(-p), or 2h(-p) in the residue classes
p≡1 mod4, p≡7 mod8, or p≡3 mod8, respectively. All geometric supersingular
j-values lie in F_{p²}. These counting inputs are classical and recorded by
[Delfs–Galbraith, pp.1–2 and §2](https://arxiv.org/pdf/1310.7789).
The global finite Ogg classification is cited, not proved by our bounded
scan. The p=2,3 cases are cited separately, not silently fed into p≥5 code.

## 2. Canonical orbit–stabilizer reading of A, Btail, C

For a Monster 2A involution t, its centralizer H is the double cover 2.B,
not the simple group B itself. This is stated explicitly by
[Cummins–Matias, introduction p.2 and §5](https://arxiv.org/pdf/1710.01071).
Consequently the 2A conjugacy-class size is

    I = |M|/(2|B|) = 97,239,461,142,009,186,000
      = 2⁴·3⁷·5³·7⁴·11·13²·29·41·59·71.

The exact support identities are

    A     = π(H) ∩ π(I),
    Btail = π(H) \ π(I),
    C     = π(I) \ π(H).

The numerical ratio |M|/|B| is twice I, and must not be presented as the
index of the selected stabilizer. The six small-prime index exponents are
(4,7,3,4,1,2), not (5,7,3,4,1,2).

Here is the complete elementary group-action argument, already proved in
session32's `SYLOW_FIXED_POINT_PARTITION.md`. If G acts transitively with
stabilizer H, a subgroup K fixes a point exactly when K lies in a conjugate
of H. Write a=v_p(|G|), b=v_p(|H|). If b=a, a Sylow subgroup of H is a Sylow
subgroup of G, so every Sylow p-subgroup fixes some point by conjugacy. If
0<b<a, Cauchy's theorem supplies an order-p subgroup fixing a point, but
no full Sylow p-subgroup can fit in a point stabilizer. If b=0, no order-p
subgroup fixes any point. Applied here:

| Shell | Regime on the 2A involution class |
|---|---|
| Btail | every Sylow p-subgroup fixes some 2A involution |
| A | some order-p subgroup does, but no full Sylow p-subgroup does |
| C | no order-p subgroup fixes any 2A involution |

The middle row does **not** assert that every order-p element fixes a point.

## 3. Exact shared mechanism: Frobenius exchange blocks and integral gluing

Let S_p be the finite set of supersingular j-values and F its p-Frobenius.
Since F²=1 on S_p, it consists of f fixed points and g transposed pairs,
where g=(s-f)/2. On the integral permutation lattice L=Z[S_p], this gives

    L ≅ Z_trivial^f ⊕ Z[C₂]^g.

This is an integral decomposition into singleton and exchange blocks. For
one exchanged pair with basis e₁,e₂, the even and odd eigenlattices are
Z(e₁+e₂) and Z(e₁−e₂). Their sum consists precisely of integer coordinate
pairs having the same parity; its determinant/index is 2. Singleton blocks
contribute index 1. Therefore, writing L±=ker(F∓1),

    rank L+ = f+g,     rank L− = g,
    L/(L+⊕L−) ≅ (Z/2)^g,
    [L:L+⊕L−] = 2^g.

Thus the Ogg condition is exactly **vanishing Frobenius-odd sector**, and
its integral even/odd decomposition has no residual index. The rational
splitting alone would miss the index-two contribution of each exchange
pair. This is an elementary consequence, not a new moonshine theorem.

The genuine level refinement is also established. For
X=X₀(N)/⟨specified Atkin–Lehner involutions⟩ and p∤N,
Aricheta proves Q_p(X)/2=g(X^p)−g(X), where Q counts non-F_p supersingular
points. His Table 1 and Appendix A give

    rationality primes of X₀⁺(2): 3,5,7,11,13,17,19,23,31,47;
    rationality primes of X₀*(6): 5,7,11,13.

Here * means **all** Atkin–Lehner involutions, not just w₆. Adjoining the
bad-level primes {2} or {2,3} is bookkeeping and gives π(B) or A; it does
not assert rationality at excluded characteristics.
[Aricheta, Theorem 1.1/3.1 and Appendix A](https://sigma-journal.com/2019/007/sigma19-007.pdf).

The already-established session32 computation gives exactly one exchanged
Frobenius pair on X₀⁺(2) at each p∈C. We rechecked all thirteen Ogg primes
p≥5 by direct F_{p²} enumeration. **Completeness in this test universe:**
point counting produces an E_j/F_p with trace zero for each geometric
supersingular j, and its Frobenius satisfies π_p²=[−p]. Thus every cyclic
2-subgroup is stable over F_{p²} (indeed all 2-torsion points are fixed for
odd p). Every corresponding X₀(2) coarse point and its u-coordinate is
therefore F_{p²}-rational. Noncuspidal points have finite nonzero u, so
enumerating all u∈F_{p²}× is complete, not merely a search restricted to
rational j. The argument is about the chosen trace-zero models, not an
assertion that arbitrary models make every isogeny rational. Thus the
level-two lattice has odd rank 1 and eigensum index 2 at all four C-primes. Its nonrational pair is the roots
of the following polynomial in s=u+4096/u, with
j=(u+256)³/u²:

| p | Irreducible quadratic in s | Parallel-edge j endpoints |
|---|---|---|
| 29 | s²−9s−13 | {2,25} |
| 41 | s²+8s−8 | {3,32} |
| 59 | s²−10s−6 | {28,48} |
| 71 | s²−22s−18 | {41,66} |

All coefficients are modulo p. Frobenius fixes these rational vertices
and exchanges two parallel level-structure edges. This is more precise
than identifying prime labels with graph vertices.

**Boundary for a metagraph comparison.** Every Frobenius permutation lattice
has rank L+≥rank L− and has no isolated integral sign block. Consequently
any independently established metagraph lattice with odd rank greater than
even rank cannot be isomorphic to such a permutation lattice. An exchange
block comparison can be exact while a whole-module identification fails.
No natural isomorphism, tournament action, or moonshine explanation for
metagraph cycles is supplied here. Crucially, “odd rank zero iff eigensum
index one” is special to this **permutation-lattice** family: on
Z_trivial ⊕ Z_sign, the odd rank is 1 but the eigensum index is already 1.
An arbitrary involution lattice can have split sign blocks with no gluing
penalty; index one alone is not an Ogg analogue there.

## 4. PRIMEGAME: four data registers, then six control registers

The printed original Conway list has data registers R={2,3,5,7} and control
registers T={11,13,17,19,23,29}. Its last three fractions are 15/2,1/7,55/1;
the alternative 15/14 variant must not be substituted without naming it.
[Conway, Theorem 1](https://gwern.net/doc/cs/computable/1987-conway.pdf).
The exact overlap table is

| Register set | A | Btail | C |
|---|---|---|---|
| R | {2,3,5,7} | empty | empty |
| T | {11,13} | {17,19,23} | {29} |

In particular the six control primes are not the first six Ogg primes.
The 11/13 reciprocal toggle uses the two primes in A beyond R. These are
valid descriptions of the printed encoding, not a demonstrated group map.

**Prime-renaming conjugacy (inherited session126 theorem).** Let an ordered
FRACTRAN list and seed use a finite prime set P. For any injection σ from P
to the primes, extend σ multiplicatively to Φ on supported positive
rationals. Then x is an integer iff Φ(x) is an integer, and
Φ(f_i x)=Φ(f_i)Φ(x). Hence every applicability guard, its priority, and the
chosen transition are preserved. The entire trajectories are conjugate.
If the output prime is fixed, pure output powers and their exponent
sequence are unchanged. This changes program and state together; it is
not an automorphism of the unmodified printed rack.

Concrete hostile: replace the control prime 29 by 67. Only 29/33 and 77/29
become 67/33 and 77/67. The four-plus-six role split and all instruction
traces survive, but 67 is not Ogg and its residue modulo 8 is 3 rather than
29's 5. The checker verifies all guards and conjugacy for 10,000 updates,
including actual use of the changed register. The universal claim follows
from the proof, not the finite run. Therefore register count or control
parity alone cannot force Monster support or an arithmetic mod-8 label law.

## 5. Inheritance, verification, and limits

The A/Btail/C support and level filtration, actual 2A index, Sylow regimes,
level-two parallel-edge pairs, and prime-renaming theorem were already in
the repository. They are recovered and primary-source checked here, not
claimed as newly discovered. `SOURCE_PINS.json` identifies seven exact
repository blobs at commit `41a9d96041ca6f970158450972c5295695387e53` and
records primary-source locators. The integration is the explicit
Frobenius-even/odd and integral-index formulation, plus the supplied-D
genus fingerprint and focused 29→67 hostile.

Reproduce with Python 3 standard library:

    python check.py > results.json
    python -O check.py > results.optimized.json
    cmp results.json results.optimized.json

45,129 explicit checks pass identically normally and optimized: class-number
and genus evaluation for primes 5≤p<1000; independent elliptic-curve point
counts for 5≤p≤199; level-two enumeration at the thirteen Ogg primes ≥5;
group-order/index/set arithmetic; and 10,000 PRIMEGAME relabelled updates.
No assertions, network access, downloaded executable, external package,
or absolute file path is needed by the checker. Its genus and supersingular
interpretations depend on the precisely cited classical formulas; finite
checking does not prove global classification. No publication performed.
