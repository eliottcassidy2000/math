# Two shifted divisor charts and a primitive Pythagorean cut

These are literal classes and inverse constraints, not a
classification of the full spectrum. Write (u,R),(v,S) for
the Hamilton and two-cover reads of physical arms A,B.
In C3[A,B,1] put d=R−u, e=S−v. The inherited read is

    H = de + 2uv.

A 7-carried right arm is one with 7 dividing both v and S;
write w=v/7, so a prescribed target H=7p reads
p=d(e/7)+2uw. The exact three-read ×7 replacement
from the [transverse packet](../transverse_lattice_prime_bits_20261117/README.md)
matches **three** independent cover values and lasts two
common wraps. The charts here only solve a scalar
factorization with two marked gaps; they neither imply
that full match nor follow from it.

## 1. A local mod-3 obstruction

The odd-pack polynomial has the parity of the vertex
order, so R≡(−1)^|A|u mod 3, and similarly for B.
If u or v is zero mod 3 then H is zero there.
Otherwise d is zero mod 3 at even order and u mod 3 at
odd order; likewise e is zero or v. Thus a prime quotient
p>3 forces at least one even arm order and 3 dividing de.
A p²qr divisor shape with all three factors prime to 3
cannot even enter this chart.

## 2. Repeated ordinary prime: de=3·7·ell²

Take a prime ell≥11 distinct from 3,7. A 7-carried v is
at least 35 (7 and 21 forbidden). Cutting full Hamilton
paths injectively into two-slot covers gives R≥(|A|+1)u,
and similarly S≥(|B|+1)v. Since no tournament of order
at most four has 35 paths, B has order at least 5 and
e≥5v≥175. Thus e=7,21 cannot occur. The only potential
factor pairs (d,e) are

    (1,21ell²), (3,7ell²), (ell,21ell), (3ell,7ell).

The inverse at a prescribed p is: take ell with
3ell²<p, factor K=(p−3ell²)/2=uw if this is an integer,
choose one of the four pairs, and test the **physical**
profiles (u,u+d) and (7K/u,7K/u+e). For d=1 the first
arm is the singleton TT_1. For d=3 it must be TT_2:
u≤3, and u=3 would violate d≥|A|u. This is an exact
smaller-profile inverse, not existence of a divisor with
the required joint reads.

For any arm of order at least two, its two-cover read
is divisible by 4: complementary subsets pair the
terms, all odd, with an even number of pairs. Thus
u≡−d and w≡−e/7 mod 4 when A is not the singleton.
Then p≡3+2·3=1 mod 8. The d=1 singleton instead has
w≡1 mod4 and p≡5 mod8. Each potential direction has
a literal prime witness:

| d | ell | quotient p |
|---:|---:|---:|
| 1 | 47 | 7069 |
| 3 | 47 | 6833 |
| ell | 4877 | 71,521,273 |
| 3ell | 1181 | 4,199,473 |

For the last row the arm profiles are (217,3760)
and (245,8512), giving gaps (3543,8267). The two
central complementary divisors cannot simply swap
while the H values stay on their own arms. The proposed
first read R'=217+8267=8484 would be zero mod3, but
its H=217 is not. No tournament at any order has that
new profile. Swapping **both entire arms** is allowed
and is a different operation.

This is a physical nonlift of the
[p²qr divisor-plane reflection](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md).
Ten proper divisors are seven squarefree plus the two
mixed nonsquarefree ones and one pure square. The
three proper primes already lie inside the squarefree
seven; F=U+S is an equality of counts, not disjoint
union. At p³ the analogous equality is 2=1+1, and at
p there are no proper nontrivial divisors. Larger odd
prime powers have F>U+S. Ordinal powers of a physical
p-path block realize p^k but raise **every** P(q) to
the kth power, not to the same scalar multiple. These
are different operations from the compensated ×7 pair.

## 3. Repeated 7: de=3·7²·r

Let r≥11 be prime, with r≡1 mod4 for the complementary
quotient channel. Again e≥175; the 7-divisible
divisors not containing r (7,21,49,147) are ruled out.
This leaves d=1,3,7,21. At d=21 a non-singleton A has
u≡3 mod4, u≤21 and 3 not dividing u. Counts 7 and 21
are forbidden, and u=11,19 give d≥|A|u≥3u>21.
Thus d=21 fails. (A singleton has d=1.)
At d=7, possible u≡1 mod4 are 1 or 5; the latter
violates d≥3u. The forced left arms are therefore
TT_1,TT_2,TT_3 at d=1,3,7 respectively.

Now p=21r+2w. The congruence on w yields
p≡3r mod8 at d=1 and p≡7r mod8 at d=3,7.
With r≡1 mod4 these give both 3 and 7 modulo 8.
Literal examples:

| d | r | quotient p |
|---:|---:|---:|
| 1 | 173 | 3823 |
| 3 | 197 | 4243 |
| 7 | 181 | 3851 |

At a prescribed p, choose one of these three d and a prime r
in the interval

    5 d p/[21(5d+2)] ≤ r < p/21

of the correct residue. Set w=(p−21r)/2 and require the
actual right profile (7w,7w+147r/d). The interval follows
from e≥5v; it guarantees that the scalar is smaller
(w≤p/(5d+2)). The unproved step is choosing r so that
this **joint profile** exists. The two charts reach all
four odd classes modulo eight via examples, not all
primes within each class.

With both arms transitive there is an independent
all-order wall: d and e are Mersenne numbers of residues
0,1,3 modulo 7. Their products miss −2. Thus
C3[TT_a,TT_b,1] never reads even a single multiple of 7.
One genuinely nontransitive arm matters.

## 4. The primitive Pythagorean lane

If d=s², e=t² for relatively prime odd s,t, the shifted
count is (st)². The three numbers
((t²−s²)/2, st, (t²+s²)/2) are a primitive
Pythagorean triple (for t>s). Conversely a primitive
triple with even leg has c±leg the two odd square gaps
in this chart. To choose primitive gap roots for a fixed
shifted square W², the full power of each prime in W
must go wholly to one root: there are 2^omega(W)
labelled assignments. Thus for W=ell²qr the middle
exponent-one divisor layer is forbidden by coprimality,
even though it exists in the divisor plane.

The checked physical example has

    (u,R)=(319,5080),   (v,S)=(679,14840),
    d=69², e=119²,     (even leg,odd leg,c)=(4700,8211,9461).

It gives H=7·9,693,389 with prime quotient.
The right arm carries 7. More generally when both arms
are nonsingletons and the right gap is a square
carrying 7, the shifted quotient de/7≡7 mod8.
Then since u w≡de/7 mod4, the p residue is 5 mod8.
This is complementary to some of the nonsingleton
repeated-ordinary-prime lanes, not an all-residue
Pythagorean selector. The literal run-length
(3k,4k,5k) family in the [blue-gate analysis](README.md)
never hits the 7-channel at all.

[replay.py](replay.py) verifies every listed small
arm by literal subset DP and every quotient by exact
trial division. Chan–Pak still requires a universal
source-choice lemma, not just nonempty infinite
families apart from the known exclusions.
