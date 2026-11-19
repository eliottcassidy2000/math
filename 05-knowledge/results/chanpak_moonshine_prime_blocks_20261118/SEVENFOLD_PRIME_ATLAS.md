# Prime targets: exact small atlas and five complementary module families

## What is still missing?

There are three different questions: does p occur, does 7p occur,
and can they be coupled by adding one vertex with all old arcs unchanged?
A miss of one of the templates below does not answer the first two.

The corrected independently reviewed S203 prefix already certifies
every odd count up to 4,294,967,295 except 7 and 21. Consequently all
odd primes p other than 3 and 7 with p <=613,566,756 have both required
witnesses (possibly unrelated). The inherited prime-frontier packet
certifies the next 1,000 odd counts beyond that cap as well, so its
joint-prefix bound is p<=613,567,042. Even the first prime beyond that
bound, 613,567,049, is not a gap: the actual inherited F1 construction
with addresses (73,37888599) yields 4,294,969,343 = 7·613,567,049 on
35 vertices, computed by the exact join/path-load replay. In fact the
same inherited F1/F2 constructions yield explicit upper witnesses for
the first 16 consecutive primes past that named bound, through
613,567,369. Fifteen use F1; the F1 miss at upper count 4,294,970,827
uses F2 left rows [28,1,2,6,14] and address 134983805. All sixteen
rows and addresses are in frontier_guard_witnesses.json. Their lower
primes are already certified by S203. The repository also reports a
longer single-run prefix, distinct in audit status. These are named
finite guarantees, not a claim to know the first unresolved prime.
Source: https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md .

## The five prime-bearing forms tested here

Write (h,P)=(H(M),P_2(M)) for an actual nonempty module M. The rules below
physically realize 7p for p=a h+bP; some also realize p as a marked
vertex deletion. The all-order module formulas do not assert that
infinitely many outputs are prime.

| ID | (a,b) | What the construction gives | Prime inputs <=4620 |
|---|---|---|---|
| 1 | (13,18) | exact single-vertex extension, original matched quotient | 157,641,1217,4253 |
| 2 | (25,284) | inherited lower and upper hosts, not deletion | 593,2297 |
| 3 | (9,79) | inherited strong upper reader | 167,641 |
| 4 | (27,127) | inherited strong upper reader | 281,4091 |
| 5 | (5,34) | new exact single-vertex extension | 73,277,1093,4357 |

ID 2 is the inherited PRIME_RAIL.md: its lower visit coefficients
are (593,284) and its upper coefficients are seven times that.
IDs 3 and 4 are the two physical right readers in
chanpak_contacts/FERMAT_TARGET.md; their formulas for arbitrary
modules are 7(9h+79P) and 7(27h+127P). No same-host parent is asserted
for IDs 3 and 4. Their p values remain separately certified in the
finite range.

ID 5 has child outgoing row masks
[268,173,256,100,367,389,295,349,10]; delete vertex 7 and substitute
the module at vertex 2, retaining labels. Its exact parent visit
coefficients are (73,34,0,0) and child coefficients
(511,238,0,0,0). The argument of UNIFORM_SEVENFOLD_EXTENSION.md
sums the counts of nonempty ordered path covers over the marked
visit slots. Q_1=h and Q_2=P−2h give p=5h+34P and 7p as claimed.
The standalone replay independently checks all coefficients,
literal deletions and path counts. Singleton and TT3 modules give
73→511 and 277→1939 respectively. Both inputs miss ID 1 exactly:
the singleton gives 49 there, and other modules would require
h<=p/49 and h≡61p mod72. These are incompatible at 73 and at 277.

A second new matched quotient is [226,260,473,99,491,6,34,106,233],
delete vertex 1, mark vertex 4. It gives p=3h+48P and its sevenfold
extension; this is a composite-only family, because p is a multiple
of 3 and at least 99. It is therefore not one of the five prime forms.

## Exact finite scope

Since P_2(M)>=2^|M| and h>=1, every module for these five forms with
p<=4620 has at most seven vertices (even the smallest b is 18:
18·2^8+13=4621, and the other forms have larger eight-vertex minima).
The atlas partitions all labelled tournaments by permutation
orbits in orders 1–7, and computes all induced path counts with
subset DP; summing products over complementary subsets gives P_2.
It recovers class counts 1,1,2,4,12,56,456. There are 620 admissible
prime inputs (odd primes excluding 3 and 7) through 4620. Thirteen
distinct inputs hit these five forms; 607 do not. The exact list
of all hits with concrete module masks is printed by the atlas;
the CSV lists all 620 with their matched form IDs, marking an
empty match as a bank miss, never as a spectrum gap. All 620
are already known to have unrelated witnesses where needed.

Compile and run:
  g++ -O3 -std=c++17 sevenfold_prime_atlas.cpp -o /tmp/primes
  /tmp/primes sevenfold_prime_status.csv

Run python3 matched_readers_replay.py for the three literal
one-vertex rules (including the composite-only rule).

## All-order classification of the parameterized families

Removing singleton strong components of M leaves a core word B.
In an ordinal sum, H multiplies, and P_2 also multiplies: split
each of the two paths at each seam, or concatenate corresponding
slots even if empty. Thus a fixed nonempty B gives the exact ray

    p_j = a H(B) + b P_2(B) · 2^j,   j>=0.

For the empty B the ray is a+b·2^j, j>=1. Conversely every nonempty
module appears on one such ray, by removing its singleton strong
components; inserting the j singletons back at the end realizes
the same profile. If gcd(aH(B),bP_2(B))>1, that ray has no prime
terms. Others are prime-candidate families, not prime infinitude
claims. In particular the old and new all-transitive rays
18·2^m+13 and 34·2^m+5 are disjoint: for both exponents at least 3
they lie in different classes mod16, and substituting m=1,2 on
either side gives no equality. The new transitive ray has prime
candidates only at odd m, since even m forces a factor 3.
Nontransitive M adds different rays and is essential to include,
as the 641 and 4253 examples show.

A universal classification of the remaining primes by these or
other physical families is still open. The finite misses show
why a bank of proportional (h,P_2) readouts needs a complementary
selector, not just longer modules on the same five forms.
