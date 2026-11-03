# Incoming work: exact overlaps and additional interfaces

**CITED current repository results; PROVED elementary corollary below.**
The six source files are pinned at commit f139791ef3 in
[INCOMING_SOURCES.json](INCOMING_SOURCES.json). Their large carrier certificates
were not rerun here. The machine producer reviewed the overlap independently;
the ROOT/lattice producer separately checked the norm corollary. These additions
do not change the four independently audited proof bodies or their receipts.

A later incoming pair-scan note and its exact-definition script are pinned
separately in [LATE_SOURCES.json](LATE_SOURCES.json) at50f4e1fa03. The note is
marked pre-audit, so its finite sweeps are not adopted as audited evidence.
Its elementary one-way inclusion is checked explicitly below.

## The first-prime scan has a precise Goldbach boundary

For even m>=6, let P(m) be the first prime congruent to either1 or-1 modulo m.
Following the incoming [pair-scan definition](../../../04-computation/experiments/circulant_fractran_goldbach_20261101/residue_pairs_posthoc.py),
retain an unordered pair {r,m-r}, with odd3<=r<=m-3, exactly when **both**
residues occur among odd primes strictly below P(m). Call the result R(m),
and let G(m) be the actual odd-prime partitions of m. Then G(m) is contained
in R(m). Every member of a Goldbach pair is at most m-3, while P(m)>=m-1,
so both members appear before the stopping prime.

If m-1 or m+1 is prime, P(m) is that first prime neighbour. In either case
the admitted nontrivial residues are exactly the odd primes at most m-3:
when P=m+1, m-1 is composite. Hence R(m)=G(m). This checks the stated
one-way theorem directly. Merely observing one end is insufficient, and
the converse equality=>prime neighbour is still open in the incoming note.
Its finite strict-containment observations are not promoted to an all-size
statement here. Beyond the first wrap, a prime's residue need not itself
be prime; this is the missing realization guard in the proposed scan.

## Three triangular charts, with different predicates

Session146's [prime290 note](../session146_carrier_primegame_20261103/prime290/README.md)
independently obtains the four CRT masses112,108,84,81 and strengthens the
inherited MRB asymptotic. Its [general wheel theorem](../session146_carrier_primegame_20261103/prime290/WHEEL_GENERAL.md)
classifies all +2 unit paths modulo2pq. At290 there are112 vertices,81 edges
and31 components, with vertex-length histogram {1:2,2:2,3:2,4:25}.
Our rank45 weighted kernel is a complementary operation on the full pair
array. The81 counts a difference-two slice; the56 unit-pair orbits count the
sum290 slice after exchanging its ports. Neither count substitutes for the
other or for the ten actual prime pairs.

The incoming [finite compiler](../session146_carrier_primegame_20261103/machines/README.md)
uses odd labels (2i+3,2j+3), 0<=i<=j<N. Their dictionary to this packet is

    Goldbach chart: (h,d)=(i+j+3,j-i),
    Conway trial chart: (n,d_trial)=(i+j+2,i+1).

The first chart tests two prime labels. The second tests (i+1)|(j+1).
Also i<=j keeps only d_trial<=n/2, whereas the native Conway trial starts
at n-1. This is an exact coordinate correspondence with different guards,
not an identification of unordered candidate enumeration with the ordered
trial dynamics.

Session146 supplies an actual finite FRACTRAN/TM/strict-CA family for each
(N,T), and obtains its primality vector from inherited PRIMEGAME phases.
Our M4 instead conjugates each literal step of the fixed Conway list on its
full typed integer carrier. Its one-token value-plane locality obstruction
does not apply to Session146's tape/head CA encoding. The consume/restore
bit test and our v2/priority guards are compatible examples of information
that a legal next instruction must retain.

## Two quadratic maps using145 states or residues

The incoming odd-numerator map F(a)=(a^2-29)/4 descends modulo290 to145
states; it has eight3-cycles and the stated finite basin census. Our trace
map D(t)=t^2-2 instead acts on56 primitive145 trace states, in four14-cycles.
The rational -29/16 cycle passes through numerator5, a nonunit modulo290.
Thus neither state count nor the occurrence of145 identifies these systems.
The incoming bounded-lift theorem applies exactly to x=a/4 with odd |a|<=7;
its divergent lift9->13->35->299 is a hostile to inferring boundedness from
a finite modular cycle.

## An exact consequence for28,145 and290

Write C_q(t)=q^3+1-t^3+3qt. In the incoming
[ARITHMETIC.md, A1](../chart_clique_norm_20261103/ARITHMETIC.md), integer
collisions of C_q are classified by Eisenstein norms:

    C_q(u)=C_q(t), u!=t  iff  u^2+ut+t^2=3q.

The explicit map N(b+c omega)=q -> (t,u)=(b+c,b-2c) supplies the collision;
the third root is -t-u. Combining that theorem with our E6 gives:

**Corollary (PROVED).** C_145 and C_290 are injective on the integers.
C_28 has exactly two nonsingleton fibres, both of size three:

    C_28^{-1}(22113)={-10,2,8},
    C_28^{-1}(21793)={-8,-2,10}.

Indeed145=5*29 and290=2*5*29 contain inert primes to odd exponent, so
neither is an integral Eisenstein norm. For28=2^2*7 the inert exponent is
even and the split exponent is one. Thus A(28)=sum_(d|28)chi3(d)=2 and,
since28 is not square, A1 gives exactly two triple fibres and no two-point
fibres. The norm representation N(6+2 omega)=28 gives (t,u)=(8,2), hence
the first displayed triple; negation gives the other. Direct substitution
gives the displayed values. Every other nonempty integer fibre is a singleton.

This gives28 a different role from the degree28 of the unweighted sieve
quotient: the map above, rather than equality of integers, is the evidence.
The modular norm on Eisenstein units modulo290 remains surjective with540
elements per fibre. Neither145 nor290 is a unit modulo290, so that statement
does not apply to their residues. A precise local/global hostile is17: it is
not an integral Eisenstein norm, yet Q(1,18)=307=17 mod290 gives a modular
unit-norm preimage. The integral/modular distinction therefore becomes an
exact injectivity/realization boundary for the cubic. It loses the exact
integer norm and height, not just a sign or an orbit label.

The independent reviewer corrected two draft fibre-value arithmetic slips;
the triples and classification were unchanged. Reproduce the exact values
and a complete small collision census with:

```sh
python3 - <<'PY'
from collections import defaultdict
from math import isqrt
for q in (28,145,290):
    fibres=defaultdict(list)
    for t in range(-isqrt(4*q),isqrt(4*q)+1):
        fibres[q**3+1-t**3+3*q*t].append(t)
    print(q, sorted((v,ts) for v,ts in fibres.items() if len(ts)>1))
print('local norm:',1-18+18**2)
PY
```

## Tournament compression after the order-ten obstruction

The latest [Session146 route](../../../00-navigation/S146.md) gives full fixed,
source-oblivious attachment-response dimensions2668,13028,65908,345106 at
orders7,8,9,10. Its independently audited
[order-ten upper proof](../session146_carrier_primegame_20261103/carrier10_upper/PROOF.md)
finds2520 invisible directions and embeds them at every higher order:

    d(n)<=c(n)-2520 binom(n,10), n>=10.

Consequently the pooled-parity carrier is not minimal from ten onward.
This is a genuine new quotient to consider before proposing another scalar
compression. It does not erase joint incidence needed for source-dependent
operations, shared-source products or native complements. Our explicit
A000568 recurrence counts unlabelled objects, a different question from the
minimal dimension of these attachment observations. Its growing arithmetic
coefficient cost, native P debt and global black-core connectivity stay open.
