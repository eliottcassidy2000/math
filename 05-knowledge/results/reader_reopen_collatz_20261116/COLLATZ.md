# The two low digits give a signed descent with an exact bypass of 7

This is arithmetic, not a premise about convergence of ordinary
positive Collatz orbits. Here the input odd integer h is the
desired Hamiltonian path count itself, not a doubled chart label.

For h=4P+1 take the + sign, and for h=4P+3 take the − sign.
The two natural readers then satisfy

| odd input | selected numerator | after two halvings |
|--:|--:|--:|
| 4P+1 | 3h+1 | 3P+1 |
| 4P+3 | 3h−1 | 3P+2 |

The two low digits have **adjacent exact outputs on a shared
high address before further division**. Keeping only positive
Collatz loses this feature: for 4P+3 its first two shortcut
steps give 9P+8, not 3P+2. The two signs here are deliberate.
After completing the cycle to 4(P+1)+1, the signed readout is
3P+4, skipping the multiple of three 3P+3. The actual odd output
is the odd part of the table entry; different extra powers of two
can spoil adjacency. Thus the readout is a useful guard, not an
asserted graph inverse.

Define D(h) to be the odd part of 3h+σ(h), where σ is +1 on
h≡1 mod4 and −1 on h≡3 mod4. For h>1, D(h)<h. It cannot be
21 (indeed it cannot be divisible by 3). A direct descent can,
however, land on 7. Then h=(7·2^e−σ)/3, e≥2, with σ=+1 on
even e and −1 on odd e.

**A uniform detour.** For e≥6 switch the first sign to −σ and
halve once, then twice at each of the next two moves with sign
−σ. The successive odd values are

    u = 7·2^(e−1)−σ,
    w = 21·2^(e−3)−σ,
    z = 63·2^(e−5)−σ.

At each step the stated power of two is the exact valuation:
the first numerator has valuation one, and the next two
have valuation two as e≥6. The last is smaller than h,
since 27h=32z+23σ. The early exceptional cases also bypass 7:
9→13→5, 19→29→11, 37→55→41→31, 75→113→85→1, always
using a suitable sign and then the full odd part. In the last
step 3·85+1=256. All bypasses avoid both forbidden values.
So the signed arithmetic has a strictly decreasing selector
after replacing the immediate 7-landing by a three-step return.
No Collatz conjecture is invoked.

## Every bypass endpoint is already physically attained

For σ=−1 and e≥7, put j=e−5. The one-apex spine word
1001 0^j 100 gives 63·2^j+1: its two fixed contexts have
left address 9 and right reverse-zero address 6. The standard
binary tape update yields the claimed intercept 1 and slope 63.
For σ=+1 with e≥8, take the root-8 tournament with outgoing rows

    (4,5,352,263,271,347,411,319,3).

Its rooted profile is (h,g,R,S)=(181,97,161,188). The
[near-sink reader](BANK.md) therefore gives 252·2^ell−1,
so set ell=e−7≥1 to realize 63·2^(e−5)−1 exactly.
At e=6 the bypass lands at 125, an ordinal product of three
count-5 tournaments. The four shallow bypass endpoints
5,11,31,1 are the one-apex words 100,1100,101010 and a singleton.
The replay verifies both symbolic lanes and their finite controls.

## What would suffice for a proof

The odd signed-inverse branches in this chart have the exact
form (2^e d−σ)/3 with the required residue and sign, for e≥2.
If they can be realized physically for admissible decreasing
parents, except the 7-landing, and if the bypass family
z→(32z+23σ)/27 can be physically lifted (with the small boundary
cases handled explicitly), strong induction from 1 proves
Chan–Pak. Indeed every target except 7 and 21 then has a smaller
realized parent and an exact witness-producing lift.

This separates the jobs precisely: the arithmetic descent is
total on admissible targets, and the exceptional parents are
already real tournaments. The missing piece is an inverse
**construction**, not a convergence argument. The 38-count
physical bank tests reader switching at exact absolute targets
but is sparse; it does not yet transport arbitrary signed
inverse parents or realize the bypass targets uniformly.
