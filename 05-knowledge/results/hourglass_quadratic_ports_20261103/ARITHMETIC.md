# Prime-pair triangles, endpoint clocks and the dimension sequence

**PROVED elementary finite-source laws / CITED Dirichlet input / inherited
machine compilation.** No Goldbach or twin-prime conjecture is resolved.

## P1. Two different addition questions

If two primes sum to a prime, exactly one summand is2: two odd summands
give an even integer larger than2. Thus these triples are exactly
(2,p,p+2), where p and p+2 are odd twin primes. In contrast, every prime
partition of an even M>=6 uses two odd primes. The exceptional even row4
has2+2. The additive step+2 connects consecutive positions in each parity
ladder; it does not preserve primality.

Represent odd labels by a=2i+3,b=2j+3. The region

    0<=i<=j, i+j<=(M-6)/2

contains exactly the unordered odd pairs with sum at most an even M.
Primality selects some cells; their row sum is2(i+j+3). Equivalently use
center c=(a+b)/2 and half-gap h=(b-a)/2. Then h>=0,c>=h+3,
c+h is odd, and the sum is2c. This is a triangular integer region with
a parity sublattice and arithmetic holes, not intrinsically a tournament.
Independent +2 steps on a and b form the four-state commuting diamond in
COMPLETIONS.md H1. In(c,h) coordinates they are(1,-1) and(1,+1).
For h>=1 the whole square remains in the unordered triangular region;
at h=0 a port exchange is needed and the two middle states merge. A wheel
or primality filter can also remove intermediate cells, so a surviving
endpoint alone does not certify both guarded paths through the square.

## P2. What the first +/-1 prime really certifies

For even M>=4 let q_M be the first prime congruent to1 or-1 modulo M.
Its existence follows from Dirichlet's theorem; see the
[translation of the original proof](https://arxiv.org/abs/0808.1408).
Necessarily q_M>=M-1: no smaller positive prime belongs to those classes.

Consequently, by the time an increasing prime stream reaches q_M, it has
already emitted every odd prime that can occur in any even-sum row<=M.
The largest possible summand in such a row is M-3. This proves a finite
coverage statement about the triangular search region, not that every row
has a selected cell. The endpoint event is a sufficient stopping clock;
waiting that long is unnecessary for many finite tasks.

The exact all-row update is particularly simple. Let P(z) be the sum of
z^p over the odd primes seen so far. The unordered pair polynomial is

    G(z)=(P(z)^2+P(z^2))/2.

When a new odd prime p arrives, update with the OLD P:

    G <- G+z^p P+z^(2p),  P <- P+z^p.                  (1)

The diagonal term retains p+p once. Each coefficient counts actual pairs,
and the rule also works for every formal Boolean subset of odd labels.
Working modulo z^(M+1) is exact for all target rows<=M. Once p>M-3, no
future update affects those rows, since all its old partners are at least3.
The bound is sharp for arbitrary source subsets: a new label M-3 can still
create the pair3+(M-3). At M=6 its diagonal does the same.

This is a source-preserving saturation theorem. The coefficient update is
the elementary convolution identity; no claim to invent convolution is made.
It keeps multiplicity, source primality bits and the diagonal correction.

Some exact clocks, independently recomputed here:

|M|q_M|unordered prime pairs summing to M|
|---:|---:|---:|
|10|11|2|
|20|19|2|
|40|41|3|
|78|79|7|
|290|1451|10|

The290 clock and its local pair wheel are inherited from Session146;
they are controls, not discoveries of this packet.

The source heights matter. At M=26 the first endpoint prime is53, but the
earlier primes5 and47 produce residue pair{5,21}. Its residues sum to26
although21 is composite; the actual primes sum to52. Generally, if
p=M u+r and q=M v+s with r+s=M, then p+q=M(u+v+1). The target row M requires
the additional guard u+v=0. A residue pair without this carry/height guard
does not represent an actual partition of M. Incoming mac-mini's
`quadratic_2916_odd_pair_triangles_n4_20261102.md` independently studies this
residue-pair excess; its finite generalization remains explicitly unproved.
Our all-row saturation law retains the original prime labels throughout.

## P3. A complete periodic-mask carrier for every odd squarefree wheel

Let R be any positive odd squarefree modulus. On the odd-label lattice
retain the predicate

    W_R(i,j)=[gcd(2i+3,R)=gcd(2j+3,R)=1].

Under independently advancing either named index and observing W_R, the
exact minimal deterministic carrier has **R^2 states**, the residues
(i mod R,j mod R). This includes all repeated future translations.
For R=1 there is one state.

Proof of minimality: a one-coordinate unit mask has no nonzero translation
period modulo R. For any p|R its projection omits exactly one residue;
a nonzero translation modulo p would force that proper nonempty subset
to be invariant under a full p-cycle. CRT transfers this to R. The support
of all future two-coordinate answers from a state is the Cartesian product
of the two translated unit masks. Its factors are nonempty. Equality of
two such products forces equality of each factor, hence equality of both
residues. Thus every pair of distinct states is distinguished by a future
pair of translations. Restricting to nonnegative advances changes nothing,
since each modular translation has a nonnegative representative.

This minimality statement uses the full ordered odd-pair lattice. Imposing
a finite height ceiling or the unordered inequality i<=j changes the
operation domain and adds a boundary/order guard; R^2 alone does not encode it.

For a fixed half-gap h, the admissible centers modulo R number exactly

    product_(p|R) (p-1 if p|h, otherwise p-2).           (2)

This follows by excluding the one or two residues +/-h at each prime.
Every factor is positive. The row has exact period R by the same projection
argument. This extends the inherited two-prime formula to all finite wheels;
the CRT mechanism itself is classical.

The common parity bit in center/gap coordinates is not two independent
bits: c+h must be odd, and a,b are already odd. Using (i,j) removes that
redundancy exactly. It does not remove the two independent odd residue
phases. Exchanging i,j loses which port the next independent advance acts on.

Wheel admissibility is not primality. A prime dividing R is an exceptional
valid label that the unit mask rejects; restore these exceptions separately.
Composites coprime to R still pass. No fixed periodic mask can equal the
prime indicator: choose a prime p beyond the putative period d, then
p and p(1+d) have the same residue and different primality.

Endpoint clocks need more than the wheel. Moduli10 and20 have the same
odd unit wheel R=5 but different +/-1 clocks. More generally, on odd labels
n=2i+3, the union n=+/-1 modulo M has exact period M/2 for even M>=6:
it is two adjacent residue classes modulo M/2, which have trivial translation
stabilizer when M/2>=3. At M=4 every odd label qualifies, so the period is1.
Prime-power precision cannot be discarded merely because coprimality can.

## P4. What transfers to PRIMEGAME, Turing machines and CA

Session146's `machines/README.md` already gives a terminating fixed-(N,T)
FRACTRAN pair decider preserving its input bits, an explicit one-tape TM,
a shift-equivariant radius-one CA, and a guarded PRIMEGAME stream-to-input
map. We inherit those exact machines. This packet supplies a different
observer: update all finite pair-sum rows by(1), and stop when the chosen
stream horizon is saturated. A macro update need not be one primitive
FRACTRAN, TM or CA step; the stream clock and execution clock remain distinct.

NATIVE.md instead compiles the eight-state quadratic graph into a separate
finite FRACTRAN program and a radius-zero CA. The two machines share the
idea of a legal source encoding, not their transition rules or complexity.
Neither is identified with the unbounded PRIMEGAME program by a shared
prime label or an equal number of states.

## P5. The exact meaning of10,20,39

These values are inherited minimum dimensions of a rational linear carrier
for separate(3,3,5) clique-count means under fixed independent attachments.
They are not prime counts. The existing dimension formula gives a useful
new presentation of its birth law. For s>=3 and n>=3,

    D_s(n+1)-D_s(n)
      =2n+sum_(j=2)^(s-2) binom(n,j)+[n+1=s],           (3)

where binomial coefficients outside range are zero. Pascal's identity
accounts for new proper clique flags; the last term is the first terminal
aggregate at size s. In particular the increments10 then19 produce
10,20,39. From n>=5,

    D_5(n)=4+2 binom(n,2)+binom(n,3)+binom(n,4),

so its fifth finite difference vanishes. The early growth is not a doubling
law. For n>=s the degree is s-1. This recurrence counts required observation
coordinates; it neither counts primes in residue classes nor proves a
Ramsey bound. The shared tool is to state which next operation can inspect
the missing coordinate, then prove the resulting closure or obstruction.
