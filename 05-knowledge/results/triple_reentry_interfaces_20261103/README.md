# Ternary return forests and nonuniform interfaces over every modulus

2026-11-03. **PROVED** for the finite-bank arithmetic compiler and the
six-root scalar forest; **FINITE-EXACT** for the interface minima and literal
seed tournaments; **OPEN** for universal source reentry and the assertion
that 7 and 21 are the only missing positive odd raw Hamiltonian-path counts.
No external novelty claim is made. The independent audit is a second
algorithmic path run by the same agent, not an external referee.

The useful new split is between **describing all interface responses** and
**finding one interface satisfying several arithmetic requirements at once**.
Quadratic interpolation solves the first. It does not make a common CRT
witness sparse. The six-root address theorem meanwhile turns the vague
ternary-return aspiration into a concrete sufficient closure target.

- [RETURN_FOREST](RETURN_FOREST.md): exact roots `1,23,25,63,65,67`, literal
  witnesses, and the source-return hypothesis that remains unproved.
- [INCOMING](INCOMING.md): newer repeated scalar returns, universal guard entry,
  and the measured normalization price on the six actual seeds.
- [TRIPLE_ATLAS](TRIPLE_ATLAS.md): older triples, exact maps, lost coordinates,
  and the transfer that is useful for this problem.
- [RESEARCH_CYCLES](RESEARCH_CYCLES.md): inheritance, hostile tests, metanalysis,
  and next decisive experiments.

## P4. Prime-support compiler for any finite named interface bank

Fix an actual tournament T and a finite family A of legal, named attachments.
For `s in A`, write `c_s>0` for the raw Hamiltonian-path count of the resulting
tournament. Define the prime-support hyperedge

    P_s = {p prime : p divides c_s}.

For every positive modulus m, all interfaces fail the unit test precisely
when `supp(m)` is a transversal of these hyperedges:

    [for every s, gcd(c_s,m)>1]
       iff [for every s, supp(m) intersects P_s].             (1)

Consequently the inclusion-minimal transversals L give all minimal
squarefree blocking moduli `product_(p in L) p`. An interface is a unit
modulo m iff **none** of these minimal blocking moduli divides m.

Proof: expand the definition of gcd in (1), then shrink a finite hitting
set to an inclusion-minimal one. Primes outside the union of the P_s are
irrelevant. If some c_s is 1, its hyperedge is empty and there is no blocking
modulus. If A is empty, the empty transversal gives blocking modulus 1.
All counts here are odd, so the prime 2 never obstructs. Prime exponents
are irrelevant to this predicate; exact residues and valuations need them.

For a budget b, retain only interfaces of cost at most b and compute its
blocker antichain `B_b`. This gives an exact **interface cost profile** for
all moduli simultaneously:

    minimum cost of a unit interface
      = min{b : no d in B_b divides m},                       (2)

with infinity if no such b exists. Keep the named cut beside each count
so that a successful arithmetic test produces an actual interface.
If two attachment constraints share a port, intersect their allowed named
cuts; intersecting only their sets of count values can invent a witness.

This is a compiler for a finite, authenticated bank. Coprimality itself
is not claimed necessary or sufficient for the missing tournament return.
It is a sharp test of what simultaneous arithmetic conditions require.

## P5. Small determining probes need not contain a unit witness

For a one-vertex birth, S is the set of old vertices pointing **into** the
new vertex. The inherited boundary formula is an integer quadratic:

    R(S)=H(T)+sum_(i in S) a_i-sum_({i,j} subset S) b_ij.

Thus values at `|S|<=2` determine every R(S), over the integers and after
reduction modulo any modulus. For a prime p, if R is nonzero at some cut
then one such small probe is nonzero modulo p: otherwise all its Möbius
coefficients vanish. For a prime power the unit test is still just the
same test modulo p.

For a composite modulus, these small witnesses for different primes can
be different cuts. No common unit follows, even if the gcd of all probe
values is 1. The following fixed-source examples make that boundary exact.

### Fixed incoming baseline

The strong six-vertex source has outgoing rows

    (32,21,1,7,13,30),          H=27.

Modulo `1335=3*5*89`, every cut of size at most two fails the unit test.
But `S={0,1,2}` gives `R(S)=131`, a unit. The complete minimal blockers for
the small bank are

    3*5*89,       3*11*13*23*89.

These classify that bank for **every modulus**, not just 1335. The full
64-cut cube has precisely three minimal blockers. With

    B=3*71*79*89*101*109*131,

they are

    35B,          715B,          3289B.                       (3)

All these statements are independently reconstructed in `audit.json`.

**Polarity matters.** This example has a unit cut with five incoming
vertices, hence only one outgoing exception to the other uniform baseline.
It does not prove that three exceptions are intrinsically necessary when
both baselines are allowed.

### Either uniform baseline

Use cost `min(|S|,n-|S|)`. For source rows

    (32,9,3,5,15,30),          H=27,

put

    M=72,320,955 = 3*5*7*71*89*109.

All cuts of cost at most two fail. The balanced cut `S={1,2,4}` has count
79, coprime to M. The small-bank minimal blockers are exactly

    3*5*7*71*89*109,
    3*5*11*71*89*109,
    3*11*13*23*71*89*109.

The full-bank blockers are again (3). This also shows that the full image
of the response does not determine its cost profile: the two six-vertex
examples have the same count image but different placements of values
among the named cuts. The receipt checks that equality explicitly.

### Three exceptions are not a universal bound either

For the strong eight-vertex source

    (244,61,232,1,236,8,170,42),          H=197,

the receipt supplies a 55-prime squarefree modulus M. Every cut with
`min(|S|,8-|S|)<=3` has nonunit count. The balanced mask 212, namely
`{2,4,6,7}`, has count 617 coprime to M. This is a literal counterexample
to a universal three-exception claim, including both baseline orientations.
The large modulus is printed in `receipt.json`, with its minimal blocker
certificate; it is not a claim of global minimality at order eight.
The source was found at the first draw of a deterministic random search
(seed 163164, 28-bit code 163817466); its frozen literal rows, not a
probabilistic assertion, are the certificate.

## Exhaustive finite boundary

The producer scans all labelled tournaments through order five and all
22,320 strong labelled tournaments at order six. A separate algorithm
partitions the universes into isomorphism orbits, counts new-vertex paths
by direct permutations, and recomputes minima using a different hitting-set
algorithm. The minima below range over **all positive moduli**, because
factor supports give a finite complete search, not a truncated modulus scan.

| Source universe | Cut cost | Smallest modulus where some cut is a unit but no cost<=2 cut is |
|---|---|---:|
|All sources, orders 1–4|fixed incoming baseline|none|
|All order-five sources|fixed incoming baseline|255|
|Strong order-five sources|fixed incoming baseline|none|
|Strong order-six sources|fixed incoming baseline|1335|
|Strong order-six sources|either uniform baseline|72,320,955|

For order five the smallest witness is source `(4,1,2,7,15)`, H=3;
mask19 gives13 modulo255. Its low blockers are `{3,5,17}` and its full
blocker is `{3,5,13,17,23}`. Both-baseline cost at most two already includes
every cut through order five, so no such obstruction can occur there.

Completeness of the modulus search: for a chosen prospective unit value v,
remove `supp(v)` from every low-bank prime edge. An empty residual edge
makes that v impossible. Otherwise every feasible modulus contains a
transversal of the remaining finite edges, and the least product occurs at
an inclusion-minimal transversal. Enumerate v in the finite response image.
No larger unseen prime or prime power can improve a minimum.

## P6. A nonuniform bound from forbidden output values

Let a Boolean-cube response R be an integer polynomial of degree at most d.
Fix a modulus m, and let b be the number of *distinct nonunit values* in its
finite image. If some cut is a unit, a unit cut exists with at most `d*b`
incoming vertices (and hence at most that many exceptions to the nearer
uniform baseline).

Proof: form `Q(z)=product_(a nonunit image value)(z-a)`. On the Boolean cube,
`Q(R(S))` is nonzero exactly at the unit cuts. Its multilinear reduction
has degree at most d*b. A nonzero multilinear polynomial over Q has a
nonzero value on a set of size at most its degree: choose a minimal-cardinality
nonzero evaluation; its top finite-difference coefficient equals that
value. Otherwise every coefficient of permitted degree would vanish.
For b=0 the empty cut is already a unit. For one birth, d<=2.

This elementary bound can be loose; it pays for the actual forbidden-value
alphabet rather than assuming a universal width of two or three. It assumes
that a unit exists and supplies no universal source-return theorem.

## P7. Parameter-dependent interfaces have a joint clock for every modulus

Suppose a finite family of named legal interfaces has counts
`c_s(k)=b_s+sum_j a_sj lambda_j^k`, with integer coefficients and fixed
integer bases `lambda_j>=2`. Retain the actual parameter k, names and costs.
For a modulus `m=product p^e`, define

    K=max ceil(e/v_p(lambda_j)) over p dividing lambda_j,
    L=lcm ord_(p^e)(lambda_j) over p not dividing lambda_j,

with empty maximum0 and empty lcm1. Then the entire count vector modulo m,
and hence its minimum-unit-cost profile, is periodic with period dividing
L for k>=K. Proof: modulo each prime power a nonunit base has vanished
by K and a unit base repeats with its multiplicative order; apply CRT to
the **same** k and interface labels.

For unit decisions alone, replacing m by rad(m) improves the certified
preperiod to at most1 and replaces every prime-power order by its prime
order. Exact residues, divisions and rooted observer values generally
need the larger clock. This gives a finite procedure for each modulus,
not one fixed clock uniform over all moduli. A periodic selector still
must use the actual source with its original k; changing k freely can
change the integer count being proved.

The incoming ready-return observers contain both2^k and3^k. Their extra
seven-adic phase bit is therefore relevant data, even where one term is
invisible modulo7. The [incoming rooted controller](../rooted_carry_closure_20261103/ROOTED_RETURN.md)
also gives a concrete positive use of shared selection: either fixed port
fails infinitely often, while selecting a at even k and b at odd k always
passes the stated primitive-return guards. These are precisely the sorts
of named joint constraints that separate modulus tests can lose.

## A concrete F=U+S bridge, with its boundary

The three exponent promotions of 1335 are

    4005=3^2*5*89,   6675=3*5^2*89,   118815=3*5*89^2.

Each is in the inherited balanced family `p^2qr`: with proper nontrivial
divisors, `(F,S,U)=(10,7,3)`. All four numbers, including 1335 itself,
have the identical unit-interface predicate by (1). But 1335 has
`(F,S,U)=(6,6,3)` and is not balanced. Thus the map to prime support
preserves unit compatibility while losing the valuation coordinate that
selects the balanced divisor family. This is an exact map, not an
identification of divisor balance with tournament return.

## Audit, corrections, and reproduction

Run from the repository root:

```sh
python3 -B 05-knowledge/results/triple_reentry_interfaces_20261103/verify.py
python3 -B 05-knowledge/results/triple_reentry_interfaces_20261103/seed_lifts.py
python3 -B 05-knowledge/results/triple_reentry_interfaces_20261103/verify.py --census
python3 -B 05-knowledge/results/triple_reentry_interfaces_20261103/independent_audit.py
python3 -O -B 05-knowledge/results/triple_reentry_interfaces_20261103/independent_audit.py
```

Use `--write` only to deliberately regenerate a producer receipt. The
independent audit does not import producer code. It verifies all 416 literal
cuts, all source-orbit minima through the declared order-six universe,
actual strongness of the unit outputs, all six root tournaments, and exact
inverse/replay addresses for every allowed odd target through 100001.
The address check is a finite control for P1, not tournament coverage.
P7 also has exact mixed-base clock controls for every modulus1..128,
including the singular preperiod failure at72. The incoming extension
checks every rooted slack at all six literal seeds and independently counts
the best-port one-ear and two-ear graphs.

The survey repaired THM-1425's path-order shift:
`I(P_m,2)=J_(m+2)`, not `J_(m+1)`. It also replaced THM-202's incomplete
historical direct proof with the noncircular consequence of the all-order
H=21 exclusion, retaining the old argument as incomplete provenance.
See the dated entry in `01-canon/MISTAKES.md`. `SOURCE_PINS.json` identifies
inherited versions; `RESEARCH_CYCLES.md` names the remaining obligations.
