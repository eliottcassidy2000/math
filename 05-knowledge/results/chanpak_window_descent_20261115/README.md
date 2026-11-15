# Exact seams, growing windows, and the first missing count

The common principle is **descent on the retained interface, not on
the readout alone**. A useful operation must solve the *absolute*
integer discrepancy and preserve whatever the next operation really
reads. Forgetting either condition explains why a construction can
have unbounded height, full modular coverage, or an exact \(+2\) at one
target and still not prove [Chan–Pak](https://arxiv.org/html/2308.10214v2).

## Main advance: the windows now really grow

[Two sliding cores](GROWING_WINDOWS.md) give arbitrarily long
intervals of consecutive odd counts in **strong** tournaments of
the same order and with the same full boundary data. A 13-count
spine can slide repeatedly by 10 in the actual apex count and a
29-count spine by 6. Each returns the same interface. Gluing the
two gives half-steps 145 and 39, so bounded integer carries fill
a genuine interval. The first displayed long interval has 183
consecutive odds at a common order; the width grows at least
linearly in the order. The accompanying absolute decoder gives
the actual two slide positions, not just residues. It still
does not reach an arbitrary prescribed target. Trying every
transitive suffix on every pump height remains extremely thin,
by a proved density bound, so recentering must use a genuinely
different exterior reservoir.

## More results

1. **[Exact CRT amplifier](CRT_AMPLIFIER.md), [general quotient fibers](CRT_WINDOWS.md).**
   For coprime ordinary spine counts \(h_i\), put \(P=\prod h_i\).
   Fix the actual baseline \(N_0\) and divide an odd target difference
   by two. Modulo each \(h_i\) there is one forced local digit residue.
   After choosing those residues, the remaining integer carry must
   lie in the sum of the local quotient sets. This is necessary and
   sufficient for that physical bank, not merely a residue test.
   At uniform local half-width \(2h_i-2\) the output inherits both the
   matching interface and the same relative width threshold. There
   is no such amplification with repeated prime factors in one
   fixed independent-switch fiber. Complete windows still need moving
   **physical baselines**, not an existential offset.
2. **[Physical tests](ONE_SIDED_STITCH.md).** A count-3 block cannot
   supply a two-sided digit, but only its right data are needed when
   it is the first piece. Face it toward a count-5 bank matched only
   on its left. At every height \(k\ge1\) this gives five consecutive
   odds \(36765\cdot2^k+\{265,267,269,271,273\}\), and separately every
   odd residue modulo 30 at that same height. The next offset 275 is
   missing from this very bank although its residue occurs. This is
   a controlled experiment separating a residue selector, an absolute
   window, and a usable continuation. [The full two-sided 67/69
   switch](CRT_WINDOWS.md) is possible at multiplier 5; see also the
   [triangle reconstruction theorem](TRIANGLE_NAF.md) and the sharp
   [all-order width bound for multiplier 5](FIVE_CORE_WIDTH.md).
   Independent full DP also gives a [5/9 Bézout seam](COPRIME_SWITCH.md).
3. **[The bounded-ratio exposure theorem](EXPOSED_INTERFACE.md).**
   If a donor of order \(n\) and count \(h\) is put into any host
   as a single uniform module, either the output is a product or
   it is at least \((n+2)h\). The host can change with every target:
   the obstruction is pointwise. A recurrence that reuses a large
   donor at bounded ratio must therefore expose part of it to the
   host or rebuild it, if it is to produce new multiplicatively
   indecomposable targets.
4. **[First missing count](PRIMITIVE_LANES.md) and [the 7/21 seams](SMALL_EXCLUSIONS.md).**
   The earlier factor argument reduces the least additional hole to
   a prime or seven times a prime. Finite dyadic fan banks miss
   infinitely many paired \(p,7p\) even under ordinal products.
   Meanwhile the two known exclusions can be checked in only
   16 and 128 seam orientations after Busch's overlapping-core
   lemma. The same decomposition is being used in two directions:
   exclusion restricts all possible small responses; construction
   must guarantee at least one response at the required excess.
5. **[Periodicity route and its limit](REGULARITY_BOUNDARY.md).**
   Given the earlier physical congruence theorem, eventual
   periodicity of the *exact* spectrum would force cofiniteness.
   Cobham makes exact recognition in two independent bases a
   possible route. But even an explicit simple profinitely full
   one-apex palette is nonautomatic in every base. Two separate
   modular controllers are not two exact recognizers.

## How the other problems sharpen the strategy

* **Frobenius/Apéry.** The two-block lower endpoint is precisely
  the Frobenius conductor in normalized half-offsets. The numbers
  3 and 5 have Frobenius hole 7, but this does **not** explain the
  tournament hole 7: the endpoints are shifted by the literal
  baseline, and multiplier 3 does not even admit the required
  two-sided switch. What transfers is the method of keeping each
  residue's least and greatest *integer representatives*. Bounds
  on both yield a genuine finite window, unlike a periodic sieve.
* **Berggren and square-sum absorbers.** In the
  [same-matrices/two-shells result](../square_sum_berggren_absorber_bijection_20261019.md),
  the norm equation and signs certify which inverse branch stays
  positive and why the height drops. The tree itself does not
  give completeness. Here the analogous certificate is a literal
  lower donor with the required seam inequalities. A branch that
  lowers the count but leaves the needed interface has no more
  force than a Berggren inverse with a negative leg. Turning the
  scalar induction into a statement about *available marked
  donors* is the counterpart of preserving the norm shell.
* **Zeckendorf and Collatz.** The signed nonadjacent expansion in
  [TRIANGLE_NAF](TRIANGLE_NAF.md) is a literal sparse-digit
  uniqueness theorem, not a numerical resemblance: two endpoint
  readings locate every mixed triangle. That rigidity obstructs
  widening. Conversely, [earlier binary-address work](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md)
  shows that even a completely prescribed low-bit cylinder can
  have the full count residue image. A Collatz-style address
  does not determine the needed exact parent merely from the
  target residue. The carry test is the missing positive analogue.
* **Lonely runner and chromatic gluing.** The older
  [two-carry endpoint result](../../../01-canon/theorems/THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go.md)
  keeps compatibility of two actual endpoint allocations rather
  than correlating two scalar totals. The
  [chromatic/selector interface work](../chanpak_crossing_scales_20261115/CHROMATIC_INTERFACE.md)
  likewise keeps the full compatibility at a seam. Here there is
  a useful asymmetry: ask for both sides only for interior blocks;
  loosen each outer boundary to the data the actual context reads.
  This change, not an increase in the size of one fixed stock of
  gadgets, produced the five-count experiment.

## A concrete next proof attempt

Work in the first-missing prime/prime-times-seven lanes. For an
absolute proposed target, pick a scale and a candidate ordinal
layout, and compute its baseline **before** accepting its
congruences. Choose the one-sided end blocks only for the forced
residues. The complementary interior block must supply *the
specific quotient carry* at that residue, or supply a certified
interval covering it. Because the target fixes only one residue
per block, it is unnecessarily strong to demand complete palettes
at every interface. Search locally for actual quotients first;
use the general worst-endpoint formula only on families with
enough redundant representatives to enlarge a whole window.

To turn this into infinite descent, prove that for every large
target in either lane some such layout has its remaining donor
target strictly smaller, **and** the induction supplies the marked
interface at which it must be realized. If a step stays at the same
count, charge it against a separate bounded nonnegative interface
rank rather than assume it eventually ends. A complete *scalar*
induction hypothesis does not say anything about interface supply.
Alternatively a fixed portable interval amplifier plus a family
of baselines with successive gaps at most the widened interval
would bypass individual donor descent. Neither supply theorem
is proved here. The contribution is an exact test, a quantitative
amplifier, and physical diagnostics that specify what each route
has to achieve.

## Reproduction

The companion replays use literal tournament paths for the small
certificates and exact integer calculations for the carry tests;
no infinite range conclusion is inferred from a finite sample.
