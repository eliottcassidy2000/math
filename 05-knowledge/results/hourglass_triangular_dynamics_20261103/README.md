# Diagonal hourglasses and complete triangular dynamics

**PROVED scoped elementary theorems / FINITE-EXACT controls / CITED inputs.**
Same-agent proof audit with separately implemented computational checks;
not an independent researcher review or proof-assistant verification.

The user's clarification selects the **diagonal-layer hourglass** of
HYP-3053, not a partial tournament or a magic-square configuration.

1. [Complete triangular dynamics](DYNAMICS.md): for every integer k, classify
   every bounded integer orbit of b->b(b+1)/2-k. Equivalently classify all
   rational preperiodic points of every x^2-c/16 with c odd. The count is

       B(k)=4([k triangular]+[k-2 triangular])
            +2[k in {1,2,5,6}]+8[k=4].

   This yields an exact generating function and seven finite graph types.
   The known eight-point example at c=29 is an inherited special case.
2. [Hourglass and Boolean operations](HOURGLASS.md): construct an integral
   rectangle/hourglass cycle basis for all adjacent-layer graphs. All
   moduli work without extra torsion. The quadratic eight-state graph has
   an explicit closed two-line Boolean rule and a source/clock-preserving
   FRACTRAN encoding. Its 8->4->3 image ladder has different fibres from
   tournament isomorphism and converse merging.
3. [All-modulus prime wheels](PRIME_WHEELS.md): classify every maximal run
   in the +2 unit wheel at every even modulus via a second-difference
   formula. Joint-offset counts have exact CRT recurrences; future ports
   distinguish equal counts. Prime arrival and actual sum pairs retain
   their separate height and source clocks.

These results do not settle A000568 enumeration, the general tournament
Hamilton spectrum, Collatz, Goldbach, or LRC(14). The eight-state quotient
is observer-specific; it does not compress arbitrary named-port attachments.

## Inheritance and current portfolio

* Closest mechanisms: [THM-4139](../../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
  HYP-3053's layer potentials, and Session146's exact modulo-290 wheel.
* Hostiles: the actual 81/27 pair, the divergent lift9->13->35->299,
  and a word satisfying every rectangle but violating the hourglass.
* Corrected near miss: equal quotient cardinalities do not identify fibres;
  a finite residue orbit does not bound its integer lift. The all-order
  small-denominator argument supplies an actual magnitude restriction.
* Least-used relevant sidecar: the chosen trace/cycle in the September
  affine-reversal work; here the original signed numerator and elapsed
  time distinguish otherwise identical future folds.

| Lane / live concept | Operation | Preserved target / loss / decisive test |
|---|---|---|
| Anchor: diagonal-layer potentials | Integrate a line word | Exact labels up to one global constant; explicit hourglass syndrome |
| Niche: triangular integer map | Restrict to bounded candidates | All bounded orbits; magnitude cannot be replaced by residue period |
| Wildcard: prime unit wheel | Append named offsets | Joint CRT mask; equal counts split after one added offset |
| Tournament presentations | Forget path and labels | Class, but not the chosen quadratic update; fibre-size hostile |
| Machines and clocks | Encode one control prime | Exact bounded-state step, source and time; multi-control hostile |
| Mean attachment carriers | Fixed named attachments | Incoming minimum10 differs from pooled parity; 81/27 remains a control |

Connection contracts and failure boundaries are explicit in the linked
proofs. The relevant method card is
[controlled forgetting and unlabeled quotients](../../../00-navigation/META-PATTERNS.md#controlled-forgetting-and-unlabeled-quotients-require-a-sidecar).
It applies here without requiring a new meta-pattern card.

## Incoming work and priority

At commit13e4842a46, [the parallel hourglass/port packet](../hourglass_quadratic_ports_20261103/README.md)
already supplies native re-mark/complement/re-mark implementation of the
eight-state graph, a strong-completion recurrence for partial orientations,
and prime-stream saturation. Its hourglass is a different declared partial
orientation. The [parallel quadratic triangle note](../quadratic_2916_odd_pair_triangles_n4_20261102.md)
also already records the fibre obstruction. We credit those overlaps;
the additions here are the all-parameter classification and generating
function, integral diagonal-layer basis, explicit two-line polynomial
rule, and arbitrary-even-modulus run formula.

Incoming Session146 has also reduced the exact order-ten mean carrier to
345106 from the pooled-parity upper bound347626 (2520 invisible directions),
using the [upper certificate](../session146_carrier_primegame_20261103/carrier10_upper/PROOF.md)
and separate lower certificate; order nine has minimum65908. These are inherited scoped results,
not consequences of the quadratic construction; see the incoming packet's
source routes. The all-size exact minimum remains a separate obligation.

Reproduction and audit: [AUDIT.md](AUDIT.md). Continuation: [HANDOFF.md](HANDOFF.md).
