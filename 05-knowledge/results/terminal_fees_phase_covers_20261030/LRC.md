# Full phase covers for an arbitrary attached speed

**PROVED:** an exact parametric decision rule for every positive integer attached
speed, with a fixed finite core and radius. **FINITE-EXACT:** 50,176 literal-comb
comparisons and native witnesses, plus large-input and hostile controls.
**OPEN:** universal LRC(14), a uniform inexpensive construction of these covers
for arbitrary cores, and a proof that the relevant covers always accept.

## Inheritance and the changed target

THM-2186 (`01-canon/theorems/THM-2186-exact-octagon-needle-profile-and-mobius-h-drift.md`)
already proves eventual rational laws on residue classes at isolated toric
maxima. THM-4448 (`01-canon/theorems/THM-4448-lrc14-general-shore-attachment-and-decoder-pair-cones.md`)
already provides a general protected-arc attachment mechanism. These are stronger
than a new claim of eventual coverage for one specially constructed family.

The present target is different: decide the **prescribed actual speed**, at a
prescribed radius, using the **whole safe set**, including its point components.
It complements the source-preserving recursive families in
[the preceding packet](../return_kernel_core_transplants_20261030/LRC.md).
It never replaces the new input by a chosen congruent speed or an asymptotic one.

The canonical near miss is repaired THM-724
(`01-canon/theorems/THM-724-single-killer-covering-min-rigidity.md`): the local top germ of (1,6)
cannot be continued down to the predicted radius 3/8 after adding7; the true
maximum is4/13. Also, independent successful arrivals of3 and12 into the core
(1,2) do not give simultaneous coverage at1/4. Its true joint maximum is3/13.
Our compiler checks both failures exactly.

## Exact attachment theorem

Let C be a finite set of positive integer speeds, 0<r<=1/2 rational, and write
its complete closed safe set in the **original time coordinate** as

    G_C(r) = {t in [0,1]: ||ct|| >= r for every c in C}
           = union_i [L_i,R_i].

Components may be points. All endpoints are rational. For a phase phi in[0,1)
define its forward distance to the next safe phase by

    g_r(phi) = 0                    if r <= phi <= 1-r,
               r-phi               if phi < r,
               1+r-phi             if phi > 1-r.

For every positive integer v,

    G_(C union {v})(r) is nonempty
      iff there is i with v(R_i-L_i) >= g_r({v L_i}).       (L1)

In that case an exact native witness is

    t = L_i + g_r({v L_i})/v.                              (L2)

**Proof.** Starting at L_i, the unwrapped new phase increases at speed v.
If already safe, L_i works. Otherwise its first arrival in any translate of
[r,1-r] is at exactly the distance g_r shown above. It lies before or at R_i
precisely when(L1) holds. If it lies later, there was no earlier safe point
in that entire component. Taking the union proves both directions. Closed
endpoints explain the non-strict inequality. This also proves the assertion
for a zero-width component: only a zero phase gap can pass. QED.

This is a decision theorem, not just a sufficient endpoint test. It does not
assume any particular component contains the maximizer of the extended core.

## Exact residue thresholds, without precomputing a huge residue table

Put Q=lcm_i denominator(L_i), with Q=1 for an empty safe set. For a residue a
modulo Q, the phases {vL_i} are fixed. Give component i the threshold

    0                                      if g_r({aL_i})=0,
    g_r({aL_i})/(R_i-L_i)                   if R_i>L_i,
    infinity                               otherwise.

Let T_C,r(a) be the minimum, infinity for an empty collection. Then exactly

    v congruent to a mod Q:
    G_(C union {v})(r) nonempty iff v >= T_C,r(a).           (L3)

This is an exact residue-by-residue conductor, including all small inputs;
there is no unspecified eventual threshold. The map is upward closed within
each residue, not necessarily across residues: at radius1/3 for core(1,2),
speed4 passes and speed6 fails. A point-only core can have entire forbidden
residue classes. If any component has positive width W, its threshold is at
most2r/W, recovering the older long-arc sufficient test.

The implementation evaluates only the residue of the actual input. Computing
every one of the Q entries could be wasteful; Q is not presumed small.

## Algorithm and complexity boundary

Once the k complete core intervals are available, (L1) needs O(k) rational
arithmetic operations and returns a rational witness. It does not enumerate
the O(v) teeth of the new speed. Bit complexity still depends on the bit lengths
of v and the endpoints; this is not a constant-time big-integer claim.
The preprocessing cost, number of core components, and denominator growth
remain explicit. For a second attachment, recompute the actual joint safe set:
the separate Boolean verdicts for two attachments cannot be composed.

The test file uses an independently organized literal interval intersection as
the reference and a separate pair-clock maximum calculation on the two hostile
examples. Universe: all196 primitive triples from1..12, radii1/5,1/4,1/3,1/2,
and every attached speed1..64. The large core26*(1,..,12), radius1/14, has416
components and phase period10,090,080. Exact witnesses are produced for339,
10^100+39 and10^1000+123. No floating point is used.

Reproduce:

    python3 05-knowledge/results/terminal_fees_phase_covers_20261030/phase_cover.py

The result is [phase_cover.json](phase_cover.json). A huge *specified* input is
now inexpensive after core compilation; infinitely many *arbitrary* cores are
still the missing quantifier. That distinction prevents an algorithmic speedup
from being misreported as universal LRC coverage.

## Next questions with decisive tests

1. Can one compress the complete interval cover while retaining(L1) for every
   future attached speed? Any proposed compression must distinguish the point
   obstruction above and the separate-versus-joint arrivals hostile.
2. Can a twelve-speed core carry a small set of weighted intervals whose
   residue thresholds accept every admissible thirteenth speed at1/14? This is
   a sufficient certificate search, not an assumed universal property.
3. Can a source-dependent bound on the rejected residues be proved without
   enumerating Q? Fourier averages alone do not certify the actual residue;
   its contribution must remain positively weighted.
