# Validity gate and reproduction scope

**Same-agent proof review plus alternate exact computations. No independent
referee, Lean proof, or external priority claim.** The packet's proved
statements are J1-J6 in [PORTS.md](PORTS.md) and L1-L6 in [LRC.md](LRC.md);
each local ID includes its file to avoid legacy-ID ambiguity.

## Mathematical gate

- Inputs to J1 are nonempty tournaments and ordered old-port tuples, with
  repetitions allowed. Decoding uses degrees, arc directions and a terminal
  triangle contraction; it does not choose a vertex equivariantly inside C3.
- J1 is a pointed-class decoder. Its possible C3 automorphism kernel is
  stated and checked. The initial stronger expectation of a universally
  automorphism-exact register would be false already at a one-vertex core.
- J2 assumes arbitrary position input and ordinary class observations.
  Port indices are operation inputs, not free hidden information. The state
  domain is literal Hamilton presentations. A fixed successor menu is weaker.
- J3 uses actual C3 substitutions, possibly within previous insertions.
  Strong connectivity makes the tags strong modules; the strict inequalities
  p>n and q>p+n-2 guarantee their unique module sizes. The last order is
  4n+(n mod2), including the n=2 edge case. The two internal C3 group factors
  are retained. The context-size bound is not claimed optimal.
- J4 distinguishes minimum positive attachment cost, minimum strong-atom
  order and inherited minimum extra port arity. There is no universal
  ordering of all possible Boolean selectors.
- J5 counts completions with old/new roles distinguished and r named new
  vertices. It does not count ordinary output classes after forgetting
  the marked embedded core. Rigid cores make the formula (1+z)^(nr);
  rigidity does not make the core relations equal.
- L1 includes isolated points and delta=1/2; endpoints of zero-length
  components contribute twice, through their two outward sides. Rates
  use minimum reciprocals, not maximum reciprocals.
- L2 retains rational threshold and velocity height. The event denominators
  are at most2M; the chosen epsilon stays below the next possible lower
  event. No positive fixed-resolution gap independent of height follows.
- L3 changes the threshold before evaluating the Gibbs observable. Its
  iff criterion covers weak witnesses at the original threshold; it does
  not contradict the old strict-branch warning at an unchanged threshold.
  The full joint depth distribution is required.
- L4 reconstructs the multiset of active owner sets, not their component
  pairing or phase. Ordinary Boolean Mobius inversion is used explicitly.
- L5 compares local multivariate functions after identifying the first
  speed-coordinate (speed1 on one side, speed2 on the other). It supplies
  a genuinely new speed7 on both sides and keeps the threshold fixed.
- L6 closes additions at a fixed threshold on the isolated-witness
  stratum. Deletions and threshold changes can recreate bulk and are
  outside its scope. Reflection swaps the two owner maxima.

These are elementary proofs and bounded controls, not a claimed solution
of LRC(14), a universal compact tournament classifier, or a new asymptotic
enumeration algorithm. The complete order-ten bank is inherited and is
not regenerated in this packet.

## Failed probes retained

1. **Bare dent decoder.** Denting the source of TT3 and denting a vertex
   of C3 gives the same strong tournament of order4. A temporary sink
   repairs the ambiguity, with a possible C3 kernel explicitly retained.
2. **First path-count draft.** Counting paths ending at the new dent as
   H(T) forgot that its predecessor cannot be the selected p. The
   one-vertex core fails first. Subtract h_p^end(T); J4 proves the repaired
   identity. This was an exploratory formula, never a promoted theorem.
3. **A scalar positive germ as a continuation state.** Equal (L,D) does
   not preserve a later runner test. The stronger hostile L5 equates the
   entire local multivariate volume germ and still fails after speed7.
   The repaired tight carrier retains phases and the relevant owner data.
4. **Spectrum as relation memory.** Both rigid seven-cores have identical
   marked attachment spectra at every arity. Their selected extensions
   differ. The strongest survivor is a symmetry/action statistic, with
   actual readout operators still needed.

No demonstrated canon error was introduced or left live. Earlier singleton
minimality text is now explicitly historical and routes to the subsequent
minimal-blue theorem and this packet's different enriched menus.

## Exact computational universes

All checks use integers and fractions; neither floating point nor supplied
sequence values decide any new result.

| Checker | Universe and independent/alternate path | Checks |
|---|---|---:|
| check_ports.py | All Hamilton presentations through7 for one port, through6 for ordered pairs including repeats; deterministic cyclic relabeling; actual group search through4; dent identity versus endpoint path DP through5 | 1374379 |
| check_strong_tags.py | All Hamilton presentations and distinct ordered ports through4; four selected order5 presentations and two pairs each; actual module search; automorphism factor through3 | 958 |
| check_lrc_germs.py | Primitive subsets of1..10 of sizes1..4; thresholds drawn from1/(n+1),1/3,1/4,1/2; endpoint intervals versus event-depth distribution and Mobius inversion | 31533 |

The LRC receipt separates1243 base cases,2994 directional tests,
9240 tight-stratum additions with speeds1..15, and429 Gibbs cases.
Large controls include AP13 and the height339 row. Dilation controls,
delta=1/2, empty/bulk/isolated sets, and full signed local perturbations of
the phase hostile are explicit.

Both the original and native-complement seven-pair are retained.
For the81/27 control, the inherited constructor recomputes the shared67
cyclic triples and462915 Hamilton paths; this packet recomputes both
automorphism actions, the176/512 attachment-orbit counts, and selected
ordered-register decodings. It uses the old exact graph-isomorphism
routine; it is an alternate check, not a separately implemented GI engine.

## Reproduction and source pins

From the repository root:

    python -B 05-knowledge/results/joint_port_germs_20261028/reproduce.py

The runner first checks LF-normalized SHA-256 pins, then copies only the
required checkers and their inherited oracle into a temporary directory
with the same relative layout. It runs each checker normally and under
python -O, compares every generated receipt with the saved one, and prints
a compact result. Checks use explicit exceptions, so optimization cannot
silently erase them. The original workspace receipts are not rewritten.

The manifest pins this packet, the imported oracle, and relevant theorem
and bridge sources. These are snapshot dependencies, not promises that
future concurrent edits will leave the source hashes unchanged. A changed
dependency requires review before resealing. The runner itself is pinned.

The signed/even-graph and Borwein literature was checked against the primary
sources linked in [BRIDGES.md](BRIDGES.md), after reading CORE-PAPERS.
The corrected THM-4835 and THM-4837 lineages take precedence over older drafts.
Concurrent Session90 and Session91 work is credited as inherited or parallel
work, not silently claimed as this session's discovery.

The later [incoming integration](INCOMING.md) records Session92's separately
implemented49-check audit of four attachment counts and the uniform phase
example, replayed normally and under python -O. Its narrower scope does not
upgrade the general proofs here to independently refereed theorems.
