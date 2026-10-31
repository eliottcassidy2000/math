# Phase fibres unlock exponential square-sum families

**Status: PROVED constructions and scoped fan theorems, with independent computation paths and same-agent proof review. No independent subagent audit or literature-priority claim. The universal tournament Hamiltonian-path spectrum remains OPEN.**

This is the successor to the frozen [packet d](../square_sum_spectrum_overnight_20261031d/README.md), inherited at `c63a347286`. It strengthens the construction rather than re-proving ordinary eventual square-sum Hamiltonicity. All tournament counts below are raw directed Hamiltonian-path counts H, not presentation counts H/Aut.

## Inheritance and live board

The closest proved mechanisms are packet d's readable 25-child tree and its complete local parity-preserving switch bank, the Session124 A/C paired and terminal interfaces, and THM-4748's one-apex fan formula. The canonical hostile gives duplicate label295 and missing293 when opposite offsets receive different source phases. The corrected near miss is to require opposite offsets to use the **same path**, when the necessary ownership invariant is only the **same label-position parity**. The least-used sidecar is this full phase function across different choices, not merely common parity within a single pair.

| Lane | Object and operation | Needed coordinate | Decisive result |
|---|---|---|---|
| Anchor | Square-sum signed copies in one phase fibre | Every label's positional parity | All49 signed offsets become independent; exponential cycle families |
| Anchor | Tournament completeness through changing words | Actual apex word and raw H | General midpoint fans, with an exact sparse-deficit obstruction |
| Niche | Protected native two-layer construction | Named ports, separated tails, integer atom clock | All2401 choices decode from the undirected cycle |
| Wildcard | Palindromic distance energy | Diameter and reciprocal binary mask | A strict carry law makes the energy injective at fixed diameter |

The shared research move has a precise scope: an interface can forget interior order only after its phase/clock requirements have been fixed. In square-sum composition that restriction enables a full product of choices. In the fan grammar a fixed diameter makes the weighted distance statistic lossless. Neither fact makes a scalar tournament count universally composable.

## 1. A uniform two-choice supply in one phase fibre

A weak-B pair of index n consists of paths on [1,n] and [1,n+1], both starting1,3, ending8 at even length and10 at odd length, with identical position parity for every common label. For a path P, its phase is the function assigning each label its zero-based position modulo2.

**Theorem 1.** At every n>=112 there are two weak-B pairs S(n,0), S(n,1) such that each corresponding member is different, while its phase is identical across the two choices. Choice0 is the incoming canonical A-only supply.

**Finite base.** For every112<=n<=2811, `phase_switches.json` gives an internal odd reversal or an exchange of two even-length blocks in each canonical member. These permutations fix the first two and last positions and move every label an even distance. The new boundary edges are explicitly square sums. The independently reconstructed permutations in `audit_phase_switches.py` check all5400 paths and7,894,800 vertices:5367 reversals and33 exchanges. The two selected moves make **each** changed member different. This is a complete finite base, not an extrapolation from sample indices.

**Induction.** For n>=2812 put m=floor((n-12)/25)>=112 and d=n-25m in12..36. Apply the fixed incoming A paired template to S(m,b), using the same choice b in every offset. Atom lengths, orientations, first two labels and endpoints are independent of b. By induction each source label also has the same phase, so its lifted value and output position parity are unchanged. Thus both resulting pairs are valid and have the same phase functions across choices. Every output member contains a full readable source tail, with a known head if split. That recovers one source member, hence b. The two members are therefore individually distinguishable. The source index strictly decreases, completing the induction.

Using A/B operator choices at large indices would **not** prove this theorem: those choices have different phases. The inherited295/293 hostile is replayed in `audit_signed_interfaces.json` and remains valid.

## 2. All signed offsets can branch independently

Let w be a finite word over A,C, where Q_A=25 and Q_C=49. Define

    R(w) = product of its radices,      R(empty)=1,
    T(w) = (225 R(w)-1)/2.

Write P_w for the incoming common prefix: P_empty=(1,3); A prepends(1,3) to the alternating lift25p_j+(-1)^j12, and C prepends(1,3,6,19) to49p_j+(-1)^j4, with j starting at1.

**Theorem 2.** For every n>=T(w), there are2^R(w) weak-B pairs with common prefix P_w. Projection to either member is injective. For each member index, **all choices in this family have one fixed phase function**.

**Proof.** The empty word is Theorem1. For w=s v with radix q, put h=(q-1)/2, m=floor((n-h)/q), d=n-qm. The threshold identity q T(v)+h=T(w) gives m>=T(v). Choose a child independently at each of the q signed offsets c=-h,...,h, and use its member selected by the fixed paired template.

For each source label x present in a selected member, its parity is the same for every choice. Therefore the lifted label qx+(-1)^j c is exactly the same label that the original common-source template assigns to that offset. All bands, low singleton labels and top partial-band ownership remain unchanged. This is stronger than relying only on a plus/minus pair: every signed offset now has the same reference phase. Connectors use only the prescribed endpoints and prefix. Each atom's length is fixed; within it the output phase is source phase plus a fixed shift. The two paired catalogs have matching shifts at every head, tail and singleton. Thus validity and the whole-family phase invariant propagate.

For injectivity, fix either member, n and w. Parse its known atoms and their lengths. At every signed offset, undo orientation, recover x=floor((v+h)/q) from each lifted label v, and restore head1 if its tail was split. Induction decodes the chosen child from that one member. Every signed offset is readable, so the number of distinct choices is

    (2^R(v))^q = 2^R(w).

All148 A/C paired states and148 terminal states are checked in `audit_signed_interfaces.py`, including symbolic phase shifts, concrete independently varied signed inputs, and decoders. The induction uses these complete interfaces; it does not infer an all-size result from the test cube. `phase_trees.py` separately tests32-choice cubes on five channels including both signs, mixed words AC/CA, and2401 decoded leaves in CC.

### Terminal cycles and exponential growth

For nonempty w, replacing the outer paired table by its terminal table leaves every signed tail readable. All cycles start1,3, so that edge fixes their orientation and rotation as undirected labelled cycles. Consequently, at every N>=T(w), there are at least2^R(w) distinct square-sum Hamiltonian cycles with **the same** prefix P_w.

Examples valid for every N above the indicated threshold:

| Word | Threshold | Lower bound |
|---|---:|---:|
| A |2812|2^25|
| C |5512|2^49|
| AA |70312|2^625|
| AC or CA |137812|2^1225|
| CC |270112|2^2401|

For all N>=2812, choose the largest power25^a<=(2N+1)/225. Then

    log_2(number of cycles) >= 25^a > (2N+1)/5625.

This is a uniform **exponential** lower bound. It replaces packet d's weaker stretched-exponential bound; it is not an asymptotic formula or upper bound.

There is also an asymptotic improvement using mixed radices:

    liminf_(N->infinity) log_2(number of cycles)/N >= 2/225.

To prove it, set Y=(2N+1)/225. The largest25^a49^b<=Y has ratio to Y tending1. Indeed log49/log25 is irrational (otherwise unique prime factorization would give7^u=5^v for positive integers). For any epsilon>0 a finite set of nonnegative multiples of this irrational rotation is an epsilon-net of the unit circle. For all sufficiently large log_25Y, choose one of these bounded b values with fractional remainder in[0,epsilon), and choose nonnegative a=floor(log_25Y-b log_25 49). Then1<=Y/(25^a49^b)<25^epsilon. Let epsilon tend0. Every such product is an actual schedule w, proving the bound. This elementary density argument claims no effective optimal cutoff.

## 3. The native path and its marked forest retain all2401 channels

The two inherited protected radix49 layers accept weak pairs and keep their named2--7 ports, orientations and cut positions. Apply Theorem2 independently at **each signed inner offset and signed outer offset**. The common phase across choices makes these2401 inputs compatible. Their prefixes and endpoints remain fixed, hence every old connector and prescribed edge remains valid.

For every word w and every

    N >= 2401 T(w)+1200 = (225*2401*R(w)-1)/2,

there are at least2^(2401 R(w)) distinct cycles retaining

* the native path(1,8,17,19,6,10,15,21,4,12,24);
* the seed482--3362--6242;
* the2303 mutually disjoint prefix rays
  (2401 p_j+(-1)^(j+1)(49u+v)), where p=P_w, u=-24,...,24 except0,2, and v=-24,...,24.

As before, cutting an edge can destroy an intact ray without destroying its complete ordered source tail. Each outer offset retains both fragments of one selected A member. Reassemble those, then read every inner tail, restore its head, and decode the child. This recovers all2401 histories. The complete196 inner states and392 outer states are audited in `audit_native_signed_readability.py`.

`native_phase_cycle.py` constructs full cycles at270112 and270113, verifies every edge and prescribed mark, **erases their starting point and direction**, and restores the integer atom clock from the canonical native edge1--8 before decoding all2401 choices. A separate control changes only the zero/zero bit and obtains a different canonical undirected cycle. Thus the count is not an artifact of rotations or reversal.

In particular every N>=270112 has at least2^2401 such protected cycles. The same mixed-radix argument gives asymptotic lower rate2/225 for this protected class as well: the fixed multiplier2401 cancels against its threshold. The forest is prescribed by this grammar. Arbitrary forests and an infinite coherent Hamiltonian ordering remain separate OPEN problems.

## 4. Tournament progress and boundary

[CENTERED-FANS.md](CENTERED-FANS.md) proves a larger family of fixed deficits below powers of two, an exact outer-pair recurrence, and a strict carry theorem for palindromic distance energies. It also proves that these fixed deficits have density zero: merely extending the Mersenne fan cannot supply all odd targets. The first excluded deficit in this grammar is9, not a globally forbidden Hamiltonian count.

The inherited raw-count result remains: all positive odds through99,999,999 except7,21 are realized;7,21 are globally absent by THM-1370 / `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. Completeness beyond the finite range is OPEN. No new unbounded acceptance assertion is made for the seven-parent bank.

## Reproduction and review

Run in this packet directory:

```sh
python3 phase_switches.py
python3 audit_phase_switches.py
python3 phase_trees.py
python3 audit_signed_interfaces.py
python3 audit_native_signed_readability.py
python3 native_phase_cycle.py
python3 centered_fans.py
python3 verify_packet.py
```

Each consumer has a matching normal/-O receipt. `VERIFICATION.json` records controls and scope; `MANIFEST.json` pins inputs and files. Producer and auditor use separate permutation reconstructions; the fan uses actual adjacency counts and small subset DP in addition to its formula. The proof review is by the same agent, with these independent computation paths, not a second-person audit. A pre-release variable-shadowing error in the native decoder was rejected by Python's type check, repaired, and fully rerun; no erroneous result was published.
