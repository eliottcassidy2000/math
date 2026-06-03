---
id: HYP-2089
status: SYNTHESIS — strong-lens reformulation + verified generic-strong split + tight=regular-rotational identification
source: opus-2026-06-02-S566
related:
  - HYP-2086
  - HYP-2080
  - THM-381
---

# HYP-2089: the LRC strong lens — observer escapes to a source, encircled by a (regular) strong block

- **Reformulation (THM-381):** LRC ⟺ the observer can always be lifted to the UNIQUE TOP SCC (sole source) of the (observer+runners) condensation at some time. As t varies the observer's SCC-membership oscillates (swallowed into a strong block whenever a runner is within 1/n); LRC = it escapes to source-hood; a counterexample = permanently entangled.
- **Verified (`lrc_strong_lens_s566.py`, 152 n=14 configs):** at a lonely time #SCC∈{1,13} only — STRONG (#SCC=1) 139 (91%, generic), transitive (#SCC=13) 13 (= exactly the semicircle/clustered ones), intermediate 0. Loneliness is GENERICALLY via a single strong block encircling the observer.
- **Moon's theorem:** strong ⟺ Hamiltonian cycle ⟺ runners cyclically SURROUND the observer; lonely-via-strong = observer in a gap of a Hamiltonian encirclement.
- **Tight = REGULAR strong block:** AP/V* are lonely via the regular rotational (doubly-regular, out-deg (m-1)/2) tournament = the roots-of-unity encirclement = the dual-Burnside FIX-side self-converse object (S565). Generic strong configs use IRREGULAR strong blocks. So the worry-set = the maximally-regular encirclement.
- **Unifies:** measure (S564) / orbit-resonance (S563) / Burnside fix-orbit (S565) / strong all point to the same object: the regular rotational encirclement.
- **PROOF SHAPE:** show the observer always reaches sole-source = the runners can always encircle/semicircle-cluster it leaving a ≥2/n gap; the only doubtful sets are forced into the regular rotational encirclement (tight family); the core = can a regular rotational encirclement always leave the observer a gap (= the danger arcs never completely cover, in strong-tournament language).
- HONEST: THM-381 & #SCC∈{1,m} known; NEW = the verified generic-strong/transitive split + tight=regular-rotational=fix-side identification + strong-language proof target; not a proof.

**See:** `07-reflections/lrc-through-the-strong-lens-observer-source-and-regular-encirclement-s566.md`, `04-computation/lrc_strong_lens_s566.py` (+.out); THM-381, S525/HYP-2000 (#SCC∈{1,m}), HYP-2086 (dual Burnside), HYP-2080 (resonance/orbit).
