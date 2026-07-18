---
id: THM-989
title: Scale-twenty-two Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + REFEREED FINITE-EXACT — all 984 scalar survivors have at most one owner-local feasible projection; independent algebraic-CRT/Python-set and literal-CRT/C++ sorted-vector implementations agree, with optimized, unoptimized, and sanitized referee builds byte-identical
source: codex-2026-07-17-S66 scale-twenty-two continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_two_hamming_six_frontier_scout_codex_c22.py
  - 05-knowledge/results/lrc13_scale_twenty_two_hamming_six_frontier_scout_codex_c22.out
  - 04-computation/lrc13_scale_twenty_two_hamming_six_referee_codex_c22.cpp
  - 05-knowledge/results/lrc13_scale_twenty_two_hamming_six_referee_codex_c22.out
---

# THM-989 — scale twenty-two has at least five impossible owners

This namespace and the companion computation filename are reserved for the
next exact primitive proper AP-centred common-scale Hamming-six face.

For `c=22`, the effective orders are `1,2,11,22`, with twenty-two literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 3,249 divisor
words and 100,975,500 literal state words per support, hence

```text
924*100,975,500 = 93,301,362,000
```

unquotiented labelled state contexts.  A from-scratch algebraic-CRT probe
leaves 984 scalar contexts on 180 supports across eight multiplicity profiles.
Exact owner-local set-union reachability gives

```text
0 feasible owners: 792 contexts,
1 feasible owner : 192 contexts.
```

The 5,904 owner rows have maximum-union histogram

```text
16:864, 17:1584, 18:2784, 19:480, 22:192.
```

Thus every candidate has at least five empty owner projections.  This is
terminal for the common-scale Hamming-six face: a global unit word would
project to a feasible word at every owner.

## Frozen primary, independent referee, and carrier audit

The standard-library Python primary solves CRT representatives algebraically,
cross-checks every representative by literal search, and checks every mask
cardinality against an independent one-period formula.  It traverses all
3,002,076 support/order rows without quotienting, hashes every scalar survivor
and every sorted owner-local reachable-mask bank, and reproduces the following
stronger telemetry:

```text
scalar support orbits under multiplication: 15, all of size 12;
distinct capacity vectors: 484;
distinct owner max-union vectors: 127;
owner-reachable-bank SHA-256:
b881af1af73fe6dde434d92ca0598bac79828881ad32feff87d711fa448a0eef;
forward/reverse layer-bank SHA-256:
b08be64149acffb006205465bbcb5825b06d4e2ffb1ece561506f0e1746b8baa.
```

Normal and `python -O` runs reproduce the frozen output byte-for-byte.  Frozen
artifact SHA-256 values are

```text
primary source  94818bd64b0492156d7bf66b4fc96a679c8bbb112664ed1311cde7cd227b5914
primary output  27b640da323604976b9172a6ffd9fc32e4c8b18c9f89227a420c0af566e40be4
```

The independent C++ referee does not translate that implementation.  It finds
every CRT representative by literal search, derives hereditary admissibility
from the two prime-power coverage conditions and separately audits all six
leave-one-out lcms, and realizes each owner-local reachable bank as a sorted
vector of integer masks.  Its optimized, unoptimized, and ASan+UBSan builds
produce byte-identical output.  The hardened proof-bearing serializations all
agree with the primary exactly:

```text
CRT-base bank          fe217f797c08702c8f607d3a936321fb7ffb6c6e73770a0097cf71c597297793
mask table             54587d940a12b70601943dbe7505d4797363395a2579ea6f3e09583db5a01282
hereditary order bank  f7c0254d8ac9108d318f4a9a21d0d2e5b244be91087b22c162bb563956e9b474
weighted grammar bank  b4ec74d190864c2a050409126bacfd79fb0fac97ce2952628b17341b1718c4dd
scalar survivor bank   29067f69b228b9956239b27a43af9bc72e8c141acfb47587f536dc557cebb1de
capacity bank          5f10732eda4cd0dcf9fe2eb0166e4191673774d40fae03ec4993d143caa3528f
owner-reachable bank   b881af1af73fe6dde434d92ca0598bac79828881ad32feff87d711fa448a0eef
forward/reverse layers b08be64149acffb006205465bbcb5825b06d4e2ffb1ece561506f0e1746b8baa
```

The multiplicity ledger and 411 tournament score vectors also cross-match at
`5e1220b6...a80d` and `fdf170ef...aa0b`.  The referee independently checks
that reversing all six provider insertions leaves every final reachable bank
unchanged and hashes both intermediate layer sequences.

The frozen referee artifact hashes are

```text
referee source  fc50462397411f76ee03bfd707c8ee9986dbce1c69257fbfa62c7d77b50bb40d
referee output  f93a14158283fd0bbf91d62d7cd9d8b7eb89653e8ac517fcb0f16db6f59aec13
```

Owner obligations are the terminal-faithful vertices: their exact summary is
`(feasible,max-union,capacity,reachable-count,maximum-mask-count)`.
Lexicographic comparison with coordinate tie breaks gives a transitive
tournament in all 984 rows (score set `0,...,5`, zero directed triangles, six
SCCs, one Hamiltonian path).  That tournament is telemetry only; it forgets
the absolute 22-sheet threshold and the number of empty owners.  Provider,
divisor, residue, isolated-sheet, and wall-event vertices lose shared-unit
incidence earlier still.

The independent replay closes this exact scale-twenty-two face.  The theorem
does not cover scale 24 or higher (scale 23 is prime and already excluded by
THM-983), H5 ramification, non-AP/deep sheets, or global sporadic emptiness.
