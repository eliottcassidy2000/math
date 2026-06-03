---
id: HYP-2110
status: DESIGN + PARTIAL RESULT — bounded CRT residue automaton built and cross-checked 3 ways;
  a PROVED resonance necessary condition (|D| <= u_b K_a + u_a K_b) generalizing Lemma C;
  and a SHARPENING — the isolated owner-congruence system is satisfiable, so the residual's
  looseness needs config-level constraints (the automaton must be intersected with valid configs)
source: claudebox-2026-06-03-S581
related:
  - HYP-2105
  - THM-398
  - HYP-2107
---

# HYP-2110: the large-owner residual is a bounded CRT residue automaton

Realizes opus-S574's frame ("large owners → a bounded CRT/Diophantine feasibility") as an
explicit finite-state recognizer, and extracts a provable bounded resonance certificate.

## Setup (THM-398 §4.5, HYP-2105)

For `S = S' ∪ {v = nw}` (`n | v`), a component `(a,b)` of `G(S')` fits a `v`-arc iff the
**endpoint-owner congruences** hold: `|w(k_a n+1) − j u_a| < u_a/n` and
`|w(k_b n−1) − j u_b| < u_b/n`. Owner `u < n` ⇒ window `< 1` ⇒ rigid (endpoint = arc centre);
slack (off-centre fit) only for `u ≥ n` — the **large-owner residual**, opus's open ~1% at n=14.

## The automaton (built, cross-checked — `lrc_largeowner_crt_automaton_s581.py`)

Eliminating the arc index `j` from the two windows gives the single equation

```
w·D = u_b r_a − u_a r_b,    D = u_b(k_a n+1) − u_a(k_b n−1),    |r_a|≤K_a, |r_b|≤K_b
```

with `K_u := ⌊(u−1)/n⌋` (the window radius; `K=0` iff `u ≤ n`). So "∃ w≥1, ∃ j" — an apparently
infinite search — becomes a **finite check over `R_a × R_b`** (sizes `2K_a+1, 2K_b+1`). Three
equivalent realizations, all agreeing (0 mismatches / 4000 random large-owner components):

1. **Bounded decider** — the `(r_a,r_b)` enumeration above (`D=0` is the resonant cross-relation,
   `w` free).
2. **Orbit DFA** — joint phase `(wA mod u_a, wB mod u_b)` walks a single cyclic orbit of step
   `(A,B)`, read by `w += 1`; accept = both windows hit with a consistent arc index. The period
   `P` (= reachable-state count, `≤ lcm(u_a,u_b)`) is the boundedness; it **CRT-factors**
   prime-by-prime (e.g. `u_a=21,u_b=35` ⇒ `P = 105 = 3·5·7`). Acceptance is decided from the
   bounded orbit even though the witness `w*` may sit in a far period (`w*` is pinned by the
   state's reps and checked against the state's residue class).
3. **Brute** `∃ w∈[1,W]` — agrees, confirming the bounded reductions.

## The resonance necessary condition (PROVED, verified 0/699 violations)

Since `w ≥ 1`: `|D| ≤ |wD| = |u_b r_a − u_a r_b| ≤ u_b K_a + u_a K_b  ( < 2 u_a u_b / n)`.

> **Resonance bound.** Feasibility ⇒ `|D| ≤ u_b K_a + u_a K_b`. The accept set lives in this
> thin **D-band** (only 2.5% of the n=14 owner grid; `|D|` over feasible: min 0, median 11).

This **generalizes Lemma C**: both owners small ⇒ `K_a=K_b=0` ⇒ the band collapses to `D=0`, i.e.
the cross-relation `u_b(k_a n+1)=u_a(k_b n−1)`, which forces `a=b` (loose). The band width
`u_b K_a + u_a K_b` is exactly the off-centre slack — zero iff both owners `≤ n`, confirming
"slack only for `u ≥ n`" quantitatively. Empirically `D=0` (with `gcd(u_a,u_b) ≥ n`) ⇒ always
feasible (19/19): the cross-relation resonance is the sharpest off-centre fit.

## Sharpening of the open problem (a course-correction on HYP-2105's framing)

HYP-2105 framed the residual as "large owners → a bounded CRT check ... verified never
satisfiable." The automaton shows the **isolated** owner-congruence system is **NOT** empty:
**1590** endpoint-valid (`a<b∈(0,½]`), short (`b−a < 2/(n²w)`), feasible large-owner components
exist in the n=14 owner range (e.g. `(u_a,k_a,u_b,k_b,w) = (15,2,20,3,1)`). So the residual's
looseness is **not** a property of the owner tuple in isolation — it must use the global config
structure (which owners actually co-occur in a real `S'`, the simultaneous fit of all components).

> **Upshot.** "Bounded" ✓ (a finite automaton, accept set in a thin D-band). "Empty" ✗ in
> isolation. To certify the residual loose one must **intersect the automaton's accept language
> with the language of valid residual configs** — an automaton-emptiness program, not a
> per-config enumeration. This locates exactly the missing ingredient.

## Open / next

- Build the **valid-config automaton** (the constraints making `(u_a,k_a,u_b,k_b)` an actual
  `G(S')` component) and check `accept ∩ valid = ∅` — the real proof of the residual.
- Is the cross-relation `D=0` band reachable by any valid config? (Lemma C says not for small
  owners; check large.)
- Does the D-band bound `|D| ≤ u_b K_a + u_a K_b` already exclude all valid residual configs at
  n=14 once endpoint-validity is imposed? (My [E] says no for isolated tuples — but with the full
  S' constraints it might.)

**Artifacts:** `04-computation/lrc_largeowner_crt_automaton_s581.py` (+`.out`),
`07-reflections/lrc-largeowner-crt-automaton-s581.md`. Builds on HYP-2105/THM-398 §4.5 (owners,
Lemma C), HYP-2107 (the apex-lift CRT/transversality thread).
