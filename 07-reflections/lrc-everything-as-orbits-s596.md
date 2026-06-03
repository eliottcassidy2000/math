---
source: claudebox-2026-06-03-S596
status: REFLECTION — everything as orbits; a taxonomy of rigidity types as orbit types, unified at
  the 2-adic seam where the apex wears every hat.
tags: [rigidity, orbits, orbit-stabilizer, doubling-map, reflection, apex, 2-adic, equidistribution,
  spectral-rigidity, Burnside, LRC, perspectives]
---

# Everything as orbits: the rigidity zoo is one orbit-stabilizer ledger

**Prompt (human):** consider types of rigidity besides global and local; see everything as orbits;
poke around; design hypotheses freely.

## The instruction is the method

"See everything as orbits" is not a metaphor here — it is the literal organizing principle, because
every structure in this project carries a group acting on it, and every quantity we have called
"rigid" turns out to be a statement about an orbit. Once you adopt the lens, the rigidity zoo —
pinning, automorphism, cascade, spectral, dynamical — stops being a list of analogies and becomes a
single ledger: `|orbit| · |stabilizer| = |G|`. Local rigidity is the small-orbit (fixed-point) end;
global rigidity is the large-orbit end; and everything else is a different group `G`.

## The six rigidities, as orbit types

- **Pinning (local).** A degree of freedom forced to one value (THM-398's owner pin, the tight
  witness). As an orbit: a **fixed point**, a size-1 orbit. The stabilizer is everything.
- **Cascade (global).** A constraint that propagates to isomorphic copies (Lemma C; the `±`-pair).
  As an orbit: a **free orbit**, where the constraint must hold on every point at once.
- **Automorphism / stabilizer.** Trivial-`Aut` = generic (HYP-2130). As an orbit: a **large orbit /
  small stabilizer**. Its dual, **large stabilizer**, is the symmetric/tight config. These two are
  literally the two factors of one product — *conservation of rigidity is orbit-stabilizer*.
- **Dynamical / orbit-closure.** The runners' geodesic on the torus: dense (Q-independent,
  equidistributed, circuit-free, Lemma A) or confined to a subtorus (relations, Lemma B). As an
  orbit: the **closure** of the flow orbit; its dimension is `k − rank(Λ)`. Rigidity here = how much
  the orbit fails to fill space.
- **Spectral / determination.** Score-determined vs cycle-determined `H`. As an orbit: whether the
  **fiber** of an invariant map is a single orbit (the invariant "rigidly" names the class) or splits.
- **Phase / doubling.** The binary phase dynamics `x↦2x mod n` (HYP-2117). As an orbit: the
  cyclotomic-coset decomposition of `⟨2⟩` on `Z/n` — one maximal-mixing cycle (2 primitive root) or
  a 2-adic collapse (even `n`).

Six rooms, one architecture: pick the group, read the orbit.

## The apex wears every hat (the 2-adic seam)

The most striking thing the orbit lens reveals is that the apex — which I have met as the
reflection's fixed point (S593), the field's zero-divisor (HYP-2063), the doubling map's collapse
(HYP-2117), and the local pin (THM-398) — is *the same point* in all four, and it appears under
*exactly one condition*: `n` even. The verification is clean and total:

```
n even  ⟺  ⟨−1⟩ gains the second fixed point n/2 (the apex)
        ⟺  ⟨2⟩ degenerates (2-adic collapse)
        ⟺  the hard LRC frontier (2·prime: n=10,14,22; deeper: 2^a)
```

This is why the single corrector dies exactly at the apex and nowhere else: the apex is the unique
point that is simultaneously fixed by the reflection, collapsed by the doubling, and a zero-divisor
in the ring. Every method keyed to one point trips on it. And the cure's *shape* is forced by the
orbit: a fixed point (the apex alone) is handled by a one-point tool; the `±`-orbit around it needs
a pair tool (the pair-sum sieve); a richer symmetry needs higher arity. **The arity of the cure
equals the size of the orbit it must clear.** That is the whole single→multi corrector story,
re-derived from orbit-stabilizer.

## Where the frontier lives, in one classification

The orbit lens hands us a clean map of difficulty across `n`:

- **odd, 2 a primitive root** (5, 11, 13, 19): one maximal-mixing orbit — the most ergodic phase
  dynamics, the fractal-regular easy case.
- **odd, `⟨2⟩` splits** (7, 17, 23): several cyclotomic cosets — more orbits to resolve.
- **`2·prime`** (10, 14, 22): the apex plus a single 2-collapse — the named hard frontier.
- **`2^a · m`**: a deeper 2-adic seam, the deep water.

The Lonely Runner is hard exactly where the multiplicative phase action is *least ergodic* — where
its orbits collapse rather than mix. Hardness is non-ergodicity; ease is mixing.

## Designing freely: the conjectures the lens suggests

Three feel most alive. **(C2)** a difficulty index `R(n) = f(v₂(n), cyclotomic-coset profile)` —
the prediction that LRC effort is read off the orbit type of `⟨2⟩ mod n`, with `2^a` and `2·prime`
the extremes. **(C3)** minimal sieve arity = the largest witness-symmetry orbit — so richer-symmetry
`n` (like `18 = 2·3²`, where the unit group is bigger) demand strictly higher-arity sieves than the
pair-sum that suffices for `2·prime`. **(C1)** conservation of rigidity: difficulty `∝ |stabilizer
of the witness|`, the orbit-stabilizer dual of genericity. Each is testable; each says the same
thing in a different group's language.

## The transcending pattern

Rigidity was never one property; it was the shadow that any group action casts, read at a point
(fixed = local), over an orbit (cascade = global), through the stabilizer (symmetry), in the closure
(dynamics), or in the fibers (spectrum). The Lonely Runner is the place where one ring's
multiplicative action — doubling on `Z/n` — has its ergodicity broken by an even modulus, and the
broken orbit (the apex) is the obstruction. To prove the conjecture is to find, for each `n`, a tool
whose arity dominates the largest orbit the broken action can produce. Everything is orbits; the
difficulty is exactly the orbits that refuse to mix.

**Artifacts:** `04-computation/lrc_rigidity_as_orbit_type_s596.py` (+`.out`); new **HYP-2140**.
Builds on HYP-2135 (local/global), HYP-2130 (perspectives), HYP-2117 (doubling/IFS), HYP-2063
(apex), HYP-2120 (orbit-closure), HYP-2085 (time-Burnside).
