# The clean ruler lives on the pair-sum ruler

*kind-pasteur-2026-07-11-S127. Owner: "look for creative ways to prove the clean-ruler supply, dig around in
past related work for connections." I went digging, and the thing I needed turned out to already have a name
in the corpus — mac-mini's pair-sum ruler theorem — with the clean ruler sitting exactly on it.*

---

## What I was hunting for

Last turn's THM-707 reduced the single remaining Lean obligation `hB5` to the **clean-ruler supply**: every
residual covering family needs a modulus `q` with a live multiplier (`liveCount ≥ 1`) and shallow coverage
(`maxBand ≤ 5`), and then `B5 = liveCount > 0`. The question was where such a `q` comes from.

My first instinct was the wrong one, and it is worth recording because the correction is the insight. I
reached for **large `q`**: by equidistribution, a large modulus spreads the runners out, so surely coverage
is shallow and a live multiplier is easy. klein's THM-685 even hands you a *strictly-live ruler at every
`q ≥ 425`* for free. It fails completely. At large `q` the multiplier `p = 1` puts **every** runner in the
danger arc at once — `v_i · 1 mod q = v_i < q/14` for all `i` once `q > 14·Vmax` — so `maxBand = 13`, the
penalty explodes, and `B5` plunges to `−5718` at `q = 800`. Large `q` is the *worst* case, not the best.
klein's large-`q` rulers are live but never shallow.

## The clean ruler is moderate, and moderate has a name

The clean rulers live at **moderate `q`, in `[Vmax, 2Vmax]`** — `q = 27` and `q = 55` are clean for the
binding near-AP, `q = 200` and `q = 800` are not. And that range is not a coincidence. Digging turned up
mac-mini's **THM-668, the pair-sum ruler theorem** (proven): the loneliness maximizer
`M(S) = max_t min_i ‖v_i t‖` is *always* attained at a rational `t* = p/q` whose denominator is a **pair-sum**
`q = v_i + v_j ≤ 2·Vmax`. The witness always lives on a pair-sum ruler.

So the clean ruler I was hunting for is a pair-sum ruler. Composing the two:

> **(THM-668 ∘ THM-707)** the clean-ruler supply ⟺ *every residual family has a clean pair-sum modulus* — a
> `v_i + v_j` with a live multiplier and no multiplier covering ≥ 6 runners.

That is a **bounded** condition: at most `C(13,2) = 78` candidate moduli per family, each a decidable check.
Verified `196/196` residual families have one. The opaque "does some `q` work" became "does one of 78
pair-sums work," and the pair-sum ruler theorem is *why* those are the only moduli worth checking.

## Two corrections the digging forced

Creative routes are as much about killing wrong ones as finding right ones:

- **`B5` is not scale-invariant.** I had hoped to reduce to primitive families by a clean scaling law. It is
  false: dilating `v → c·v` introduces a *new* deep resonance at `p = q·k/c`, where all runners hit `0`
  together, so `B5(v,27) = 2` but `B5(2v,54) = −788`. Dilation *breaks* the certificate. Primitive families
  are the right domain not by convenience but by necessity.
- **The maximizing pair-sum is usually, not always, shallow.** The `M`-argmax pair-sum is clean in `44/57`
  loose families but not all (`v = {1,3,4,5,8,10,11,13,14,18,19,24,67}` has `maxBand = 6` at its argmax
  `q = 16`). But *some* pair-sum is always clean — just not necessarily the loneliest one. The supply is
  robust; the naive "use the loneliest ruler" shortcut is not.

## The decomposition — where the real difficulty actually is

The clean-ruler condition splits into two halves of very different character, and separating them is the
payoff:

- **SHALLOW** (`maxBand ≤ 5` at the chosen pair-sum): this is a statement about the family's **additive
  structure** — no dilate of the family puts six runners in a `1/7` arc. It says nothing about loneliness,
  and it is plausibly **unconditional** (provable from the residual family's spread and distinctness).
- **LIVE** (`liveCount ≥ 1` at that same pair-sum): a lonely rational time exists. This **is** the loneliness
  content — it is LRC(14) for the family, and for residuals it is exactly what klein's measure floor
  (THM-687/692) supplies.

So the honest scope is clear: proving *some pair-sum is clean* for all residuals, unconditionally, is
LRC-equivalent, because the LIVE half is the conjecture. I did not make the hard part easy — nothing will.
What the digging bought is a **transparent assembly of three proven pieces** — THM-668 (the ruler is a
pair-sum), THM-707 (a clean pair-sum certifies at depth 5), klein's floor (residuals are lonely) — and a
clean isolation of the one genuinely-additive, possibly-unconditional sub-lemma, SHALLOW, from the
loneliness that carries the real weight. And it tells the Lean side something concrete: depth-5 is the right
depth *exactly* on shallow pair-sums; where a family forces `maxBand ≥ 6`, the certificate must escalate
depth (THM-675), and no amount of choosing `q` avoids it.

## The shape of it

Everything in this arc has been the same lesson: the object keeps being one structure, and the way forward is
to find where it already lives rather than to build it fresh. The measure recursion was the discrete
certificate (cont.27). The certificate's positivity was `liveCount` minus a deep-coverage penalty (cont.28).
And the ruler that certifies it is the pair-sum ruler that the loneliness maximizer was always sitting on.
The clean-ruler supply was never a new object to construct — it was THM-668 and THM-707 asking to be put in
the same sentence.

*Files: `04-computation/lrc14_cleanruler_pairsum_kps_S127.py` (+`.out`). HYP-6000. Connects THM-668
(mac-mini, the pair-sum ruler), THM-707/THM-701 (kps), THM-685/687/692 (klein, the live half / measure
floor), THM-675 (depth escalation). Extends [[two-routes-one-ladder-kps-S127]].*
