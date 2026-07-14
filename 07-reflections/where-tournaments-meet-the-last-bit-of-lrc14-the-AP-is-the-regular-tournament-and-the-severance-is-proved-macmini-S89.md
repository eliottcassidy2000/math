# Where the tournament meets the last bit of LRC(14): the AP is the regular tournament, the meeting is four-times proved, and the severance is proved too

*mac-mini-2026-07-13-S89. Owner asked: synthesize the recent closure work, and see how tournaments
relate to the last bit of LRC(14). This reflection does both. The honest headline: the two domains
of this project — parity in tournaments, and the lonely runner — **meet** at the last bit on one
Z₂-folded staircase, and that meeting is architecturally exact and proved in four places; but every
attempt to carry a **proof step** across the bridge is either refuted (the variance transfer is
miscalibrated) or an open program. The tournament's gift to the last bit is the **concept** — the
AP-rigidity IS the regular-tournament rigidity — not the metric bound.*

---

## 1. The last bit, precisely (synthesis of the S128 / S270–290 cascade)

As of today the covering case of LRC(14) has been reduced to a single metric statement and is being
closed stratum by stratum by a fast-moving disc_v / peeling machinery:

- **klein THM-731** (the certificate): peel any runner `v`; `L = (6/7)|G'_{~v}| − ε_v`,
  `|ε_v|² ≤ (6/49)·disc_v`, so `L ≥ (6/7)|G'_{~v}| − √((6/49)disc_v)` — rigorous, the metric door.
- **klein/kps THM-732** (the arithmetic): `disc_v` is a **generalized Dedekind sum** over the good-set
  edges (`disc_v = (1/2v²)Σ_{e,e'} σ_eσ_{e'} B₂({v(e−e')})`), so `disc_v ≤ r²/(3v²)` and per-family
  certificates become **exact ℚ proofs** of `L>0`.
- **kps THM-733/734** (the near-AP region): every 13-set with **≥11 speeds in {1,…,14}** is
  1/14-lonely — all 364 bodies, 245 994 exact-ℚ pairs; tight census = exactly `{AP, GW-doubling}`.
- **opus THM-735** (the compressed tower) and **kps THM-735** (bonferroni simultaneous multi-peel):
  the coherent-pack / clustered stratum, closed by sharing one clock across the pack.
- **THM-736** (mac-mini, this arc): the closed-form deep-well instance, `disc₁₈₂` via the Farey/
  three-gap arc count `N=φ(13)=12`, `|G'|=H₁₂/91`, reduces to `H₁₂>√2`.

**What remains — the last bit** (three routes name the same wall, klein-S290 / opus-S271 / kps-S128):
the **non-isolated, gcd-incoherent clustered residual** — covering families with no coherent
sub-pack (klein's `{1,90..101}`). Equivalently, klein-S290's exact form for `S={1}∪C`:
`L = |G(C)|·(1 − conc/7)`, `conc := 14|G(C)∩[0,1/14)| / |G(C)|`, so **`L>0 ⟺ conc<7`**, with the
**AP `{1..13}` the unique `conc=7` tight case**. Equivalently again: the covering-min
`M ≥ 14/183 > 1/14` (a uniform gap below the AP).

**My census this session (HYP-6545).** The `conc` forbidden band is real but subtle: covering
families top out not at a bounded drop-`x` set but at **far near-AP families** — `{1..11,13,84}`
reaches `conc = 6.436`, closer to the AP wall 7 than any `{1..14}\{x}` (drop-6 = 6.177). And
`conc<7` **uniformly ⟺ the covering-min** (`M>1/14` with a gap). So klein-S290 is honest that this
is a *restatement, not a reduction*: the band edge is approached, and pinning the uniform gap **is**
the last bit. This is the "multiple non-coinciding axes" phenomenon (HYP-4061) once more — the
`conc`-closest family (near-AP residue) is not the `M`-closest (deep well).

---

## 2. The four proved meeting points

The two domains genuinely touch — these are proved identities, not analogies:

1. **One partition function.** The danger count `X(t)=#{i:‖v_i t‖<1/14}` plays the tournament
   adjacency/score role; `Q(w)=E[w^X]` gives the tournament Hamiltonian count at `Q(2)` and the LRC
   loneliness at `Q(0)=L`. The AP extremizes both (S82) — the regular circulant is the joint
   extremal.
2. **One parity in the 2nd moment.** Tournament: `Var(H)=(n!/4^{n-1})(W(n)−n!)`, `CV(H)²~2/n`,
   with **odd shared runs cancelling by orientation parity, even surviving** (THM-589, PROVED). LRC:
   the covering floor's exact **2-adic descent — even survives, odd peels** (THM-580). THM-589's own
   text calls itself "the `S_n`-side mirror of the LRC 2-adic descent."
3. **The AP is the rotational (Paley-adjacent) tournament.** The oriented danger circulant
   `C_m(1..j)` is a rotational tournament `R_m` with **odd Hamiltonian-path count** (Rédei) and
   `H(R_m)=I(Ω,2)=` the OCF; the covering-min is the **fractional-chromatic / independent-set**
   problem on that circulant (HYP-3733). The apex-7 odd core shares the Paley/QR Cayley spectral
   gap. *Tournaments and runners meet at the circulant `Z₇`.*
4. **The last-bit object is the tournament's odd-graph shadow.** `disc_v` is a Dedekind sum
   (THM-732); the covering-min **margin `= −12 s(n,Φ₆)/n²`** is the `ι`-odd **cusp form `f₁₄`**
   (genus `X₀(14)=1`) (HYP-3768); and under the merged-metagraph Z₂-fold the tournament **odd
   (blue) graph** pairs with exactly this `ι`-odd sawtooth (HYP-3813, Bridge 3). The margin lives on
   the odd side of both folds.

---

## 3. The proved severance (why the naive bridge is dead)

The flagship transfer — THM-589's clean tournament variance → the LRC gatekeeper `CV(N_R)²`
(THM-579) — **cannot be a proof step**, and this is proved, not merely stalled:

- `CV(N_R)²` is **unbounded** over 14-free `R` (sup ≥ 8.74, grows as `m_R→0`; HYP-3554), and worse,
  **miscalibrated**: the covering family with the *lowest* actual floor (`R'=0.315`, drop-7) *passes*
  the `CV²<6` gatekeeper, while gatekeeper-*failing* families have higher floors (opus-S58). "A
  uniform floor proof cannot repair the CV bound — it must **replace** the proxy."
- **The mechanism has no tournament analog.** `CV(H)` is bounded because `S_n` acts transitively (no
  vanishing fiber); `CV(N_R)` blows up because dense `R` drives `m_R→0`. The tournament model is
  *too clean* — its symmetry hides exactly the pathology (apex-7: the same 7 that is the band's
  Fourier-zero, the `Var(N_R)` sheet-folder, and the floor-minimizer) that the LRC last bit turns
  on.

So every attempted **transfer of a proof step** — THM-589 → `CV(N_R)²`, the Γ₀(N) congruence lift
(HYP-3553, Han–Lee), `Var(W)≤c·R₂` (LEM-007, the "real open mile," core inequality open at the
barely-covers wall), OCF = covering-min — carries an honest *SYNTHESIS / OPEN / not-a-proof* status.
The bridges are coherent; the load-bearing analytic transfer is unbuilt.

---

## 4. What the tournament actually gives the last bit

Not the bound — the **structure**. The last bit is an **AP-rigidity**, and the tournament explains
*why the AP is isolated*: it is the **regular circulant**, the free `(Z/13)*` orbit, the `χ=2` orbit
(s581o: "χ separates regular tournaments beyond vertex-transitivity; the LRC extremal is the χ=2
orbit"). Regularity is the ground rung of a **quantized ladder with forbidden bands** (HYP-4306:
"regularity ⟹ divisibility ⟹ quantization," proved ≥9 times in the repo). Every LRC magic number is
a rung: `1/13` (AP, rung 0), `2/25` (mediant, rung 1), the covering-min `14/183`. My `conc`
forbidden band `(6.436, 7)` is the **newest surface of this same quantization** — covering is pushed
off the AP ground-rung `conc=7` by a full band, exactly as it is pushed off `M=1/13` to `M≤14/183`.

The tournament tells you the AP *must* stand alone (a regular free-orbit cannot be approximated by a
covering set without paying a full rung); it does not tell you the *width* of the gap. The width is
metric — and it is closed on the LRC side, family by family, by the Dedekind/disc_v certificate.

---

## 5. The one unbuilt transfer that would matter

If any tournament fact could become a proof step for the last bit, it is meeting point **4**. The
margin is the `ι`-odd cusp form `f₁₄`, and the tournament **odd-graph** is the same Z₂-fold. The
Joukowski bridge supplies the analytic scaffold: the tournament `I(Ω,x)` is **real-rooted**
(Chudnovsky–Seymour, claw-free), the LRC miss-count PGF is **circle-rooted** (Lee–Yang), and
`w=z+R²/z` carries one locus to the other, with the apex-7 ideal = the **7th cyclotomic**
`1+z+⋯+z⁶` (Joukowski image = the de Moivre angles). **The unbuilt theorem:** that the tournament
odd-graph's regularity-rigidity (regular = extremal, the χ=2 orbit) *controls* the cusp-form/Dedekind
obstruction — i.e. that "not the regular tournament" forces "margin bounded below." That would turn
the tournament's structural gift into the metric bound. It is the single bridge worth trying to
build, and it is not built.

---

## Coda

The recent work closed everything except the incoherent-clustered residual, and localized that to
one metric fact: covering ⟹ `conc<7` uniformly (⟺ the covering-min). On this fact the two halves of
the project stand face to face. They meet — four proved identities, one Z₂-folded staircase,
`14=2·7`, apex-7 = Paley = Klein quartic — and they are severed, provably, at the variance transfer,
because the tournament's transitive symmetry cannot see the runner's `m_R→0` pathology. The last bit
will be finished on the runner side, in Dedekind sums and disc_v certificates; the tournament's role
is to say, with the authority of the regular-orbit rigidity, *why there had to be a last bit, and why
the AP is the thing it defends.* If the odd-graph ↔ cusp-form fold (§5) can be made quantitative,
the tournament finishes it too. Until then: coherent meeting, proved severance, one bound still owed.

---

*Cross-links: the closure cascade — THM-731/732/733/734/735/736; the last-bit residual — klein-S290
(`conc`), opus-S271 (dilation-shadow), kps-S128 (census); the four meeting points — THM-589+THM-580
(parity), HYP-3733 (rotational tournament = covering-min), THM-732+HYP-3768+HYP-3813 (Dedekind =
odd-graph shadow), S82 (joint `Q(w)` extremal); the severance — HYP-3554 + opus-S58 (miscalibration),
LEM-007 (`Var(W)≤c·R₂` open); the quantization frame — HYP-4306; the analytic scaffold — the
Joukowski / apex-7 reflection, s581o (χ=2 orbit). My census: HYP-6545,
`04-computation/lrc14_conc_forbidden_band_macmini_S89.py`, `lrc14_conc_band_uniform_macmini_S89.py`.*
