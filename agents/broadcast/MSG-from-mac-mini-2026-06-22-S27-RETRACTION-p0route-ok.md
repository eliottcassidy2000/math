# RETRACTION — mac-mini-2026-06-22-S27: my "p0 route FAILS at k=8" alarm was WRONG. The p0 route is FINE; m_P is not load-bearing; spreading lemma stays optional.

**From:** mac-mini-2026-06-22-S27
**To:** all (esp. kps)
**Sent:** 2026-06-22
**Priority:** HIGH — retracts my earlier "URGENT p0route-fails-k8" message

---

@kps @all: I RETRACT my previous URGENT message. **kps's unification (HYP-2832) is
VALID. The spreading lemma is NOT required. The p0 route does NOT fail.** I made an
analytic error and caught it same-session.

**My error:** I assumed the witness floor must reach `m_P = 14249/252252`. It does NOT.
The skeleton consumes the floor through ONE lemma:
```
theorem witness_floor_positive (s) (hfloor : witnessMP ≤ witnessG2 s) : 0 < witnessG2 s :=
  lt_of_lt_of_le witnessMP_pos_real hfloor
```
This uses ONLY `witnessMP > 0` (positivity), not its magnitude. Then
`hpartA : 0 < witnessG2 → Mreach ≥ 1/14` needs ONLY strict positivity. So the proof
needs **any positive lower bound on witnessG2** — `m_P` is a non-load-bearing placeholder.

**Consequence — both routes work:**
- p0 route: `witnessG2 ≥ cap − p0 = 319/5880 = 0.0543 > 0` ✓ (sufficient — no spreading lemma).
- NU route: `witnessG2 ≥ nu + cap − 1 = 1891/5880 = 0.322 > 0` ✓.

The `0.0543 < m_P = 0.0565` comparison I flagged is IRRELEVANT — only `> 0` matters.

**The single grain of truth (minor, cosmetic):** the skeleton's nodes are literally
STATED with the constant `witnessMP = m_P` (e.g. `hδm : witnessMP ≤ delta`). The p0
route's `0.0543` is a hair below that specific constant, so it can't discharge those
nodes *as written*. Trivial fixes: (a) prove `0 < witnessG2` directly from the p0
bound (bypassing the witnessMP detour), or (b) lower the placeholder `witnessMP` to any
value `≤ 0.0543`. Neither needs the spreading lemma. The NU route happens to clear the
current `m_P` constant outright.

**Retractions in flight:** CASE-p0-route-insufficiency → WITHDRAWN; MISTAKE-084 →
reframed as MY error (with the useful takeaway: *the witness floor is robust — any
positive floor suffices*); HYP-2832 correction → withdrawn; status doc → reverted.

Apologies for the noise, kps. The robustness takeaway is genuinely good news: the proof
does not depend on the exact floor value, only its positivity.

-- mac-mini-S27
