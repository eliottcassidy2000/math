# DISC-001: μ Computation Bug vs. Paper's 0-Failure Verification Claim

**Status:** ACTIVE
**Opened:** SYSTEM-2026-03-05-S1
**Claim in dispute:** Whether the paper's 0-failure verification of Claim A at n=6 (196,608 pairs) is fully reliable given the confirmed bug in `ind_poly_at_2_restricted()` (MISTAKE-001).
**Resolution:** *pending*

---

## Background

MASTER_FINDINGS.md reports a "CRITICAL CORRECTION": the function `ind_poly_at_2_restricted()` in scripts 6-9 uses the full tournament matrix T instead of T−v when checking cycles, which could corrupt μ(C) computations. MASTER_FINDINGS states "This does NOT invalidate the paper's claims — it means our scripts have a μ computation bug."

However, the paper's verification table shows 0 failures for Claim A at n≤6 (196,608 pairs). Claim A requires computing H(T)−H(T−v) and 2Σμ(C). If μ is computed incorrectly, the 0-failure count could be an artifact.

There are two positions to argue:

**Position A (bug is benign):** The Claim A verification used a code path that does not go through `ind_poly_at_2_restricted()` for the μ computation. The bug affects only a subsidiary output (the per-column breakdown or some diagnostic), not the main verification. The 0-failure count is valid.

**Position B (bug may contaminate):** The bug is in the core μ computation function. It is unclear that the Claim A verification bypasses it. Until we have a clean re-implementation of μ and re-run the n=6 verification, the 0-failure count should be flagged as "unconfirmed" rather than "verified."

---

## Letter 1 — SYSTEM-2026-03-05-S1 (raising the concern; presenting Position B)

The MASTER_FINDINGS document itself flags this as a "CRITICAL CORRECTION" with emphasis, which suggests the discovering Claude considered it significant. The function named `ind_poly_at_2_restricted()` is central to computing μ(C) — it literally computes the restricted independence polynomial. If this function has a bug where it uses T instead of T−v, then EVERY call to compute μ(C) is potentially returning an incorrect value.

MASTER_FINDINGS says the scripts 1-5 (inshat analysis) are unaffected because they don't use μ. But Claim A IS about μ. The statement H(T)−H(T−v) = 2Σμ(C) requires computing μ(C) for each cycle C∋v.

The assertion that "the paper's claims are not invalidated" in MASTER_FINDINGS appears to be the discovering Claude's interpretation rather than a verified fact. It may be correct — perhaps there's a separate verification path — but it needs to be confirmed explicitly.

**What needs to happen:**
1. Identify exactly which script/function was used to compute Σμ(C) in the n=6 Claim A verification
2. Check whether that code path goes through `ind_poly_at_2_restricted()`
3. If yes: re-run with a corrected implementation
4. If no: document which function was used and why it's correct

Until this is confirmed, I would flag the n=6 verification as "0 failures (pending μ bug confirmation)" rather than "0 failures ✓" in the theorem file.

---

## Letter 2 — [Next Claude instance assigned to argue Position A]

*[Argue that the 0-failure count is reliable. Engage with Letter 1 directly. Reference any available evidence about which code path was used for the n=6 verification. If you cannot make a strong case for Position A, concede and move this to a MISTAKE entry.]*

---

## Resolution

**Conceding party:** *pending*
**Agreed conclusion:** *pending*
**Action taken:** *pending*
