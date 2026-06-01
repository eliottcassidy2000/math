---
id: HYP-2004
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S526
related:
  - HYP-2003
  - HYP-1998
  - THM-381
---

# HYP-2004: LRC = harmonic covering ( |SAFE| = (1-2/n)^{n-1} + resonance corrections ); n=3 PROVED

**Methodology (S525 covering + S522/S523 roots of unity).** B_s={t:||st||<1/n}
open, measure 2/n; SAFE = complement of union; LRC(n) <=> SAFE != empty. Fourier
(safe indicator g, g_0=1-2/n, g_k=-sin(2pi k/n)/(pi k)):
> |SAFE| = sum_{ (k_i): sum k_i s_i = 0 } prod g_{k_i} = (1-2/n)^{n-1} + resonance corrections.
The main term (1-2/n)^{n-1} = opus-S524's independence value ((6/7)^13=0.1348 at n=14).

**PROVED — LRC(n=3).** Two runners a<b, gcd=1. m=2 => only 2-term resonances; the
mod-3 character is the Legendre symbol chi. Closed form (verified vs numerics):
> |B_a ∩ B_b| = 4/9 + (2/9) chi(a)chi(b)/(ab),  i.e.  |SAFE| = 1/9 + (2/9)chi(a)chi(b)/(ab).
Since chi(a)chi(b) in {-1,0,1} and ab>=2: |SAFE| >= 1/9 - (2/9)/(ab) >= 0, equality
iff {a,b}={1,2}. So a lonely time always exists ({1,2}: t=1/3 boundary). LRC(n=3)
holds, {1,2} the unique tight case (the AP / regular polygon). THM-worthy.

**Why n>=4 is open.** m>=3 => |SAFE| is an m-fold intersection whose resonance sum
has 3-term and higher resonances sum k_i s_i = 0; these do not factor to one
character (= opus's multi-way correlation). Main term stays positive (-> e^-2);
the proof needs a SIGN/lower bound on the higher-resonance correction.

**VERIFIED (`lrc_small_n_covering_proof_s526.py`):** n=4..7 exact scans, 0 LRC
failures; AP always tight; AP NOT unique tight at n=5,6 ({1,3,4,7}; {1,3,4,5,9}).

**OPEN / next:**
- (A) Bound the 3-term resonance correction for n=4 (the first case past pairwise);
  if the n=3 character trick has an n=4 analogue (mod-4 odd-harmonic character),
  n=4 closes.
- (B) General lower bound: resonance correction >= -(1-2/n)^{n-1} for all speed
  systems = LRC(n). Identify the worst resonance configuration (the AP wall).
- (C) Refine S525's "wall-only => AP": there are several margin-1/n tight sets.

**Files:** `04-computation/lrc_small_n_covering_proof_s526.py` (+.out). Reflection:
`07-reflections/lrc-small-n-covering-proof-s526.md`.
