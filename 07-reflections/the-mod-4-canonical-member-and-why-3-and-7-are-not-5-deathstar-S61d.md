# The mod-4 canonical member: why odd switching classes have a name, why 3 and 7 are not 5, and where the parity actually lives

**CREDIT / COLLISION (MISTAKE-199, added post-hoc):** the central theorem here — that
Babai–Cameron Remark 7.4 is 0 at **every** odd n, via the all-even canonical member at n≡1
mod 4 and the all-**odd** canonical member at n≡3 mod 4 — was derived FIRST and IDENTICALLY by
**klein-S338** (the score-parity law + n≡1 case), **opus-S409 THM-1460**, and **kp-S128c116
THM-1465** (the σ-equivariance argument: σ fixes a member iff it fixes a vector of the coset,
= a union of cycles; the empty union at n≡1, a single odd cycle at n≡3). They pushed before my
checkpoint; my "THM-1465" was a rediscovery and is CEDED (file deleted). This reflection is
therefore **convergent synthesis, not a primary claim**, on §2–§5. My one genuinely distinct
contribution is **§1, the 3/8-mass-at-K8 confirmation** of boxeph's HYP-8295 handoff
(E[eps over S₈]=3/8), which the Babai–Cameron trio did not touch. Credit accordingly.
**death-star-2026-07-20-S61d** (confirming boxeph HYP-8295; Babai–Cameron part convergent w/ opus/kp/klein; owner: run the n=8 census and confirm the 3/8 mass
at K8; each odd-sized tournament corresponds to a natural number; odd-valued functions ↔
tournaments ↔ even graphs/even functions; Babai–Cameron Remark 7.4 is 0 at n≡1 mod 4 via the
unique even member; 3 and 7 are alike mod 4 while 5 resonates with 1 and 9). All of it is one
theorem (THM-1465) plus one confirmed prediction, and together they say: **the parity of C(n,2)
gives every odd-order switching class a canonical name, and forces every parity phenomenon into
the even-n case.**

## 1. The confirmed handoff: 3/8 mass at K8

boxeph (HYP-8295) predicted that the DFGPR eps-parity statistic, quarter-quantized for small n,
first leaves the quarter lattice at n=8 with mass 3/8. Confirmed exactly:
E[eps over S_n] for n=2..10 is **0, 0, 1/2, 1/2, 0, 0, 3/8, 3/8, 0**, with
eps(σ) = (−1)^{Σ over even-length cycles of (length/2)}. At n=8, E[eps] = **3/8**, so
f+ = (1+3/8)/2 = **11/16** — boxeph's K8 numbers exactly. The complete-core class K_n has
automorphism group S_n, so this S_n-mean *is* the K8 class mass; the census confirms the S_m
computation. The quarter law breaks at 8 because the Wallis sequence (1/2)_k/k! = C(2k,k)/4^k
governing it first takes a non-{0,1/2} value there — the repo's two-sheeted (1−x)^{−1/2}
fiber-fraction constant, now running the switching-class parity statistics (boxeph's fourth
domain).

## 2. The canonical member (THM-1465): every odd-order class has a name

The parity of the out-degree sequence gives each odd-order switching class a **unique canonical
member**:
- **n ≡ 1 (mod 4)**: a unique member with **all even** out-degrees;
- **n ≡ 3 (mod 4)**: a unique member with **all odd** out-degrees.

The mechanism is pure F₂ linear algebra: for odd n, switching moves the out-degree parity
vector p(T) by the sum-zero subspace, bijectively, so exactly one member hits the all-even
vector 𝟎 (possible iff C(n,2) even iff n≡1) or the all-odd vector 𝟙 (possible iff C(n,2) odd
iff n≡3). C(n,2) = n(n−1)/2 has parity 0,0,1,1 as n ≡ 0,1,2,3 (mod 4) — that single arithmetic
fact runs the whole show.

The canonical member is fixed by every automorphism (relabeling preserves "all even"/"all
odd"), so **Aut(class) = Aut(canonical member)** and **no odd-order switching class has a
fixed-point-free automorphism.** This settles Babai–Cameron's Remark 7.4 ("we cannot do this")
for **all** odd n — their n≡1 mod 4 case is the all-even member; the n≡3 case, which they did
not state, is carried by the all-**odd** member. Verified: 0 of 2 (n=3), 0 of 64 (n=5). For
**even** n there is no canonical member and the count is positive — 8/8 at n=4, 640/1024 at
n=6 — which is exactly the enumeration they left open, and exactly where the eps-parity mass
of §1 lives. **Odd n is organized (count 0); even n carries the parity (mass 3/8, count > 0).**

## 3. Why 3 and 7 are not 5 — the Paley dichotomy is the canonical-member dichotomy

The owner's observation that 3 and 7 pattern together while 5 goes with 1 and 9 is the mod-4
split made structural:
- **n ≡ 3 (mod 4): 3, 7, 11.** C(n,2) **odd**; canonical member **all-odd** (near-regular);
  this is exactly where **Paley tournaments** exist (q ≡ 3 mod 4: the 3-cycle, Paley₇, …). The
  tournament side.
- **n ≡ 1 (mod 4): 5, 9, 13.** C(n,2) **even**; canonical member **all-even**; this is where
  **Paley graphs** (self-complementary) exist (q ≡ 1 mod 4). The graph side.

So 3 and 7 are alike because both are **tournament primes** — odd C(n,2), all-odd canonical
member, a Paley tournament — while 5 is a **graph prime**, all-even member, a self-complementary
Paley graph. The QR sign χ(−1) = (−1)^{(q−1)/2} that decides "tournament vs graph" is the same
(−1)^{C(q,2)/…} parity that decides "all-odd vs all-even canonical member." One Legendre symbol,
two faces. And the exhaustiveness gap {7, 21}: **7 ≡ 3 mod 4** sits on the tournament/all-odd
side, the same side as 3 — consistent with 7 being the anomalous H-value.

## 4. Each odd tournament is a natural number — now with a canonical address

S60 gave the correspondence natural numbers ↔ tournaments (the tiling hypercube = integers in
binary; H = the multiplicative norm). THM-1465 supplies the missing piece: at odd n the
switching class — the base-path-free object of S61 — has a **canonical representative**, the
all-even (n≡1) or all-odd (n≡3) member. So "each odd-sized tournament corresponds to one of the
natural numbers" becomes precise: the odd-order switching classes are named by their canonical
members, and the canonical member is the *fixed point* every symmetry must respect — the
tournament-theoretic analogue of a number's canonical decimal form. The even-n classes, having
no canonical member, are the "ambiguous" numbers where the parity character (eps, 3/8) measures
the ambiguity.

## 5. The even/odd braid, completed

Three "evens" now sit in one picture, and they are genuinely different (S61b's lesson kept):
- **even function** (H is complement-invariant, S61b) — the *symmetry* parity;
- **even graph** (cycle space, opus THM-1430 = two-graphs, odd-n bijection) — the *degree*
  parity;
- **all-even out-degree** (THM-1465 canonical member, n≡1 mod 4) — the *score* parity.

H, the odd-valued even function, projects to a constant on the switching quotient (S61b); the
switching quotient is named, at odd n, by the canonical-member score parity (THM-1465); and the
switching classes themselves are the even graphs at odd n (opus). The odd invariant (Rédei H),
the even object (even graph / cycle space), and the even function (complement symmetry) are the
three faces, and C(n,2) mod 4 is the hinge that turns them.

## 6. The high-leverage question and the proof clarifications

**Clarification 1 (Babai–Cameron).** Their Remark 7.4 is 0 for *all* odd n, not only n≡1 mod 4;
the all-odd member supplies the missing half. The open enumeration is purely an **even-n**
problem, and THM-1465(3) locates it precisely.

**Clarification 2 (the eps proof target).** boxeph's DFGPR equinumerosity / the eps character
is exactly the obstruction to a canonical member at even n. The **high-leverage question**: is
the eps-parity statistic at even n computable as the *deficiency* of the odd-n canonical-member
bijection — i.e., does the failure of "unique all-even/all-odd member" at even n have a
generating function that is the even-n restriction of the (1−x)^{−1/2} Wallis series boxeph
found? If so, the 3/8 at n=8 is the first term where the canonical-member obstruction is
genuinely two-dimensional, and the DFGPR proof reduces to counting canonical-member deficiencies
— a concrete, finite, mod-4-graded target. *(Honest negative, this session: eps is NOT a
standard sign — eps ≠ sgn(σ), eps ≠ sgn(σ on edges), eps ≠ (−1)^{#even cycles} at n=4,5,6, and
E[sgn_edges] is the trivial 0/1 while E[eps] gives 3/8. So eps is genuinely boxeph's Wallis-mean
character and the canonical-member/deficiency bridge, if real, is not a cheap sign identity — it
remains the open target.)*

**Clarification 3 (Paley/exhaustiveness).** The tournament-vs-graph split at n≡3 vs n≡1 predicts
that the H-spectrum anomalies (the {7,21} gap, 7 ≡ 3 mod 4) and the LRC/QR resonances should be
sharpest on the n≡3 (tournament, all-odd, Paley) side. Worth testing directly: is the {7,21}
gap a mod-4 (all-odd-member) phenomenon?

## 7. Honesty and credit

§1 confirmed exactly (E[eps over S_n], the K8 = S_8 mean); credit boxeph HYP-8295 for the
prediction and the eps/Wallis framework. §2 THM-1465 is proved and verified n=3,4,5,6; it
extends Babai–Cameron Remark 7.4 (their n≡1 mod 4 statement) to all odd n. §3–§5 assemble the
mod-4 / Paley / even-odd structure — the arithmetic facts are exact, the "Paley dichotomy =
canonical-member dichotomy" is a genuine structural match (same Legendre parity), and the
natural-number framing (§4) builds on S60. §6 poses targets, not results; Clarifications 2 and 3
are conjectural. Credit: opus THM-1430 (even graphs = two-graphs, odd-n), S61b (the two-evens),
S60 (numbers ↔ tournaments), boxeph HYP-8295 (eps, 3/8), kp THM-1445/1455 (OCF truncation, the
odd/even skew-Pfaffian thread).

## Cross-links
THM-1465 · Babai–Cameron Remark 7.4 · boxeph HYP-8295 (eps, 3/8 at K8) · opus THM-1430 · S61b
odd/even reflection · S60 two-arithmetics · Paley/QR (q≡3 mod 4 tournament thread) · the {7,21}
exhaustiveness gap · CONJ-001.
