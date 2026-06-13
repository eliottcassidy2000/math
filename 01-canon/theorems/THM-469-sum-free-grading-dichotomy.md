# THM-469: the sum-free grading dichotomy — why the seam is 2-adic (parity closure, leading digits, and the Schur arity)

**Status:** PROVED (parts A–B by hand, short proofs below) + COMPUTED/VERIFIED
(parts C–E: Glucose3 + complete verifiers, witnesses independently re-verified;
script `04-computation/erdos592_parity_closure_kp0611.py`, output
`05-knowledge/results/erdos592_parity_closure_kp0611.out`).
**Source:** kind-pasteur-2026-06-11-S1 (HYP-2390/2391). Answers THM-464 D's sharp
open note and THM-465 C's caveat; corrects THM-464 D's interpretive slogan.

**Setting.** The (sign, v_p) feature algebra on pairs of [t]^n (per coordinate:
('=',) for zero gap, else (sign, v_p(|gap|))), the b-ary subgrid game Q_b(n,t)
(THM-464 A), and the feature-quotient SAT of THM-465: a feature rule is a set T of
classes, blue iff feature ∈ T; it must be triangle-free and hit every b-ary subgrid.

---

## A. The unit dichotomy (PROVED)

Let p be prime, L_v = {g ≥ 1 : v_p(g) = v}, and for u ∈ (Z/p)^× let
L_{v,u} = {g : v_p(g) = v, g/p^v ≡ u (mod p)} (the leading-digit refinement).

**Lemma A1.** L_v is sum-free for every v ⟺ p = 2.
Proof. x = p^v a, y = p^v b, p ∤ ab ⟹ x+y = p^v(a+b), and v_p(x+y) = v ⟺
p ∤ (a+b). For p = 2: a, b odd ⟹ a+b even — escape is FORCED (parity closure:
the unit group of Z/2 is trivial, so units must cancel). For p ≥ 3: a = b = 1
gives a+b = 2 ≢ 0, so (p^v, p^v, 2·p^v) is a Schur triple inside L_v at every
scale v. ∎

**Lemma A2.** For EVERY prime p and every (v,u), the refined class L_{v,u} is
sum-free; indeed x, y ∈ L_{v,u} ⟹ x+y ∈ L_{v,2u} when p is odd (2u ≢ u since
u ≢ 0), and v_p(x+y) ≥ v+1 when p = 2.
Uniform statement: **the leading-digit grading is sum-free for every p, and p = 2
is the unique prime where it coincides with the bare valuation grading** (the unit
group (Z/2)^× is trivial). ∎
(Range-verified to 4096 for p = 2, 3, 5; violations for bare p = 3, 5 start at
(1,1,2).)

## B. The mono-collapse theorem (PROVED): odd-p algebras die by arithmetic progressions

**Theorem B1.** Fix n ≥ 1, b ≥ 2, t ≥ 3 and an odd prime p ≥ b. The (sign, v_p)
feature algebra is feature-UNSAT for the b-ary game on [t]^n.
Proof. Call a class a *cone class* if every entry is ('=',) or (±, 0). (i) For
each direction d ∈ {−1,0,1}^n \ {0} (lex-leading sign +), the 3-term AP
x, x+d, x+2d fits in [t]^n (t ≥ 3) and its three pair-features all equal the cone
class f_d, since v_p(1) = v_p(2) = 0 for p ≥ 3. The triangle clause of a mono
triple (f,f,f) is the UNIT clause ¬f: every cone class is forced OFF. (ii) The
standard b-ary subgrid on values {0,…,b−1} at every node has all per-coordinate
gaps of magnitude ≤ b−1 < p, hence v_p = 0 when nonzero: its hitting clause is a
disjunction of cone classes only. (i)+(ii) contradict. ∎

**Theorem B2 (parity closure: p = 2 has NO mono triples, ever).** In the
(sign, v₂) quotient on [t]^n, no triple (f,f,f) is realizable as a triangle, at
any t, n. Proof. The lex-leading coordinate of f is (+, v); the three gaps there
are a, b, a+b > 0; v₂(a) = v₂(b) = v forces v₂(a+b) ≥ v+1 ≠ v (Lemma A1). ∎
(Census: 0 mono classes at t = 3, 4, 5 vs 13 for p = 3 from t = 3 — the 13 =
(3³−1)/2 cone classes.)

So the answer to THM-464 D's sharp open note: **the algebraic property the v₂
classes have and the v₃ classes lack is sum-freeness of the gap grading** —
forced by triviality of (Z/2)^× ("odd+odd=even"); for odd p the unit escape
1+1=2 plants a Schur triple (d, d, 2d) inside every class at every scale, and the
game detects it through the cheapest possible objects, the 3-term APs.

## C. The mechanism is exactly the UNSAT core (COMPUTED)

At (n,t,b,p) = (3,4,2,3) the minimized UNSAT core of the feature instance is
**14 clauses**: the hitting clause of the {0,1}³ subgrid (a disjunction of the 13
cone classes) + the 13 AP unit clauses, one per cone class. Nothing else is
needed — the zero-CEGAR death of THM-465 C is Theorem B1 verbatim.

## D. The branching control (COMPUTED, settles HYP-2390): the seam follows the SCHUR ARITY, not the branching

b = 3 game at n = 3: (sign,v₂) SAT at t = 4 (32 edges), t = 5 (240), t = 6 (1910),
each independently re-verified; (sign,v₃) feature-UNSAT at t = 4 (and B1 proves it
for all t ≥ 3). The v₂/v₃ asymmetry persists IDENTICALLY at ternary branching.
**Correction to THM-464 D's slogan**: the seam does not follow the branching
through the algebra — it follows the arity of edge composition. Triangles close
over 2-term gap sums (the doubling map d ↦ 2d); the "2" in 2-adic is the "2" of
pairs-of-edges, not the "2" of binary subgrids. (THM-464's data are untouched;
only the interpretive note is superseded.)

## E. The leading-digit rescue (COMPUTED, settles HYP-2391)

(sign, v₃, leading digit): SAT at (3,4) (474 edges) and (3,5) (1625 edges), both
independently re-verified — the refined 3-adic algebra carries witnesses exactly
where bare (sign,v₃) is feature-UNSAT, as Lemma A2 predicts. Control: bare
(sign, v₅) is feature-UNSAT at (3,6) (lazy=0). Quotient-identity note: at t = 4
the algebras F3LD, F2-jet, and full gap discrimination coincide (gaps {1,2,3} all
separated), explaining identical instance statistics; they diverge from t = 5,
where F3LD's gap partition {1,4}|{2}|{3} and F2's {1,3}|{2}|{4} are BOTH sum-free
partitions of {1,…,4} — the partition-shape isomorphism trap of THM-465 C
dissolves: what matters is the additive structure, not the partition shape.

## Honesty

- Sum-freeness of the grading defeats the AP/mono collapse — the dominant
  small-t obstruction — but is NOT proved sufficient for per-size SAT in general;
  F3LD beyond t = 5 and its uniformity are untested (the t-uniform question for
  ANY algebra is THM-470's business).
- B1 needs p ≥ b only for the cone-purity of the standard subgrid clause; for
  p < b the statement is unproved here (no computed instance either).
- All feature-UNSAT verdicts have THM-465 semantics: the algebra is too coarse;
  they say nothing negative about Q_b(n,t) itself.

**Cross-refs:** THM-446 (dyadic rungs = sum-free classes in the Sidon/Schur
ladder — the same grading recurs), THM-453 F4 (composition-freeness as the graded
Schur condition), THM-464/465, HYP-2390/2391 (both CONFIRMED → this file).
