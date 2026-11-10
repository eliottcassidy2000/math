# Dyadic cones, the Pillai chain law and bounded-response closures

**Session:** kind-pasteur-2026-11-09-S20. **Canon:** THM-5109, HYP-9656.

**Status.**
- PROVED (new): the window bound and the census of THM-5109 (C), the chain law (D) and the Pillai rigidity (E).
- PROVED (elementary): (F), a reformulation of known decoders.
- KNOWN and restated:
  - (A)-(B) are this machine's own THM-4879 (MISTAKE-750);
  - the cone mechanism of (C) is Session236 and this machine's THM-5085 (3) (MISTAKE-751).
- CONJECTURE with saturation evidence: HYP-9656.
- FINITE-EXACT: the censuses and the pre-registered tests P1-P6, all CONFIRMED.
- OPEN: Chan-Pak (binary and ternary). Nothing here claims it.

**Independent audit:** PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md)). No proof was wrong and every number
reproduced. Five MAJORs (priority, scope, two misquotations, a literally false HYP statement) and nine MINORs, all
applied; MISTAKE-751.

Proofs are in [PROOF.md](PROOF.md); the predictions and their outcomes are in [PREREGISTRATION.md](PREREGISTRATION.md).

## Owner prompt (2026-11-09)

> spend another session on chan-pak and come up with creative infinite constructions that can go beyond our attempts
> at finite constructions for proofs

## The idea and the answer in one paragraph

The even-graph proofs (THM-5040, THM-5055) work by *cones*: a move that adds 2^n to the count, independently of the
structure. Tournaments have cones too. When a root response is an exact power of two, B_v(T) = 2^j (a *dyadic port*),
the near-join gives H(T) + 2^K for every K > j. This mechanism is already in Session236 and in this machine's S19.

This session adds four things.
- **Window bound.** Cones reach only half of each dyadic window. Since B >= h, every cone value lies in
  (2^K, 3 * 2^(K-1)], so odd values whose binary expansion begins 11 are never cone values.
- **Product family.** The strong tournaments of order 3-8 whose source response Q is a power of two have
  H = 5, 19, 217, 223, 599. So 5^a 19^b 217^c 223^d 599^f + 2^K are Hamiltonian-path counts.
- **Chain law.** A cone step lands on a dyadic port at the exit inside the attached transitive gadget exactly when
  h + 3^a = 2^t, or a variant involving the root end-count e holds. This is Pillai's equation for the primes 2 and 3.
  Ellison's bound leaves only the near-misses 2 - 1 = 4 - 3 = 1, 32 - 27 = 5 and 256 - 243 = 13. The realisable steps
  are 1 -> 9 and one e-relation family (33 -> 289), and no two of them chain. Re-landings at old vertices do occur;
  whether those chain is OPEN.
- **Bounded responses.** Self-feeding one-arc closures whose root responses have odd part at most r0 appear to have
  finitely many values, yet each is complete on a growing prefix (HYP-9656).

Read together, a cone-style induction needs both moves outside the cone window and root responses whose complexity
grows with the target.

## Inheritance (read first)

- **Near-joins** h + 2^k B: HYP-9453 and THM-4865. HYP-9453 also has B >= h and proposes the response-interval
  induction.
- **The dyadic cone itself.**
  - Session236 `offsets/PROOF.md` section 2: rays 2^(m+5) + 5 and 2^(m+4) + 9 from the seeds (5, 0, 32) and
    (9, 2, 16). It also proves R = 1 iff T is transitive and rooted at its sink.
  - This machine's THM-5085 (3): "pure dyadic rays" (S19).
  - I missed both; MISTAKE-751.
- **Two-path cover count** Q(X) = sum_S H(X[S]) H(X[V - S]):
  - Session227 `multiroot` section 5, with (5, 32) and Q multiplicative;
  - Session228 and Session229, with (19, 256), the product law (`dynamic` eq. (6)) and the source-response formula
    (`overlap` eq. (1));
  - THM-5107 (Session242), with 4 | Q.
- **The one-arc join transfer**: SEAM.md and THM-4875 (B). The descent lemma: THM-4875 (C).
- **The transitive gadget table** with its 3^(z-u-1) entry, and the plane S = 0 (= the e = 0 sheet) with its covering
  densities: THM-4879 and HYP-9471, this machine's own S11. Missed before writing; MISTAKE-750.
- **Fixed-port tails** and their decoders: THM-4888, THM-5085, THM-5096, Session197 ONE_CORE_BINARY_INVERSE.
  - Session206's (ln 2) X mass is the quantitative precursor of (F)'s "critical" reading.
- **Root changes:** HYP-9480 and HYP-9494.
- **Strong fibres and floors:** THM-4810.
- **The r0 = infinity baseline:** THM-5070 (4) (S18). The same order-<= 5 gadget library, with no odd-part bound,
  attains every odd value up to 2 * 10^5 except 7 and 21.

## Results

**(A)-(B) Inherited tools.**
- **Source.** THM-4879 (A)-(B), this machine's S11, in SEAM's converse convention. Restated in the direct convention
  in PROOF.md sections 1-2.
- **The response row.** In the one-arc join of (T, v) with (G, u, z), R' = B_z(G) h - J_1 e + (J_2 + K) R, and
  B_z(G) >= 1 always. So the new response always contains a copy of the old value.
- **The closed transitive table.**
  - Exits at or after the entry have R' = 2^(m-q) h'.
  - Exits before the entry have R' = 2^(m-q) h + 3^(p-q-1) 2^(m-p+q) R.
  - Re-checked here: h' on 4,440 cases, e' on 3,360 and R' on 3,768 (audit: 7,956).

**(C) Dyadic cones.**
- **Mechanism (KNOWN).** If B_v(T) = 2^j, then H(T) + 2^K is a count for every K >= j + 1.
- **New: the window bound.** Every cone value lies in (2^K, 3 * 2^(K-1)]. So no cone reaches 7, 13, 15, 25-31, 49-63,
  ...
- **New: the census.** Exactly 11 strong classes of order 3-8 have Q a power of two: H = 5 (2^5), 19 (2^8), 217, 223
  (2^12) and 599 (2^13). The values 217, 223 and 599 are new.
- **The product family.** 5^a 19^b 217^c 223^d 599^f + 2^K is a count whenever K > 5a + 8b + 12c + 12d + 13f.
- **How common.** Dyadic ports occur for 69 of the 148 attained odd values <= 300 at order <= 8.

**(D) The Pillai chain law** (new port at the exit inside the attached gadget). A cone step from a dyadic port
(h, e, 2^j) lands on a dyadic exit iff:
- q >= 2, p = 2q + j and h + 3^(q+j-1) = 2^t; or
- q = 1 and 2^(p-1) h - e + 3^(p-2) 2^(j+1) = 2^s.

**(E) Rigidity** (same scope).
- **(E1)** The only solutions of 1 <= 2^x - 3^y <= 2^y are (1, 0), (2, 1), (5, 3), (8, 5). This uses Ellison's
  |2^x - 3^y| > 2^x e^(-x/10) (Waldschmidt section 2.2.2) and a finite check.
- **(E2)** The realisable steps are:
  - 1 -> 9 (two base triples);
  - the family e = 2^(p-1) h - d_p 2^(j+1), (p, d_p) in {(3, 1), (5, 5), (7, 13)}, realised by 33 -> 289 from
    (33, 4, 64).

  The near-miss values 5 and 13 carry no suitable port (THM-4810's singleton fibres).
- **(E3)** No two such steps chain. This holds even for abstract triples (PROOF.md), and was cross-checked on all
  abstract triples with j <= 14 (audit).
- **Outside the scope.** One-step re-landings at old vertices exist (1 -> 5, (1, 0, 4) -> 9, (33, 0, 128) -> 289).
  None chains in the audit's search (outputs of order <= 15). Re-rooted chains are OPEN.

**(F) Fixed-port sheets are exactly critical** (elementary reformulation of known decoders). On a post-bridge sheet
h_u = c, the fixed-port letters S^j B tile the odd lattice (h, R) off the line R = -c. Every point with R > 0 has
exactly one predecessor, and the total mass sum_j 2^(-(j+1)) is 1.

**HYP-9656: bounded-response closures.**
- **Setup.** Take triples with the odd part of B at most r0, closed under one-arc gadgets (exit in the gadget) and
  sinks, starting from one vertex.
- **(A) Finite value set** (conjectural, OPEN). The triple set is always infinite through sinks; the conjecture is that
  the value set is finite. Saturation evidence:
  - r0 = 1: 197 values with gadgets of order <= 7, unchanged up to 10^13 (P1, audit);
  - r0 = 63 with order <= 6, and r0 = 1023, 2047 with order <= 5 (P2, P3, P6);
  - audit stress runs up to r0 = 4095.

  Saturation needs a wide margin. The audit's order-6, r0 = 255 run has 2 values above 10^9.
- **(B) Complete prefix.** The closure contains every odd n < N(r0) except 7 and 21:

  | r0 | 3 | 7 | 15 | 31 | 63 | 127 | 255 | 511 | 1023 | 2047 | 4095 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | N(r0) | 13 | 35 | 35 | 67 | 127 | 431 | 2,303 | 5,153 | 20,921 | 55,343 | 145,433 |

  N(2047) was pre-registered (P5); N(4095) is the audit's. Growth is about r0^1.75 so far, possibly slowing.
- **N(r0) -> infinity.** This is equivalent to completeness of the full one-arc closure of a point (every odd n except
  7 and 21), and it would prove Chan-Pak.

## What this means for a proof

- **Proved.**
  - Cones reach only half of each dyadic window.
  - Cone steps whose new port is the exit inside a transitive gadget never chain.
- **Conditional on HYP-9656 (A).** A self-feeding one-arc construction must let the odd parts of its root responses
  grow with the target. This does not concern fixed-port tails or finite banks: they keep the root in the core, their
  responses are unbounded, and their coverage is HYP-9480 (C).
- **Remaining routes.**
  - HYP-9453's response intervals, which must come from non-dyadic responses;
  - the full one-arc orbit (HYP-9640), whose odd-part-bounded slices are the closures studied here.
- **Quantitative target.** Prove N(r0) -> infinity, with the data suggesting N(r0) is about r0^1.75.

**Note on the pre-registration (audit minor 6).** PREREGISTRATION "Already seen" says "No persistence through the
standard near-join: 15 cases tested; none has a dyadic port afterwards". That is true only of those 15 cases. P4's own
33 -> 289 is a near-join (TT_3 with p = m = 3) landing on a dyadic exit, and re-landings at old vertices also occur.
The pre-registration text is left unedited.

## Files

| file | purpose |
|---|---|
| `resp.c` | H, Q (two-path covers), all root responses B_v and end counts e_v of a tournament |
| `tl.py` | Python helpers: one-arc joins, transitive tournaments, batch calls to `resp` |
| `gadget.c`, `gadget2.c` | 3 x 3 one-arc gadget matrices: all labelled gadgets of order <= 5, all class gadgets of order <= 7 |
| `check_tt.py`, `check_tt_p1.py` | Theorem B against brute force (p >= 2 and p = 1) |
| `check_e_identity.py` | e_v = O - E on 468 random ports |
| `h1_list.py` | (E1): the four Pillai near-misses, and a numerical check of Ellison's inequality for 28 <= x <= 60 |
| `e_classify.py` | (E2)/(E3) over abstract triples: no classification violation, no two-step chain |
| `p4_chains.py` | P4: realised dyadic bases of order <= 8 |
| `d_finite.py` | a first, coarser two-step enumeration; its two arithmetic solutions violate B >= h (PROOF (E2)) |
| `dyadic_census.py`, `q_dyadic_strong.py` | the dyadic census and the strong classes with Q a power of two |
| `closure.py`, `rclosure3.c`, `rclosure4.c` | bounded-response closures (P1-P3, P5-P6, complete prefixes) |
| `out/` | every quoted output, and the SHA-256 of the two regenerable gadget tables |
| `audit/` | the independent audit: report, own programs and outputs, `run_all.sh` |

**Reproduce.**
1. Generate the classes with S19's recipe (`classgen 8 cls`).
2. Build the gadget tables: `gcc -O2 -o gadget2 gadget2.c && ./gadget2 7 <class dir> > gadgets7c.txt`, and
   `./gadget 5 > gadgets5.txt`. Their SHA-256 are in `out/gadget_tables.sha256`.
3. Run a closure: `./rclosure4 gadgets5.txt 5 1000000000 1023 26`.
4. The Python checks need `resp.exe` built from `resp.c` in the same directory.
