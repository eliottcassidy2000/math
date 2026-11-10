# Dyadic cones, the Pillai chain law and finite bounded-response closures

**Session:** kind-pasteur-2026-11-09-S20. **Canon:** THM-5109, HYP-9656.

**Status.**
- PROVED (new): THM-5109 (C)-(F). Ellison's 1971 bound and THM-4810's strong fibres are cited inputs.
- KNOWN and restated: THM-5109 (A)-(B), which are this machine's own THM-4879 (A)-(B) (S11, converse convention). The
  pre-audit draft typed them as new; MISTAKE-750.
- FINITE-EXACT: the censuses and the pre-registered tests P1-P6, all CONFIRMED.
- OPEN: Chan-Pak (binary and ternary). Nothing here claims it.

Proofs are in [PROOF.md](PROOF.md); the predictions and their outcomes are in [PREREGISTRATION.md](PREREGISTRATION.md).

## Owner prompt (2026-11-09)

> spend another session on chan-pak and come up with creative infinite constructions that can go beyond our attempts
> at finite constructions for proofs

## The idea and the answer in one paragraph

The even-graph proofs (THM-5040, THM-5055) work by *cones*: a move that adds 2^n to the count, independently of the
structure. Tournaments do have cones. When a root response is an exact power of two, B_v(T) = 2^j (a *dyadic port*),
the near-join gives H(T) + 2^K for **every** K > j. Explicit infinite families follow, such as 5^a 19^b + 2^K. A cone
induction, however, needs the cone to be repeatable: the new tournament must again carry a dyadic port. For transitive
one-arc gadgets this happens exactly when h + 3^a = 2^t (or a variant with the root end-count e). The equation is
Pillai's, for the primes 2 and 3. Ellison's bound leaves only the four near-misses 2 - 1 = 4 - 3 = 1, 32 - 27 = 5 and
256 - 243 = 13, and no two cone steps ever chain (THM-5109 (E)). Generalising, every self-feeding one-arc closure whose
root responses have odd part at most r0 is finite. Yet it is complete on an initial segment that grows like a power of
r0 (HYP-9656). A construction that reaches N therefore needs root responses with odd part growing with N. Bounded
banks and pure cones cannot do it; that is the precise sense in which a proof must go beyond finite constructions.

## Inheritance (read first)

- **Near-joins** h + 2^k B: HYP-9453 and THM-4865. HYP-9453 also has B >= h and proposes the response-interval
  induction.
- **The one-arc join transfer** and its closed entries: SEAM.md and THM-4875 (B). The descent lemma: THM-4875 (C).
- **The transitive gadget table** with its 3^(z-u-1) entry, and the plane S = 0 with its covering densities: THM-4879
  and HYP-9471, this machine's own S11. I missed them before writing; MISTAKE-750.
- **Fixed-port tails** and their decoders: THM-4888, THM-5085, THM-5096, THM-5091.
- **Root changes:** HYP-9480 and HYP-9494.
- **Strong fibres** and floors: THM-4810.
- **The two-path cover count** Q(X) = sum_S H(X[S]) H(X[V - S]): Session229's proper-graft row, and THM-5107
  (Session242, pushed during this session), which also proves 4 | Q.
- **Concurrent sessions.** Session237 (THM-5090..5093) and Session238 (THM-5094..5096) worked on Chan-Pak today on
  other prompts. Session242 (THM-5107, THM-5108) works with grafts and the same Q. No one had the dyadic-port lane.

## Results

**(A)-(B) Inherited tools.**
- **Source.** THM-4879 (A)-(B), this machine's S11, in SEAM's converse convention. Restated in the direct convention
  in PROOF.md sections 1-2.
- **The response row.** In the one-arc join of (T, v) with (G, u, z), R' = B_z(G) h - J_1 e + (J_2 + K) R, and
  B_z(G) >= 1 always. So the new response always contains a copy of the old value.
- **The closed transitive table.**
  - Exits at or after the entry have R' = 2^(m-q) h'.
  - Exits before the entry have R' = 2^(m-q) h + 3^(p-q-1) 2^(m-p+q) R.
  - Re-checked on 4,440 brute-force cases.
- **The S = 0 plane.** THM-4879's S is e here, so its plane S = 0 is the e = 0 sheet. Its covering-density analysis,
  including the sharp cone condition of at most 0.27 at order 8, is HYP-9471, and is not repeated.

**(C) Dyadic cones and certified infinite families.**
- **Cones.** If B_v(T) = 2^j then H(T) + 2^K is a Hamiltonian-path count for every K >= j + 1.
- **Source ports.** The source port over X has response Q(X). Exactly 11 strong classes of order <= 8 have Q a power
  of two: H = 5 (Q = 2^5), 19 (2^8), 217, 223 (2^12) and 599 (2^13).
- **Products.** Since Q is multiplicative, 5^a 19^b 217^c 223^d 599^f + 2^K is a Hamiltonian-path count whenever
  K > 5a + 8b + 12c + 12d + 13f.
- **How common.** Dyadic ports exist for 69 of the 148 attained odd values <= 300 at order <= 8.

**(D) The Pillai chain law.** A cone step from a dyadic port lands on a dyadic port inside the gadget iff:
- q >= 2, p = 2q + j and h + 3^(q+j-1) = 2^t; or
- q = 1 and 2^(p-1) h - e + 3^(p-2) 2^(j+1) = 2^s.

**(E) Rigidity.**
- **(E1)** The only solutions of 1 <= 2^x - 3^y <= 2^y are (1, 0), (2, 1), (5, 3) and (8, 5), with values 1, 1, 5,
  13. This uses Ellison's |2^x - 3^y| > 2^x e^(-x/10) and a finite check.
- **(E2)** The dyadic cone steps are exactly:
  - 1 -> 9 (two realisations);
  - a family needing e = 2^(p-1) h - d_p 2^(j+1), with (p, d_p) in {(3, 1), (5, 5), (7, 13)}. The realised example is
    33 -> 289, from (h, e, B) = (33, 4, 64).

  The values 5 and 13 are excluded by THM-4810's singleton strong fibres.
- **(E3)** No two dyadic cone steps chain. This holds even for abstract triples, cross-checked for j <= 17.

**(F) Fixed-port sheets are exactly critical.** On a post-bridge sheet h_u = c, the fixed-port letters S^j B tile the
odd lattice (h, R): every point has exactly one predecessor, and the total mass sum_j 2^(-(j+1)) is 1. A finite bank
has no surplus mass; mass comes only from root changes.

**HYP-9656: bounded-response closures are finite.**
- **Setup.** Take triples with odd part of B at most r0, closed under one-arc gadgets (exit in the gadget) and sinks,
  starting from one vertex.
- **(A) Finiteness.** The closure is finite:
  - r0 = 1: 197 values with gadgets of order <= 7 (P1);
  - r0 = 63: saturation with gadgets of order <= 6 (P2);
  - r0 = 1023 and 2047: saturation at 93,159 and 233,115 values (P3, P6).
- **(B) Complete prefix.** The complete prefix grows: every odd n < 13, 35, 35, 67, 127, 431, 2,303, 5,153, 20,921,
  55,343 except 7 and 21, for r0 = 3, 7, ..., 2047 (P5 out of sample). Over r0 = 63 to 2047 the growth is roughly
  r0^1.8.

## What this means for a proof

- **Pure cone inductions are impossible.** This is a theorem for transitive gadgets and evidence (P1) for all gadgets
  of order <= 7.
- **Bounded-response families are finite** (HYP-9656 (A)). This covers finite banks and fixed-port tails, whose
  dynamics is exactly critical (F).
- **A proof must carry unbounded root responses.** Two routes remain that do this:
  - HYP-9453's response intervals, which must now come from non-dyadic, unboundedly complex responses;
  - the full one-arc orbit (HYP-9640), whose odd-part-bounded slices are the finite closures here.
- **A quantitative target.** Show that the r0-closure is complete up to some N(r0) that tends to infinity.
  - This would prove Chan-Pak.
  - It is equivalent to cofiniteness of the full one-arc closure of a single vertex, a statement of the same kind as
    HYP-9640.
  - The data suggest N(r0) is about r0^1.8.

## Files

| file | purpose |
|---|---|
| `resp.c` | H, Q (two-path covers), all root responses B_v and end counts e_v of a tournament |
| `tl.py` | Python helpers: one-arc joins, transitive tournaments, batch calls to `resp` |
| `gadget.c`, `gadget2.c` | 3 x 3 one-arc gadget matrices: all labelled gadgets of order <= 5, all class gadgets of order <= 7 |
| `check_tt.py` | Theorem B against brute force (p >= 2; `out/check_tt.out`); `out/check_tt_p1.out` holds the p = 1 rows |
| `check_e_identity.py` | e_v = O - E on 468 random ports |
| `h1_list.py` | (E1): the four Pillai near-misses, and a numerical check of Ellison's inequality for 28 <= x <= 60 |
| `e_classify.py` | (E2)/(E3) over abstract triples: no classification violation, no two-step chain |
| `p4_chains.py` | P4: realised dyadic bases of order <= 8 |
| `d_finite.py` | a first, coarser two-step enumeration; its two arithmetic solutions violate B >= h (PROOF (E2)) |
| `dyadic_census.py`, `q_dyadic_strong.py` | the dyadic census and the strong classes with Q a power of two |
| `closure.py`, `rclosure3.c`, `rclosure4.c` | bounded-response closures (P1-P3, P5-P6, complete prefixes) |
| `out/` | every quoted output, and the SHA-256 of the two regenerable gadget tables |

**Reproduce.**
1. Generate the classes with S19's recipe (`classgen 8 cls`).
2. Build the gadget tables: `gcc -O2 -o gadget2 gadget2.c && ./gadget2 7 <class dir> > gadgets7c.txt`, and
   `./gadget 5 > gadgets5.txt`. Their SHA-256 are in `out/gadget_tables.sha256`.
3. Run a closure: `./rclosure4 gadgets5.txt 5 1000000000 1023 26`.
4. The Python checks need `resp.exe` built from `resp.c` in the same directory.
