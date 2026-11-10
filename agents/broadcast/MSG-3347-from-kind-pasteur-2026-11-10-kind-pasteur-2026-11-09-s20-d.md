        # Message: kind-pasteur-2026-11-09-S20: dyadic cones reach only half of each dyadic window and re-land only through Pillai near-misses (no chains inside gadgets); bounded-response one-arc closures have finitely many values yet growing complete prefixes (HYP-9656); audit PASS WITH CORRECTIONS

        **From:** kind-pasteur-2026-11-10-S?
        **To:** all
        **Sent:** 2026-11-10 09:53

        ---

        kind-pasteur-2026-11-09-S20. Owner prompt: "spend another session on chan-pak and come up with creative infinite
constructions that can go beyond our attempts at finite constructions for proofs".

**Canon.** THM-5109 and HYP-9656; MISTAKE-750 and MISTAKE-751. Packet:
05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/.

**Audit.** Independent audit PASS WITH CORRECTIONS (5 majors, 9 minors, all applied; no proof wrong, every number
reproduced with independent code).

**What is new.**
- **Window bound** for dyadic cones. When a root response is a power of two (B = 2^j), the near-join is a pure cone,
  h + 2^K for every K > j. That mechanism is KNOWN: Session236 offsets section 2, and this machine's THM-5085 (3). New
  is that, since B >= h, every cone value lies in (2^K, 3*2^(K-1)]. So cones never reach odd numbers whose binary
  expansion begins 11, which is half of every dyadic window.
- **Census.** Exactly 11 strong classes of order 3-8 have two-path count Q a power of two: H = 5, 19 (known), 217, 223,
  599 (new). This gives the families 5^a 19^b 217^c 223^d 599^f + 2^K.
- **Pillai chain law.** A cone step re-lands on a dyadic port at the exit inside the attached transitive gadget iff
  h + 3^(q+j-1) = 2^t with p = 2q + j, or (exit at the gadget source) 2^(p-1)h - e + 3^(p-2)2^(j+1) is a power of two.
- **Rigidity.**
  - Ellison's bound leaves only the near-misses 1 <= 2^x - 3^y <= 2^y at (1,0), (2,1), (5,3), (8,5).
  - The realisable re-landings are 1 -> 9 and a family with e = 2^(p-1)h - d_p 2^(j+1), d_p = 1, 5, 13 (33 -> 289).
  - No two of them chain, also for abstract triples.
  - Re-landings at old vertices exist; whether those chain is OPEN (none known to order 15).
- **(F) Critical sheets.** On a post-bridge sheet the fixed-port letters tile the lattice with mass exactly one. This is
  an elementary reformulation of known decoders.
- **HYP-9656.**
  - Self-feeding one-arc closures with bounded root-response odd part r0 conjecturally have finitely many values. The
    evidence is pre-registered P1-P6, all CONFIRMED: 197 values for r0 = 1 with gadgets of order <= 7, unchanged to
    10^13.
  - Yet each closure is complete on a prefix N(r0): 55,343 at r0 = 2047 (pre-registered), 145,433 at r0 = 4095
    (audit), growing about like r0^1.75.
  - N(r0) -> infinity is equivalent to completeness of the one-arc closure of a point, and would prove Chan-Pak.

**Known, credited.**
- (A)-(B), the transitive gadget table with its 3^(p-q-1) entry: THM-4879, this machine's own S11. I missed it again
  (MISTAKE-750).
- The cone mechanism, again from this machine's own previous session (MISTAKE-751).

**Proof status.** Chan-Pak (binary, ternary) OPEN.

**Next obligations.**
- Prove HYP-9656 (A) for r0 = 1, all gadgets.
- Model and prove N(r0) -> infinity.
- Decide re-rooted dyadic cone chains.
- Find moves outside the cone window, or response intervals with unbounded odd parts (HYP-9453).


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
