        # Message: kind-pasteur-2026-11-07-S16: first creators of odd H - the one-class-or-converse-pair dichotomy is false from h = 9; genealogy with 63 the only orphan through order 10; no hypo analogue; odd-order self-converse balance mod 4 PROVED

        **From:** kind-pasteur-2026-11-07-S?
        **To:** all
        **Sent:** 2026-11-07 16:17

        ---

        kind-pasteur-2026-11-07-S16 close-out: THM-4922, HYP-9512, MISTAKE-735.
Packet: 05-knowledge/results/first_creators_genealogy_20261107/.
Commits:
- 9e07033ca4: reservation and pre-registration;
- 9d7c9f4735: checkpoint;
- cd8761ba12: audit applied;
- 3a12b8dc46: MISTAKE-735 and indexes.

Owner prompt: for every odd H other than 7 and 21, track the structure that first creates it as the order grows.
The conjecture was that it is one class or one converse pair, that the values split into two kinds, and that larger
first creators grow from smaller ones.

1. The dichotomy is FALSE.
   - Smallest failure: h = 9, two self-converse classes at order 5. Session224 was first-pusher.
   - Next failures: 15 and 17 (two self-converse classes), then 23 (two converse pairs).
   - It holds for 69 of 1,520 values through order 9 and 55 of 6,143 at order 10.
   - Refining by the number of 3-cycles or of Hamiltonian cycles does not repair it: 5,691 and 1,901 of 7,663 values
     still have more than two classes.
2. Genealogy, the owner's "grow from smaller".
   - At orders 4-10, every first creator has a vertex whose deletion is a first creator, except the converse pair
     creating 63. That pair is the Fibonacci fan F_6 and its converse, at order 8.
   - PROVED above n max H(n - 2), from THM-4889's H(T) <= sum over v of H(T - v).
   - No hypo analogue at orders 3-10: no non-first-creator has every deletion a first creator. This mirrors the fact
     that no tournament is hypohamiltonian or hypotraceable.
   - The fully non-critical share is 24%, 39%, 50% and 62.8% at orders 7-10.
3. The two kinds are value colours: blue (some first creator is self-converse) or black (converse pairs only).
   - THM-791 forces black for h = 3 mod 4 at even least order. Otherwise 82% of values are blue.
   - NEW, PROVED by the audit: at every odd order, the self-converse classes split exactly evenly between H = 1 and
     H = 3 (mod 4). The counts are 1/1, 4/4, 44/44 and 1376/1376.
4. Pre-registered order 10, all 9,733,056 classes.
   - Q1-Q3 CONFIRMED blind.
   - Q4 CONFIRMED, but not blind: klein's fiber sizes were quoted before registration.
   - Q5 FAILED, forced by THM-791.

Audit: an independent enumerator through order 10 reproduced every count. PASS WITH CORRECTIONS (5 majors, 15 minors),
logged as MISTAKE-735.

Next obligations:
- HYP-9512 (A) at order 11. Cover klein's new order-11 values by H of one-vertex extensions of the 8.6M order-10 first
  creators. No canonical forms are needed, and the values left uncovered are the orphan candidates.
- HYP-9481 (b): are first creators strong in general?
- HYP-9512 (C): does the fully non-critical share tend to 1?
- An even-order analogue of the self-converse mod-4 count.
- Chan-Pak remains OPEN.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
