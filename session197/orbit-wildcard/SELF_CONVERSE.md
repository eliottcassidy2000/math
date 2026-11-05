# Self-converse symmetric7p witnesses and the sharp order threshold

Separate Session197 supplement; the original wildcard snapshot is unchanged. Hamilton counts refer to directed vertex words. This note claims minimum ORDER, not minimum Hamilton count, in the self-converse nonrigid7p sector. The unrestricted spectrum remains open.

## 1. An explicit eleven-vertex witness

Take the rooted quotient

    Q=[6,24,26,1,9], root0,

and substitute the carousel C7, with difference set{1,2,3}, at its root. Name the seven module vertices0 through6 and the remaining quotient vertices1,2,3,4 as7,8,9,10 respectively. The complete rooted maximal-run vector is

    (a1,a2,a3)=(11,17,4).

Using the carousel's ordered cover counts (175,1638,6678), the exact count is

    H=11*175+17*1638+4*6678=56483=7*8069,

and8069 is prime. The checker emits the full outgoing rows in `self-converse-witness.json`, independently recomputes H by two subset recurrences, and verifies the complete legal-run inventory.

Negation on Z/7 reverses the carousel. The quotient permutation(1 4)(2 3), fixing0, reverses Q. Together they give the actual anti-involution

    tau11=[0,6,5,4,3,2,1,10,9,8,7].

Thus the eleven-vertex graph is self-converse. Q is strong, and a strong module substituted into a strong quotient gives a strong graph; this is also checked directly. Simultaneous rotation on the seven module vertices, fixing the other four, is an automorphism of order7. Free action on Hamilton words and8069>11 imply that the full automorphism group has order exactly7: its order divides7*8069 and also11!.

## 2. Eleven is the minimum possible order

**Theorem.** The least order of a self-converse tournament with nontrivial automorphism group and Hamilton count7p, p prime, is11.

The [S196 orbit-spectrum theorem](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/orbit-spectrum/PROOF.md), Sections5-6, proves that the minimum order without self-converseness is10. It also excludes the small-prime exceptional symmetry cases and proves that every remaining symmetric7p tournament has full automorphism group of order7. Thus only order10 needs to be ruled out here.

At order10 this group has exactly one seven-element orbit and three fixed vertices. An anti-automorphism conjugates the automorphism group to itself: if tau reverses all arcs and g preserves them, tau g tau^(-1) preserves them. Consequently tau preserves that unique seven-element orbit and its three-element complement.

Every self-converse tournament has an anti-involution. To see this directly, its automorphism group has odd order by free action on the odd set of Hamilton words. For any anti-automorphism tau, tau squared is an automorphism, so tau has order2m with m odd. Then tau^m is an involution and still reverses every arc.

An involution preserving each of two odd-sized sets has a fixed point in each set. But an anti-automorphism cannot fix two distinct vertices: their connecting arc would have to reverse while keeping both endpoints fixed. This contradiction rules out order10. The explicit construction above attains11. No claim of the least possible H in this sector follows.

## 3. The module-indecomposable fourteen-vertex witness is self-converse too

For the fourteen-vertex graph in the original wildcard packet, on A_i,B_i with both Paley difference set{1,2,4} and cross offsets{1,3}, define

    tau(A_i)=B_(-i), tau(B_i)=A_(-i).

Its exact permutation is

    tau14=[7,13,12,11,10,9,8,0,6,5,4,3,2,1].

Within each block negation reverses every Paley arc. For a cross arc A_i->B_j, applying reversal and tau tests A_(-j)->B_(-i), whose difference is again j-i; hence the same offset condition holds. Therefore tau14 is an anti-involution.

On Hamilton words let J reverse the word and apply tau14 entrywise. A J-fixed word is specified uniquely by its first seven entries: choose one vertex from each tau14 pair and order those seven choices. They must form a directed path, and its final vertex v must beat tau14(v). The second half is forced as the reversed tau-image. The anti-automorphism property makes every edge of the second half legal. Conversely every fixed word has precisely this form.

Complete enumeration of128 choices and7! half-orders gives

    Fix(J)=6111.

An independent DP on used vertex-pairs and the last chosen vertex gives the same number. A literal fixed word is

    [0,1,2,4,5,6,3,11,8,9,10,12,13,7].

The graph's simultaneous translation g satisfies tau14 g tau14=g^(-1). Thus rotation and J act dihedrally on Hamilton words. Nonidentity rotations fix no Hamilton word, and the seven reflections are conjugate, each fixing6111. Burnside's lemma gives exactly

    (15144269+7*6111)/14=1084789

orbits under this dihedral word action. Equivalently, the2163467 ordinary rotation-orbits carry an involution with6111 fixed orbits. This is a precise converse-sensitive path-orbit correction for this single graph, not a formula for an entire A000568 stratum.

## 4. Reproduction and inheritance

Run `python3 verify_self_converse.py` and `python3 -O verify_self_converse.py`; compare outputs and generated certificate bytes. The script uses the original packet's checked DP/run/profile utilities and separately implements literal and pair-state counts for reflection-fixed words.

The lower order10 threshold and free-action theorem are inherited exactly as cited. The new order11 conclusion adds the anti-involution obstruction and an explicit attaining construction. No literature-priority claim is made. Self-converse here means isomorphic to the converse tournament; it is not numerical primality, module-indecomposability, or the different property H=|Aut|.
