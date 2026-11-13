# Exact pairs, hidden sites, and a prime refill that lands on 7p

kind-pasteur-2026-11-13-S26 · [THM-5151](../../../01-canon/theorems/THM-5151-exact-strong-decomposable-slot-pairs-census-r-factor-sites-one-apex-b-and-a-prime-pair.md) · [HYP-9678](../../hypotheses/HYP-9678-exact-pair-members-are-uniquely-hamiltonian-and-spent-odd-cores-hide-no-sites.md) · **Chan–Pak remains open.**

The owner asked us to turn Session 310's triangle–pentagon trade into an unbounded changing-module return for
prime and 7p targets. This session first mapped every exact pair that exists in small orders. It then worked out
what such pairs can and cannot supply to a return.

## What was found

**The complete stock is tiny and uniquely Hamiltonian** ([PROOF](PROOF.md) §3).

- Through order 9, exactly 7 strong classes share their full slot polynomial with a decomposable tournament: H = 15
  at order 6, and H = 45 and 95 at order 8. All 7 have exactly one Hamiltonian cycle.
- Session 310's component budget forces even order and a two-component partner (§1).
- Inside the uniquely Hamiltonian class, generated to order 12 with Session 197's Fibonacci counts, order 10 adds
  exactly three more profiles: H = 153, 219, and the prime **317**.

**Session 310's B is a one-apex tournament** (§4). B is the transitive tournament on seven vertices plus one apex
that beats the two lowest. Every cycle runs through the apex along an increasing chain, so the closed form

    Q = t^(m+1) + 2 Σ binom(gap, 2l-2) t^(m+1-2l)

counts its 10 triangles, 30 pentagons and 7 heptagons directly. The trade partner A is not one-apex, and neither
is any other partner block except C3.

**A neutral rebase creates sites only by R-factor recombination** (§2). Write Q_T = t^β ∏ R_S over the strong
components. Equal Q means equal β and equal R-product, so the number of singleton components changes by
(o − o') + 2(e − e'). Through order 9 the hidden capacity is exactly s + d over all 198,946 classes. Here s counts
singleton components and d counts the even strong components that have a detachable partner. So an odd core never
hides a site, and an even component hides at most one. In particular, spent cores cannot be refilled
(Propositions 3 and 3′).

**The prime pair breaks the Π clock's divisibility by 5** (§5). Every value of Session 310's B^r clock is divisible
by 5. The order-10 pair with H = 317 instead gives

    H_k = 60 P_M(3) 9280^r 2^(r−k) 3^k + 317^r H(M).

With M = TT3 and r = 2, the neutral refill sits at the **prime 558,047,332,489**. Two literal one-arc transfers then
end at **1,255,606,372,489 = 7 × 179,372,338,927**, a 7p target. A Rabinowitsch-type criterion (Proposition 7)
explains the 7: since 3/2 is a primitive root mod 7, any run of six consecutive transfers meets a multiple of 7
unless 7 divides the intercept.

## What this means for a return rule

**Proved:** call cap(Q_N) the site budget of the right module. It is the largest number of trailing singletons
any tournament with the same slot polynomial can have. It does not change under neutral rebases and drops by
exactly one per transfer, so a clock fed only by neutral refills stops after cap(Q_(N_0)) transfers (PROOF
Proposition 8).

Under HYP-9678 (C), checked through order 9, a neutral rebase of a right module can only cash in detachable even
components that are already present. After a run the right module consists of odd cores, which hide nothing.

An unbounded return must therefore alternate two kinds of move:

- a **non-neutral** move that inserts a fresh detachable even component (Session 309's migration or recharge, or a
  change of frame);
- a **neutral** rebase that splits it.

That matches the Lonely Runner "collision-backed handoff", where the refill is installed before the last site is
spent ([HISTORY_TANGENTS](HISTORY_TANGENTS.md) B). Session 309's fixed-bank theorem then says the detachable stock
must keep growing. This is why HYP-9678 (B), infinitely many detachable pairs, matters.

## Tangents

[HISTORY_TANGENTS](HISTORY_TANGENTS.md) gathers the three history sweeps.

- **2647.** Session 310's prime frame value 2647 is also the largest value no uniquely Hamiltonian tournament
  attains. 231, the frame value of the order-6 pair, is also such a miss. No mechanism is known, so this is
  recorded as a coincidence.
- **Heegner numbers.** 7, 43, 67 and 163 recur among misses of restricted constructions.
- **Lonely Runner.** LRC proved a fixed-bank no-go four months before Session 309, and escaped it with
  core-dependent adaptive clocks.
- **Hypergraphs.** THM-4803 says the cyclic-triple 3-uniform hypergraph determines Q, so an exact pair is a
  connected and a disconnected hypergraph with the same profile.

## Pre-registered searches

[PREREGISTRATION](PREREGISTRATION.md) was pushed at b1b011e4d2 before these ran. Both searches cover only strong
order-11 tournaments with at most 5 backward arcs. The same program at orders 7 and 9 with at most 4 backward arcs
returns exactly the known partners.

- **S1** (order-12 detachable partners, predicted yes): **FAILED, no match.**
  - 1,291 distinct target profiles Q_X/t from the 2,584 uniquely Hamiltonian order-12 classes.
  - 3,847,591 candidates, of which 1,243,347 strong; 312,197 were fully profiled; 18 min 42 s.
  - Not pre-registered: a follow-up over "biclique" backward-arc patterns L × U (L entirely below U, |L|, |U| <= 5;
    7,469 tournaments) also found nothing. Session 310's partner A has this shape ({5,6} × {0,2}); the order-9 H =
    317 partners do not.
  - Order 12 is therefore undecided. Detachable partners there, if any, need at least 6 backward arcs and are not
    bicliques.
- **S2** (order-11 odd merges, predicted no): **CONFIRMED (no match).**
  - 655 distinct target profiles Q(Y5 => Y7)/t.
  - 3,847,591 candidate tournaments, of which 1,243,347 strong; 189,829 passed the triangle and H filters and were
    profiled.
  - No strong order-11 tournament with at most 5 backward arcs is an odd merge of a 5-core and a 7-core.

## Replay

Run `python verify_s26.py` (add `--full` for the order-10 Fibonacci count). It is standard-library only and passes
under `python -O`. `scripts/` replays the complete census from the class lists; see its README.
