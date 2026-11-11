# Arbitrarily many fresh return/birth steps inside one full-tail controller

**PROVED construction, conditional only on the explicitly inherited kernel
and profile-trade identities.** Local initial profiles are FINITE-EXACT.
This gives any prescribed finite number of actual reusable opportunities,
not a regenerative finite carrier and not a Chan–Pak proof.

## 1. The added idea: keep unused modules as physical fuel

Let A,B be the Session261 six-vertex complete-profile trade, with
(H,Q2,Q3)=(17,224,1449). Here Qj permits empty slots in an ordered spanning
j-path cover. Equality of the complete cover profile means an untouched A
module can be replaced by B in ANY marked surrounding host, even after
other modules have changed. Thus earlier births do not consume later
modules' exchange guards. This is stronger than having one collision at
one initial state, but it consumes a fresh module each time.

Replace A's vertex1 by a transitive pair, or B's vertex0 by a transitive
pair, to obtain A+ and B+. Literal initial moment triples are

    A+: (35,676,5649),       B+: (33,612,5163).          (1)

The rows and direct exact subset counts are in `fuelled_receipt.json`.
A still-unused slot can therefore be reselected physically as A or B
without changing ANY current host response, and then take one of two
different actual births. Do not assert a new A+/B+ trade: their profiles
are different and that slot is now spent.

## 2. A controller that preserves the full outer tail guard at every stage

Use the inherited two-port quotient Q=(198,196,136,211,167,207,20,64),
ports5,7 and exterior marked root0. Put a word of k>=1 slots in each port.
Initially every slot is A or B. Perform the same birth choice in matched
slots on both sides, permitting arbitrary permutation of the slots between
sides. After a A-births and b B-births, a+b<=k, each side has

    h=17^(k-a-b)35^a33^b,
    q=224^(k-a-b)676^a612^b,
    r=1449^(k-a-b)5649^a5163^b.                        (2)

Proof of the moments: every ordered j-path cover of an ordinal sum is
exactly a choice of ordered j-path covers in its factors, concatenating
each slot separately. Hence Qj multiplies, including empty portions.
The quantities needed by the inherited kernel are

    H=h, C2=q-2h, C3=r-3q+3h.

Both sides have the same triple throughout, so their normalized C2/H and
C3/H agree. This is precisely the inherited guard for preservation of
EVERY positive tail and the exceptional bare-apex head under outer-port
exchange. The physical outer-port exchange therefore remains available at
every stage, not only before the first birth.

Moreover each unused slot is still a module of the whole host. Applying
the inherited full-profile trade to it preserves the current full tail,
regardless of how many other slots have already grown. This gives k
successive exact-return/choice/birth stages in one initially supplied
carrier. A,B reconfiguration within either side can be performed separately;
the synchronized births are the steps that retain the outer exchange guard.

## 3. Explicit common response and initial cases

Inserting (2) into the inherited formulas gives, at EVERY stage,

    K=9h^2-4hq+37q^2,
    U0=-8h^2+19hq+10q^2-4hr+62qr,
    E0=7hq.                                           (3)

The rooted source has actual H=D=K and W=0. Its common positive tail and
bare-apex head are respectively

    K+2^t U0, t>=1,
    K+U0-E0=h^2+8hq+47q^2-4hr+62qr.                   (4)

The negative coefficients in (3) arise only from the change from Cj to Qj;
the inherited Cj formula has nonnegative coefficients and literal graphs
are emitted. No virtual subtraction is used as a construction.

For k=1 the full starting state and the two possible spent states are:

    initial: K1843881, U20596980, E26656,
             bare22414205, tail(t1)43037841;
    A birth: K16824497, U240979528, E165620,
             bare257638405, tail(t1)498783553;
    B birth: K13787145, U199343808, E141372,
             bare212989581, tail(t1)412474761.

These are exact substitutions into a proved inherited kernel, not a
claim of large-order independent recounting. The initial kernel itself is
checked literally: H149,U0=430,E0=14. Each new local atom in (1) is counted
literally. This keeps verification small while supplying actual starting
cases for the arbitrary-depth formula.

## 4. Reuse boundary and next constructive gate

This realizes arbitrarily long FINITE controlled sequences by prepaying k
slots. A source after a+b births has 12k+6+2(a+b) vertices; exactly k-a-b
matched slots remain. It does not produce an infinite trajectory from one
fixed initial graph, nor give a digit recurrence hitting arbitrary target N.
The outer full-tail response survives arbitrary exchanges throughout, while
births intentionally change it according to (2)-(4).

There is a precise new design target: replenish one spent slot with a fresh
A/B trade without changing the required interface, or find a target coding
in the three multiplicative moment coordinates. Neither follows from having
unbounded pre-supplied fuel. The construction nevertheless separates a true
multi-use return supply from unsupported extrapolation of a one-use trade.

Run local checks:

    python3 -B 05-knowledge/results/session262_interface_supply_20261111/constructor/fuelled_controller.py

Emit a literal three-slot carrier after the A birth, B birth, and an unused
slot, with its common tail:

    python3 -B 05-knowledge/results/session262_interface_supply_20261111/constructor/fuelled_controller.py --states 1,2,0 --tail 2

States0,1,2 mean unused, A-birth, B-birth. The two port words are emitted in
opposite orders; this makes the permutation aspect physical, rather than
only exchanging indistinguishable scalar labels. No count-level coverage
claim is attached to this restricted grammar.

## 5. Actual queue rearming at a fixed exposed position

The word can do more than keep fresh slots scattered through its interior.
For every j, Qj of an ordinal word is the product of its factors' Qj and
is therefore invariant under any factor permutation. For nonempty modules,

    Qj = sum_(i=1)^j binomial(j,i) Ci,

so triangular inversion shows that the COMPLETE nonempty cover profile is
also permutation-invariant. Replacing a port word by any factor permutation
therefore preserves every marked host context, even after its factors have
grown. This imports no bounded-three-moment inference.

Expose the first block as a service position. Grow it in both sides using
the selected A/B birth, then rotate that spent block to the end of its word.
The rotation is an exact full-profile return at the NEW post-birth count;
a fresh block now occupies the same exposed first-block position. If k>=2,
this is a genuine rearm after the first use. Repeat until all k blocks have
been used. A and B have different internal root labels, so each service
step retains the corresponding literal root map rather than identifying
vertex1 with vertex0 without a rebase.

The rotation physically reorients the interblock arcs according to the new
ordinal order; it is not a claim that mere relabelling transports the menu.
All modules remain unselected and root0 stays exterior. This gives an exact
queue of k reusable guarded operations, with a transparent remaining-fuel
coordinate. Its failure boundary is exhaustion, not failure at the second
step. Infinite regeneration remains the separate missing operation.
