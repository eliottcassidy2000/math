# Retaining a chain mark gives an exact indefinitely reusable state

**PROVED by the universal endpoint-cover coefficient identity below;
FINITE-EXACT literal controls.** A particular port retention rule closes
on two integer coordinates, preserves their gcd, and has a unique exact
inverse. The H17 preparation enters its primitive domain. More generally,
every source satisfying the exact prefix-repair guard enters it without a
clock restriction. This is a constructive recursive grammar, not coverage
of all odd Hamilton counts or all positive coprime coordinate pairs.

## 1. The physical mark that must not be forgotten

For an arbitrary rooted source `(T,r)` with profile z=`(h,s,q,F)`, perform
k≥2 births at r. The k children form a homogeneous transitive module.
Let w be the last-born child, hence the first vertex in that module's
directed order. **Retain r and w**. Add y at r and retain y,w, discarding r
only from future eligibility. Add z at y and retain y,w, discarding the
newborn z from eligibility. Now choose the already retained w for future
births. Call this macro `L_k`.

Its full unmarked graph is exactly the graph after `R^k N R`; the retained
ports differ. The earlier balanced return retains y,z. This macro retains
y,w. It does not recover a forgotten vertex or assert that three ports
have been maintained with a two-port controller. Two live names suffice:

1. Until the last of the k births, keep r and the newest child.
2. At the last such birth, keep r,w.
3. At y's birth, replace r by y, keeping w.
4. At z's birth, discard z, keeping y,w.

In the six-command convention of
[LIVE_PORTS L2](../../adaptive_state_induction_20261103/LIVE_PORTS.md),
after the root-child interface `(r,w)` the last two commands are **02**.
Command0 births at slot0 and replaces it by y; command2 births at y but
discards the newborn. The new active root w is slot1. For the next macro,
birth at that available slot, retain its newest child, and forget the
previous filler y. No nonlocal source selection is hidden in iteration.

With fixed controller-slot names, the complete literal words are
`1^k02` when the active root occupies slot0 and `3^k45` when it occupies
slot1. The active slot alternates after each macro. Thus repeated macros
need no implicit slot renaming. The receipt retains the input slots,
exact command word, output slots, and new active slot for the literal
entry and renewals.

The chain remains a homogeneous TT_k after the last two births: r,y,z
all beat it; it beats every old nonroot vertex. The chosen w has indegree
exactly3, regardless of the old root's indegree or k. All k+2 births
preserve strongness when the input is a nontrivial strong tournament
(order at least3). The vacuously strong singleton is excluded from this
strongness assertion: its first birth gives TT2. There are no temporary
vertices and no old-arc changes inside `L_k`.

## 2. A universal finite coefficient identity

First take just one chain child w and then the two births above. At w the
resulting profile, for every actual rooted source, is

    h* = h−s−2q+5F,
    s* = F−q,
    q* = h+2F,
    F* = 4h−s−q+10F.                                (7)

In particular `F*=h*+3q*−s*` and the new defect is exactly `2q*`.
This small map is a recovered operation, not a new claim about the size
of the rooted grammar. In the old physical coordinates
`(h−s,q−s,s,F−h−q+s)`, its row-action matrix is

    [[5,2,1,6],
     [3,2,0,4],
     [3,3,0,6],
     [4,1,1,4]].

It is the word003 representative in the 28-map alphabet of
[ROOTED_RESETS R2](../../two_port_return_induction_20261103/ROOTED_RESETS.md)
and its `rooted_resets.json` certificate. That representative forgets y
instead of z before the final seal; the pointed parent at w is the same.
The new use is the R-delay closed chart, its inverse and invariant content,
and its connection to exact-target neutral preparations. The direct proof
below also verifies the precise retained-chain interpretation.

The identities in (7) are not an interpolation claim. Here is a finite
coefficient proof, directly using the inherited endpoint-cover bijection
from [ROOTED_RESETS R1](../../two_port_return_induction_20261103/ROOTED_RESETS.md).

Temporarily retain only the old root r. Write `*` for a forgotten endpoint
and e for the number of old Hamilton paths ending at r. For an old source
of order at least two, the only relevant endpoint-cover signatures and
their exact weights are

| signature | weight |
|---|---|
| `(*,*)` | h−s−e |
| `(*,r)` | e |
| `(r,*)` | s |
| `{(*,*),(r,r)}` | q |
| `{(*,*),(*,r)}` | F−e−q |

There may be at most one component ending at a forgotten vertex: future
newborns beat every forgotten vertex, so such a component can never be
continued. For a birth at a retained root, choose at most one component
ending there to enter the newborn and at most one component starting at a
different vertex to leave it. Do not choose the same component on both
sides, since that would close a cycle. Merge through the newborn, or make
it a singleton, then apply the declared forgetting map. Deleting the
newborn reverses this construction, so each child cover is counted once.

Apply this kernel three times: create w retaining `(r,w)`; create y at r
and replace r by y; create z at y and forget z. Read at w the four exact
observations: one-component covers; one-component covers starting at w;
two-component covers having singleton w; and covers of one or two
components having a component ending at w. In formal variable order
`(h,s,q,F,e)`, the coefficient rows are exactly

    (1,−1,−2, 5,0),
    (0, 0,−1, 1,0),
    (1, 0, 0, 2,0),
    (4,−1,−1,10,0).

The e column vanishes identically. `chain_reuse.py` implements this finite
symbolic calculation independently, with integer coefficient vectors; it
does not derive these rows from sample tournaments. Together with the
kernel's cover bijection this proves (7) for all sources of order≥2.
The singleton boundary follows by adjoining a forgotten global sink before
starting, exactly as in ROOTED_RESETS R1: all four input observations and
all resulting observations at the live mark are unchanged.

For the full macro `L_k`, first perform k−1 R births, then apply (7).
Put

    a=h−s, d=F−q, v=2^(k−1), c=a+3q−d.

Since k≥2, the pre-interface rooted state is
`(a+q+(v−1)d,q,q,vd+q)`. Substitution into (7) gives

    s* = vd,
    q* = 3vd+c,
    h* = q*+3s*,
    F* = 4q*+2s*.                                  (8)

Thus every `L_k` lands on the same algebraic two-coordinate chart. Its
positive domain in (9) requires d>0; arbitrary weak sources can have d=0,
whereas the strong sources and admissible prefix entries used here have
d>0. Deleting w is
also a direct check of q*: it reduces the R-chain length by one, so its
Hamilton count is the balanced-return count at k−1. The k≥2 restriction
keeps that deletion on the positive-chain formula; the one-child identity
(7) remains valid separately.

## 3. An exact closed primitive grammar

Let an **actual** rooted source lie on the chart

    C(q,s)=(q+3s, s, q, 4q+2s),    q odd>0, s even>0.       (9)

Its defect is `2q`; its rooted content is `gcd(q,s)`. Formula (8) proves
that, for every k≥2, the physical macro `L_k` returns to this chart with

    s' = v(3q+2s),
    q' = q+3s',        v=2^(k−1).                    (10)

The matrix on `(q,s)` is

    [[9v+1,6v],
     [3v,  2v]],         determinant 2v=2^k.

The exact content identity is stronger than determinant intuition:

    gcd(q',s') = gcd(q,v(3q+2s)) = gcd(q,s),          (11)

because q is odd and v is a power of two. Positivity, q odd, s even, and
strongness are retained. Therefore any actual primitive seed in this
chart admits **every finite or infinite choice of k≥2**, with primitive
rooted state at each finite stage. Every fresh root-child seal has content1
by ROOTED_RESETS R3. No repeated neutral preparation or source lookup is
needed after entering (9).

There is also a universal primitive balanced exit at every stage. Since
delta=2q and `d=3q+2s`, the inherited `R^ell N R`, ell≥1, has content
`gcd(q,d)=gcd(q,s)=1`. If additionally `3∤s`, a C3 cycle always enters
the usual three-divisibility guard: its balanced pair has

    t=2^ell(3q+2s)−q,
    b=3·2^ell(3q+2s)+q,
    b−2t ≡ 2^(ell+1)s ≠0 (mod3).

The predicate `3∤s` is preserved by (10), giving a convenient one-cycle
bridge. It is not required for a bridge to exist: the general finite
connector below handles every primitive balanced exit.

**General finite connector, as a corollary of inherited R and C3.** For
any actual coprime odd balanced pair `(b,t)`, R itself stays balanced and
primitive, with `(b,t)→(2b−t,t)`. It retains the physical old root and
its new child; it does not rename a forgotten mark. C3 can introduce a
common divisor only at prime3, and is primitive exactly when
`b−2t≠0 mod3`. If this condition fails, b and t are both nonzero modulo3;
one R makes the new b zero modulo3 while t stays a unit. C3 is then legal
and gives the usual guard `3|t'`, `3∤b'`. Hence **every primitive balanced
source reaches that guard using at most one R and one C3**. The inherited
mod11 bridge then enters the eight-state C3/C5 controller: apply C5 unless
its ratio is7, and otherwise apply C3 followed by C5.

Conversely, an actual primitive balanced source has d=b−t,c=3t in (8),
so direct `L_k` entry has content `gcd(3,b−t)`. If b−t is a unit modulo3,
it enters the primitive chart immediately for every k≥2. Otherwise b=t
is nonzero modulo3; C3 is legal and gives the usual guard, from which
every `L_k` is primitive. Thus **every primitive balanced source enters
the chain chart with at most one preliminary C3**.

All cases are the four projective phases `r=t/b` in P1(F3):

| r | bridge to usual C3 guard | bridge to primitive chain |
|---|---|---|
| 0 | C3 | L_k |
| 1 | C3 | C3 then L_k |
| 2 | R then C3 | L_k |
| infinity | C3 | L_k |

At r=0 the guard already holds, so the displayed C3 is optional. The
script checks these rules on all eight nonzero residue pairs, without
claiming that every formal pair is an actual source. The nearest inherited
entry proofs, RESET_JUNCTION J4 and
[ROOTED_RETURN RC2](../../rooted_carry_closure_20261103/ROOTED_RETURN.md),
give source-specific bridges; the table records this general corollary of
their inherited physical maps. These connectors change the full count.
They do not preserve unit readiness or identify the two grammars' images.

## 4. Entry from the H17 preparation and from exact prefix repairs

For the prepared H17 state `(17,4,9,28)`, take k=2. The actual graph has
ten vertices, with retained-chain root profile

    C(135,38) = (249,38,135,616).                     (12)

The two ports are the N newborn and the newest R-chain child. Both names
appear in the receipt; the root has indegree3. Here gcd(135,38)=1 and
`3∤38`, so (10) can be iterated indefinitely and all balanced exits above
are permitted. One literal `L_2` renewal gives H5907 on fourteen vertices;
an `L_3` renewal from the same seed gives H11679 on fifteen vertices.

The [one-arc H9 repair](../history/ONE_ARC_H9.md) gives a smaller scalar
entry at the same order. Its repaired root state `(9,4,5,10)` first uses
one N birth, giving `(15,5,9,19)` with the newborn retained. Now a=d=10
and c=27, so every `L_k`, k≥2, enters primitively. At k=2 this is
`C(87,20)=(147,20,87,388)` on ten vertices, with full literal rows and
fixed-slot replay recorded. Direct `L_k` from the H9 repaired source,
without that N, has content5; the preparatory N is material.

More generally, a unit-ready source `(h,q−1,q,h+1)` after one N has
a=d=h+1 and c=3h. It enters primitively by every `L_k` exactly when
`3∤(h+1)`. In particular this holds for the two-width ready sources in
[H9_INFLATION](H9_INFLATION.md) when both widths are odd, giving h≡0
modulo3. Every primitive balanced source `(b,t,t,b)` already satisfying
`3|t`, `3∤b` also enters by every `L_k`: here d=b−t,c=3t, so their gcd
is1. This supplies the reverse interface from the accepted C3/C5 domain
back into the retained-chain grammar.

More generally the one-position W_n preparation in [PROOF.md](PROOF.md)
has L=`2^(n−2)`, d=4L+3 and c=5L+1. Every k≥2 therefore reaches

    s*=2^(k−1)(4L+3),
    q*=3s*+5L+1,
    gcd(q*,s*)=gcd(5L+1,4L+3).

This gcd divides11, by `5(4L+3)−4(5L+1)=11`. It equals11 precisely
when L≡2 modulo11, equivalently n≡3 modulo10, since 2 has order10
modulo11. Hence every n≥3 **except n≡3 mod10** enters the primitive
chart at every k≥2, with no accepted-clock restriction. At the excluded
sizes it enters the same actual chart with permanent content11; (11)
prevents this particular grammar from repairing that content. The earlier
balanced-root clock theorem still applies at those sizes.

The wider exact-target entry follows from the separately proved
[prefix repair](../EXACT_REPAIR.md). For its repaired state

    z=(15h,0,5h,21f),
    h odd, f even, gcd(5h,21f)=1,

one has `d=21f−5h` and `c=35h−21f`, both odd. Formula (8) gives
primitive chart entry for every k≥2 because

    gcd(c,d)=gcd(30h,d)=1.

Indeed d is odd; its residues modulo3 and5 are nonzero under the stated
guard; and `gcd(h,d)=gcd(h,21f)=1`. The physical exact-prefix construction
supplies the actual source, so no arbitrary positive vector is declared
realizable. It also has `3∤d`, hence `3∤s*` and the accepted balanced exits
remain available indefinitely. This joins a count-preserving exact-source
repair to a closed, content-preserving sequence of further constructions.

## 5. Unique inverse addresses and the coverage boundary

Every output of (10) uniquely determines its previous state and k. From
an output `(Q,S)` recover

    q = Q−3S,
    v = 2^v2(S),       k=v2(S)+1,
    s = (S/v−3q)/2.                                  (13)

Here `3q+2s` is odd, so the power of two in S identifies v exactly. A
positive-chart predecessor exists exactly when q is positive and the last
numerator is a positive multiple of4; this ensures s is positive and even.
It is still necessary to check that the predecessor belongs to the
chosen actual seed language. Repeated inversion strictly decreases H,
since every forward step has `s'=v(3q+2s)>s`.

Every exact-prefix entry in §4 is outside this positive-even-s image:
its hypothetical inverse would have q=c=`35h−21f<0`, since the two
empty-slot contributions already give f≥2h. Independently, its formal
`s=42f−55h` is odd because h is odd. Either obstruction supplies a
definite seed boundary; the negative q already settles these actual inputs.

For the H17 seed `(135,38)`, inversion would give q=21 and s=−22; it is
outside the positive image. Therefore the grammar from that seed is a
free rooted tree of finite words in k≥2: each profile has one last letter
and one parent, and no word returns to the seed. Equations (13) give a
terminating membership test for exactly this seed language: invert until
reaching `(135,38)` or a failed positivity/evenness condition. The stored
letters replay the actual graphs and port transport. This is uniqueness
of **rooted profiles and addresses**, not of full Hamilton counts or
unmarked tournament classes.

The growth law makes the remaining boundary explicit:

    H' = 4v H+(14v+1)q.

It is not an exact ternary-digit recurrence. The primitive chart contains
formally positive points that need not be realized by these constructors;
the inverse and graph program retain the source-image condition. Also
q'≡q modulo3, so the H17 seed has only counts divisible by3 while it stays
inside this chart. Its balanced exit followed by C3 can leave that residue
restriction. Indefinite physical reuse solves a genuine interface problem;
it does not prove the full odd spectrum.

## 6. Reproduction

    python3 -B 05-knowledge/results/h17_exact_repairs_20261112/constructors/chain_reuse.py

The script emits the universal formal five-signature coefficient
certificate, not a fit; independently evaluates literal endpoint-DP
states for the ten-vertex H17 entry and two renewals through order15;
checks retained names and root indegree; and records a longer address
using the proved two-coordinate recurrence, explicitly marked formula-only.
Normal and optimized Python retain every check via `require`.
