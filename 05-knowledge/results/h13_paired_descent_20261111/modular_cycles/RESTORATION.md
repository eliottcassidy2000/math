# Native gap compilation and paired count restoration

**PROVED IN NOTE:** a J-based difference compiler and the exact equations
for two proposed infinite resynchronization families. **FINITE-EXACT:**
the least nonempty O/Z resynchronization of the two211/213 sources and
the declared larger address-box experiment. **OPEN:** a repeatable
all-target restoration mechanism. This continues [PROOF.md](PROOF.md).

## 1. A root change turns a pure-H difference into a physical singleton

Use the actual operation J: append a vertex y beating only the retained
root r, then retain y as the new root. Its state map is

    J(H,E,U,d)=(H+U-E,H-E,H+U,H).                     (1)

The first coordinate is the inherited Z count. Deleting y gives H.
Paths ending at y are the old Hamilton paths ending anywhere except r,
giving H−E. A weighted root-start path either consists of y alone,
contributing H, or continues y→r, contributing U. This proves (1)
directly. The accompanying [dynamics proof](../dynamics/PROOF.md) develops
the full positive-coordinate order for O/Z/J.

Suppose two actual rooted tournaments have equal E,U,d and H difference
k. Their state difference is `(k,0,0,0)`. After one common J it is

    (k,k,k,k)=k*(1,1,1,1).                           (2)

But `(1,1,1,1)` is exactly the rooted singleton's state, with H(empty)=1.
All three native maps are linear on these four coordinates. Therefore:

**Theorem1 (coupled gap compiler).** For every common O/Z/J word W,

    state(W J T+)−state(W J T−)=k*state(W K1).        (3)

In particular the actual Hamilton-count gap is `k H(W K1)`. Every graph
in (3) is built by the same explicitly prescribed sequence of births and
root changes. The theorem does not assert that the W(K1) family is the
entire tournament spectrum.

For the original31/33 H13 pair, J followed by OZ gives actual counts
221/227, because OZ builds a directed triangle from K1 and hence has
H3. The gap is6. Ordinary ordinal-product multiplication already gives
gaps `2H(T)` from any unit pair and any T, so **gap existence alone is not
new**. The extra content of (3) is the same native word, retained source
names/root trajectory, and its full state difference. For comparison,
ordinal multiplication of31/33 by C3 would give93/99, a different baseline.
Equation(3) does not grant independent addition of its gap to a chosen H.

The hidden endpoint pair73/73 from the new construction illustrates a
different phase. Its difference is `(0,-2,0,0)`. O erases that difference;
Z restores a pure H difference2; J produces `(2,2,0,0)`, retaining an
endpoint difference. Thus root changes enlarge the actual defect language.
The first attempt to hide the H13 unit difference after an arbitrary old
word fails already at address B1: the complete paired one-birth search
has no output difference `(0,-2,0,0)` there. J, rather than an assumed
uniform replay of the asymmetric birth, supplies the positive extension.

## 2. Why the gap compiler does not automatically cancel the reset cost

At the211/213 same-count replacement, common address C lowers both
counts by

    R(C)=90C-12(g(C)-1).                             (4)

A standalone application of the k=2 gap compiler produces `2H(W K1)`.
That half-gap is always odd. But `R(C)/2=45C-6(g(C)-1)` has the parity
of C. For every positive even C it therefore cannot equal `H(W K1)`.
This is one exact boundary of that particular correction proposal.
For odd C the parity gate passes, but a matching gap would still not
place its lower endpoint at the required current count. No independent
gap-addition operation has been supplied.

Instead the next probe retained the actual source coupling and allowed
two different continuation words, one on each side of the same-count
replacement.

## 3. Actual paired resynchronization after a nonempty continuation

For the lower old/new211 sources, O/Z addresses C,D give respectively

    A(C)=13g(C)+334C+198,
    B(D)=25g(D)+244D+186.                            (5)

Their companion upper sources always give exactly2 more. Thus one equality
`A(C)=B(D)` restores **both** target counts simultaneously, with actual
programs and source identities. This is not the addition of an unrelated
gap certificate.

**FINITE-EXACT leastness.** The least common value with C,D positive is

    C=79, D=102,
    A(79)=B(102)=26649.                              (6)

Both upper continuations have26651. Their lower rooted states are

    old: (26649,13,42843,13),
    new: (26649,25,33157,25).                        (7)

The new source again reduces U, now by9686. Applying Theorem2 of
[PROOF.md](PROOF.md) gives another valid all-address source comparison,
with a12,b9686. The actual original tail addresses are636 in the H13
family and205 in the25 family: `636=4+8*79`, `205=1+2*102`.

The leastness assertion is bounded by the target itself, not merely by
the chosen search box. Since g>=1, any equality at H<=26649 has

    C<=(26649-211)/334=79,
    D<=floor((26649-211)/244)=108.

Those complete finite ranges are included in the replay. The claim is
only leastness for the two specified normalized O/Z source languages.

The larger declared box `1<=C,D<=32767` contains227 matched address
pairs. At152 of them the new U is no larger; at75 it is larger. Thus
resynchronization by itself does not guarantee a cost decrease. The
actual states, not just the common scalar H, decide which direction is
usable. Every match is retained in `synchronization.json` for subsequent
native source work.

## 4. Two explicit attempts to make the matches repeat

The first test appends all-one blocks in ordinary binary notation:

    C_k=2^k(C+1)-1, D_k=2^k(D+1)-1.

Because g is unchanged by appending1, equality in (5) for every k is
equivalent to

    C+1=122t, D+1=167t,
    13g(122t-1)-25g(167t-1)=78.                    (8)

The second test appends all-zero blocks. Set

    f(n)=n+1-(n mod2)-g(n).

The inherited recurrence gives `f(2n)=f(n)`. For k>=1,

    A(2^k C)=347 C2^k+13(1-f(C))+198,
    B(2^k D)=269 D2^k+25(1-f(D))+186.

Therefore equality for all k>=1 is equivalent to

    C=269t, D=347t,
    13f(269t)=25f(347t).                             (9)

These are exact reductions of two specific infinite-family proposals.
No odd t<=200000 satisfies (8), and no odd t<=500000 satisfies (9),
in the recorded finite tests. Restricting to odd t loses no primitive
solution, since doubling t simply extends the same run and leaves the
g/f conditions unchanged. Neither test proves a global no-go.

The positive residue is a clean coupled equation for future investigation:
instead of demanding unmarked H coverage, seek a native paired solution
of (8), (9), or a mixed periodic-word analogue, and carry both actual
states through the continuation. The matched states in (6)-(7) already
give a genuine resynchronization. What remains missing is a law that
repeats or selects such matches at arbitrary targets.

## 5. Reproduction

```bash
python3 -B 05-knowledge/results/h13_paired_descent_20261111/modular_cycles/probe.py
python3 -B 05-knowledge/results/h13_paired_descent_20261111/modular_cycles/synchronize.py
```

The first script now checks all3280 O/Z/J words of lengths0..7 in (3),
with four literal J-prefix graph controls, in addition to the previous
new-source tests. The second records the complete specified address box,
the target-dependent leastness bounds, and both finite repetition probes.
These are constructive correlation tests, not a new spectrum coverage
range and not a revalidation of the older modular literature.
