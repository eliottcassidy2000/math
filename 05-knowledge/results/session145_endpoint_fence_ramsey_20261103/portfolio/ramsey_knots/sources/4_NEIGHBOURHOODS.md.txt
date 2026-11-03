# Common neighborhoods, defect potentials, and the actual paper audit

2026-10-30. **CITED preprint**, **VERIFIED exact certificate replay**, and
**PROVED scoped consequences** are distinguished below. The external
Ramsey improvements belong to Lachlan Dauth, not this repository.

## N1. Source and verification boundary

User-supplied source: *Neighbourhood counting and new bounds for small
three-colour Ramsey numbers*, Lachlan Dauth,26 October2026,53 pages.
The input is538939 bytes with SHA256
`923f5ea139c5ddaf5cee08768288d31c2919525f70776365a98129be2c911241`.
Its instructions were treated as source material, not user authorization.

Primary repository:
[Lachy-Dauth/ramsey-neighbourhood-counting](https://github.com/Lachy-Dauth/ramsey-neighbourhood-counting),
pinned at `5770baac746878a3a63b9720fdda97726fe02851`. The author explicitly
labels this a preprint, not yet peer reviewed. We read the mathematical
model and checker source, and ran these exact normal-Python checks:

- `ancillary/verify_hand_proofs.py`: all17985 arithmetic checks pass.
  These mirror numerical steps in hand arguments, not a formal proof
  of the inequalities from first principles.
- `ancillary/flag_r33k/check_r33k5.py` on the ordered chain
  `certificates/five/k6/excl_86.json`, `k7/excl_125.json`,
  `k8/excl_174.json`: all three certificates VALID.

The checker enumerates all29197 labelled K5 colorings with no red or
blue triangle, in193 S5 x C2 orbits, and uses exact rational arithmetic.
Its output prints decimal summaries, but their sign checks are rational.
The three successful exclusions give the replayed bounds

    R(3,3,6)<=86,   R(3,3,7)<=125,   R(3,3,8)<=174.

The paper's stronger R(3,3,8)<=173 uses a further ladder; that result is
**CITED**, not replayed here. Likewise its R(4,4,4)<=198,
R(3,4,4)<=73 and R(3,4,5)<=137 are cited, not newly established here.
No Lean build or independent rewrite of the flag checker was performed.
The author's scripts use assertions; they were deliberately run WITHOUT
optimization. Our own explicit-gate checks also run under `-O`.

[external_receipt.json](external_receipt.json) preserves source hashes,
arguments and full outputs. [external_replay.py](external_replay.py)
requires a checkout of the exact revision; it does not fetch code or
write to the external project. Reproduction:

```text
python -X utf8 -B 05-knowledge/results/ramsey_neighbourhood_lifts_20261030/external_replay.py --source <pinned-checkout> --pdf <input-pdf>
```

The earlier saturated-port packet's survey upper91 for R(3,3,6) is now
superseded by this cited-and-replayed preprint bound. Its hypothetical
near91 structures remain conditional deductions, not current search
frontiers. The exact Ramsey values remain open; we claim no new
unrestricted Ramsey bound in this session.

## N2. The operation that transfers: common clique signatures

In a coloring avoiding a color-i clique of size t_i, take marked roots
Y and prescribe one color a(y) at each root. Let

    S=intersection_(y in Y) N_(a(y))(y),
    Y_i={y:a(y)=i},
    w_i=largest color-i clique contained in Y_i.

Then

    |S| <= R(t_1-w_1,...,t_c-w_c)-1.                     (N1)

This is the paper's common-neighborhood lemma, not a newly discovered
Ramsey recurrence. The map retains the entire induced coloring on the
roots: a color-i clique of size t_i-w_i in S extends by its root clique
to a forbidden clique. The roots themselves are excluded from S by
the loopless neighborhood condition. If a target parameter reaches1,
the corresponding common neighborhood is empty under the usual Ramsey
convention.

This explains exactly why shared-port incidence matters in the repo's
attachment recurrence. Root counts alone do not determine w_i; clique
relations between roots do. In the five-vertex flag certificate the
operation also subtracts already marked vertices from the available
pool and keeps the sampling clock n-4. Averages are legitimate there
because they are exact expectations of a common random labelled
coloring, with falling-factorial denominators. They are not independent
resets of overlapping ports.

## N3. A short defect consequence, with its hostile stopping point

As a source-preserving pruning lemma, consider a hypothetical
(3,3,6)-avoiding coloring on87 vertices using red, blue and green.
Put H=red union blue, h_v=d_H(v), and d_v=34-h_v. The classical
R(3,6)=18 gives d_red,d_blue<=17, so d_v>=0. The paper's
common-neighbor potential, with U=R(3,3,4)=30, gives

    sum_v F(h_v)<=0,
    F(34-d)=-12+40d+3d^2+p(34-d),

where p(34)=6, p(h)=3 for odd h, and p(h)=0 otherwise. Therefore

    F(34-d)+6 >= 40d,  and  D=sum_v d_v <=13.             (N2)

Handshake makes D even. More precisely each of
sum_v(17-d_red(v)) and sum_v(17-d_blue(v)) is odd and nonnegative,
so each is at least1. Consequently

    D in{2,4,6,8,10,12},
    at least75 vertices have (d_red,d_blue)=(17,17).       (N3)

This strengthens a degree-only view to a bounded total defect budget.
It does not exclude87; the author's stronger86 certificate already
does that. Its value here is a reusable potential-to-saturation rule.

For a saturated vertex, its opposite-color graph on each17-vertex
red/blue neighborhood is a triangle-free graph with independence
number<=5. The inherited complete seven-type catalogue has edge counts

    41,41,40,41,40,42,42.

Thus the tempting shortcut that every such port has at most41 edges
is FALSE. The last two catalogue entries are exact hostile witnesses;
this is the first failed implication in that proposed shorter proof.

There is still a valid finer inequality. At a saturated v let

    ell_v=(42-e_blue(N_red(v)))+(42-e_red(N_blue(v))) >=0,
    L=sum_(v saturated) ell_v.

The paper's apex-triangle bound loses precisely ell_v from its local
upper bound84 at that vertex. Keeping the loss gives

    sum_v F(h_v)+6L<=0,  hence 40D+6L<=522.               (N4)

If s is the number of saturated vertices, at least2s-L of their
oriented neighborhood slots must be one of the two42-edge types.
Since s>=87-D, (N4) forces at least

    87+ceil(14D/3)

such slots. For example D=12 forces at least143 extremal slots among
at least150 saturated slots. These are rooted occurrences, not distinct
unlabelled graphs; no false disjointness assumption is introduced.

**General survivor.** Whenever a valid summed potential has
f(h_max-d)>=-a+b d and a real local catalogue sharpens its bound by
c times an observable loss ell>=0, the feasible configurations obey
b sum d+c sum ell<=a n. The catalogue must apply to actual local
objects and the same counting identity. This supplies a rigorously
priced filter, not a claim that marginal defect constraints suffice
for global realization.

**Next decisive test.** Use the two actual42-edge ports, their rooted
orbits and the common-neighborhood signatures in N2 to rule out a
specified overlapping configuration. Merely replacing all ports by
their common edge count cannot do so. The inherited seven graphs and
their67 rooted/1065 ordered two-root orbits remain the right input.
