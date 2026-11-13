# A neutral H17 preparation unlocks an all-length primitive return

**PROVED in note; finite literal controls and independent source-state
review.** The new connection moves a physical root **inside** the incoming
resource carrier, then uses the inherited root-child return. It changes a
specific continuation capability at the exact target17. It does not claim
new realization of17, a new unmarked six-cycle, or Chan–Pak completion.

## 1. The named six-state carrier supplies the neutral preparation

Use the marked resource-and-cut controller in
[cycle_opening_six_return, dynamics S3–S3a](../../cycle_opening_six_return_20261112/dynamics/PROOF.md),
with A=B=singletons and resource radius m=1. Keep the physical names

    A1=0, A2=1, B1=2, B2=3, X=4, Y=5.

For phase (r,s), the four resource vertices have transitive order

    A1 ... Ar B1 ... Bs | A(r+1) ... A2 B(s+1) ... B2,

and the exterior is `X→Y→module→X`. Retain **root A1**, which lies inside
the resource module. The autonomous phase update is

    (r,s) → (2−s,r+s−1).

Two moves give the following actual strong tournaments, with all six
physical names unchanged:

| phase | resource order | outgoing rows | state at A1 |
|---|---|---|---|
| (2,1) | A1,A2,B1 \| B2 | (30,28,24,16,32,15) | (17,8,9,18) |
| (1,2) | A1,B1,B2 \| A2 | (30,16,26,18,32,15) | (17,8,9,18) |
| (0,2) | B1,B2 \| A1,A2 | (18,16,27,19,32,15) | (17,2,9,44) |

The state convention is `(h,s,q,F)`: h is the full Hamilton-path count,
s counts paths starting at the root, q is the deletion count, and F counts
ordered spanning covers (P,Q), with P nonempty ending at the root and Q
possibly empty. The table is an exact six-vertex path/cover count; the
script emits the full marked six-cycle and the actual rows.

Every phase has H17 because the resource module remains a transitive
four-tournament, so its slot polynomial stays x^4. The wrapper has
`H=1+2^4=17`. It is strong because its quotient is a directed triangle.
The incoming all-context neutrality theorem applies to roots **outside**
the module: in particular the physical root Y keeps state `(17,1,1,17)`.
It does not freeze A1's rooted response. Moving A1 to a different position
inside the module is precisely the nonuniform marked interface used here.
The six-period assertion remains about the named resource-and-cut state;
several phases have the same unmarked tournament.

The operational change is not initial primitivity: both displayed A1
states have coordinate gcd1. Instead the birth-return defect

    delta=F−h−q+s

changes from0 to20. The graph count17 and root name0 are retained during
the two preparation moves. Its indegree changes from1 to3. In particular
this operation does **not** preserve the incoming paper's separate
indegree-one readiness defect `q−s`, which changes from1 to7. The two
defect definitions are not interchangeable.

## 2. The prepared return is primitive at every positive length

Use the actual commands from
[ROOTED_RESETS](../../two_port_return_induction_20261103/ROOTED_RESETS.md):
add x with `root→x→all other old vertices`; R retains the old root and N
marks x. Their state maps are

    R(h,s,q,F)=(h+F−s,q,q,2F−q),
    N(h,s,q,F)=(h+F−s,h−s,h,h+F).

Starting from the prepared state `(17,2,9,44)`, for **every k≥1** the
actual word `R^k N R` returns

    (b_k,t_k,t_k,b_k),
    b_k=105·2^k+7,
    t_k=35·2^k−11,
    b_k−3t_k=40,
    gcd(b_k,t_k)=1.                                  (1)

Proof: put a=h−s=15, d=F−q=35 and delta=d−a=20. Direct iteration, as in
[RESET_JUNCTION J3](../../coverage_machine_induction_20261103/RESET_JUNCTION.md), gives

    t=a+q+(2^k−1)d,
    b=a+3q+(3·2^k−1)d.

Substitution proves the formulas. The integer t is odd and is−1 modulo5,
so `gcd(b,t)=gcd(40,t)=1`. Equivalently the inherited defect guard is
`gcd(20,35·2^k+9)=1`, since that argument is odd and4 modulo5.

All outputs are strong and have k+8 vertices. Every root-child birth
preserves strongness from a strong source. The same physical root A1 is
kept during preparation and all initial R commands; the command N then
deliberately moves the active root to its newborn. The first counts are
217,427,847. All have the form `7(15·2^k+1)`; no prime-infinitude claim
is made for their quotient.

The existing final root-child pair is primitive as well as the returned
rooted profile. Immediately before the final R its coordinate gcd is
`gcd(t_k,9,35·2^k+9)=gcd(t_k,9,20)=1`. The exact root-child divisor
theorem therefore gives content1 for that pair. Independently, (1) gives
content1 for every fresh root-child seal of the returned rooted graph.

The boundary k=0 is excluded deliberately. The actual word N R gives
`(105,17,17,105)`, not the formula obtained by setting k=0 in(1).

## 3. What the neutral preparation repairs

Without the two resource moves, start instead at `(17,8,9,18)`. Its
a=d=9 and delta=0. The same word `R^k N R`, k≥1, returns

    t=9(2^k+1),       b=3t.                           (2)

The returned rooted profile has gcd t>1, while the existing final
root-child pair has content9. Thus initial coordinate gcd1 alone did not
protect this chosen continuation from a singular factor trap. The neutral
H17 preparation changes precisely that outcome, yielding(1) at every
positive length.

At k=2, (1) gives `(b,t)=(427,129)`, already satisfying gcd(b,t)=1,
`3|t` and `3∤b`. It therefore enters the inherited repeatable C3/C5
balanced controller in
[NATIVE_SELECTION section7](../../source_selection_closure_20261103/NATIVE_SELECTION.md).
More generally every positive even k has that guard. This is a concrete
exact-target preparation followed by a reusable physical continuation,
not merely a matrix orbit without named source graphs.

For the eight-state branching controller, take the same explicit entry
bridge as before: apply C5 if `t/b` is not7 in P1(F11); if it is7, apply
C3 first (the ratio becomes3), then C5. The resulting ratio is0. At k=2
the original ratio is7, so that example uses the latter branch. For odd k,
t≡2 mod3 and a direct C3 return would introduce a factor3; no automatic
entry by that branch is asserted for those k.

The neutral six-state phase moves are proved neutral **before** these
selective births. After a birth singles out A1, the original four-vertex
module is no longer uniformly seen; one cannot blindly reapply the old
phase-neutrality theorem. Repeatability of the completed primitive source
comes from the separately proved balanced controller, not from pretending
that the resource move remains invisible under nonuniform attachments.

## 4. Reproduction and the exact remaining gate

Run

    python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/constructors/six_return.py

The receipt contains all six marked resource phases, the two fixed-root
source states, full literal prepared returns for k=1,2,3, the corresponding
unprepared states, and the k=0 hostile boundary. The source table was also
independently recounted by permutation/cut enumeration of the720 six-vertex
words, without importing the endpoint-DP implementation.

This resolves a genuine interface question raised by the incoming neutral
controller: an **internal** physical mark can turn its constant-target
phase freedom into a useful rooted defect change. What remains missing
for Chan–Pak is a target-directed supply of such prepared sources, or an
operation that performs this repair at the needed arbitrary full count
while retaining its next exact target branch. No such universal selector
is inferred from the H17 family.
