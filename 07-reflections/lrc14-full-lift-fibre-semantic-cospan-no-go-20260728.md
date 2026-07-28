# The full physical lift SCC has no semantic cospan

**Status: FINITE-EXACT WHOLE-OPEN-CYLINDER NO-GO + INDEPENDENTLY
HOSTILE-AUDITED.**  This note refines THM-2707 on its exact fixed packet
skeleton.  It excludes no scalar row and does not prove LRC(14).

## 1. The tempting composition

THM-2707 produces an unusually large positive object above the THM-2694
terminal vertex:

```text
3346 packet addresses on one common open base interval I;
eleven active residue parts;
complete directed eleven-partite physical-lift graph;
one recurrent SCC with 10,177,910 edges.
```

The obvious next test is not another support cycle.  THM-2305 requires a
literal prescribed-clock cospan

```text
E_j at q  ->  Q_(j,sigma) at D^(2j)q,        D(q)={13q}.  (1)
```

If two residue parts of the THM-2707 simplex contained positive instances of
`(1)`, its complete multipartite graph would provide a recurrent semantic-
source subgraph.  The exact answer is no, and the reason is a one-line
congruence rather than a sparse search failure.

## 2. Exhaustive whole-cylinder classification

Write the THM-2707 packet addresses as

```text
q_n(x')={13x'+7n/13^6},                 x' in I.          (2)
```

The audit reconstructs all `3346` addresses and partitions every open
physical cylinder `q_n(I)` at every source and prescribed-endpoint danger
wall.  No cylinder contains an internal semantic wall.  Hence its midpoint
signature is its whole-open-cylinder signature.

The exact source census is

```text
exclusive E_1:                         0;
exclusive E_2:                       265;
exclusive E_3:                         0;
A_0 with zero dangerous blockers:   3081.                (3)
```

The `265` source cylinders occur in every active residue part: `24` in each
part except residue `9`, which has `25`.  They split as

```text
241 nonzero-lift targets + 24 same-residue vertices.      (4)
```

Thus excluding only the old terminal residue would not repair the semantic
failure.

## 3. The address-independent failed bit

Every source in `(3)` is `E_2`, whose prescribed clock is four.  Since
`c_2=13^3`, equation `(2)` gives

```text
c_2 D^4(q_n)
 =13^7(z+7n/13^6)
 =13^7 z+91n                         (mod 1).              (5)
```

The address term is integral.  The endpoint source-blocker phase is therefore
the same on all eleven residue parts.  Exactly,

```text
{13^7 z}=1598741/1599416,
centered phase=-675/1599416.                              (6)
```

The effective radius of the full open cylinder after the `13^7` multiplier
is `1/1599416`.  Its remaining strict danger-band slack is

```text
1/14-675/1599416-1/1599416=12/169>0.                     (7)
```

Consequently the old source blocker `c_2` remains dangerous throughout every
one of the `265` prescribed endpoints.  But every terminal word
`Q_(2,sigma)` requires `c_2` to be safe.  This is the first universal failed
predicate.

The boundary is sharp.  Exactly `16` nonzero-lift cylinders satisfy every
other predicate of the pure word `Q_(2,{1})`: endpoint guard and five ordinary
speeds are safe, `c_1` is dangerous, and `c_3` is safe.  Only the single
`c_2`-safe bit fails.  The no-go is therefore not caused by diffuse packet
damage.

## 4. What the result changes

The obstruction separates two objects that had looked nearly composable:

```text
THM-2707:
  abundant physical support recurrence
  + common interval
  + carry/root/unit sidecars
  - no prescribed source-to-word cospan;

half/C221 full source fibres:
  abundant strict E_3 -> Q_(3,{1,2}) cospans
  - attachment to the previously declared reciprocal packet cycles.
```

The second line is audited in the companion source-fibre census: for the half
carrier, the dynamic present packet removes the `E_3` rows; in the selected
`C_221` subgenerator, the forced carry/root slice removes them.  The two
failures are complementary.  The old large SCC has the sidecars but not the
semantic cospan; the unrestricted transverse fibres have the cospan but not
the sidecar-compatible reciprocal cycle.

The next constructive operations are therefore tightly typed:

1. re-root the half carrier's present packet at an `E_3` source;
2. enlarge the `C_221` carrier from forced carry/root `6` to an all-carry,
   variable-root microphase transport; or
3. build a new coloured transit/private-row grammar, rather than attempting
   to terminalize the unchanged THM-2707 skeleton.

This audit is distinct from reserved THM-2712.  THM-2712 targets the frozen
following-atom semantic support and its same-residue congruence lock.  The
present result tests the literal THM-2305 cospan on every original THM-2707
packet cylinder.

## 5. Reproduction and scope

Run

```bash
python3 04-computation/lrc14_full_lift_fibre_semantic_cospan_no_go_20260728.py
python3 -O 04-computation/lrc14_full_lift_fibre_semantic_cospan_no_go_20260728.py
```

Both modes byte-match
`05-knowledge/results/lrc14_full_lift_fibre_semantic_cospan_no_go_20260728.out`.
The companion contains no Python `assert` nodes and explicitly partitions
every open cylinder at all semantic walls.

SHA-256:

```text
script  e0f7de84b7d294e2e86ea4cdd6b1ba746da260467e4f5e5a28dc238299d0acd0
output  e1d32978ed175d0c9600fe926100c1952817f5f88083a4e23f89e82858b497b0
```

The universe is exactly the THM-2707 fixed packet skeleton over its inherited
common interval.  No other packet grammar, changed present bank, variable
root carrier, semantic endpoint current, scalar-row exclusion, or LRC(14)
conclusion is claimed.
