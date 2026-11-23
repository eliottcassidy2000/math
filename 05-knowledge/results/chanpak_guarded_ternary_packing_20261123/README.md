# Guard-restoring ternary moves: a pure unit switch and the golden packing tape

**PROVED scoped / FINITE-EXACT.** This extends the
[five-measurement apex interface](../chanpak_six_five_boundary_atlas_20261123/APEX-TRANSFER.md).
Write \((D,A,B,Q,H)\) for deletion mass, the two scaled positive
endpoint reads, the two-slot deletion read, and the Hamilton count.
Guards are ordinal singleton \(O\) first (apex beats it), \(I\)
last (it beats the apex). The target is a physical recipe that
has a legal input, carries enough state to execute a digit, and
returns to the same ready class. These requirements are separate.

## 1. The exact move that works

[The unit switch](UNIT-SWITCH.md) gives one actual arc reversal
in a six-vertex interior, leaving the apex star and both guards
unchanged:
\[
 (29,898,964,1280,751)\ \longleftrightarrow\
 (29,898,964,1280,753).
\]
The height is the only changed coordinate, not just the only
changed residue. Under any common ordinal singleton collars the
unit difference persists exactly. The marked original guard
vertices can be kept at the two extremes by naming the new
same-contact singleton in the old guard's intermediate
position; the old arcs among the retained vertices are unchanged.
All these hosts are strong.

In particular the guarded height \(h=12693\) is attained. Using
the same interior switch in two explicit contexts gives
\[
 38077\ \longleftrightarrow\ 38079,\qquad
 38079\ \longleftrightarrow\ 38081,\qquad 38079=3h.
\]
Both literal signed exits preserve their own four auxiliary reads,
but the two centers do **not** have the same continuation state.
This is a two-context fork, not three mutually accessible heights
from one chosen center. The switch also gives unbounded rails
of both signs around exact multiples of 3. Independent parent
witnesses are supplied at the stated starts and at 12693; a
uniform witness for every quotient along those rails is not proved.

### Where the earlier 11-rays and direct forks fit

The incoming [direct physical fork](../chanpak_ternary_buckets_20261122/DIRECT-FORK.md)
already supplies a much stronger *same-host* scalar statement on
\(h=21\cdot2^s-13\): three literal arcs from the same center
give \(h,3h-2,3h+2\). Its deletion masses for the parent,
center, lower and upper are \(3,9,5,9\), so the lower exit is
not a pure height switch on the five-read interface with that apex.
The [two 11-rays](../chanpak_ternary_buckets_20261122/ELEVEN-RAYS.md)
likewise give actual arcs \(h\to3h\to3h+2\) at
\(h=11w2^{r-1}+7\) for \(w=3,5\); all-arc scanning of those
prescribed hosts leaves just the two exceptional lower exits at
\(h=117,139\). The 11 there is an affine coefficient of a physical
count chart, not the number of packing states or a golden norm.
Our new control is orthogonal: it gives height-only unit switches
and exact return of the continuation ledger under common contexts,
but not a common actual tripling center with both such switches.

## 2. Why multiplication must be a separate operation

[The composition theorem](PRODUCT.md) gives exact guarded
concatenation, sharing only the apex. In particular it satisfies
\(H(WV)\ge H(W)H(V)\) and
\(H(WV)\equiv H(W)H(V)\pmod4\).
Writing \(C=OI\) for the triangle, the right and left rules are
\[
 H(WC)=3H(W)+4E_O(W),\qquad H(CW)=3H(W)+4E_I(W).
\]
The \(E\)'s count complement-weighted even paths whose endpoints
are both in the indicated contact class. So neither sign
\(\pm2\) can be produced by that move. A stronger obstruction
in the [switch note](UNIT-SWITCH.md) rules out every fixed
mass-three external graft that restores the original kind of
guards while keeping the old word contiguous: its height is
at least \(3H+6D\). Replacements and internal flips are not
covered by this graft obstruction.

There is a physical neutral digit. Inherit the
[strong ternary module upgrade](../session302_ternary_ray_birth_20261113/README.md)
and keep both marked singleton guards. Replace an interior
uniform-contact transitive pair \(TT_2\) by a uniform-contact
directed triangle, with all ordinal contacts to other blocks
unchanged. The modules have \((H,Q_2)=(1,4),(3,12)\), hence
the same normalized middle transfer \(d_+=4\) and \(d_-=0\).
Under uniform contact \(a=d-1,b=0\) or the reverse, and \(c=0\).
The two normalized transfer matrices are identical, so this
replacement multiplies **all five reads** by three. For example,
\(O\,TT_2\,I\) with the pair of type \(O\) goes from
\((1,14,8,16,9)\) to \((3,42,24,48,27)\).

**A finite resource theorem even without using that specific pair.**
No infinite ternary transport can keep applying the homothety
\((D,A,B,Q,H)\mapsto3(D,A,B,Q,H)\). If the deletion has \(n\)
vertices, split each of its \(D\) Hamilton paths at each of its
\(n-1\) cuts. This injects into ordered two-nonempty-slot
covers, so \(Q\ge(n+1)D\) including the two empty-slot choices.
A fixed ratio \(Q/D\) bounds \(n\le\lfloor Q/D\rfloor-1\),
and \(D\le n!\). This forbids indefinite homothetic tripling,
whether or not one can redesign the entire host. An actual
renewable multiplier must let at least this normalized coordinate
grow, or move to a different interface.

If one drops a guard, [the signed endpoint theorem](PRODUCT.md)
replaces the vanished signed sum by a negative ordinary endpoint
count. Then a genuine \(3H+2\) construction works for infinitely
many strong attained inputs. Its deletion family has \(H=105\) at
nine vertices (the guarded theorem is not being applied there):
the apex input at that member has \(H=1255\) and the output has
\(H=3767\). Re-arming the missing singleton without changing
those exact heights remains a different problem.

## 3. What 5, 11 and 105 actually contribute

- **5:** depth-normalize the five outputs \(19,61,63,65,23\).
  The five residues are exactly all of \(\mathbf F_5\), and the
  [packing tape](PACKING.md) implements each centered ternary
  digit by exactly four legal toggles, returns to a singleton,
  and never empties. Its four clock phases are necessary under
  neutral no-motion and bijective labelling. The golden
  eigenvalues exchange modes between phases; after two fair
  digits the nonconstant component is multiplied by \(-1/9\).
- **11:** this tape lives in the eleven independent-set states
  of the five-cycle, the older
  [golden guarded reader](../session91_packing_golden_carriers_20261029/README.md).
  The singleton golden read \(4^i\) gives the squares mod 11;
  the intermediate pairs give the nonsquares; the empty state
  is separate. This is an exact modular controller, not the
  [literal 11-square packing octic](../square11_octic_langlands_collatz_20261006.md).
  The useful import from packing is the requirement to carry
  legality at intermediate moves, not to equate the geometry's
  \(S_8\) field with the quadratic golden field.
- **105:** the old [golden bridge contraction](../chanpak_55_105_golden_20261110/golden/NOTE.md)
  has prime support \(3\cdot5\cdot7\); its raw zero mod 5
  cannot intertwine a bijective five-state digit update unless
  the observation is constant. The normalized remainder is
  needed. The [older 105 source split](../chanpak_55_105_golden_20261110/README.md)
  likewise shows why restoring the scalar alone is not
  restoring ready state. Its targets were 317 and 319 from
  105, i.e. the different digits \(3h+2,3h+4\).
  The packing note gives a separate all-order obstruction
  to both *centered* signs for a fixed marked source in the
  one-reversed-cross-contact architecture.

The next construction needs to combine a **nonhomothetic exact
tripler** with the pure unit switch while producing a center at
the prescribed \(3h\). The finite packing tape can plan phases,
and the signed-count identity makes the boundary error explicit;
neither chooses the absolute center. Chan–Pak remains open.

Run **python3 switch-replay.py**, **python3 product-replay.py** and
**python3 packing-replay.py** for independent literal and finite
controls. The same unit switch was also found by a parallel
search in a differently labelled interior.
