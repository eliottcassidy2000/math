# Six three-cusp scout words admit a nonabelian braid-group quotient

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED for the literal
group presentation; underlying geometric paths remain HEURISTIC.** The exact statement below concerns
the six literal words, not an already certified affine complement.
Agreement of the numerical scouts at 256 and 512 subdivisions is not
a rational path certificate.

## 1. Exact input and the first failed implication

Use chronological Hurwitz action

    H_i(...,u,v,...)=(...,uvu^-1,u,...),
    H_i^-1(...,u,v,...)=(...,v,v^-1uv,...).

The six input words on `(a,b,c,d)` are

| Name | Chronological word |
| --- | --- |
| cusp_plus | `[-1,-2,1,1,1,2,1]` |
| cusp_minus | `[2,1,3,3,3,-1,-2]` |
| cusp_two | `[2,1,1,1,-2]` |
| node0 | `[2,-1,2,2,1,-2]` |
| node1 | `[2,1,2,3,3,-2,-1,-2]` |
| node2 | `[2,3,2,2,-3,-2]` |

Let `T` be the group on `a,b,c,d` with relations saying that each
of these six automorphisms fixes the entire generator tuple. Then
**`T` surjects onto `B_3=<p,q | pqp=qpq>`**. In particular fixedness
does not force `a=b=c=d`, and does not force cyclicity. This is the
decisive stopping witness for that proposed group-theoretic shortcut.

The closest mechanism is the successful five-word arbitrary-group
elimination used by the infinity-eleven lane. The changed object is
the full common-access word system. Its source and target are exact
group presentations; the map is the quotient substitution in
Section 3. It preserves every stated relation and generator image.
It loses the actual geometric path certification, any additional
relations of the affine complement, and all retained/deleted-sheet
data needed by a Keller passport. None of those is reconstructed
from this quotient.

## 2. Complete prefix calculation

Write `e=bcb^-1`. Fixedness of a Hurwitz square on a pair is
equivalent to commutation; fixedness of its cube is equivalent
to the braid relation on that pair. A word of the form
`P,core,P^-1` fixes the original tuple exactly when its core fixes
the tuple after applying the chronological prefix `P`.

The relevant prefix pairs are

| Word | Prefix | Core pair | Exact relation |
| --- | --- | --- | --- |
| cusp_plus | `[-1,-2]` | `(b,c)` | `bcb=cbc` |
| cusp_minus | `[2,1]` | `(b,d)` | `bdb=dbd` |
| cusp_two | `[2]` | `(a,e)` | `aea=eae` |
| node0 | `[2,-1]` | `(e^-1ae,b)` | `[e^-1ae,b]=1` |
| node1 | `[2,1,2]` | `(a,d)` | `[a,d]=1` |
| node2 | `[2,3]` | `(e,bdb^-1)` | `[c,d]=1` |

Every displayed inverse is the literal reversed prefix with signs
changed. No braid reordering or choice of peripheral transporter
is being assumed. In the last row the commutator is conjugated
by `b^-1` to give the relation shown.

## 3. The exact infinite nonabelian quotient

In `B_3` set

    a=c=d=q,             b=p.

The first two cusp relations are the defining braid relation.
Put `e=pqp^-1`. The same relation gives

    e=q^-1pq,             e^-1 q e=p.

For the second identity, the braid relation also gives
`p^-1qp=qpq^-1`; hence
`e^-1qe=pq^-1(p^-1qp)qp^-1=p`.

Conjugating the braid relation on `(q,p)` by `q^-1` proves the
relation on `(q,e)`, paying the third cusp. The first node becomes
`[p,p]=1`; the other two become `[q,q]=1`. All six constraints
therefore hold. The image contains both `p` and `q`, so this is
a surjection, not just a nontrivial representation of one relator.

For a separate literal nonabelian and infinite control, use

    p=[[1,1],[0,1]],       q=[[1,0],[-1,1]]

in `SL_2(Z)`. They satisfy `pqp=qpq`, do not commute, and
`p^n=[[1,n],[0,1]]` has infinite order. The source checks every
full input word directly on this matrix tuple, independently of
the prefix reduction.

There is also a transitive `S_3` quotient with `p=(12)`, `q=(23)`.
This does not assert a degree-three Keller cover; the actual
fibre counts and classical low-degree exclusions are separate
conditions.

The substitution `a=c=d` itself is not forced by the words.
Another exact control in `S_4` is

    a=d=(34),       b=(23),       c=(12).

Here `e=(13)`, and the first node compares the disjoint
transpositions `(14)` and `(23)`. The other node relations also
hold, while `a!=c`. Adjacent transpositions generate a transitive
`S_4`. Thus the proof claims a quotient of `T`, not a complete
identification of its presentation with `B_3`.

This same control refutes a more subtle proposed replacement for
cyclicity: generation of the global group by one of the three
declared cusp subgroups. Their literal access pairs are
`((23),(12))`, `((23),(34))`, and `((34),(13))` in this quotient.
Each generates an `S_3` fixing one label, while all four generators
generate `S_4`. Hence no one of those cusp pairs generates `T`.
An actual-complement use of cusp-local generation would require
additional geometric relations beyond this word bank.

## 4. Verification and precise stopping point

The standalone exact source checks every prefix decomposition,
the full six Hurwitz actions in both permutation controls and
the integer matrix control, and the infinite-order formula.
No finite group census is used to infer an arbitrary-group
statement.

The next geometric step must retain exact shared access paths and
the actual sheet passport, or find genuinely additional relations.
Certifying these six words alone cannot justify the rejected
cyclicity conclusion. No claim about the actual three-cusp
complement, a full presentation for it, or a Keller realization
has been promoted here.

Reproduction, from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_three_cusp_group.py
python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_group.py
```

Both runs and the frozen output agree byte for byte: **115 always-active
gates**, 360 output bytes. The finite universe comprises the six named
words, their six literal prefix decompositions, the displayed matrix
tuple and its powers `1..20`, the named `S_3` and `S_4` tuples, and the
three declared cusp subgroups in the `S_4` control. The proof of infinite
order is the symbolic matrix-power identity, not extrapolation from the
twenty controls. Group closures have orders `6`, `24`, and `6` for each
of the three proper cusp subgroups, respectively.

- [Exact source](../../04-computation/planar_jc48_sep08_three_cusp_group.py):
  SHA256 `0b38b41600521ac417d23ad6399ff2c6fedc40e75b35030a8fed67a14fa84f4b`.
- [Frozen output](planar_jc48_sep08_three_cusp_group.out):
  SHA256 `cbb8780d4447aa68a30c374392d4ac22a49a9e1c2deb12fa003f14481f0583d2`.
- Semantic digest:
  `1beec165faf7fcfc193ac52e673375e0816f54f2592913a5091eb1c3c4c92f64`.

The source and output are frozen; independent analytic review of the
group theorem is pending. The geometric status remains **HEURISTIC**
regardless of that algebraic audit.

The [independent root audit](planar_jc48_sep08_three_cusp_group_audit.md)
passes the full group/source argument and both115-gate frozen replays.
This promotion does not certify the underlying geometric paths.
