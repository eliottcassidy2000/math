# LRC Pressure DAGs S491

The question was whether "pressure searches returning DAGs" is a concept worth
separating.  Yes: in the current LRC pressure setup, a DAG is not merely the
absence of a cyclic obstruction.  It is an ordered dependency certificate.

## Pressure Relation

For runners `{0} union V` at time `t`, define deletion relief:

```text
relief_i(j) = nearest_distance_i(after deleting j) - nearest_distance_i.
```

For an unordered pair `{i,j}`, orient:

```text
j -> i  if  relief_i(j) > relief_j(i).
```

This means `j` is the more irreplaceable blocker of `i`.  Ties stay missing in
the strict pressure graph.  If a complete tournament is needed, the repo's base
Hamiltonian path supplies the tie completion, but the strict graph is the right
object for peeling.

## DAG Versus SCC

The old phrasing was: if pressure has no SCC, the row is pressure-peelable.
S491 sharpens that.

A strict pressure SCC is the disproof-like signal.  It says every runner in
the component blocks another runner in a closed dependency loop.  That loop
would need arithmetic labels and should be killed using the THM-365 endpoint
cycle style.

A strict pressure DAG says the dependence can be topologically sorted.  That
sort is a proof object:

```text
source layers: runners whose blocking pressure is upstream
sink layers:   runners that can be peeled from the dependency order
```

The proof target is to pair these layers with endpoint-private rows.

## Computation

`lrc_pressure_dag_s491.py` audited:

```text
n14 initial
n14 d=7
n14 d=14
n18 initial
n18 d=3
n18 d=9
n18 d=18
```

Across bounded pressure search windows, every sampled pressure graph was a
DAG:

```text
case          times  cyclic  max_scc  max_tri  max_chain
n14 initial      42       0        1        0          1
n14 d=7         113       0        1        0          3
n14 d=14        183       0        1        0          4
n18 initial      40       0        1        0          1
n18 d=3         100       0        1        0          3
n18 d=9         233       0        1        0          3
n18 d=18        425       0        1        0          4
```

The source/sink layers are concrete.  For example:

```text
n14 d=7 at t=3053/25872
source: {1,14,49} {7,35,56,77,91} {0,84}
sink:   {0,35,56,84,91} {7,14,77} {1,49}

n18 d=18 at t=8681/114048
source: {1,36,90,162,270} {18,54,126,180,198,288,306} {0}
sink:   {0,18,54,126,180,288,306} {36,90,162,198,270} {1}
```

These layers are small, structured, and exact.  They should be stored in future
pressure searches alongside gap/debt numbers.

## Interpretation

The current hard LRC rows keep falling into the DAG regime:

```text
n14 d=7:  gap/th=5/924,  unprotected=84,  product=5/11
n14 d=14: gap/th=5/1848, unprotected=168, product=5/11
n18 d=9:  gap/th=1/176,  unprotected=176, product=1
n18 d=18: gap/th=1/352,  unprotected=352, product=1
```

Thus the main signal is not "small gap plus no cycle."  It is:

```text
small gap + exported endpoint debt + pressure DAG.
```

That package points toward proof, not disproof.  The counterexample-like
object would need to defeat endpoint-private peeling and also produce a
labelled pressure SCC.

## Next Work

1. Add pressure DAG layer vectors to the endpoint/IP ledger.
2. Match sink layers against endpoint-private rows: a sink with a private
   endpoint should be an immediate branch deletion.
3. Search bounded perturbations specifically for `pressure_largest_scc > 1`;
   scalar gap alone should not drive the search priority.
4. If a pressure SCC appears, label each arc by `(src_speed,dst_speed,t,margin)`
   and compare it to THM-365 endpoint-cycle arithmetic.
