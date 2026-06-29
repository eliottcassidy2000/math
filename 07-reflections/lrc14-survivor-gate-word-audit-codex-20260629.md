# LRC14 Survivor-Gate Word Audit Reflection

HYP-3438 executes the HYP-3436 local-to-global hook.  The useful object is not
the total survivor measure and not a raw count of mixed even-safe components.
It is the local gate word inside each mixed component:

```text
bad-core block / survivor gap / bad-core block
+ endpoint labels
+ adjacent minimal B0/B1 owner covers
+ parent even wall
+ branch mask
+ legal sidecar route
```

The exact audit on the same `135` primitive covering rows as HYP-3436 finds
`6228` mixed even-safe components and `8702` survivor gates inside them.  The
branch masks split symmetrically:

```text
both=1064, branch0=3819, branch1=3819
```

Most gates are edge gates (`3515` left-bad edge and `3515` right-bad edge);
`1672` are two-sided.  Every parent endpoint pair is `E|E`, so the survivor
gate is always inside a finite even-wall chamber.  The tight AP-with-84 row is
especially clean: its four gates are one-gap edge words with endpoint labels
`B1:7|E:84`, `E:84|B1:5`, `B0:5|E:84`, and `E:84|B0:7`.

## Proof Pull

The next lemma should be a gate-word theorem.  A survivor-mass argument is
legal only if it reconstructs or routes:

```text
adjacent bad-core covers
endpoint labels
branch mask
parent even wall
```

The exits are finite: direct branch relocation, corridor-fence,
endpoint-spine/wall certificate, owner-current exception, two-adic loss,
overlap-cut bridge, or signed-SPEC.  Raw survivor measure and harmonic budgets
remain negative controls.

## Assumption Challenge

Runners, gaps, fixed sections, section boundaries, wall-crossing events,
residues, cover arcs, endpoint walls, branch-bad owners, cover-pair deltas,
mixed even components, survivor gaps, owner-current labels, and proof
obligations all remain candidate vertices.  The selected quotient preserves
the two-adic relocation predicate while making explicit what scalarization
destroys.
