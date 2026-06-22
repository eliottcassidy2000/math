# Goddyn-Wong sporadic tilers as accelerated AP tails

S106 started from the new THM-560 signal: difference-closed exact tilers are
exactly AP dilates.  That makes the Goddyn-Wong row necessarily
non-difference-closed, but not unstructured.

The right structural lens is Tao's restatement of the Goddyn-Wong construction:
start with `{1,...,n}` and replace a velocity `v` by `2v` when every integer in
`[n-v+1, 2n-2v+1]` has a common factor with `v`.  This is a Jacobsthal-style
nonunit interval.  For LRC14, `n=13`, `v=12`, and `[2,3]` is nonunit modulo
`12`, giving `{1,...,11,13,24}`.

This explains the repeated tail pattern:

```text
n=7:  {1,2,3,4,5,7,12}
n=13: {1,...,11,13,24}
n=19: {1,...,17,19,36}
...
```

The infinite subfamily is `n ≡ 1 mod 6`, with `v=n-1`.  The broader family is
formed by tail speeds whose obstruction window is entirely non-coprime modulo
`v`; exact audit also verifies the multi-acceleration row
`n=73`, tail `...,69,71,73,140,144`.

The important LRC14 interpretation: `24` is not a first-order blocker at the
isolated grid witnesses.  The q-grid blocker pairs for `{1..11,13,24}` are
still AP-like.  The accelerated speed repairs the petals opened by deleting
`12`.  So Goddyn-Wong is a second-order petal repair, not a new equal-spacing
tiler.

Proof consequence: classify exact tilers by two branches before invoking
generic cap/residual machinery:

```text
1. difference-closed / self-similar branch -> THM-560 -> AP dilates;
2. accelerated-tail branch -> Jacobsthal nonunit windows -> Goddyn-Wong atoms.
```

Everything outside these branches should have a positive safe interval.  For
LRC14 specifically, this suggests treating `{1..11,13,24}` as a named finite
atom in the low-depth AP/Freiman atlas feeding HYP-2890, not as an arbitrary
wide counterpressure row.

External sources checked:

- https://terrytao.wordpress.com/2017/01/10/some-remarks-on-the-lonely-runner-conjecture/
- https://arxiv.org/pdf/2409.20160
