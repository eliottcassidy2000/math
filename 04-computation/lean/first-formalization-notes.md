# First Lean Formalization Notes

Session focus: get oriented enough to choose the first Lean target without
overbuilding a mathlib project.

## Recommended First Target

Start with the local vocabulary needed for THM-001 (Redei's theorem):

- A finite tournament on a type `V` as an orientation predicate
  `edge : V -> V -> Prop`.
- Irreflexive edge relation: no `edge v v`.
- Completeness/antisymmetry: for distinct `u v`, exactly one of
  `edge u v` and `edge v u`.
- A directed Hamiltonian path as a list with no duplicates, covering all
  vertices, and satisfying consecutive edge constraints.

This is intentionally below the ambition level of proving Redei. The first
useful Lean milestone is: "the transitive tournament on a finite linear order
has exactly one increasing Hamiltonian path" or, even smaller, "a one-vertex
tournament has one Hamiltonian path."

## Why Not Start With OCF

THM-002 depends on conflict graphs of all directed odd cycles and an
independence polynomial. That is a better second-stage formalization once the
basic tournament/path vocabulary is stable.

## Suggested Lean File Shape

Use a scratch file first, not a tracked theorem file:

```lean
structure Tournament (V : Type u) where
  edge : V -> V -> Prop
  irrefl : forall v, not (edge v v)
  total : forall u v, u <> v -> edge u v \/ edge v u
  antisymm : forall u v, edge u v -> not (edge v u)
```

Once mathlib is introduced, prefer `Finset`, `List.Nodup`, and `Fintype` rather
than inventing custom finite-set machinery.

## Environment Status

The Windows machine has `elan` installed under `%USERPROFILE%\.elan` with Lean
4.30.0 and Lake available. The first toolchain archive is large (516 MB on
Windows), so the bootstrap uses `aria2c` when available. Future cloud
environments should use the Linux bootstrap path where persisted `/root/.elan`
and `/root/.cache` avoid repeated downloads.

The scratch file `.scratch/lean/tournament_smoke.lean` compiled successfully
with:

```powershell
lean .scratch\lean\tournament_smoke.lean
```

Small Lean lesson from the first pass: use `u ≠ v` or `Ne u v` for inequality;
`<>` is not Lean syntax.

## First Lemmas Proved in Scratch

The scratch file now proves three tiny foundation lemmas:

- `tournament_no_loop`: no vertex has an edge to itself.
- `tournament_oriented_pair`: distinct vertices have at least one directed edge.
- `tournament_no_reverse`: an edge `u -> v` rules out `v -> u`.

It also constructs the two-vertex tournament with `left -> right` and proves:

- `two_left_beats_right`
- `two_right_does_not_beat_left`

## Suggested Next Lean Session

Define a directed path as a `List V` whose consecutive pairs satisfy `T.edge`.
Then prove the base cases for Hamiltonian path counting:

1. A one-vertex tournament has exactly one Hamiltonian path.
2. The two-vertex tournament above has exactly one Hamiltonian path.
3. Reversal sends Hamiltonian paths of `T` to Hamiltonian paths of `T^op`.

That third item is the small formal bridge used in THM-002's proof chain:
`ham(T^op) = ham(T)`.
