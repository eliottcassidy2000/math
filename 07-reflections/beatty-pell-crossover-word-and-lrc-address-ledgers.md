# Beatty-Pell Crossover Word And LRC Address Ledgers

codex-2026-06-13

The useful correction from this session is that "Sturmian/Beatty-like" was
too blurry.  The triangular crossover has three layers:

```text
1. Beatty address:       d_m = floor((sqrt(2)-1)m + sqrt(2)/4)
2. Pell carry:           r_m = m^2 - 2m*d_m - d_m^2
3. Visible side token:   state(B_m.L) state(B_m.R)
```

The first-difference word of `d_m` is genuinely Sturmian.  The visible token
word is not.  It is a six-interval rotation coding with two zero-density Pell
atoms.  That is a tiny but important distinction: the proof-bearing coordinate
is one layer below the thing we first see.

## The Pattern

The old floor-sqrt classifier said:

```text
given a B side [u,v], compare it to floor(sqrt(u)), the A midpoint, and the
square boundary.
```

The new normal form says the same thing without ad hoc interval comparison:

```text
address d  -> where the triangular row sits among square shells
carry r    -> how far the row is from the midpoint/boundary walls
token      -> the visible coarse state after forgetting d and r
```

This is exactly the repo's recurring repair grammar:

```text
scalar/product collision -> attach address coordinate -> recover the word.
```

The Pell atoms are the moments where a coarse inclusion statement tightens into
an endpoint equality.  They are not noise; they are the wall labels.

## LRC Transfer

For LRC14, the analogous mistake would be to call a denominator "blocked" and
stop.  HYP-2241, HYP-2443, and HYP-2444 already show that this mixes rows:
AP/Vstar/2AP are endpoint wall atoms, but neighboring carry rows open.

A better LRC14 row should remember something like:

```text
(q, shell class, unit quotient class, divisor fiber, owner support,
 carry residue, endpoint atom, deletion/opening target).
```

The triangular toy says why this is not bookkeeping excess.  Without `r_m`,
the crossover word looks like a vague Beatty phenomenon.  With `r_m`, each
side state is an exact inequality and the rare Pell walls are classified.

## Practical Next Computation

Build a Q27 analogue of the script:

```text
for each candidate row S and each q in Q27:
  record blocked twists
  record minimum owner support tau_q(S)
  record shell-27 unit class
  record divisor fiber d|14
  record carry/opening/deletion class from HYP-2241/HYP-2443
```

Then factor the visible status word:

```text
strict/open/wall token
owner-private deletion token
fibered q=91 rescue token
pure shell failure token
```

The target is not a pretty slope like `sqrt(2)`.  The target is the same
normal-form move: find the address coordinate where the mixed projection
becomes an exact finite set of inequalities plus rare boundary atoms.

## Caution

This triangular model is exact, but it does not prove LRC14.  Its value is as a
toy that prevents a bad slogan.  "The word is Sturmian" is false for the visible
two-side token.  "A hidden address increment is Sturmian, and the visible word
is its carry-decorated projection" is the proof-friendly version.
