---
id: THM-445
name: lrc-momentum-twistor-is-the-discrete-log
status: PROVED (prime shell) + VERIFIED (exact, n=4..14; linearization all shells incl. ramified 27)
date: 2026-06-08
session: claudebox-2026-06-08-S718
depends_on:
  - HYP-2370  # LRC dual conformal symmetry (S717): multiplier (Z/m)*, inversion = special conformal
  - THM-420   # transversal core: bad multipliers = {+-v_i^-1}; non-transversal => dodge
  - THM-415   # signed-LRC half-system = homometry
  - THM-403   # cyclotomic witness orbit
provisional_id: true   # THM counter contested; renumber at PR
---

# THM-445: The LRC momentum twistor is the discrete logarithm; it linearizes dual conformal symmetry

Hodges' momentum twistors make the (hidden) dual conformal symmetry of N=4 amplitudes LINEAR. The LRC
dual conformal symmetry (HYP-2370) is the MULTIPLIER group `(Z/m)*` acting by `v -> a v` on residues at
shell `m`, with the INVERSION `v -> v^{-1}` as special-conformal generator. A multiplicative action is
linearized by a logarithm, so the LRC momentum twistor is the **discrete log**.

## The twistor map

Fix a primitive root `g` mod `m` (exists for `m` prime or a prime power). Define
```
   ell : (Z/m)* -> Z/phi(m),     v |-> dlog_g(v).
```
Then (VERIFIED exactly for `m = 7,11,13,19,23,27`):
- the **multiplier** `v -> a v` becomes a **translation** `ell -> ell + ell(a)`;
- the **inversion** (special conformal) `v -> v^{-1}` becomes **negation** `ell -> -ell`;
- `ell(-1) = phi(m)/2 =: c` (the unique order-2 element), so `v -> -v` is the shift `ell -> ell + c`.

The hidden dual conformal symmetry is now linear (additive). This is the LRC analog of momentum twistors.

## The dodge criterion in twistor coordinates (EXACT, all shells)

For a config with unit residues `R = {v_i}` and log-set `L = ell(R)`, a multiplier `a` is BAD for runner
`i` iff `a v_i in {+1,-1}`, i.e. `ell(a) in {-ell(v_i), c - ell(v_i)}`. So the bad-multiplier log-set is
```
   B = (-L) cup (c - L)   subset Z/phi(m),     a union of TWO TRANSLATES of -L.
```
A good multiplier (THM-420 dodge) exists **iff `B != Z/phi(m)`**, i.e. `|B| < phi(m)`. VERIFIED exact
against direct search: `600/600` at every shell `m = 7,11,13,19,23,27` (including the ramified `27`). The
multiplicative dodge problem becomes **additive covering on the cyclic group `Z/phi(m)`**.

## The transversal core = a half-system tiling, EXACTLY when `2n-1` is prime (PROVED)

`|B| = phi(m)` (no dodge) iff `(-L)` and `(c-L)` are disjoint, i.e. `L` contains no pair `{l, l+c}`. In
residues `l' = l + c <=> v' = v g^c = -v`, so:
```
   no dodge  <=>  R has no +-pair {v,-v}  <=>  R is a +-transversal      (recovers THM-420),
```
and, since `|L| = n-1` and there are `phi(m)/2` pairs `{x, x+c}`:
```
   no dodge  <=>  L is a HALF-SYSTEM: a transversal of {0,c}, one element of each pair  <=> {-L, c-L} TILES Z/phi(m).
```
This tiling characterization requires `|L| = n-1 = phi(m)/2`, i.e. `phi(2n-1) = 2n-2`, which holds **iff
`2n-1` is prime**. VERIFIED: "no-dodge <=> half-system" holds `600/600` for prime shells
`m = 7,11,13,19,23`, and **breaks** (`325/600`) at the ramified shell `m = 27` (`n=14`), where
`phi(27) = 18 < 26 = 2(n-1)`.

## The twistor explains the prime/ramified dichotomy and n=14's two heads

- **`2n-1` prime:** `phi = 2(n-1)`, sizes match, the transversal core = unit half-systems; THM-420 is the
  clean statement. The frontier `n = 15,19,21,22` (with `2n-1` prime) lives here.
- **`2n-1` ramified (e.g. `27 = 3^3`, `n=14`):** `phi(2n-1) < 2(n-1)`, no unit half-system exists, so the
  unit residues ALWAYS dodge at shell `m`; the difficulty migrates to the **non-unit residues**
  (`3 | v`), which are OFF the twistor `(Z/m)*` — exactly the S643 ramified-fiber / divisor story
  (S710). The momentum twistor cleanly separates the dodgeable unit part from the ramified off-twistor
  part — the two heads of `n=14`.

## Twistor dimension

`(Z/m)*` is cyclic (1D "P^1-like" twistor) iff `m in {p^k, 2p^k}`. For `m = 2n-1` with several prime
factors the twistor is HIGHER-dimensional (one log coordinate per cyclic factor) — a structural measure
of arithmetic complexity of the shell.

## Significance

The discrete-log twistor (i) makes the LRC dual conformal symmetry manifest and linear, (ii) converts the
THM-420 dodge into additive covering / tiling on `Z/phi(m)`, (iii) recovers THM-420 and reveals the
transversal core as a half-system (linking THM-415), and (iv) turns the prime/ramified dichotomy into the
size-matching `phi(2n-1) = 2n-2`. It does not close the conjecture (the half-system core stays open at
prime shells), but it gives the cleanest coordinate for it and pins n=14's hardness to off-twistor
(ramified) residues.
