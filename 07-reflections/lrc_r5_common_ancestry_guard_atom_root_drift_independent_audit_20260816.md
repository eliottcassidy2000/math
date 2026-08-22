# Independent audit of the common-ancestry guard-atom/root-drift package

**Status: ACCEPT, SCOPED FINITE-EXACT.**  The candidate correctly realizes
the actual 39 guard atoms as a lawful refinement of the THM-2471/2594
common-base Boolean service.  The common gauge has exactly `48/52` active
support, and retaining its absolute offset gives full drift-Fourier support
for every nonzero owner character.  This is a support/spectral-capacity
result only.  It proves no equality with the frozen U_full endpoint weights,
physical current, grouped coefficient, scalar-row exclusion, or LRC(14).

The independent verifier is
`04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_independent_audit_20260816.py`,
with matching output
`05-knowledge/results/lrc_r5_common_ancestry_guard_atom_root_drift_independent_audit_20260816.out`.
It hash-pins but never imports or edits the candidate implementation.

## 1. Inheritance and exact object

The closest proved mechanism is THM-2471 equation (31): any rational Boolean
partition may be evaluated separately on the linked current node
`X_(u,a)` and source node `Y_(q,e')` before finite-sheet marginalization.
THM-2594 supplies the proved `r=5` instance and its structural `u=q` zero.
THM-3514 supplies the rational 39-atom partition

```text
I_(a,C),  a in F_13,
C in {L=[0,1), M=[1,6), R=[6,7)}.
```

THM-3518 supplies the owner-phase warning: the absolute left label must be
retained before character contraction.  It does not identify ancestry
masses with endpoint response weights.

The canonical hostile is the opposite drift orientation `s=d`, which has
the same support count.  Therefore support cannot choose the sign; the typed
labels and owner phase remain load-bearing.

## 2. Reconstruction independent of the candidate buckets

The verifier imports only the hash-pinned THM-2594 stage constructor.  It
then independently:

1. rebuilds `E`, `Q`, and `f=1_Q P^2 1_E`;
2. intersects their weighted interval pieces with the contiguous 39-atom
   partition using a forward two-pointer sweep;
3. multiplies the atom indicators into `f` and `E` **before** applying
   `P_(13^5)`;
4. verifies pointwise recombination of all 39 folded atom profiles;
5. extracts all thirteen collision-root windows; and
6. integrates every atom/root pair with a separately written step-profile
   product integrator.

This reproduces the candidate's independent pre-fold census:

```text
input pieces (f,E):       120234, 57072
nonempty atom groups:     21, 20
active atom pairs:        420/1521.
```

The current-origin revision also exports the unmarginalized
`39 x 39 x 13` atom-pair/common-offset tensor for a downstream experiment.
The audit constructs that tensor inside its own atom/root loop and reproduces
the candidate's serialization digest exactly:

```text
e4d0d4fa674e1f54496e613f7a3e1f057af033fa8992322a5f414ea176e1c3d4.
```

Because the split occurs before `P_(13^5)`, the labels are exactly
`P_omega(X_(u,a))` and `P_nu(Y_(q,e'))`.  The audit never identifies `X`
and `Y` as one circle point and never inserts an endpoint diagonal.

## 3. Gauge, mass, and support laws

With address drift `d=b-a`, root drift `s=u-q`, and offsets

```text
c_L=a-u,  c_R=b-q,
```

the verifier exhausts all `13^4=28,561` label quadruples and checks

```text
c_L=c_R  iff  a-u=b-q  iff  s=-d.
```

The offset predicate is imposed before root, sheet, or atom marginalization.
Summing the refined tensor reproduces the whole-table root marginal exactly.
The result has

```text
root support:        12/13
s=0 mass:            0, atom by atom
total service:       48602521488933856/1996669193081744181
                   = 169 I_5.
```

The common-gauge support is `72/117`.  On the active `{L,R}^2` carrier it
is exactly

```text
{L,R}^2 x F_13^*,  hence 48/52.
```

All four missing cells have `d=0`, which forces `s=0`.  The opposite
orientation `s=d` also has `72/117` support, confirming that the sign is a
typing statement rather than a census consequence.  All thirteen gauge
defect classes occur, so the common-gauge slice is not the whole tensor.

## 4. Exact K4 and zero-owner spectrum

The four integer rows `LL,LR,RL,RR` have rational rank four.  Their ranks
are also four modulo each of

```text
547, 911, 1093, 2003, 2549.
```

Each row and each Walsh row has raw support `12/13`.  The exact zero-mode
ledger is

```text
W_0 !=0,  W_L=0,  W_R=0,  W_X!=0.
```

For primitive drift frequencies, the independent verifier does not rely on
a lucky prime.  It collects the integer coefficients of the cyclic
polynomial `P(x)` in `Z[x]/(x^13-1)`.  At a primitive thirteenth root,
`P` vanishes exactly when its thirteen cyclic coefficients are all equal,
that is, when its representative is a multiple of `Phi_13`.  This gives

```text
k=0 Walsh drift support = (13,12,12,13)
```

exactly over `Q(zeta_13)`.

## 5. Retaining the common offset

The audit retains

```text
c=a-u=b-q
```

inside each active chamber/drift bucket.  For every owner character `k` and
drift character `h`, it forms the exact integer cyclic polynomial

```text
sum_(c,d) W(c,d) x^(-kc-hd).
```

The same `Phi_13` criterion decides all `13*4*13=676` coefficients.  The
complete result is

```text
k=0:      (13,12,12,13)
every k!=0: (13,13,13,13).
```

The five split primes independently reproduce this decision.  Writing a
five-bit mask for the primes at which a coordinate survives gives

```text
mask 0:   2 coordinates   (the exact k=0 Haar zero modes)
mask 30:  4 coordinates
mask 31:  670 coordinates.
```

Thus every exact nonzero coordinate survives at least one certified split
prime, while the only zero masks are the two proved zeros.

For each owner frequency, the verifier also chooses a rank-four minor in a
split field, recomputes its determinant as a polynomial in
`Z[x]/(x^13-1)`, and checks that the determinant is not a multiple of
`Phi_13`.  All thirteen owner-frequency chamber matrices therefore have
exact rank four over `Q(zeta_13)`, not merely rank four in one reduction.

## 6. Sign and denominator audit

- Replacing `c` by `-c` permutes nonzero owner characters; replacing `d` by
  `-d` permutes drift characters.  Full-support and rank claims are invariant,
  while the typed equality remains `s=-d` in the declared labels.
- The service has one common denominator
  `13^2 * 13^10 * T_DEN`.  It is a unit at all five split primes.  Support
  and rank are computed on integer numerators, so no denominator was silently
  dropped coordinatewise.
- A nonzero split reduction is only used as a certificate of nonvanishing.
  Exact cyclotomic reduction independently prevents a false conclusion from
  a bad embedding or an unlucky prime.

## 7. What is and is not closed

The package closes a genuine former capacity question:

```text
lawful common-base guard ancestry exists;
the common gauge realizes all 48 primitive active buckets;
its retained offset restores every missing nonzero-owner spectral mode.
```

It does **not** close the endpoint-weight or current problem.  The source
entries are THM-2471 service masses.  They are not equal to the THM-3514
`q_H-q5` bucket values, with its phase and `13^-3` normalizer.  The result
therefore supplies no frozen bridge value, tree determinant, physical
chronology/current, all-unit grouped coefficient, row exclusion, or LRC(14)
conclusion.

While this audit was running, commit `b5950b68e` used the newly exported
pair tensor in a separate finite address-weighted endpoint pullback.  That
downstream package is not audited here and does not alter this verdict.  Even
on its own stated terms it distinguishes the auxiliary finite pullback from
a same-integral temporal/Fubini current.  The remaining frontier is therefore
weight-and-chronology sensitive, not another support census.

Concurrent sidecar `77ae45486` independently isolates the same zero-fibre
boundary: the ancestry relation has no `s=0` mass, whereas the Cartesian
endpoint bank has four nonzero `d=0` corners.  Its off-diagonal transplant
therefore supports the typed `s=-d` reading but, consistently, still leaves
the connection/current debt open.

## 8. Reproduction

Run from the repository root:

```text
python -B 04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_independent_audit_20260816.py
```

The current-origin candidate/output hashes are

```text
83f1fa49ac4d02e21a1d76fed169d101715a6620342714ed05b9172ae967a730
7727571722f69a9c59af0183c50c8bae0b7944d6acb7bbb1b3c0e9bd02d54565.
```

Candidate and independent verifier were each replayed in normal and `-O`
mode; all four runs reproduced their respective stored outputs exactly, and
each normal/optimized pair agreed.  The independent script/output hashes are

```text
ec16f5124c5b83a10337fba6046d08572251d0776285b7e250a86539a26ddf97
ebe5acbe2e19465ed10f1380eaa2bd19975849d72be64151514cb0087e645d4f.
```

The pinned independent semantic digest is

```text
ad81e945207703956f5d6ec300430d562c9f98a9ec8788011119a4857ab34e01.
```
