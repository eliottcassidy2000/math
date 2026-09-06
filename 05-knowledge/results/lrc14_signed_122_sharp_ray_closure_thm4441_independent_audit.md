# Independent referee: signed `(1,2,2)` all-height closure

## Verdict

**PASS after two wording repairs; no mathematical constant, cutoff, family,
or equality-case correction is needed.**  The theorem is safe at status

```text
PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT +
INDEPENDENTLY AUDITED.
```

It removes only the signed `(1,2,2)` residual from THM-4437.  It does not
close `(1,1,1)`, `(1,1,2)`, chart entry, synchronization, or LRC(14).

## Wording repairs before promotion

1. In Section 3, replace “choose the primitive integral covector `eta` on
   `L` such that `w cross eta=u`” by the typed statement:

   ```text
   Since primitive w makes w cross Z^3 -> ker_Z(w) surjective, choose
   eta in Z^3 with w cross eta=u, and define phi(C)=eta dot C on L.
   Then ker(phi)=Z u.
   ```

   The bound is `|phi(C)|=|u dot e|<15/14`.  This is the argument already
   intended, but the current sentence conflates a covector on `L` with an
   ambient integral vector.

2. Call the supplied C++ check “independent of the ray formulas” or an
   “orthogonal literal-engine replay,” not a fully independent implementation:
   it deliberately imports the already audited THM-4434 literal sheet engine.
   The new clean-room Python referee below supplies the genuinely independent
   raw-carrier implementation.

An optional useful endpoint sentence is:

```text
No admissible address can lie on a roof: if 3 does not divide k, then
14 k |u_i| is nonzero modulo 3 whereas 3(w_j+w_k) is zero modulo 3.
```

Thus the strict and weak roof sets happen to agree on eligible ray addresses,
although the proof correctly retains strict `<`, `R_<(T)`, and zeroing at the
cutoff.

## Exact audit results

- Exhaustive normalized sign vectors modulo global reversal: `12`; exactly
  four survive positivity and strict sorting, with the claimed F1--F4
  equations.  The only algebraic proper-cone overlap is F2/F3 at multiples of
  `(2,5,6)`, excluded by the ternary-unit condition.
- H170 typed head: `1,951` rows, split `280/559/744/368`.
- Every H170 row was recomputed by brute-force `(C_1,C_2)` enumeration under
  all three strict roofs, without the producer's modular solver.  Every raw
  carrier set equals exactly the claimed `Z u` ray with `3`-deleted addresses.
- H611 independent ray extension: `24,876` rows, split
  `3553/7103/9483/4737`; no leader changes and no `6/77` failures.
- Row-by-row H611 C++ literal/ray comparison: `124,380` value checks,
  `5,279,340` native contacts, exact agreement in all three projections,
  physical mass, and the factor-three contact multiplicity.  O0/O3 agree.
- All displayed selected projection integrals and the physical identity
  `integral f_phys=(7/8)r^2` were reconstructed from exact rational
  breakpoints.  Every claimed switch has equality at its stated point, and
  the omitted profile is never lower.
- The exact first strict tail cutoffs reproduce as network
  `61/65/51/171` and physical `51/46/50/151`; hence H170 is sufficient.
- Sharp family leaders and complete equality loci reproduce:

  ```text
  network: F1 5/77 {(1,10,22)}
           F2 51/770 {(1,11,20)}
           F3 46/665 {(2,19,20)}
           F4 3/49 {(10,14,17),(13,14,20)}
  physical: F1 5/77 {(1,10,22)}
            F2 51/770 {(1,11,20)}
            F3 173/2660 {(2,19,20)}
            F4 731/12740 {(13,14,20)}
  ```

  Therefore the global sharp constants are exactly `46/665`, uniquely at
  `(2,19,20)`, and `51/770`, uniquely at `(1,11,20)`.  Both are strictly
  below `6/77`.

## Reproduction

```powershell
python -B 04-computation/lrc14_signed_122_cleanroom_referee_thm4441.py
g++ -std=c++20 -O0 04-computation/lrc14_signed_122_literal_referee_thm4441.cpp -o cross_O0.exe
g++ -std=c++20 -O3 -DNDEBUG 04-computation/lrc14_signed_122_literal_referee_thm4441.cpp -o cross_O3.exe
./cross_O0.exe
./cross_O3.exe
```

The clean-room run reports `91,681` explicit checks.  The supplied producer
also passes under ordinary and `python -O` execution with byte-identical text
and `204,876` explicit checks.  Fresh O0/O3 builds of its H170 literal wrapper
agree, and the H611 optimized extension reproduces every leader.
