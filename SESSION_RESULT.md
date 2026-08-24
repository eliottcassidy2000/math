# Session Result

## Task

Continue the LRC(14) mathematical frontier for a long session, treating
incoming work and older decoder/tournament ideas as signals to extend rather
than merely summarize.

## Outcome

LRC(14) remains **OPEN**. The session proved and independently audited two
new reduction theorems in the exact rank-eleven `11+2` branch.

1. **THM-4003 — two-boundary body gain and scale-two component erosion.**
   If `U` is the body maximum, `Q=91^6`, and `H=min(U,Q)`, then

   ```text
   U lambda_+(u)>=1/42+1/(84(H-1)).
   ```

   This combines the two endpoint owners of a cited deep body component with
   THM-3818's internal pair-height cap. In scale two, THM-3995's two separate
   retained cores then force

   ```text
   3t(2H-1)<=8(H-1)(U-1).
   ```

   For `U<=Q` and odd `t`, this is exactly
   `U>=floor(3t/4)+2+1_(t=1 mod 4)`. A separate top-balance selector, with
   `V` the second-largest body speed, closes all 17 conditional certificate
   types for `U/V>=11` and scale two already for `U/V>=1001/189`.

   Direct wall grids prove the four oriented owner-residue formulas and the
   reflection identities `e_1=e_4`, `e_2=e_3`. In 62,989 old-strip cells
   through `t<=1001`, the symbolic gate closes 742, exact residues close 77
   more, and 62,170 remain. Nineteen closures require both improvements in
   the same inequality.

2. **THM-4004 — three-detuned divisor combs and the `t<U` floor.**
   The common labelled lift branches give the necessary multiplier-prime
   profile

   ```text
   ell>=5: at most 7 body coordinates divisible by ell;
   ell=3:  at most 8;
   ell=2:  at most 9, with equality only at scale one and
           reduced odd exception sum >7.
   ```

   Every exact `t<U` survivor satisfies `U>=3,208,300,859`. The literal
   component-swapped width proof cannot fire there because
   `t lambda(u)<6/7`. The prime-three displayed row is a sharp selector
   hostile, not an LRC hostile: it is lonely at `4/33`.

## Incoming signal and corrected boundaries

- Incoming THM-4002 was compared before promotion. It retains full signed
  cross-phase and closes fixed bodies in `[1,21]`; THM-4003 instead prunes
  arbitrary-body parameter and ratio regions. Neither subsumes the other.
- The discovery phrase calling `(3,7)` a minimal two-lift hostile was false.
  `(3,5)` is already hostile at `3/16`; the theorem retains only the correct
  reduced-sum boundary.
- A uniform scale-one endpoint shrink is refuted by exact congruence
  alignment. The positive survivor is the deletion-gcd filter: every prime
  dividing `t` must miss at least two body coordinates.
- Exact decoder controls show why graph density alone stops. Among all 31,434
  primitive connected eleven-subsets of `[1,18]`, the unique minimizer of
  `U lambda_+` is `(1,2,3,4,8,9,10,11,12,13,14)` with `7/110`; a connected
  hostile has `r_+=432>U=355`. The missing sidecar is the circular owner-wall
  word, not graph connectivity.

## Verification

Primary and independent certificates:

```text
04-computation/lrc14_scale_two_component_erosion_boundary_strip_thm4003.py
04-computation/lrc14_scale_two_component_erosion_boundary_strip_independent_audit_thm4003.py
04-computation/lrc14_tltu_divisor_comb_profile_thm4004.py
04-computation/lrc14_tltu_divisor_comb_profile_independent_audit_thm4004.py
```

Normal and optimized executions agree; stored LF outputs and SHA-256 hashes
are pinned in THM-4003/4004. The direct residue audit covers 7,764,075 exact
conditions; the divisor-comb paths independently audit open walls, branch
multiplicities, primes, controls and strict rounding.

## Next frontier

The highest-value next tests are the joint four-endpoint owner assignment
under the divisor profile, the exact phase-labelled multi-component
projection in `t<U` at primes 2 and 3, cross-phase on the trimmed rather than
untrimmed scale-two cores, and compatibility between the connected decoder
edge word and the circular wall-owner word. See
`07-reflections/lrc14-component-erosion-divisor-comb-frontier-root-20260824.md`.
