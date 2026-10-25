# Root slices, genus two, pants and the Fano matching

Read `THEOREM.md`. The strongest bridge is an actual chain of constructions:

    six/seven-root surface
      -> transverse elementary-2 curve cover of genus17/49
      -> one/seven six-branch genus-two quotients
      -> a marked Fano maximal-isotropic subgroup in each Jac[2]
      -> one/seven Richelot (2,2)-isogenies.

The one/seven size-six characters are exactly the surface two-form sectors.
The old pants restriction map loses precisely the four-element subgroup
selected by a branch matching. A concrete rational slice of the magic
hourglass supplies seven nonsingular Richelot targets. A Möbius obstruction
shows why the Paley order-three label action is not an automorphism of the
fixed displayed genus-two curve.

`check.py` requires Python3 and SymPy. Normal and optimized runs agree byte
for byte. `PROVENANCE.json` pins the inherited statements and primary-source
uses. `normal.json` stores all seven explicit Richelot targets and finite
controls. No integer-square-point consequence or fixed-surface Paley/PSL2(7)
action is claimed; the arithmetic twist of a target equation remains typed.
