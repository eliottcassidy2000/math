# LRC14 Private-Firewall Nerode Reflection

HYP-3513 changes the shape of the HYP-3491 follow-up.  I expected the
private-label firewall to need a new component-label sidecar.  The finite
Nerode audit says the existing HYP-3474 colored axes already decide the
private-firewall bit: `C`, `F`, `N`, and `T` are pure for both the boolean
status and the stronger `nondead/nonprivate/private` status.

The route split is the real obstruction.  No subset of `K,N,T,S,F,C,M,A`
preserves HYP-3490's five-way route.  The route sidecar `R` is compact and
complete, while incidence/projection sidecars `I` and `Q` are compact
three-fiber carriers for private status only.

The proof-facing move is therefore more specific.  Prove an incidence-cut
lemma that identifies `I/Q` with dead-component/blocker-label multiplicity
one, without naming rows.  Then keep `R` as a route sidecar unless a separate
route-reconstruction theorem can be proved.  HYP-3486..HYP-3481 still own the
random031 route, and HYP-3480/HYP-3478 still own the six non-hard private
rows; HYP-3513 just makes the private predicate itself compatible with the
existing colored-gate quotient interface.
