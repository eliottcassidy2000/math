# Field-sensitive monomial jet lengths

19 October 2026. Sharp finite-word approximation formulas for the scalar-time flows of x²y and x³.

For the target (x,y+x^r), r≥4, through degree N≥r, put k=r−2 and D=N−2.

- Real scalar fields: the minimum nonzero pole count is max(k,D−1). The shortest word has length 2k+2 when D=k, and 2D when D≥k+1. All upper bounds use rational times.
- Algebraically closed characteristic-zero fields: the minimum pole count is max(k,D−k+1). The shortest word has length 2k+2 for k≤D≤2k−1, and 2(D−k+1)+1 thereafter. Roots of unity provide the constructions.
- The first length separation is the degree-six jet of (x,y+x⁴): eight real-time factors versus seven factors over Q(ζ₃).
- An explicit rational shear over Q belongs to the extended-field group but not to the real-time group at any length, since its finite poles are nonreal.

The proof reduces the problem to a finite derivative functional and split-polynomial coefficient gaps. Strict real-root interlacing gives the additional real obstruction. This is a result about the specified generating pair, not arbitrary Keller maps or the Jacobian conjecture.

## Files

- `field_sensitive_monomial_jet_lengths.md`: full proofs, optimal constructions, and field scope.
- `check_field_sensitive_jets.py`: exact finite controls.
- `field_sensitive_jets_exact_controls.json`: 3011 passing checks summarized by category, with exact examples and a digest of the full ordered check labels.

Run `python3 check_field_sensitive_jets.py` with SymPy installed. General minimality is proved in the note, not inferred from the finite checks. No historical-novelty claim is made.
