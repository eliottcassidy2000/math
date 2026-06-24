## codex-S130 -- Mutated Farey Payloads Rank LRC14 Side Channels

Tested the owner's four Farey mutations (`p+q`, `p*q`, `q^p`, `p^q`) against the HYP-2930 mediant/tournament interface, identifying that these "mutated operators" serve as a hierarchy of side-channel ledgers rather than scalar theorem replacements (commit `1bb78f96`). This checkpoint stabilizes the "fraction address" ledger for the LRC14 proof tree.

### 1. The Mutated Farey Operator Scout
Conducted an exact audit of four mutated Farey operators to rank how proof quotients handle exact fraction addresses.
- **Finding:** The ordinary denominator $q$ remains the irreducible binding-scale coordinate. The mutations serve as **side-channel ledgers**:
    - **Additive Mutation ($p+q$):** Acts as the **Farey recursion ledger**. It tracks the unit-excess chain recursion (e.g., $p/(14p-1)$) where $q$ walks by $+14$.
    - **Product Mutation ($p*q$):** Acts as the **multiplicative/coimage ledger**. It carries divisor profiles, exact-period units, and product-Mobius packet data.
    - **Power Mutations ($q^p, p^q$):** Function as **magnitude-blindness stress tests**. If a tournament invariant is unaffected by these mutations, it is likely leaking only finite class data and missing the unbounded scale.
- **Tournament Ranking:** Established a transitive majority tournament of proof roles: $q > (p+q) > (p*q) > (p^q) > (q^p)$.

### 2. Role in Proof Tree (R_sf and Resonance)
The mutated operators address the **resonance-packet** and **residual-leak** problems by providing a structured ledger for non-scalar data:
- **Resonance Packets:** The multiplicative ledger $(p*q)$ captures the coimage/periodicity data that a raw denominator ignores, providing a carrier for resonance-graph state.
- **Residual Leak ($R_{sf}$):** The additive ledger $(p+q)$ tracks the progress of the $n+2$ recursion and $LRC(14)$ residuals without trying to collapse them into a single scalar bound.

### 3. Handoff: The Fraction Address Ledger
The proof now treats the lonely runner witness as a multi-labeled **fraction address**:
```text
(excess e, binding scale q, additive ledger p+q, multiplicative ledger p*q)
```
This address is used to compare tournament spectra without losing the metric magnitude information (the "blindness" risk identified in HYP-2926).

### 4. Integration with Unit-Distance (S129)
The "scalar-plus-side-channel" lesson from the $u(22)$ universality audit is reinforced: a single quotient is too loose. Just as unit-distance extension requires geometric cocyclicity data beyond graph traceability, LRC14 tightness requires side-channel ledgers ($p+q, p*q$) beyond the raw binding scale.

### 5. Net Impact
This checkpoint stabilizes the project's analytical "ledger" system. By formalizing the mutated Farey operators as side-channel carriers, the cluster has unified the project's Farey-geometric and multiplicative-coimage branches, providing a robust guardrail against magnitude-blindness in the final LRC14 closure.
