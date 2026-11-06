# Finite-defect polynomial universality

For each fixed j, the coefficient of n^(b+c−j) in the eventual ordinary tournament multiplicity polynomial at H=3^b5^c is Q_j(b,c)/(b!c!), where Q_j has rational coefficients and total degree at most2j. Only strong atoms of prime weight at mostj+1 can affect it. See PROOF.md for the precise domain, proof, and inherited inputs. This is not a minimum-order formula or a global count-completeness theorem.

Run the dependency-free exact controls with assertions enabled:

    python3 check.py > /tmp/session207-defect-check.out
    cmp check.out /tmp/session207-defect-check.out

An independent optional symbolic audit requires an already installed SymPy:

    python3 independent_symbolic_audit.py > /tmp/session207-defect-symbolic.out
    cmp independent_symbolic_audit.out /tmp/session207-defect-symbolic.out

SymPy is required only for that optional checker, not the proof or standard-library verifier. REVIEW.md records independent mathematical audit; PACKAGING_REVIEW.md records final artifact checks. manifest.json freezes every other file. The checks are finite controls, not interpolation proofs of the all-order statement.
