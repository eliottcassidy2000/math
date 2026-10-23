# Fifth audit: branch partition functions and the complete frequency region

Status: PROVED theorem review; FINITE-EXACT replay PASS.

The exact multivariate prefix polynomials follow from the already audited
Cartesian language theorem. Unlike products of independently summarized
shared-port marginals, these factors have a proved mutual-compatibility basis.
A single canonical-C position pinned to B removes exactly one free bit in R;
other fixed prefix choices pin their own factors. The entropy rates are unchanged
by finitely many feasible pins, while the same B pin is infeasible in D.

The complete limiting-frequency region is correctly parametrized by the
absolute densities u of changed C positions and v of changed B positions:
f_C=alpha-u, f_A=alpha+v, f_B=beta+u-v. Necessity follows from the counts;
sufficiency follows by choosing subsets of the two independently enumerated
eligible position sequences with prescribed relative densities. The forced-A
positions imply lower asymptotic A density at least alpha>1/3 even for branches
without limiting frequencies, so equal three-way frequencies are impossible.

The uniform finite-prefix law gives limiting EXPECTED prefix frequencies
(1/2, phi^-1/2, phi^-2/2) in order A,B,C. If claiming almost-sure limiting
frequencies for the induced infinite product measure, one additionally uses
the independent-bit strong law. This distinction does not affect any finite
partition function or expected-prefix calculation.

Both normal and optimized checker replays agree and PASS all374 controls.
The scope remains occurrence-dependent rows, not a single globally shared
noncanonical block dictionary. The arithmetic selector-tautology warning is
correct and does not diminish the independently fixed canonical register.

## All-branch aperiodicity strengthening

The proposed A-indicator theorem is valid for every allowed R stripe, hence
also D. Canonical A positions floor(k*phi^2) are forced to have A-indicator1;
canonical C positions one less are forced to have A-indicator0 even after
all allowed recolorings. For each modulus q, density of the irrational
rotation {k*phi^2/q} implies that both position sets meet every residue class
infinitely often. No residue class can therefore have the constant eventual
A-indicator required by an eventual period q. Thus the indicator is never
eventually periodic, regardless of whether its frequency limit is rational
or exists at all.

A rational ordinary generating function of this binary sequence would give
a deterministic linear recurrence of fixed order eventually. There are only
finitely many binary history states, so the sequence would then be eventually
periodic, a contradiction. Only the A-indicator conclusion is claimed:
other color classes can be erased by allowed choices and have rational GFs.

The finite mod144 clock was independently replayed from a substitution-generated
word, without the radical floor formula: the first5000 canonical A positions
cover95 residue classes, and complete coverage first occurs at occurrence11001.
The separately included check_residue_clock.py produces identical normal/-O
receipts. This is not a uniform quantitative bound for arbitrary moduli.
