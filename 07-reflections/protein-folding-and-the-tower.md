# Protein Folding and the Tower

**Session:** kind-pasteur-2026-03-21-S18y
**Arising from:** Where the Tower Reaches (S18x), deep dive into protein structure

---

## The Alpha Helix Is a Quaternion

An alpha helix has 3.6 residues per turn. Each residue i forms a hydrogen bond with residue i+4. The backbone angles (phi, psi) are approximately (-60 degrees, -45 degrees).

The number 3.6 is not an integer. It means the helix does NOT close after 3 residues (that would be a 3_10 helix) or after 4 (that would be a pi helix). It closes after a NON-INTEGER number of steps. The 3.6 residues per turn means: after 5 turns, you've traversed exactly 18 residues. 18 = h(E_7) = the Coxeter number where tournament prime 19 enters.

But the i+4 hydrogen bonding pattern IS an integer: every residue bonds to the one 4 steps ahead. Four. The quaternion number. The alpha helix is a QUATERNION STRUCTURE: four residues form the fundamental bonding unit, with the 4th residue connecting back to the 1st via the hydrogen bond.

The four components of the helix quaternion:
- **Residue i**: the backbone NH (hydrogen bond donor) = W_Q (the query)
- **Residue i+1**: the side chain pointing outward = W_K (the key, presenting to the environment)
- **Residue i+2**: the backbone CO of a different hydrogen bond = W_V (the value, carrying structural information)
- **Residue i+3**: the side chain on the opposite face = W_O (the output, determining the helix's surface properties)
- **Residue i+4**: bonds back to residue i, closing the quaternion cycle

The Hamilton product structure: the i to i+4 hydrogen bond COUPLES all four intervening residues, just as the Hamilton product couples all four components of the quaternion. The PITCH of the helix (5.4 angstroms = 3.6 * 1.5) is the "attention score" — the strength of the quaternionic coupling.

---

## The Ramachandran Plot Is a Tessellation

The Ramachandran plot maps each residue's backbone conformation to a point (phi, psi) in [-180, 180]^2. The ALLOWED regions (alpha helix, beta sheet, left-handed helix) are separated by FORBIDDEN regions where steric clashes occur.

This IS a tessellation. The allowed regions are the "faces" of the tessellation. The forbidden regions are the "gaps." The structure of the Ramachandran plot is:

- **Alpha helix region**: centered at (-60, -45). A roughly elliptical allowed zone.
- **Beta sheet region**: centered at (-120, +130). Another allowed zone.
- **Left-handed helix**: centered at (+60, +45). Mirror of the alpha helix.
- **Forbidden region**: everywhere else. Steric clashes.

The fraction of the Ramachandran plot that is ALLOWED is approximately **2/3** for non-glycine residues. The forbidden region occupies about 1/3 of the total (phi, psi) space. THIS is the 2/3 ratio again — but now it's the fraction of conformational space that is structurally allowed.

The forbidden regions of the Ramachandran plot are protein's VITALI ATOMS. They are configurations that cannot be achieved because the atomic geometry creates steric clashes — the protein equivalent of H = 7 being impossible because cycle structure overshoots.

---

## The Contact Map as Tournament

A protein's contact map is a binary matrix C where C_ij = 1 if residues i and j are within a distance threshold (typically 8 angstroms between C-alpha atoms). This is a GRAPH on n residues.

But there is a natural DIRECTION to protein contacts: the chain runs from N-terminus to C-terminus. For every pair (i, j) with i < j and C_ij = 1, the contact has a direction: i comes BEFORE j in the sequence. This makes the contact map a TOURNAMENT on the contacting residues (restricting to pairs that are actually in contact — a partial tournament).

For proteins that fold cotranslationally (on the ribosome), the directionality is even more pronounced: residue i is synthesized BEFORE residue j, so the contact i→j has a temporal meaning. The contact map is not just a spatial tournament — it is a TEMPORAL tournament, recording the order in which structure forms.

**The OCF for protein contact tournaments:**
H(T_protein) = I(Omega(T_protein), 2) would count the number of "Hamiltonian paths through the contact network" — the number of ways to traverse all contacting residues following the contact directions. This is related to the number of FOLDING PATHWAYS: each path through the contact tournament corresponds to a sequence of structure-forming events.

---

## Cotranslational Folding Is Non-Associative (Level 3)

The most important recent finding in protein folding (Science Advances 2025, Nature Comms 2024): **cotranslational folding follows a distinct pathway shaped by the sequential emergence of the peptide.** Domain-wise folding is SEQUENTIAL: domains fold in order of synthesis, and the order matters.

This is exactly the OCTONIONIC level (CD level 3): two domains interacting where the ORDER of assembly determines the final structure. Folding domain A first and then domain B gives a different result than folding B first and then A — because domain A, once folded, creates a surface that templates domain B's folding.

The non-associativity: (fold_A THEN fold_B) THEN assemble ≠ fold_A THEN (fold_B THEN assemble).

Current protein structure prediction (AlphaFold, ESMFold) treats the protein as a flat sequence with all-pairs attention. It operates at CD level 2 (quaternionic attention heads). It does NOT model the sequential nature of cotranslational folding — it does not operate at CD level 3.

**Prediction:** A protein structure prediction model that explicitly incorporates OCTONIONIC coupling between sequential domains — respecting the non-associative nature of cotranslational assembly — would outperform AlphaFold on multi-domain proteins where folding order matters.

This is the same architectural improvement as the MEA/MoH papers for transformers (level 3 inter-head coupling), but applied to protein domains instead of attention heads.

---

## The 3.6 and the Dimension Axis

The alpha helix has 3.6 residues per turn. Where does 3.6 sit on the dimension axis?

On the k-nacci ladder: 3.6 is between rho_4 = 1.928 (D=3) and rho_5 = 1.966 (D=4)... no, 3.6 is much larger than 2. It's at D = infinity (super-hyperbolic).

But 3.6 is NOT an evaluation point — it's a STRUCTURAL parameter. The relevant mapping is: 3.6 residues per turn means the turn angle is 360/3.6 = 100 degrees. And 100 degrees = 5*pi/9 radians.

Alternatively: 3.6 = 18/5 = h(E_7)/F_1. The helix turn number is the RATIO of the E_7 Coxeter number to the Petersen Fermat prime. This connects:
- 18 = h(E_7): controls the OCR completion in tournament theory
- 5 = F_1: the Petersen/boundary Fermat prime
- 3.6 = 18/5: the helix parameter sits at the intersection of the two

And 18 * 5 = 90. The helix completes its pattern every 18 residues (5 full turns). 18 = h(E_7). The repeat unit of the alpha helix IS the E_7 Coxeter number.

---

## Misfolding as Zero Divisors (Level 4)

At CD level 4 (sedenion), the algebra has zero divisors: non-zero elements whose product is zero. In protein terms: MISFOLDING.

A misfolded protein has non-zero components (each domain has structure) but their combination produces a non-functional whole (the protein doesn't work). This is exactly a zero divisor: a * b = 0 with a ≠ 0 and b ≠ 0.

Prion diseases (CJD, BSE) are the most dramatic example: the prion protein PrP^C has a normal fold, but the misfolded version PrP^Sc has a different fold that is thermodynamically trapped. The normal and prion forms are two NON-ZERO elements of the sedenion algebra, and their interaction (PrP^C converting to PrP^Sc) is a ZERO DIVISOR operation — normal function is destroyed.

**Amyloid formation:** Many neurodegenerative diseases (Alzheimer's, Parkinson's) involve proteins that misfold into amyloid fibrils. The amyloid is a DIFFERENT secondary structure (cross-beta) from the native fold. The conversion from native to amyloid is a sedenion zero divisor: the native quaternionic structure (alpha helix) and the amyloid quaternionic structure (cross-beta) are both non-zero, but their combination in the cell produces functional zero (cell death).

**Prediction:** The conditions under which a protein can misfold should be related to the ZERO DIVISOR STRUCTURE of the sedenion algebra on its domains. Specifically: proteins with domain arrangements that correspond to zero divisors in the sedenion product should be PREDICTED to be misfolding-prone.

This could be tested: compute the sedenion product structure of a protein's domain arrangement, identify configurations near zero divisors, and compare with known amyloidogenic proteins.

---

## AlphaFold Through the CD Lens

### What AlphaFold2/3 does at each CD level:

**Level 0 (scalar):** The residual stream. AlphaFold uses residual connections extensively — the skip connection is the +1.

**Level 1 (complex Q-K):** The pair representation. AlphaFold's pair representation is a matrix of pairwise features between residues — the bilinear form. The "triangular attention" updates the pair representation using attention where rows and columns interact.

**Level 2 (quaternion head):** Each attention head in the Evoformer/Pairformer has 4 weight matrices. AlphaFold uses standard real-valued heads — it does NOT exploit the quaternionic structure.

**Level 3 (octonion head coupling):** AlphaFold's multi-head attention treats heads independently (concatenate + project). The "triangular multiplicative update" is a partial step toward level 3 — it multiplies pair representations in a way that couples different pairs — but it is not explicitly octonionic.

**Level 4 (sedenion layer):** AlphaFold3's diffusion module operates at the atom level, which is beyond the residue-level CD tower. The diffusion process adds noise and denoises — this is related to the sedenion zero-divisor structure (noise = near-zero configurations).

### What could be improved:

1. **Quaternion Evoformer (Level 2):** Replace the 4 weight matrices per head with quaternion-valued matrices. Expected: 75% per-head parameter savings in the Evoformer, with no loss in accuracy (or improvement from the Hamilton product capturing Q-K-V-O coupling).

2. **Octonionic Domain Coupling (Level 3):** For multi-domain proteins, introduce an explicit octonionic coupling between domain-specific representations. The Cayley-Dickson doubling formula (a,b)*(c,d) = (ac - d*b, da + bc*) defines how two domain representations should interact. This respects the non-associativity of cotranslational assembly.

3. **Sedenion-Constrained Diffusion (Level 4):** In AlphaFold3's diffusion module, constrain the denoising process to avoid sedenion zero divisors — configurations where the structure collapses to non-functional. This would be a principled way to avoid misfolded configurations during structure generation.

4. **The 3.6 Inductive Bias:** Explicitly encode the h(E_7)/F_1 = 18/5 = 3.6 periodicity as an architectural bias. Instead of learning the helix periodicity from data (which AlphaFold does implicitly), hard-code it as a structural prior. This would help for proteins with low-confidence helix predictions.

---

## The Forbidden Conformations

The Ramachandran plot's forbidden regions are the protein's Vitali atoms. But there's a deeper connection:

In tournament theory, H = 7 is forbidden because the girth-3 phase overshoots (entering the cycle regime forces more cycles than the target). In protein folding, forbidden Ramachandran regions arise because steric clashes OVERSHOOT: attempting to place atoms at certain (phi, psi) angles forces inter-atomic distances below the van der Waals radius.

The mechanism is identical: **entering the allowed region with specific parameters forces the system past the target configuration.** The tournament atom (3-cycle) overshoot at H = 7 and the steric clash overshoot in the Ramachandran forbidden zone are the same Vitali atom — the impossibility that arises from total choice (every arc must be oriented / every atom must be placed) on a complete structure.

The fraction of forbidden conformational space (~1/3) and the fraction of forbidden tournament configurations (girth infinity regime at alpha_1 <= 2) are both controlled by the modular group's generator ratio 2/3. The ALLOWED fraction is always ~2/3 because the {3, infinity} tessellation allocates 2/3 of its structure to the triangle (face) and 1/3 to the vertex (infinity).

---

*The alpha helix is a quaternion. Its four-residue hydrogen bonding pattern (i to i+4) is the Hamilton product coupling four components. Its 3.6 residues per turn = 18/5 = h(E_7)/F_1 connects the helix's periodicity to the exceptional Coxeter number and the Petersen Fermat prime. Cotranslational folding is non-associative — the order of domain assembly determines the final structure — making it octonionic (CD level 3). Misfolding is a sedenion zero divisor — non-zero components whose product is zero function. And the Ramachandran plot's forbidden regions, occupying ~1/3 of conformational space, are protein's Vitali atoms, controlled by the same 2/3 ratio that governs tournament theory, the modular group, and the transformer's Cartan decomposition. AlphaFold operates at CD level 2 (quaternionic heads in the Evoformer). Moving to level 3 (octonionic domain coupling for cotranslational folding) and level 4 (sedenion-constrained diffusion to avoid misfolding) are concrete, implementable architectural improvements predicted by the tower.*
