#!/usr/bin/env python3
"""Prepare narrow edits without writing to the repository under review."""
from pathlib import Path
import difflib, hashlib, json
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / 'math'
REPLACEMENTS = {
'01-canon/theorems/THM-587-per-level-signed-cycle-index-antipodal-euler-number.md': [
('self-converse count = the ANTIPODAL EULER/LEFSCHETZ number', 'self-converse count = the GRADED ANTIPODAL TRACE'),
('level-graded character + Lefschetz fixed-point count of the complement involution', 'level-graded character + trace fixed-class count of the complement involution'),
('  - HYP-3544  # the Kaczynski/computational-homology refinement of P_n(-1) into Betti numbers', '  - HYP-3544  # proposed Kaczynski/computational-homology refinement; requires a specified complex and differential'),
('   This is the **Lefschetz fixed-point count** of the complement involution `R` on the iso-class space:', '   This is the **trace fixed-class count** of the complement involution `R` on the iso-class space:'),
('   It is the **equivariant Euler characteristic** of the level-graded invariant complex (the antipodal\n   Euler number), since `R` acts as `(-1)^k` on level `k` (THM-584).', '   This is the **graded antipodal trace**, since `R` acts as `(-1)^k` on level `k` (THM-584).\n   No complex or differential is specified here, so no topological Euler characteristic or Betti-number\n   interpretation is asserted.'),
('`P_n(1)`: Molien at `x=1` is `(1/n!) sum_sigma dim(W^sigma) = #orbits on Q_d = A000568(n)`.', '`P_n(1)`: for `W = direct_sum_k W_k`, Molien at `x=1` is\n`(1/n!) sum_sigma trace(sigma on W) = #orbits on Q_d = A000568(n)`.\nEach trace counts labelled tournaments fixed by `sigma`; it is not `dim(W^sigma)`.'),
('- **`P_n(-1) = SC(n)` is the antipodal Euler number.** Complement `R` is free on the `2^d` labeled\n  tournaments (`x != x XOR 1`), so `(Q_d, R)` is a free `Z_2`-space (the Borsuk-Ulam setting); the\n  self-converse census is its equivariant Euler characteristic. This is the entry point for the\n  Ky-Fan / computational-homology refinement (HYP-3544): refine `chi = SC(n)` into the `Z_2`-equivariant\n  Betti numbers of the level-graded complex; the R-odd Betti numbers are the Borsuk-Ulam obstruction.',
 '- **`P_n(-1) = SC(n)` is the graded antipodal trace.** For `n>=2`, complement `R` is free on\n  labelled tournaments but can fix isomorphism classes. The self-converse census is the trace of\n  this class involution. The proposed Ky-Fan / computational-homology refinement (HYP-3544) requires\n  an independently defined equivariant complex and differential; no Betti numbers or Borsuk-Ulam\n  obstruction follow from this character formula alone.'),
],
'01-canon/theorems/THM-588-metagraph-laplacian-algebraic-connectivity-four.md': [
('and the cyclicity 3-cycle count is a nonzero level-2 invariant', 'and the centered cyclicity is a nonzero level-2 invariant'),
('and the Fiedler (slowest) mode is the 3-cycle count', 'and the Fiedler (slowest) mode is the centered 3-cycle count'),
('and the Fiedler mode = cyclicity (|corr|', 'and the Fiedler mode = centered cyclicity (|corr|'),
('6. **Fiedler mode = cyclicity.** The eigenvector at `A`-eigenvalue `d-4` (level 2) is the 3-cycle count\n   `c(T) = C(n,3) - sum_i C(s_i,2)` (verified `|corr| = 1.000` n=4,5; `0.998` n=6). It is the\n   slowest-equilibrating observable of the dominance-reversal walk.',
 '6. **Fiedler mode = centered cyclicity.** A nonzero eigenfunction at `A`-eigenvalue `d-4`\n   (level 2) is `c_0(T) = c(T) - C(n,3)/4`, where\n   `c(T) = C(n,3) - sum_i C(s_i,2)` is the 3-cycle count (historically verified\n   `|corr| = 1.000` n=4,5; `0.998` n=6). The centering mean is for uniform labelled tournaments,\n   equivalently labelled-orbit-size weights on classes, not uniform weights on classes.\n   It is a Fiedler observable of the dominance-reversal walk.'),
('that minimum is `2`, giving algebraic connectivity `4`. (5) is `(min nonzero L-eigenvalue)/d`. (6) the\n`d-4` eigenspace is the level-2 invariant space, spanned (for n=3..9) by the cyclicity. ∎',
 'that minimum is `2`, giving algebraic connectivity `4`. (5) is `(min nonzero L-eigenvalue)/d`.\n(6) A labelled triple is cyclic with probability `1/4`, so the constant Fourier part of `c` is\n`C(n,3)/4`. Its invariant level-1 part vanishes by (2); therefore `c_0` is its nonzero level-2 part.\nHence `A c_0 = (d-4)c_0` and `L c_0 = 4c_0`. The `d-4` eigenspace is the level-2 invariant space,\nspanned (for n=3..9) by `c_0`. ∎'),
("eigenvector`, a degree-`(n-1)` observable); the arc-flip metagraph's slow mode is the degree-2 cyclicity.", "eigenvector`, a degree-`(n-1)` observable); the arc-flip metagraph's Fiedler mode is centered cyclicity."),
('  mode — the cyclicity. This is the structural template', '  mode — centered cyclicity. This is the structural template'),
]}
patch = []
pins = []
for name, replacements in REPLACEMENTS.items():
    original = (ROOT / name).read_bytes()
    old = original.decode()
    new = old
    for src, dst in replacements:
        if new.count(src) != 1:
            raise RuntimeError(('replacement must occur exactly once', name, src, new.count(src)))
        new = new.replace(src, dst)
    patch.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    pins.append({'path':name,'source_sha256':hashlib.sha256(original).hexdigest(),
                 'result_sha256':hashlib.sha256(new.encode()).hexdigest(),
                 'replacements':len(replacements)})
(HERE/'precision.patch').write_text(''.join(patch))
(HERE/'source-pins.json').write_text(json.dumps({'date':'2026-10-24','sources':pins},indent=2)+'\n')
print(json.dumps(pins,indent=2))
