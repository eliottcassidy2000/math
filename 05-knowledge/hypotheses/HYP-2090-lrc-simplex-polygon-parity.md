---
id: HYP-2090
status: SYNTHESIS — tournaments simplex(mesh)+polygon(dihedral); regular/dihedral=odd-m (verified); LRC even-n tightness = observer-marked regular n-gon
source: opus-2026-06-02-S567
related:
  - HYP-2086
  - HYP-2089
  - HYP-2045
---

# HYP-2090: tournaments as simplex+polygon; the LRC worry-set's parity-dihedral structure

- **Two faces (trinity):** tournament = SIMPLEX (oriented K_n, whole mesh, A000568) AND POLYGON (permutohedron/zonotope; ROUND = reconstructible from the cyclic gap-necklace 'outside', A000016 = S565 orbit side). Round 'acts like its outside'; general 'acts like the whole mesh'.
- **Dihedral face = odd m only (VERIFIED):** a regular tournament (out-deg (m-1)/2, Aut=C_m, +reversal=dihedral D_m) exists iff m odd. Even m has NO regular tournament. 'Every other n.'
- **LRC parity dichotomy:** the loneliest config is the geometric regular n-gon (D_n) for all n, but the runner TOURNAMENT (m=n-1) has the dihedral/regular face iff m odd iff n EVEN. So EVEN n (doubled primes 10,14,22) = polygon-shaped worry-set (apex (q,q) reflection, S547); ODD n = simplex/irregular. The repo's even-fold/(q,q)/mod2×mod7/dual-Burnside-fix all live on the even-n polygon face.
- **Precision:** the AP block is the n-gon MINUS the observer vertex (near-regular {5:6,6:7} at n=14, with antipodal ties), exact symmetry = the apex reflection (S547), not full D_{n-1}; 'dihedral' lives on the full n-gon (observer+runners), marked at the observer.
- **n+2 recursion = the parity carrier (Mode B, n->n-2 Cayley-Dickson):** the dihedral face recurs along n+2; regular tournaments = m=3,5,7,... (D_3,D_5,...); LRC polygon-tight worry = even-n sub-sequence with doubled primes 2q as the prime-clean rungs.
- HONEST: trinity + regular<=>odd-m grounded; NEW = the LRC even/odd parity dichotomy + its identification with n+2/Mode-B; precise symmetry is apex-reflection not full D; not a proof.

**See:** `07-reflections/lrc-tournaments-as-simplex-and-polygon-the-parity-dihedral-face-s567.md`, `04-computation/lrc_simplex_polygon_parity_s567.py` (+.out), `07-reflections/the-polygon-simplex-staircase-trinity.md`; HYP-2086 (dual Burnside), HYP-2089 (strong lens), HYP-2045 ((q,q)/n=2q), S547, Mode B.
