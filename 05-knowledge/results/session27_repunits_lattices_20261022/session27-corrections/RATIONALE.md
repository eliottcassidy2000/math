# Historical repair rationale and application boundary

Status: PROPOSED, source-bound corrections only. No theorem IDs reserved and no repository files or Git state changed. See manifest.json for pinned source commit, blob IDs and before/after SHA-256. Two tracked Python files are absent from the sparse worktree, so their exact originals were read from the pinned HEAD, not reconstructed from the fetched previews. All six originals were inspected; none had an earlier maintained warning. The new byte-zero warnings therefore displace no warning; the near-miss shebang remains first.

The unnumbered proposed mistake entry is in MISTAKE-PENDING.md. The publisher must assign a fresh ID and reconcile any intervening changes before applying. The unified patch passed git apply --check against the originals copy. This is not permission to overwrite a changed upstream file.

## 1. Trace identity versus tournament spectrum

Target: 01-canon/theorems/THM-227-knacci-mersenne.md.
First failed implication: Tr(M_k^k)=2^k−1 was promoted from a matrix trace identity to a forbidden tournament Hamiltonian-path value. Matching the k=3 value gives no source-to-target obstruction. MISTAKE-028 already records attained 31,63,127; THM-344-n8-h63-complete-omega.md supplies 63, and THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md proves the fan constructions. The reviewed two-line proposal is retained with precise path links. Strongest survivor: all Newton-range trace identities, first-departure formula, and arithmetic Mersenne-prime statements, with no forbidden-H consequence. This patch does not adopt MISTAKE-028's own stale 'only two gaps' wording.

## 2. Failure of split-prime pattern is not universal nonmodularity

Target: 07-reflections/the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md.
First failed implication: the inert factors at p=19 (and p=23) refute a specific uniform split-rich explanation, then are claimed to exclude every modular interpretation and all stable arithmetic structure. No universal map class is defined, and prime-factor smoothness does not bound distance from special j-values. Strongest survivor: the retained numerical table, exact r(11)=1729, and failure of the proposed all-split pattern at the next tested genus-one Paley case. The header overrides the historical title and claims while preserving provenance and numerical output.

Important inherited result: THM-438-paley-cluster-integrals-are-catalan.md separately marks R(p)→e PROVED. This task has not audited that proof and must neither silently retract it nor assert globally that the R-limit is open. The correction instead says the reflection's finite table does not prove a limit and distinguishes R(p)=H2^(p−1)/p! from r(p)=H/|Aut|. No universal Paley-maximizer theorem is demonstrated by that table; this is an evidence-scope repair, not a new refutation of Paley maximality. Standard full-subset Held–Karp has exponential space and time; the preserved runtime promise has no benchmark and is explicitly uncertified, not declared impossible.

## 3. Vertex growth is not a group generator

Target: 07-reflections/the-modular-tournament.md.
First failed implication: assigning modular names S and T to reversal and vertex addition is treated as a group action. A group action sends each group element to a permutation of its carrier. Vertex addition increases order, is not onto the finite-tournament universe, and needs extra incidence choices even to be a function. It cannot fix a finite transitive tournament. Independent order-two/order-three examples do not verify a common pair satisfying S²=(ST)³=1. Strongest survivor: arithmetic identities and modular analogies as possible research inspiration; no blanket impossibility of enriched carriers or other actions. The banner intentionally does not audit all unrelated cycle-count claims in this wide historical reflection.

## 4. A formal coefficient solution is not a realizable independence polynomial

Targets: 04-computation/repunit_gap_synthesis.py and 05-knowledge/results/repunit_gap_synthesis.out.
First failed implication: the asserted vector alpha=(T−13,3,1) does not even solve the scalar equation: 1+2(T−13)+4·3+8=2T−5. Correcting its first coordinate to T−10 only repairs that equation, not realizability. An independent triple in the cycle-conflict graph has three singleton subsets and three pairs. Thus alpha_3>=1 forces alpha_1>=3 and alpha_2>=3, contradicting the original claimed vector (0,3,1) at T=13. These necessary constraints imply H>=1+2·3+4·3+8=27 but cannot prove all H>=27 attainable. Strongest survivor: 7=Phi_3(2), 21=Phi_3(4), the exact base identities, and the established omissions 7,21. Global completeness remains OPEN in THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md. Independently, the printed H=21 mod8 inference is false: it gives alpha_1+2alpha_2≡2 mod4, which does not force the claimed separate residues. Warnings preserve every old output byte and do not rewrite the experiment. Some historical bytes are non-UTF-8 (including byte 0x97); the patch and snapshots retain them exactly through surrogate-escape handling.

## 5. Preserve the recurrence, correct the geometric/unit commentary

Target: 04-computation/ramanujan_near_miss_kpc1.py.
First failed implications:
- Degree confusion: X³+Y³=Z³+W³ has homogeneous degree three, so is not a quartic equation. The raw cubic surface and associated elliptic K3 construction are different objects.
- Missing unit assumption: finiteness of Z[omega]^× does not force residual x³+y³−z³ to belong to that group. Positive witness (1,1,2) has residual −6. The intersection of the unit group with Z is {−1,1}, but only conditional on unit membership.
The adjacent factorization is correct and unchanged: N(x+y omega)=x²−xy+y², so (x+y)N(x+y omega)=x³+y³. An initial audit suspicion about its sign was rejected by this direct algebraic check.

The recurrence proof and its two-sided extension do not depend on either of the two requested commentary claims. The proposed header preserves printed historical commentary and executable statements rather than rewriting the proof. The existing independently replayed near-miss output is the proper evidence for the computational survivor, not the uncorrected coda.
