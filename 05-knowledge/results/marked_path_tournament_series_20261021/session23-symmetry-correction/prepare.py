from pathlib import Path
from difflib import unified_diff
import hashlib,json,sys
ROOT=Path(sys.argv[1]).resolve() if len(sys.argv)==2 else Path('/workspace/scratch/07d5532f2fce/math')
OUT=Path(__file__).parent
thm='01-canon/theorems/THM-260-rooted-tournament-layer-decomposition.md'
ledger='01-canon/MISTAKES.md'
before={p:(ROOT/p).read_text() for p in (thm,ledger)}
after=dict(before)
anchor='**Filed by:** kind-pasteur-2026-03-20-S6\n\n'
correction='''**Correction (2026-10-21):** The general layer denominator below now uses
factorials of the multiplicities of each nontrivial cycle length, rather than
`r!`. Unequal cycle lengths cannot be interchanged in the centralizer. The
mixed `(5,3)` layer is `128/15`, as already recorded correctly in the table;
the rooted count prefix and all displayed table entries are unchanged. This
is a scoped denominator repair, not a re-audit of unrelated interpretations.

'''
if after[thm].count(anchor)!=1:
    raise RuntimeError('theorem insertion anchor changed')
after[thm]=after[thm].replace(anchor,anchor+correction,1)
old='''For partition lambda = (k_1, k_2, ..., k_r, 1^m) with all k_i >= 3 odd:
D_lambda(n) = 2^{f(lambda,n)} * (n - m) / (prod k_i * r! * m!)

where f(lambda, n) = number of free edge-orbit pairs.
'''
new='''For partition lambda = (k_1, k_2, ..., k_r, 1^m) with all k_i >= 3 odd,
let a_l be the multiplicity of the nontrivial cycle length l. Then:
D_lambda(n) = 2^{f(lambda,n)} * (n - m) / (m! * prod_l l^{a_l} a_l!)

where f(lambda, n) = number of free edge-orbit pairs. The denominator is the
permutation centralizer size z_lambda. Its rooted Burnside term is
m * 2^{f(lambda,n)} / z_lambda, so subtracting it from n times the ordinary
term proves the displayed formula. For lambda = (5,3), f = 4 and z_lambda =
15, giving 2^4 * 8 / 15 = 128/15; the former r! denominator gave 64/15.
'''
if after[thm].count(old)!=1:
    raise RuntimeError('theorem formula anchor changed')
after[thm]=after[thm].replace(old,new,1)
anchor='''# Mistakes Log

**Purpose:** preserve failed implications, minimal witnesses, and repaired
statements so future agents do not repeat them. Search this ledger before using
a historical synthesis or promoting computational evidence.

'''
entry='''- **2026-10-21, THM-260 rooted Burnside layer — unequal cycle lengths were treated as interchangeable.** The general denominator used `r!` for all r nontrivial cycles, instead of the product of the factorials of their equal-length multiplicities. At cycle type `(5,3)`, the displayed formula therefore gave `64/15`, contradicting the theorem's own correct `128/15` table entry. The repaired denominator is `m! * product_l l^(a_l) a_l!`, where m counts fixed points and a_l counts nontrivial cycles of length l; the layer is `(n-m) * 2^f / z_lambda`. The rooted Burnside term is `m * 2^f / z_lambda`, which proves the repair by subtraction from n times the ordinary term. Correct table entries and the rooted prefix are unchanged. Missing coordinate: cycle-length multiplicity, not merely the number of cycles. [Corrected theorem](theorems/THM-260-rooted-tournament-layer-decomposition.md).

'''
if not after[ledger].startswith(anchor):
    raise RuntimeError('maintained ledger header changed; do not prepend')
after[ledger]=anchor+entry+after[ledger][len(anchor):]
# Scope guards: preserve both maintained headers, numeric table, and unrelated ring sentence.
if not after[thm].startswith('# THM-260: Layered Decomposition of the Rooted Tournament Sequence\n'):
    raise RuntimeError('theorem banner displaced')
table='| 8 | 8192/15 | 128/3 | 512/3 | 16 | 128/15 | 784 |'
if before[thm].count(table)!=1 or after[thm].count(table)!=1:
    raise RuntimeError('correct table was not preserved')
ring=before[thm][before[thm].index('The layers form a GRADED RING structure:'):]
if not after[thm].endswith(ring):
    raise RuntimeError('unrelated graded-ring text modified')
patch=''.join(''.join(unified_diff(before[p].splitlines(keepends=True),after[p].splitlines(keepends=True),fromfile='a/'+p,tofile='b/'+p)) for p in (thm,ledger))
(OUT/'THM260-correction.patch').write_text(patch)
(OUT/'correction-receipt.json').write_text(json.dumps({'status':'PREPARED; NOT APPLIED','original_sha256':{p:hashlib.sha256(before[p].encode()).hexdigest() for p in before},'patched_sha256':{p:hashlib.sha256(after[p].encode()).hexdigest() for p in after},'patch_sha256':hashlib.sha256(patch.encode()).hexdigest(),'scope_guards':['original headers preserved','correct mixed-cycle table preserved','unrelated graded-ring text preserved','no scarce ID reserved']},indent=2)+'\n')
print('Prepared exact two-file patch outside repository; no repository files modified.')
