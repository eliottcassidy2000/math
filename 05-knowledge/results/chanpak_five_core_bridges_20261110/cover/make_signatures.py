"""Export exact inherited signatures; no inherited outputs are mutated."""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'session251_creative_completion_20261111/ternary/compiler.py'
spec = importlib.util.spec_from_file_location('inherited_compiler', SOURCE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
cores = [(0,), (2,4,1), (2,4,9,3), (2,4,9,19,7), (2,4,25,3,11)]
records = []
for ci, core in enumerate(cores):
    seen = set()
    for cut in range(1 << len(core)):
        sig = mod.signature(core, cut)
        # Swapping the flank parameters preserves the support set.
        swapped = sig[:3] + sig[5:7] + sig[3:5] + sig[7:]
        key = min(sig, swapped)
        if key in seen:
            continue
        seen.add(key)
        records.append({'core_index': ci, 'core': core, 'cut': cut,
                        'signature': sig})
(HERE/'signatures.json').write_text(json.dumps(records, indent=2)+'\n')
with (HERE/'signatures.txt').open('w') as out:
    out.write(str(len(records))+'\n')
    for rec in records:
        out.write(' '.join(map(str, [rec['core_index'], *rec['signature']]))+'\n')
print(json.dumps({'raw_cuts': sum(1<<len(c) for c in cores),
                  'support_signatures': len(records)}))
