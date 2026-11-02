#!/usr/bin/env python3
"""Verify only the explicitly curated research sources, independent of cwd."""
from pathlib import Path
import hashlib,json
here=Path(__file__).resolve().parent
root=here.parent
manifest=json.loads((here/'CURATED_FILES.json').read_text())
rows=manifest['files']
paths=[row['path'] for row in rows]
if paths!=sorted(set(paths)):raise RuntimeError('manifest paths must be unique and sorted')
for row in rows:
    rel=Path(row['path'])
    path=(root/rel).resolve()
    if rel.is_absolute() or '..' in rel.parts or not path.is_relative_to(root):raise RuntimeError('unsafe manifest path')
    raw=path.read_bytes()
    if len(raw)!=row['bytes']:raise RuntimeError('byte count mismatch: '+row['path'])
    if hashlib.sha256(raw).hexdigest()!=row['sha256']:raise RuntimeError('SHA256 mismatch: '+row['path'])
expected=''.join(row['sha256']+'  '+row['path']+'\n' for row in rows)
if (here/'SHA256SUMS').read_text()!=expected:raise RuntimeError('SHA256SUMS differs from curated manifest')
pin=json.loads((here/'PRIMARY_FINAL_REVIEW.json').read_text())
primary=next(row for row in rows if row['path']==pin['reviewed_file'])
if primary['sha256']!=pin['sha256']:raise RuntimeError('primary proof review pin mismatch')
for a,b in [('exact_results.json','exact_results.optimized.json'),('mutation_results.json','mutation_results.optimized.json'),('arithmetic_slices.json','arithmetic_slices.optimized.json'),('producer_arithmetic.replayed.json','producer_arithmetic.replayed.optimized.json')]:
    if (here/a).read_bytes()!=(here/b).read_bytes():raise RuntimeError('normal/optimized mismatch: '+a)
if (here/'producer_arithmetic.replayed.json').read_bytes()!=(root/'stability/eisenstein_bernstein.json').read_bytes():raise RuntimeError('producer replay differs from source receipt')
if (here/'producer_arithmetic.replayed.optimized.json').read_bytes()!=(root/'stability/eisenstein_bernstein.optimized.json').read_bytes():raise RuntimeError('optimized producer replay differs from source receipt')
print(json.dumps({'status':'PASS','curated_files':len(rows),'source_bytes':sum(row['bytes'] for row in rows),'primary_proof_sha256':primary['sha256'],'all_source_digests_verified':True,'normal_optimized_and_producer_replays_identical':True,'scope':manifest['scope']},indent=2,sort_keys=True))
