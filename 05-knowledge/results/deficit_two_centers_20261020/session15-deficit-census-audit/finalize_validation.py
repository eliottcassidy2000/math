#!/usr/bin/env python3
"""Validate completed runs and freeze a plaintext-first audit artifact index."""
import gzip,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'session15-deficit-centers'
def need(t,m):
    if not t:raise RuntimeError(m)
def sha(data):return hashlib.sha256(data).hexdigest()
specs=json.loads((SRC/'residual_frame_specs.json').read_text())
runs={}
for file in ['independent_full_scan.jsonl','producer_certificate_replay.jsonl','sanitized_full_scan.jsonl']:
    data=[json.loads(x) for x in (ROOT/file).read_text().splitlines()]
    need(len(data)==11 and sum(x['checked'] for x in data)==2121984,'incomplete run')
    need(all(x['survivors']==0 for x in data),'survivor')
    need({(x['n'],x['q'],x['low']) for x in data}=={(x['n'],x['q'],-1 if x['low_out_orbit'] is None else x['low_out_orbit']) for x in specs},'missing family')
    runs[file]=dict(families=11,presentations=2121984,survivors=0)
certificates=[]
for spec in specs:
    n,q,low=spec['n'],spec['q'],spec['low_out_orbit'];name=f'n{n}_q{q}_low{low}'
    # Exact agreement of producer plaintext frame descriptors with audited JSON.
    ints=list(map(int,(SRC/f'frames_{name}.txt').read_text().split()))
    expected=[n,q,spec['orbits'],len(spec['free'])]
    for k,rows in enumerate(spec['bits']):expected += [spec['fixed'].get(str(k),-1)]+[z for row in rows for z in row]
    need(ints==expected,'producer descriptor and JSON differ')
    raw=(ROOT/f'independent_{name}.model.witnesses.bin').read_bytes()
    need(raw==(ROOT/f'independent_{name}.model.sanitized.bin').read_bytes(),'sanitizer changed witness stream')
    need((ROOT/f'producer_{name}.model.replayed.bin').read_bytes()==(SRC/f'frames_{name}.txt.witnesses.bin').read_bytes(),'replay stream differs')
    need(len(raw)==5*(1<<len(spec['free'])),'independent witness length')
    gz=gzip.compress(raw,compresslevel=9,mtime=0)
    filename=f'independent_{name}.witnesses.bin.gz'
    (ROOT/filename).write_bytes(gz)
    need(gzip.decompress(gz)==raw,'independent gzip broken')
    certificates.append(dict(file=filename,records=len(raw)//5,raw_sha256=sha(raw),gzip_sha256=sha(gz),gzip_bytes=len(gz)))
(ROOT/'independent_certificate_index.json').write_text(json.dumps(certificates,indent=2)+'\n')
summary=dict(verdict='APPROVED',total_presentations=2121984,runs=runs,producer_plaintext_descriptors_equal_to_json=True,producer_replay_streams_equal=True,all_sanitized_streams_equal=True,address_sanitizer='passed all presentations',undefined_behavior_sanitizer='passed all presentations',leak_sanitizer='disabled because ptrace environment is unsupported',independent_compressed_certificates=11,python_normal_and_optimized_reconstructions_equal=(ROOT/'reconstruction.log').read_bytes()==(ROOT/'reconstruction_optimized.log').read_bytes())
need(summary['python_normal_and_optimized_reconstructions_equal'],'python optimization affects reconstruction')
(ROOT/'final_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
# Pin all inspected inputs, including the full producer proof/certificates and exact consequences.
source_files=[p for p in SRC.iterdir() if p.name in ['DEFICIT_TWO_COMPLETE_CENSUS.md','PREFIX_CONNECTIVITY_CUTOFF.md','EXACT_LEAF_SERIES.md','HAND_LOCAL_CONFIGURATION_SUPPLEMENT.md','generate_residual_frames.py','probe_periodic_frames.py','scan_residual_frames.cpp','residual_frame_specs.json','residual_all_scan.jsonl','residual_sample_checks.json','witness_certificate_index.json','leaf_series_coefficients.json'] or p.name.startswith('frames_n') and (p.name.endswith('.txt') or p.name.endswith('.gz'))]
source_files += [ROOT.parent/'session15-deficit-cutoff-audit/AUDIT.md',ROOT.parent/'session14-black-leaf-reduction/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md',ROOT.parent/'session14-black-leaf-reduction/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md']
(ROOT/'AUDITED_SOURCE_SHA256SUMS').write_text(''.join(f'{sha(p.read_bytes())}  ../{p.relative_to(ROOT.parent)}\n' for p in sorted(source_files)))
package=[p for p in ROOT.iterdir() if p.suffix in ['.py','.cpp','.md','.json','.jsonl','.model','.gz'] or p.name=='AUDITED_SOURCE_SHA256SUMS']
(ROOT/'SHA256SUMS').write_text(''.join(f'{sha(p.read_bytes())}  {p.name}\n' for p in sorted(package)))
print(json.dumps(summary,indent=2))
print('AUDIT.md',sha((ROOT/'AUDIT.md').read_bytes()))
print('HAND_SUPPLEMENT_AUDIT.md',sha((ROOT/'HAND_SUPPLEMENT_AUDIT.md').read_bytes()))
print('SHA256SUMS',sha((ROOT/'SHA256SUMS').read_bytes()))
