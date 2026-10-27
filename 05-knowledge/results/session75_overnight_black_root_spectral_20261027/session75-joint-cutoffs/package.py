#!/usr/bin/env python3
"""Freeze a verified publication-ready package; excludes mutable logs/checkpoints."""
from pathlib import Path
import argparse,hashlib,json,shutil
ROOT=Path(__file__).resolve().parent

def need(ok,s):
    if not ok:raise RuntimeError(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--destination',default='publish');args=ap.parse_args();dest=ROOT/args.destination
    summary=json.loads((ROOT/'SUMMARY.json').read_text());need(summary['status']=='COMPLETE','all registered runs and verifications must complete')
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text());need(pins['commit']=='96e35642110104cc547516a63f37eebcf930f707','source commit')
    for z in pins['sources']:need(sha(ROOT/z['local_path'])==z['sha256'],'pinned source '+z['local_path'])
    names=['README.md','AUDIT.md','COVERAGE.md','SUMMARY.json','SOURCE_PINS.json','PREREGISTRATION.json','PREREGISTRATION.phase2.json','frontier.txt',
           'scan.py','literal_verify.py','mutation_controls.py','mechanisms.py','pair_check.py','boundary_invoice.py','summarize.py','package.py','run_all.sh']
    names += [str(p.relative_to(ROOT)) for d in ['sources','input'] for p in (ROOT/d).iterdir() if p.is_file()]
    for N in summary['declared_N']:
        stem=f'results/receipt_N{N}';lit=f'results/literal_receipt_N{N}';source=f'results/native_sources_N{N}.json'
        need(sha(ROOT/(stem+'.json'))==sha(ROOT/(stem+'.optimized.json')),'sweep optimized comparison')
        need(sha(ROOT/(lit+'.json'))==sha(ROOT/(lit+'.optimized.json')),'literal optimized comparison')
        x=json.loads((ROOT/(stem+'.json')).read_text());y=json.loads((ROOT/(lit+'.json')).read_text())
        need(sha(ROOT/source)==x['source_ledger_sha256']==y['source_ledger_sha256'],'complete source binding')
        need(x['unresolved_sources']==0 and x['holds_all_rectangles'],'finite test scope result')
        names += [source,stem+'.json',stem+'.optimized.json',lit+'.json',lit+'.optimized.json']
    for stem in ['mutation_receipt','N4096_clock_border_mechanism','N4096_depth_boundary_stability_hostile','resource_port_diagnostics','pair_birth_receipt_N65536','boundary_invoice_receipt']:
        a=f'results/{stem}.json';b=f'results/{stem}.optimized.json';need(sha(ROOT/a)==sha(ROOT/b),'auxiliary normal optimized comparison '+stem);names += [a,b]
    need(len(names)==len(set(names)),'unique package members')
    dest.mkdir(parents=True,exist_ok=True)
    for name in names:
        p=dest/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,p)
    manifest=dict(status='FROZEN_FINITE_EXACT_VERIFIED',source_commit=pins['commit'],declared_N=summary['declared_N'],
                  rectangles=summary['total_rectangles'],original_source_labels=summary['total_source_labels'],
                  limitations=['No universal R3 or joint cutoff theorem.','Second implementation was producer-authored, not an independent-agent audit.','Generic deletion sharpness family is not an authentic A_N family.'],
                  files={name:sha(dest/name) for name in sorted(names)})
    mp=dest/'MANIFEST.json';mp.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    for name,h in manifest['files'].items():need(sha(dest/name)==h,'post-copy verification')
    print(json.dumps(dict(destination=str(dest),files=len(names),manifest_sha256=sha(mp),status=manifest['status']),sort_keys=True))
if __name__=='__main__':main()
