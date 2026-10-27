#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import json,hashlib
ROOT=Path(__file__).resolve().parent

def need(ok,s):
    if not ok:raise RuntimeError(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    prereg=json.loads((ROOT/'PREREGISTRATION.json').read_text());extra=json.loads((ROOT/'PREREGISTRATION.phase2.json').read_text())
    declared=sorted(set(prereg['primary_N']+prereg['structured_N']+prereg['optional_extension_N']+extra['additional_N']))
    rows=[];waiting=[]
    for N in declared:
        path=ROOT/'results'/f'receipt_N{N}.json'
        if not path.exists():waiting.append(N);continue
        x=json.loads(path.read_text()); op=path.with_name(f'receipt_N{N}.optimized.json'); lit=path.with_name(f'literal_receipt_N{N}.json');lo=path.with_name(f'literal_receipt_N{N}.optimized.json')
        flags=dict(sweep_optimized_equal=op.exists() and sha(op)==sha(path),literal_exists=lit.exists(),literal_optimized_equal=lit.exists() and lo.exists() and sha(lit)==sha(lo))
        if op.exists():need(flags['sweep_optimized_equal'],'normal/-O sweep mismatch')
        if lo.exists():need(flags['literal_optimized_equal'],'normal/-O literal mismatch')
        row={k:x[k] for k in ('N','sources','S_full','rectangle_count','holds_all_rectangles','full_denominator_margin','unresolved_sources','max_exterior_state')}
        row.update(D=x['max_joint']['D'],T=x['max_joint']['T'],D_plus=x['max_joint']['D_plus'],max_depth_only=x['max_depth_only'][-1],max_clock_only=x['max_clock_only'][-1],D_plus_over_S=str(Fraction(x['max_joint']['D_plus'],x['S_full'])),verification=flags)
        rows.append(row)
    worst=max(rows,key=lambda r:Fraction(r['D_plus'],r['S_full'])) if rows else None
    out=dict(status='COMPLETE' if not waiting and all(all(r['verification'].values()) for r in rows) else 'PARTIAL_CHECKPOINT',
             source_commit=prereg['source_commit'],declared_N=declared,completed_N=[r['N'] for r in rows],pending_N=waiting,
             total_rectangles=sum(r['rectangle_count'] for r in rows),total_source_labels=sum(r['sources'] for r in rows),
             all_completed_rectangles_hold=all(r['holds_all_rectangles'] for r in rows),worst_D_plus_over_S=worst,rows=rows,
             scope='Finite exact at these N only. Joint extension is stronger than open depth-only R3; no asymptotic claim.')
    (ROOT/'SUMMARY.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    md=['# Exact finite coverage','','| N | sources | rectangles | S | max D+ | D,T | depth-only max | clock-only max |','|---:|---:|---:|---:|---:|---|---:|---:|']
    for r in rows:md.append(f"| {r['N']} | {r['sources']} | {r['rectangle_count']} | {r['S_full']} | {r['D_plus']} | {r['D']},{r['T']} | {r['max_depth_only']} | {r['max_clock_only']} |")
    md+=['',f"Status: **{out['status']}**. Completed: {len(rows)} scales, {out['total_rectangles']} breakpoint rectangles. Registered but incomplete N: {waiting}.",'', 'Normal/-O and fresh literal verification flags are explicit in SUMMARY.json. Missing verification never counts as completed verification.']
    (ROOT/'COVERAGE.md').write_text('\n'.join(md)+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
