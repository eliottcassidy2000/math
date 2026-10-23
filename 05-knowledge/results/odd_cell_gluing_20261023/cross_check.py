"""Reconcile concurrent actual order31 work by marked cases and literal cells."""
from pathlib import Path
import argparse,gzip,hashlib,json,subprocess
import check_cells as C
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
REF='d45865c332'
PATH='05-knowledge/results/endpoint_overlap_ports_20261023/receipt.json.gz'
def main():
    raw=gzip.decompress(subprocess.check_output(['git','show',REF+':'+PATH],cwd=ROOT))
    theirs=json.loads(raw);ours=json.loads((HERE/'cells.json').read_text())
    by={tuple(r['marks']):r for r in ours['nondiscrete_DDP7_records']}
    A=C.double(C.double(C.paley(7)))
    C.need(theirs['matrix_outmasks']==[sum(x<<j for j,x in enumerate(row)) for row in A],'source gauge differs')
    expected={(3,1):'triple_negative_closed',(3,3):'triple_C3',(7,1):'seven_negative_closed',(7,3):'seven_to_C3',(7,21):'seven_C21'}
    repaired=0
    for r in theirs['cases']:
        marks=tuple(r['marks']);cells=C.refine(A,marks)
        C.need(sorted(sum(1<<v for v in D) for D in cells)==sorted(r['initial_cells']),'initial literal cells differ')
        if marks not in by:
            C.need(r['class']=='already_discrete','missing nondiscrete case');continue
        o=by[marks]
        C.need(expected[o['largest_cell'],o['automorphisms']]==r['class'],'closure differs')
        C.need(o['minimum_extra_root_marks']==r['minimal_extra_actual_owners'],'minimum owner count differs')
        if o['largest_cell']==7:
            cc,_=C.motif_repair(A,cells)
            C.need(sorted(sum(1<<v for v in D) for D in cc)==sorted(r['after_overlap_cells']),'flag and subset-equality partitions differ')
            repaired+=1
    C.need(len(theirs['cases'])==930 and repaired==132,'incomplete universe')
    return {'status':'FINITE-EXACT independent incoming reconciliation','incoming_commit':REF,
            'incoming_receipt_sha256':hashlib.sha256(raw).hexdigest(),'ordered_mark_cases_agree':930,
            'small_flag_vs_subset_equality_literal_partitions_agree':132,'checks':C.CHECKS}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=main();f=HERE/'cross_check.json'
    if a.write:f.write_text(json.dumps(out,indent=2)+'\n')
    else:C.need(json.loads(f.read_text())==out,'receipt mismatch')
    print(json.dumps(out,indent=2))
