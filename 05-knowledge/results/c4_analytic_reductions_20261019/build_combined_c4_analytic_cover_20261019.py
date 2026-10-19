#!/usr/bin/env python3
"""Build exact, separately scoped analytic coverage without rerunning prior controls."""
from collections import Counter
from pathlib import Path
import gzip,hashlib,importlib.util,json
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('defects',ROOT/'verify_one_defect_coverage_20261019.py')
D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D)


def main():
    old_path=ROOT/'c4_one_defect_analytic_coverage_20261019.json.gz'
    controls=json.loads((ROOT/'one_defect_coverage_controls_20261019.json').read_text())
    D.B.check(hashlib.sha256(old_path.read_bytes()).hexdigest()==controls['coverage_sha256'],'Prior coverage hash mismatch')
    D.B.check(hashlib.sha256((ROOT/'verify_one_defect_coverage_20261019.py').read_bytes()).hexdigest()==controls['verifier_sha256'],'Prior verifier changed')
    D.B.check(hashlib.sha256(D.BASE.read_bytes()).hexdigest()==controls['dependency_sha256'],'Prior dependency changed')
    prior=json.loads(gzip.decompress(old_path.read_bytes()))
    old={r['core_index']:r for r in prior['certificates']}
    stats=Counter();records=[];catalogue=hashlib.sha256();residual=[]
    for pred in D.B.c4_cores():
        succ=D.B.successors(pred);degrees=[8-(pred[v]|succ[v]).bit_count() for v in range(9)]
        if max(degrees)>7:continue
        index=stats['degree_valid_cores'];stats['degree_valid_cores']+=1
        catalogue.update((json.dumps(pred,separators=(',',':'))+'\n').encode())
        span=[]
        for z in range(3,7):
            for u,v in ((0,z),(z,1)):
                delta=D.directed_defect(pred,succ,u,v);slack=7-degrees[z]
                if delta+slack<=1:
                    span.append({'pair':[u,v],'common_vertex':z,'core_directed_defect':delta,'exterior_slack_bound':slack})
        stats['ambient_balanced_pair_cores']+=index in old
        stats['max_span_directed_cores']+=bool(span)
        if index in old:
            record=old[index]
            D.B.check(tuple(record['predecessor_masks'])==pred,'Prior catalogue mismatch')
            record={**record,'scope':'AMBIENT_BALANCED_PAIR','certificate_type':'two_sided_defect_slack'}
        elif span:
            record={'core_index':index,'predecessor_masks':pred,'core_degrees':degrees,
                    'scope':'MAXIMUM_SPAN_MAJORITY_CONTRADICTION','certificate_type':'one_sided_defect_slack',
                    'witness':span[0]}
        else:
            residual.append(index);continue
        records.append(record)
    stats['union_analytic_cores']=len(records);stats['residual_cores']=len(residual)
    expected={'degree_valid_cores':79121,'ambient_balanced_pair_cores':49592,
              'max_span_directed_cores':54235,'union_analytic_cores':56630,'residual_cores':22491}
    D.B.check(dict(stats)==expected,'Combined coverage count mismatch')
    D.B.check(catalogue.hexdigest()==controls['degree_valid_catalogue_sha256'],'Catalogue hash mismatch')
    coverage={'scope':'Union of explicitly scoped analytic c=4 exclusions; no full c=4 theorem yet.',
              'degree_valid_catalogue_sha256':catalogue.hexdigest(),'counts':dict(stats),
              'certificates':records,'residual_core_indices':residual}
    data=gzip.compress((json.dumps(coverage,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0)
    target=ROOT/'c4_combined_analytic_coverage_20261019.json.gz';target.write_bytes(data)
    result={'status':'PASS','counts':dict(stats),'coverage_file':target.name,
            'coverage_sha256':hashlib.sha256(data).hexdigest(),'degree_valid_catalogue_sha256':catalogue.hexdigest(),
            'previous_coverage_sha256':controls['coverage_sha256'],
            'builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'combined_analytic_coverage_controls_20261019.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
