"""Verify pinned inputs and normal/-O receipts without network access."""
from pathlib import Path
import sys,subprocess,json,hashlib,platform
import sympy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
DEPENDENCIES=[
 '05-knowledge/results/session160_ordered_programs_20261103/axioms/NOTE.md',
 '05-knowledge/results/session160_ordered_programs_20261103/priority_tilings/NOTE.md',
 '05-knowledge/results/session160_ordered_programs_20261103/metagraph/NOTE.md',
 '05-knowledge/results/session160_ordered_programs_20261103/adaptive_codes/NOTE.md',
 '05-knowledge/results/apex_wedge_constructor_20261103/mechanisms.py',
 '05-knowledge/results/apex_wedge_constructor_20261103/profile.py',
 '05-knowledge/results/apex_wedge_constructor_20261103/THEORY.md',
 '05-knowledge/results/apex_wedge_constructor_20261103/PROFILE.md',
 '05-knowledge/results/apex_wedge_constructor_20261103/CARRIER.md',
 '05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md',
 '05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md',
 '05-knowledge/results/session156_spectrum_20261103/states/witnesses.json',
 '05-knowledge/results/polarized_path_calculus_20261103/CALCULUS.md',
 '05-knowledge/results/session36_information_curl_contexts_20261023/root/CURL_COCENTER_AND_LIE_BOUNDARY.md',
 '05-knowledge/results/session81_counterexample_families_20261028/LIFTING.md',
 '07-reflections/tournament-subset-decks-higher-order-faces-and-kneser-middle-duality-codex-20260821.md',
 '04-computation/tournament_subset_deck_hierarchy_scout_codex_20260821.py',
 '05-knowledge/results/tournament_subset_deck_hierarchy_scout_codex_20260821.out',
 '05-knowledge/results/black_descent_20261020c/snarks.md',
 '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
 '05-knowledge/results/mixed_returns_20261103/controllers.json',
]

def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    source_file=HERE/'sources.json'
    if '--pin-sources' in sys.argv:
        entries=[{'path':p,'sha256':digest(ROOT/p)} for p in DEPENDENCIES]
        prior=json.loads((ROOT/DEPENDENCIES[-1]).read_text())
        data={'status':'INHERITED_PINNED_NOT_RETESTED','base_commit':'34fc4d716863','hash_scheme':'sha256 of UTF8 text with LF newlines',
              'dependencies':entries,'hostile81_27':prior['hostile81_27'],
              'one_birth55_orientations':prior['one_birth55_orientations']}
        source_file.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
        print('Pinned',len(entries),'inherited source files');return
    pins=json.loads(source_file.read_text())
    for entry in pins['dependencies']:
        need(digest(ROOT/entry['path'])==entry['sha256'],'changed dependency: '+entry['path'])
    records=[]
    for name in ('normal_forms','bank','controllers','neutral_carries'):
        receipts=[];outputs=[]
        for optimized in (False,True):
            cmd=[sys.executable,'-X','utf8','-B']+(['-O'] if optimized else [])+[str(HERE/(name+'.py'))]
            result=subprocess.run(cmd,cwd=ROOT,check=True,capture_output=True,text=True,encoding='utf8')
            receipts.append(digest(HERE/(name+'.json')));outputs.append(result.stdout.strip())
        need(receipts[0]==receipts[1] and outputs[0]==outputs[1],'normal/-O discrepancy: '+name)
        records.append({'script':name+'.py','script_sha256':digest(HERE/(name+'.py')),
                        'receipt':name+'.json','sha256':receipts[0],'normal_equals_optimized':True,'output':outputs[0]})
        print(outputs[0])
    data={'status':'PASS','python':platform.python_version(),'sympy':sympy.__version__,
          'dependency_count':len(pins['dependencies']),'source_manifest_sha256':digest(source_file),'runs':records}
    (HERE/'reproduction.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: all pinned inputs; eight executions; byte-identical normal/-O receipts')
if __name__=='__main__':main()
