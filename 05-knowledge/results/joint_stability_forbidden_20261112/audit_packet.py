"""Close-out consistency, replay, source and dependency manifest.

This checks existing receipts; it does not replace their original programs.
One historical scope spelling difference is allowed only for the order9
two-block replay, after requiring identical complete mathematical payloads.
"""
from pathlib import Path
import ast,gzip,hashlib,importlib,json,re,sys
from collections import Counter
from census_nine import read,save,need
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]

def sha(path):
    raw=path.read_bytes()
    if path.suffix in ['.py','.cpp','.md','.json','.txt']:raw=raw.replace(b'\r\n',b'\n')
    return hashlib.sha256(raw).hexdigest()

def run():
    replays=[]
    for optimized in sorted(HERE.glob('*.optimized.json.gz')):
        if optimized.name.startswith('packet_audit'):continue
        normal=HERE/optimized.name.replace('.optimized','');need(normal.exists(),'normal receipt for every replay')
        a=json.loads(gzip.decompress(normal.read_bytes()));b=json.loads(gzip.decompress(optimized.read_bytes()));kind='byte-identical'
        if normal.read_bytes()!=optimized.read_bytes():
            need(normal.name=='two_block_probe.json.gz','only the documented metadata wording difference is allowed')
            need({k:v for k,v in a.items() if k!='scope'}=={k:v for k,v in b.items() if k!='scope'},'entire two-block mathematical payload identical')
            expected=a['scope'].replace('order-nine','order-9').replace('order eleven','order 11')
            need(expected==b['scope'],'exact documented scope-only wording change');kind='identical mathematical payload; historical scope spelling only'
        replays.append({'normal':normal.name,'optimized':optimized.name,'comparison':kind,'normal_sha256':sha(normal),'optimized_sha256':sha(optimized)})
    need(len(replays)>=37,'all finished normal/optimized receipts present')
    u=read('unit_ten.json.gz');a=read('audit_unit_slopes.json.gz');need(a['sources']==len(u['records'])==465023 and a['marked_interfaces']==u['strong_unit_orbits']==478684,'complete per-record independent replay')
    need(a['sum_U']==sum(m['U'] for r in u['records'] for m in r['marks']),'entire U-bank checksum')
    stable=read('stable_ten_selected.json.gz');need((stable['births'],stable['stable_states'],stable['descendants'])==(850,1253,403),'complete stable classification')
    family=read('birth_union_ten.json.gz');need(family['births']==1150 and family['counts_n1_to40'][10]==1487,'family/DAG overlap count')
    need(Counter(len(f) for r in family['records'] for f in r['faces'])==Counter({1:967,2:145,3:38,4:2}),'physical height-bound multiplicities')
    d=read('descendants_eleven.json.gz');need(d['marked_descendants']==1487 and d['edges_to_order_ten']==1597 and d['terminal_birth_count_histogram']=={'1':1486,'2':1},'complete eleven descendant classification')
    band=read('coefficient_band.json.gz');need(band['covered_odd_prefix_through']==53557 and band['lower_heads']==3954 and band['head_omissions']==[7,21],'complete finite coefficient band')
    closed=read('completed_bulk_heads.json.gz');need(closed['fully_closed_laws']==540 and closed['resolved_head_obligations']==786 and not closed['unresolved'],'complete binary-library closure')
    need(read('higher_odd_families.json.gz')['complete_odd_slope_families']==320,'higher-mode parity families')
    need(read('audit_affine_births.json.gz')['full_kernels']==1379,'every full U skeleton independently audited')
    need(read('extreme_ten.json.gz')['distinct_marked_crossing_quotients']==0,'extreme-gate scoped negative result')
    own_sources=sorted([*HERE.glob('*.py'),*HERE.glob('*.cpp')])
    for p in HERE.glob('*.py'):ast.parse(p.read_text(encoding='utf8'),filename=str(p))
    modules=['joint_nine','stable_census','birth_union','coefficient_band','complete_rays','natural_recursion','head_repair','rank_raise','reflection_crossings','extreme_ten','bulk_head_interfaces','coefficient_intersections','library_boundary']
    for name in modules:importlib.import_module(name)
    external=set()
    for m in list(sys.modules.values()):
        f=getattr(m,'__file__',None)
        if not f:continue
        p=Path(f).resolve()
        if p.is_file() and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):external.add(p)
    for rel in ['session288_twin_lengths_20261112/stable/cones.py','session251_creative_completion_20261111/ternary/compiler.py']:external.add(HERE.parent/rel)
    input_names=['blue_quotient_20261020b/n9_classes.tsv.gz','metagraph_aspect_sequences_20261111/records.jsonl.gz','five_stress_connection_atlas_20261112/cyclic_port_receipt.json','five_stress_connection_atlas_20261112/contraction_dag.json','endpoint_coverage_compiler_20261111/wheel.json','path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt']
    inputs=[HERE.parent/p for p in input_names]
    for p in inputs:need(p.is_file(),'inherited input remains available')
    for doc in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text(encoding='utf8')):
            if '://' in target or target.startswith('#'):continue
            need((doc.parent/target.split('#')[0]).exists(),'local documentation link '+target)
    def entries(paths):return [{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256_lf_for_text':sha(p)} for p in sorted(paths)]
    outputs=[p for p in HERE.glob('*.json.gz') if not p.name.startswith('packet_audit')]
    result={'scope':'Close-out consistency and source manifest for this packet. Same-agent mathematical review with independent algorithms, not independent-agent peer review. Checks do not assert global Chan-Pak completeness.','status':'PASS','replay_pairs':replays,'own_sources':entries(own_sources),'loaded_repository_dependencies':entries(external),'inherited_data':entries(inputs),'output_receipts':entries(outputs),'external_catalogue':{'url':u['source_url'],'sha256':u['source_sha256'],'ordinary_classes':u['ordinary_classes']},'optional_figure':{'path':'crossings.png','sha256':sha(HERE/'crossings.png'),'visually_inspected':True}}
    save('packet_audit.json.gz',result);print(json.dumps({'status':'PASS','replay_pairs':len(replays),'own_source_files':len(own_sources),'loaded_repository_dependencies':len(external),'inherited_data_files':len(inputs),'output_receipts':len(outputs)}),flush=True)
if __name__=='__main__':run()
