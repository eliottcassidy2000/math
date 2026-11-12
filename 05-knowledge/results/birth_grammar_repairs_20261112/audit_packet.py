"""Deterministic receipt audit and pinned input/source manifest; no implicit asserts."""
from common import *
import importlib,re,subprocess
from urllib.parse import unquote
from collections import Counter
ROOT=HERE.parents[2]
BASELINE='47de5a3f10fd'

def digest(path):
    raw=path.read_bytes()
    if path.suffix in {'.py','.md','.json','.svg'}:raw=raw.replace(b'\r\n',b'\n')
    return hashlib.sha256(raw).hexdigest()

def run():
    normals=[p for p in sorted(HERE.glob('*.json.gz')) if '.optimized.' not in p.name and p.name!='audit_packet.json.gz']
    need(len(normals)==33,'exact declared33 receipt pairs')
    pairs=[]
    for p in normals:
        q=p.with_name(p.name.replace('.json.gz','.optimized.json.gz'))
        need(q.exists() and p.read_bytes()==q.read_bytes(),'normal/optimized byte equality: '+p.name)
        d=read(p.name);need(bool(d.get('scope')),'explicit receipt universe: '+p.name)
        if 'complete' in d:need(d['complete'] is True,'completed checkpoint: '+p.name)
        pairs.append({'receipt':p.name,'sha256':digest(p),'compressed_bytes':p.stat().st_size})
    v=read('even_valuations.json.gz')
    need(v['classification']=={'bounded_constant':352,'bounded_parity':60,'unbounded':66},'all478 valuation partition')
    need(len(v['records'])==len({(r['birth_id'],r['port']) for r in v['records']})==478,'distinct inherited axes')
    need(read('repaired_rays.json.gz')['closed_families']==5,'five complete physical repairs')
    need(len(read('repaired_rays.json.gz')['literal_whole_graph_controls'])==164,'whole graph repair controls')
    need(read('mixed_family_probe.json.gz')['complete'] is True,'complete mixed discovery')
    m=read('mixed_grammar.json.gz')
    need(m['families']==46 and m['counts']=={'13':46,'15':46,'17':46},'mixed birth grammar counts')
    need(all(t['intersection_m']==0 for t in m['intersections']) and len(m['intersections'])==20,'all positive mixed charts disjoint')
    need(all(t['sign']==1 and t['Stirling_coefficients'][0]>0 and min(t['Stirling_coefficients'])>=0 for r in m['records'] for t in r['frozen_nonvanishing']),'positive frozen certificates')
    need(read('synchrony_grammar.json.gz')['families']==24,'even birth grammar count')
    b=read('birth_inventory.json.gz')
    need(b['retained_births_through14']==4042 and b['birth_counts_by_order']=={'3':1,'5':4,'6':8,'7':25,'8':60,'9':202,'10':850,'11':464,'12':858,'13':1502,'14':68},'scoped retained union')
    need(b['order13_plain_mixed_intersection']==0,'disjoint order13 routes')
    need(read('hamilton_clocks_9.json.gz')['distinct_profiles']==784832,'complete pointed profile inventory')
    need(read('repair_obstructions.json.gz')['affine_TT_barrier_families']==126 and read('repair_obstructions.json.gz')['periodic_barrier_families']==470,'distinct repair barriers')
    need(read('audit_ozj.json.gz')['periodic_root_barrier_families']==478 and read('audit_ozj.json.gz')['semigroup_states']==25,'complete modular parity barrier')
    nine=read('periodic_nine.json.gz')
    need(sum(nine['period_counts'].values())==278915 and not nine['joint_eigenvalue_one_nine_periods'],'entire norm-bounded65-period screen')
    need([r['H'] for r in nine['physical_finite_controls']]==[(9**m+1)//2 for m in range(7)] and max(len(r['rows']) for r in nine['physical_finite_controls'])<=14,'formal-family finite witnesses only')
    need(read('lifting_roots.json.gz')['root_histogram']=={'irrational_in_Z2':65},'two-mode rationality test')
    need(read('ordinal_chart_audit.json.gz')['statistics']=={'slot_identities':156,'nonempty_identities':156,'physical_words':39},'independent ordered-module controls')
    links=[]
    for p in sorted(HERE.glob('*.md')):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text(encoding='utf8')):
            target=unquote(target.strip('<>').split('#')[0])
            if not target or re.match(r'[a-zA-Z]+://',target):continue
            q=(p.parent/target).resolve();need(q.exists(),'local link: '+p.name+' -> '+target)
            links.append([p.name,target])
    # Import the inherited compiler chain, without running its producers.
    for name in ['hamilton_clocks','repaired_rays','mixed_grammar','audit_math','periodic_nine','construct','two_block_probe']:
        importlib.import_module(name)
    inherited=set()
    for module in list(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if not name:continue
        p=Path(name).resolve()
        if p.suffix=='.py' and p.is_relative_to(ROOT) and not p.is_relative_to(HERE):inherited.add(p)
    inherited.update(OLD/name for name in ['affine_births.json.gz','birth_union_ten.json.gz','stable_nine.json.gz','stable_ten_selected.json.gz','unit_ten.json.gz','coefficient_band.json.gz'])
    inherited.update(HERE.parent/name for name in ['five_stress_connection_atlas_20261112/cyclic_port_receipt.json','metagraph_aspect_sequences_20261111/records.jsonl.gz','blue_quotient_20261020b/n9_classes.tsv.gz'])
    rels=sorted(p.relative_to(ROOT).as_posix() for p in inherited)
    check=subprocess.run(['git','diff','--quiet',BASELINE,'--',*rels],cwd=ROOT)
    need(check.returncode==0,'pinned inherited inputs unchanged from baseline')
    inputs={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(inherited)}
    packet={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.suffix in {'.py','.md','.png','.svg','.json'}}
    out={'scope':'Deterministic33-pair receipt audit, decisive cross-receipt identities, explicit completion flags, local Markdown links, packet source/figure hashes, and inherited inputs pinned to the stated Git baseline. This validates the certificates, not a replacement proof of all mathematical assertions.','baseline':BASELINE,'receipt_pairs':pairs,'normal_optimized_pairs':len(pairs),'classification':v['classification'],'complete_repairs':5,'unclosed_old_families':473,'infinite_birth_families':70,'retained_births_through14':4042,'local_links_checked':len(links),'pinned_inherited_inputs':inputs,'packet_sources_and_documents':packet}
    save('audit_packet.json.gz',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['receipt_pairs','pinned_inherited_inputs','packet_sources_and_documents']}),flush=True)

if __name__=='__main__':run()
