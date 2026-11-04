"""Frozen artifact verification; replay never rewrites its expectations."""
import hashlib,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ENGINES=(
 ('algebra/check.py',[], 'algebra/receipt.json','algebra/receipt.optimized.json'),
 ('geometry/verify.py',[], 'geometry/receipt.json','geometry/receipt.optimized.json'),
 ('lists/run.py',[], 'lists/normal.json','lists/optimized.json'),
 ('lists/core_audit.py',[], 'lists/core_audit.json','lists/core_audit.optimized.json'),
 ('recovery/check.py',[], 'recovery/normal.json','recovery/optimized.json'),
 ('recovery/transfer_check.py',[], 'recovery/transfer.json','recovery/transfer.optimized.json'),
 ('response/run.py',[], 'response/normal.json','response/optimized.json'),
 ('self_converse/run.py',['--max-half','4'], 'self_converse/through9.json','self_converse/through9.optimized.json'),
 ('self_converse/bank_audit.py',[], 'self_converse/bank_audit.json','self_converse/bank_audit.optimized.json'),
 ('residual/run.py',['--limit','10000'], 'residual/prefix10000.json','residual/prefix10000.optimized.json'),
 ('residual/audit.py',[], 'residual/audit.json','residual/audit.optimized.json'),
)
def need(b,message):
    if not b:raise RuntimeError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads((HERE/path).read_text(encoding='utf8'))
def pins():
    manifest=load('MANIFEST.json')['files']
    actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
    need(set(manifest)==actual,'artifact inventory')
    for name,digest in manifest.items():need(sha(HERE/name)==digest,'artifact changed: '+name)
    source=load('SOURCE_PINS.json')['files']
    for name,digest in source.items():need(sha(ROOT/name)==digest,'inherited source changed: '+name)
    return len(manifest),len(source)
def receipt_sources(script,receipt):
    hashes=receipt['source_sha256']
    if isinstance(hashes,str):need(sha(HERE/script)==hashes,'receipt source: '+script);return
    for name,digest in hashes.items():
        candidates={p.resolve() for p in (ROOT/name,HERE/name,(HERE/script).parent/name) if p.is_file()}
        need(len(candidates)==1,'source path ambiguous/missing: '+name)
        need(sha(candidates.pop())==digest,'receipt source changed: '+name)
def main():
    need(sys.argv[1:] in ([],['--replay'],['--full-censuses']),'usage: verify.py [--replay|--full-censuses]')
    counts=pins();expected={}
    for script,args,normal,optimized in ENGINES:
        expected[script]=load(normal);need(expected[script]==load(optimized),'normal/optimized mismatch: '+script)
        receipt_sources(script,expected[script])
    full11=load('self_converse/full11.json');receipt_sources('self_converse/run.py',full11)
    need([r['classes'] for r in full11['orders']]==[2,8,88,2752,279968],'complete self-converse universe')
    for r in full11['orders']:
        need(r['classes']==r['sidecar_classes'] and not r['collisions'],'self-converse collision')
        need(r['class_collisions']==r['extra_sidecar_symmetries']==0,'extra self-converse symmetry')
        h=(r['n']-1)//2;need(r['source_labelled_mass']==2**(h*h),'self-converse labelled mass')
    need(full11['orders'][:4]==expected['self_converse/run.py']['orders'],'bounded/full family mismatch')
    need(sha(ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives.json')==full11['bank_sha256'],'self-converse source bank')
    full10=load('residual/full10.json');receipt_sources('residual/run.py',full10)
    need(full10['classes']==9733056 and full10['complete_bank'],'complete modular-sieve universe')
    need(full10['endpoint_parity_balanced']==12670,'complete endpoint-parity count')
    need(not full10['combined_residual_addresses'],'combined modular survivor')
    t=expected['response/run.py'];need(sum(x['queries'] for x in t['controls'])==8996,'actual query inventory')
    need([x['sidecar_aut'] for x in t['controls']]==[1,55],'automorphism controls')
    for x in t['controls']:
        n=x['n'];need(n==11 and x['relabelled_controls']==64,'tensor control universe')
        for name,power in [('S',1),('E',1),('M',2),('C',2),('U',3),('V',3),('D',4)]:need(len(x[name])==n**power,'full named tensor: '+name)
        need(sum(x['S'])==sum(x['E'])==sum(x['M']),'Hamilton endpoint sums')
        need(sum(x['C'])==sum(x['U'])==sum(x['V']),'two-path joint sums')
        need(sum(x['D'])%6==0,'three-path order factor')
    g=expected['geometry/verify.py'];need(g['proper_fivecolour_boundary_patterns']==42 and g['required_contacts']==14 and g['maximum_capacity']==8,'full geometry boundary')
    need(g['connected_control_counts']['5']['supported_boundary_partitions']==42,'actual interacting control')
    c=expected['lists/core_audit.py'];need(c['all_graphs']==2**21 and c['metric_filter_graphs']==8190 and c['eligible_set_equals_orbit_union'],'independent exhaustive core universe')
    need(expected['residual/audit.py']['integer_endpoint_comparisons']==2123,'parity audit universe')
    replayed=0
    if sys.argv[1:]==['--replay']:
        for script,args,_,_ in ENGINES:
            path=HERE/script
            for flags in ([],['-O']):
                out=subprocess.run([sys.executable,'-B']+flags+[str(path)]+args,cwd=path.parent,check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
                need(json.loads(out.stdout)==expected[script],'frozen replay differs: '+script);replayed+=1
                print('PASS '+script+(' -O' if flags else ''),flush=True)
    if sys.argv[1:]==['--full-censuses']:
        for lane,value in [('self_converse',full11),('residual',full10)]:
            out=subprocess.run([sys.executable,'-B',str(HERE/lane/'run.py')],check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
            need(json.loads(out.stdout)==value,'full census replay differs: '+lane);replayed+=1
    need(pins()==counts,'files changed during verification')
    print(json.dumps({'status':'PASS','artifact_files':counts[0],'inherited_source_files':counts[1],'matching_receipt_pairs':len(ENGINES),'replayed_engines':replayed,'complete_self_converse11_classes':279968,'complete_order10_sieve_classes':9733056,'expected_receipts_modified':False},sort_keys=True))
if __name__=='__main__':main()
