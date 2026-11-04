"""Frozen verification. Replay commands never update their expectations."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ENGINES=(
    ('algebra/check.py',[],'algebra/normal.json','algebra/optimized.json'),
    ('geometry/verify.py',[],'geometry/receipt.json','geometry/receipt.optimized.json'),
    ('audit/parity_check.py',[],'audit/parity.json','audit/parity.optimized.json'),
    ('audit/self_converse_check.py',[],'audit/self_converse.json','audit/self_converse.optimized.json'),
    ('audit/endpoint_scout.py',[],'audit/endpoints.json','audit/endpoints.optimized.json'),
    ('response/run.py',[],'response/tensors.json','response/tensors.optimized.json'),
    ('census/run.py',['--limit','100000'],'census/prefix100000.json','census/prefix100000.optimized.json'),
)


def need(b,message):
    if not b:raise RuntimeError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):return json.loads((HERE/path).read_text(encoding='utf8'))


def pins():
    manifest=load('MANIFEST.json')['files']
    actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file()
            and p.name!='MANIFEST.json' and '__pycache__' not in p.parts and p.suffix!='.pyc'}
    need(set(manifest)==actual,'artifact inventory mismatch')
    for path,digest in manifest.items():need(sha(HERE/path)==digest,'artifact changed: '+path)
    sources=load('SOURCE_PINS.json')['files']
    for path,digest in sources.items():need(sha(ROOT/path)==digest,'inherited source changed: '+path)
    return len(manifest),len(sources)


def source_receipts():
    for script,_,normal,_ in ENGINES:
        a=load(normal);values=a['source_sha256']
        if isinstance(values,str):need(sha(HERE/script)==values,'receipt source: '+script)
        elif script=='response/run.py':
            for name,digest in values.items():need(sha(HERE/name)==digest,'tensor source: '+name)
        else:
            for name,digest in values.items():need(sha((HERE/script).parent/name)==digest,'receipt source: '+name)
    c=load('census/full10.json')
    for name,digest in c['source_sha256'].items():need(sha(HERE/'census'/name)==digest,'full census source: '+name)


def main():
    need(sys.argv[1:] in ([],['--replay'],['--full-census']),'usage: verify.py [--replay|--full-census]')
    counts=pins();source_receipts()
    expected={script:(HERE/normal).read_text(encoding='utf8') for script,_,normal,_ in ENGINES}
    complete=(HERE/'census/full10.json').read_text(encoding='utf8')
    for script,_,_,opt in ENGINES:need(expected[script]==(HERE/opt).read_text(encoding='utf8'),'normal/optimized difference: '+script)
    c=json.loads(complete)
    need(c['complete_bank'] and c['classes']==9733056 and not c['collisions'],'complete census consequence')
    need(c['sidecar_class_collisions']==c['sidecar_automorphism_collisions']==0,'sidecar collision')
    need(c['sidecar_labelled_mass']==2**45,'labelled mass')
    old=json.loads((ROOT/'05-knowledge/results/consecutive_flags_order10_20261028/verification10.json').read_text(encoding='utf8'))
    need(c['sidecar_automorphism_histogram']==old['automorphism_order_histogram'],'symmetry histogram')
    need(c['classes_with_endpoint_ties']<=c['full_cover_canonicalizations']<=c['classes'],'sieve scope')
    p=json.loads(expected['census/run.py']);need(not p['complete_bank'] and p['classes']==100000,'bounded control scope')
    t=json.loads(expected['response/run.py']);need(len(t['controls'])==8,'tensor universe')
    need(sum(x['queries'] for x in t['controls'])==19400,'actual response query count')
    for x in t['controls']:
        n=x['n']
        for name,power in [('S',1),('E',1),('M',2),('C',2),('U',3),('V',3),('D',4)]:
            need(len(x[name])==n**power,'full tensor inventory')
        need(sum(x['S'])==sum(x['E'])==sum(x['M']),'endpoint totals')
        need(sum(x['C'])==sum(x['U'])==sum(x['V']),'joint two-path totals')
        need(sum(x['D'])%6==0,'ordered three-path normalization')
    need(t['controls'][-1]['sidecar_aut']==81,'high symmetry control')
    g=json.loads(expected['geometry/verify.py'])
    need(g['spindle_external_triple_centres']==0 and g['spindle_maximum_external_contacts']==2,'spindle contact certificate')
    need(g['bank_assemblies']==256 and g['boundary_partitions']==203 and g['blocked_q5_extensions']==0,'full boundary scope')
    need(g['support_counts']['free_angle']=={'5':200,'6':203},'genuine conditional zeros')
    replayed=0
    if sys.argv[1:]==['--replay']:
        for script,args,_,_ in ENGINES:
            path=HERE/script
            for flags in ([],['-O']):
                r=subprocess.run([sys.executable,'-B']+flags+[str(path)]+args,cwd=path.parent,
                                 check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
                need(r.stdout==expected[script],'frozen replay differs: '+script)
                replayed+=1;print('PASS '+script+(' -O' if flags else ''),flush=True)
    if sys.argv[1:]==['--full-census']:
        r=subprocess.run([sys.executable,'-B',str(HERE/'census/run.py')],cwd=HERE/'census',
                         check=True,stdout=subprocess.PIPE,text=True,encoding='utf8')
        need(r.stdout==complete,'full census frozen replay differs');replayed=1
    need(pins()==counts,'files changed during verification')
    print(json.dumps({'status':'PASS','artifact_files':counts[0],'inherited_source_files':counts[1],
                      'matching_receipt_pairs':len(ENGINES),'complete_census_classes':c['classes'],
                      'replayed_engines':replayed,'expected_receipts_modified':False},sort_keys=True))


if __name__=='__main__':main()
