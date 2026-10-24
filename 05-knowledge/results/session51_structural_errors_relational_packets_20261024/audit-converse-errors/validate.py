#!/usr/bin/env python3
"""Relocation and hostile mutation controls for the independent audit packet."""
from pathlib import Path
import hashlib,json,shutil,subprocess,tempfile
HERE=Path(__file__).resolve().parent

def need(x,msg):
    if not x:raise RuntimeError(msg)
def invoke(path, optimized=False):
    return subprocess.run(['python3']+(['-O'] if optimized else [])+[str(path/'independent.py')],cwd='/tmp',text=True,capture_output=True)
def run():
    baseline=(HERE/'normal.json').read_bytes();result={}
    need(baseline==(HERE/'optimized.json').read_bytes(),'normal optimized byte identity')
    with tempfile.TemporaryDirectory(prefix='session51-converse-audit-realcopy-') as tmp:
        target=Path(tmp)/'packet';shutil.copytree(HERE,target,ignore=shutil.ignore_patterns('__pycache__'))
        need(not target.is_symlink() and all(not p.is_symlink() for p in target.rglob('*')),'literal relocation not symlink')
        for optimized in (False,True):
            p=invoke(target,optimized);need(p.returncode==0 and p.stdout.encode()==baseline,'real relocated replay exact match')
        result['real_copy_normal_and_optimized']='PASS'
        source=target/'sources/order-ten-audit-receipt.json';source.write_bytes(source.read_bytes()+b' ')
        p=invoke(target);need(p.returncode!=0 and 'inherited source pin' in p.stderr,'tampered inherited source rejected')
        result['tampered_source_rejected']='PASS'
        shutil.copyfile(HERE/'sources/order-ten-audit-receipt.json',source)
        script=(target/'independent.py').read_text()
        wrong=script.replace('m[-j%5]>>(-i%5)','m[j]>>i')
        need(wrong!=script,'anti action mutation applied')
        (target/'independent.py').write_text(wrong)
        p=invoke(target,True)
        need(p.returncode!=0 and ('matrix converse agrees full isomorphism' in p.stderr or 'affine anti fixed set formula' in p.stderr),'missing negations in converse action rejected')
        result['wrong_converse_action_rejected_under_O']='PASS'
    j=json.loads(baseline);ordinary=dict(j['order_ten']['ordinary_spectrum']);fixed=dict(j['order_ten']['self_converse_spectrum']);merged=dict(j['order_ten']['merged_spectrum'])
    need({e:c//2 for e,c in ordinary.items()}!=merged,'halving ordinary count without fixed correction rejected')
    need(any(r['defect']==40 and not r['anti_automorphism_order'] for r in j['order_ten']['representatives']),'paired stable-cell sizes not sufficient control')
    need(j['checks']['missing t square mutation rejected']==1,'weighted pair exponent hostile executed')
    result.update({'naive_halving_rejected':'PASS','stable_cell_pairing_not_sufficient':'PASS','missing_t_square_rejected':'PASS','normal_sha256':hashlib.sha256(baseline).hexdigest(),'status':'PASS'})
    return result
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
