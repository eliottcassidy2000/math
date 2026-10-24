#!/usr/bin/env python3
"""Verify a frozen producer packet, reproduce it in a real copy, and compare audit data."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,tempfile
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def execute(root,script,optimized=False,extra=()):
    p=subprocess.run(['python3']+(['-O'] if optimized else [])+[str(root/script),*map(str,extra)],cwd='/tmp',capture_output=True,text=True)
    need(p.returncode==0,'producer replay failed '+script+' '+p.stderr)
    return p

def main():
    ap=argparse.ArgumentParser();ap.add_argument('candidate');args=ap.parse_args();source=Path(args.candidate).resolve()
    initial={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    replay=[]
    execute(source,'verify_manifest.py');execute(source,'verify_manifest.py',True)
    with tempfile.TemporaryDirectory(prefix='session51-converse-candidate-realcopy-') as tmp:
        target=Path(tmp)/'packet';shutil.copytree(source,target,ignore=shutil.ignore_patterns('__pycache__'))
        need(not target.is_symlink(),'candidate copied directory not symlink')
        for script,receipt in [('check.py','normal.json'),('matrix_check.py','matrix-normal.json'),('series_check.py','series-normal.json'),('prime_check.py','prime-normal.json'),('port_check.py','port-normal.json')]:
            for optimized in (False,True):
                out=Path(tmp)/'replay.json';execute(target,script,optimized,('--output',out))
                need(out.read_bytes()==(source/receipt).read_bytes(),'producer real-copy receipt match '+script)
                replay.append({'script':script,'optimized':optimized,'status':'PASS'})
        execute(target,'verify_manifest.py');execute(target,'verify_manifest.py',True)
        # Manifest rejects an edited theorem even though all receipts are untouched.
        theorem=target/'THEOREM.md';theorem.write_bytes(theorem.read_bytes()+b'\n')
        bad=subprocess.run(['python3','-O',str(target/'verify_manifest.py')],cwd='/tmp',capture_output=True,text=True)
        need(bad.returncode!=0 and 'SHA256 mismatch' in bad.stderr,'producer theorem mutation rejected by manifest')
    producer=json.loads((source/'normal.json').read_text());audit=json.loads((HERE/'normal.json').read_text())
    for a,b in [('ordinary_positive_histogram','ordinary_spectrum'),('self_converse_positive_histogram','self_converse_spectrum'),('merged_positive_histogram','merged_spectrum')]:
        need(producer[a]==audit['order_ten'][b],'producer independent spectrum match '+a)
    # Align each actual representative by exact graph when available; full graph
    # isomorphism establishes equality for the differently built SCC paddings.
    import independent as G
    for r in producer['representatives']:
        a=tuple(r['rows']);matches=[s for s in audit['order_ten']['representatives'] if s['defect']==r['error'] and G.maps(a,tuple(s['rows']),True)]
        need(len(matches)==1 and bool(matches[0]['anti_automorphism_order'])==r['self_converse'],'producer representative aligns uniquely to independent slice')
    ports=json.loads((source/'port-normal.json').read_text());auditports=json.loads((HERE/'ports-normal.json').read_text())
    index={tuple(r['rows']):(r['hamiltonian_paths'],r['blue_ports'],r['black_ports']) for r in auditports['records']}
    need(all(index[tuple(r['rows'])]==(r['Hamiltonian_paths'],r['blue_ports'],r['black_ports']) for r in ports['representatives']),'every producer port count independent DP match')
    final={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    need(initial==final,'producer frozen source unchanged throughout replay')
    print(json.dumps({'status':'PASS','producer_manifest_sha256':initial['MANIFEST.json'],'producer_pins':initial,'literal_real_copy_replays':replay,'producer_theorem_mutation_rejected':True,'representative_alignment':160,'port_alignment':24,'source_unchanged':True},indent=2,sort_keys=True))
if __name__=='__main__':main()
