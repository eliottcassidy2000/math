#!/usr/bin/env python3
"""Real-copy replays and adversarial validation, with checks active under -O."""
from pathlib import Path
import argparse,copy,hashlib,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
EXPECTED='cd89cad144786c67f076c235409c90e3b4de21e0a5f81926139184068f1d676d'

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def inventory(root):return {str(p.relative_to(root)):digest(p) for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}

def run(script,optimized=False,args=(),success=True):
    p=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(script),*map(str,args)],cwd='/tmp',capture_output=True,text=True)
    need((p.returncode==0)==success,'unexpected subprocess result: '+str(script)+'\n'+p.stderr[-2000:])
    return p

def compare(a,b):
    for key in ('positive_nodes','blue_ports','induced_blue_edges','induced_blue_loop_edges','induced_nonloop_blue_edges','internal_blue_ports','boundary_blue_edges','adjacency_multiplicities','induced_cycle_rank','bridge_edge_indices','articulation_node_indices','induced_degree_histogram','endpoint_transport_histogram','adjacency_multiplicity_histogram','components','nodes','internal_edges','boundary_edges'):
        need(a[key]==b[key],'independent mismatch: '+key)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('candidate',type=Path);args=ap.parse_args();candidate=args.candidate.resolve()
    need(digest(candidate/'MANIFEST.json')==EXPECTED,'exact frozen candidate manifest')
    before=inventory(candidate); audit=json.loads((HERE/'normal.json').read_text());result={'status':'PASS','candidate_manifest_sha256':EXPECTED,'candidate_inventory':before,'checks':[]}
    for opt in (False,True):run(candidate/'verify_manifest.py',opt)
    result['checks'].append('candidate strict file-set/SHA256 manifest verified normally and under -O')
    with tempfile.TemporaryDirectory(prefix='session51-blue-audit-') as td:
        temp=Path(td);real=temp/'actual-copy';shutil.copytree(HERE,real,ignore=shutil.ignore_patterns('__pycache__'))
        need(not real.is_symlink() and (real/'independent.py').resolve()!=HERE/'independent.py','actual independent relocated copy, not symlink')
        for opt in (False,True):
            p=run(real/'independent.py',opt);need(p.stdout==(HERE/'normal.json').read_text(),'independent real-copy normal/-O byte replay')
        result['checks'].append('independent actual-copy relocation from /tmp reproduces complete receipt byte-for-byte in normal and -O modes')
        relocated=temp/'candidate-copy';shutil.copytree(candidate,relocated,ignore=shutil.ignore_patterns('__pycache__'))
        need(not relocated.is_symlink(),'candidate literal copied directory')
        for opt in (False,True):
            run(relocated/'verify_manifest.py',opt)
            out=temp/('producer-optimized.json' if opt else 'producer-normal.json');run(relocated/'check.py',opt,('--output',out))
            need(out.read_bytes()==(candidate/'normal.json').read_bytes(),'candidate real-copy byte replay');compare(audit,json.loads(out.read_text()))
        result['checks'].append('candidate actual-copy normal/-O replays byte-identical and all 21 wiring/witness fields match independent audit')
        # Changed pinned source must fail in both modes, with no reliance on assert.
        source=real/'sources/converse-receipt.json';original=source.read_bytes();source.write_bytes(original+b' ')
        for opt in (False,True):need('source SHA256 mismatch' in run(real/'independent.py',opt,success=False).stderr,'source corruption diagnostic')
        source.write_bytes(original);result['checks'].append('independent modified source rejected normally and under -O')
        # Complement bug: do not flip non-path arcs.
        script=real/'independent.py';original=script.read_text();needle='else 1-matrix[i][j]';need(original.count(needle)==1,'unique complement mutation anchor')
        script.write_text(original.replace(needle,'else matrix[i][j]'))
        for opt in (False,True):need('literal C free involution' in run(script,opt,success=False).stderr,'semantic complement mutation diagnostic')
        script.write_text(original);result['checks'].append('semantic C-as-identity bug rejected normally and under -O')
        # Wrong blue predicate: count black presentations instead.
        needle='if not is_blue(matrix,path):continue';need(original.count(needle)==1,'unique blue mutation anchor')
        script.write_text(original.replace(needle,'if is_blue(matrix,path):continue'))
        for opt in (False,True):need('literal transpose-reflection exact blue' in run(script,opt,success=False).stderr,'semantic R predicate diagnostic')
        script.write_text(original);result['checks'].append('reversed blue/nonblue predicate rejected normally and under -O')
        # Candidate manuscript corruption caught by manifest, independent of scripts.
        theorem=relocated/'THEOREM.md';text=theorem.read_bytes();theorem.write_bytes(text+b'\nMutation.\n')
        for opt in (False,True):need('SHA256 mismatch: THEOREM.md' in run(relocated/'verify_manifest.py',opt,success=False).stderr,'candidate theorem mutation diagnostic')
        theorem.write_bytes(text);result['checks'].append('candidate theorem mutation rejected by strict manifest normally and under -O')
        # Receipt-level hostile controls test adjacency and actual coordinates.
        target=json.loads((candidate/'normal.json').read_text())
        for name,mutator in [('loop counted once',lambda d:d.__setitem__('internal_blue_ports',176)),('wrong actual target path',lambda d:d['internal_edges'][0]['target_path'].reverse()),('wrong boundary owner',lambda d:d['boundary_edges'][0].__setitem__('source_index',999)),('heuristic endpoint condition treated as sufficient',lambda d:d.__setitem__('internal_blue_ports',2037))]:
            changed=copy.deepcopy(target);mutator(changed)
            try:compare(audit,changed)
            except RuntimeError:result['checks'].append('receipt hostile rejected: '+name)
            else:raise RuntimeError('hostile passed: '+name)
    need(inventory(candidate)==before,'candidate source unchanged throughout audit')
    result['checks'].append('frozen producer source inventory unchanged; no producer packet edits or publication')
    result['independent_receipt_sha256']=digest(HERE/'normal.json');result['producer_receipt_sha256']=digest(candidate/'normal.json')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
