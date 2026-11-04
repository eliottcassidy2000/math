"""Pinned exact replay; benchmark times are intentionally nondeterministic."""
from pathlib import Path
import hashlib,json,subprocess,sys,platform
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PREV=HERE.parent/'return18979_ports_20261103'
SCRIPTS=('integer_census','bucket_census','exterior_search','exterior_two','audit_feedback','feedback_search','offset_fibres','next_fibres','family_certificate','source_pins','splice_mobius','first_ready','audit_first_ready')

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,d):(HERE/name).write_text(json.dumps(d,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        old=json.loads((PREV/'sources.json').read_text());paths={d['path'] for d in old['dependencies']}
        paths.update(p.relative_to(ROOT).as_posix() for p in PREV.iterdir() if p.suffix in ('.py','.md','.json') and not p.name.startswith('verification.'))
        incoming=HERE.parent/'return18979_palindrome_20261103'
        paths.update(p.relative_to(ROOT).as_posix() for p in incoming.iterdir() if p.suffix in ('.py','.cpp','.md','.json') and not p.name.startswith('verification.'))
        paths.update('05-knowledge/results/session98_tournament_growth_20261029/'+p for p in ('recurrence/verify.py','recurrence/PROJECTIVE_RESPONSE.md','README.md'))
        paths.add('01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md')
        paths.update('05-knowledge/results/'+p for p in (
            'session175_tournament_spectrum_20261104/a000568/projective.py',
            'session175_tournament_spectrum_20261104/a000568/inherited_local_scale.py',
            'session175_tournament_spectrum_20261104/a000568/README.md',
            'session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md',
            'session87_four_vertices_enumeration_20261028/fast_counts.py',
            'session87_four_vertices_enumeration_20261028/FAST_COUNTS.md',
            'session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md',
            'lca_ports_reverse_20261025/counting/README.md',
            'orbit_spectrum_session_20261103/counting/NOTE.md',
            'orbit_spectrum_session_20261103/README.md'))
        changes=[{'path':d['path'],'previous_sha256':d['sha256'],'current_sha256':digest(ROOT/d['path'])} for d in old['dependencies'] if digest(ROOT/d['path'])!=d['sha256']]
        need(all(x['path'].endswith('.md') for x in changes),'inherited computational input changed; review before pinning')
        save('sources.json',{'status':'PINNED_SOURCE_INHERITANCE','startup_commit':'9af43540fd6e','hash_scheme':'SHA256 after LF newline normalization',
            'dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],'incoming_document_changes':changes,'one_birth55_orientations':old['one_birth55_orientations']})
        save('manifest.json',{'code':[{'file':p.name,'sha256':digest(p)} for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.cpp')],
            'receipts':[{'file':k+'.json','sha256':digest(HERE/(k+'.json'))} for k in SCRIPTS]})
        print('PINNED',len(paths),'inherited files and',len(SCRIPTS),'receipts');return
    sources=json.loads((HERE/'sources.json').read_text());manifest=json.loads((HERE/'manifest.json').read_text())
    for row in sources['dependencies']:need(digest(ROOT/row['path'])==row['sha256'],'changed inherited input: '+row['path'])
    for row in manifest['code']:need(digest(HERE/row['file'])==row['sha256'],'changed code: '+row['file'])
    expected={z['file']:z['sha256'] for z in manifest['receipts']};out=[]
    only=next((a.split('=',1)[1].split(',') for a in sys.argv if a.startswith('--only=')),None)
    scripts=SCRIPTS if only is None else tuple(only)
    need(all(s in SCRIPTS for s in scripts),'declared replay scripts')
    for script in scripts:
        args=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(script+'.py'))]
        p=subprocess.run(args,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(p.returncode==0,script+' failed:\n'+p.stdout+p.stderr)
        name=script+'.json';need(digest(HERE/name)==expected[name],'receipt drift: '+script)
        out.append({'script':script,'sha256':expected[name],'stdout':p.stdout.strip()});print('PASS',script,flush=True)
    mode='optimized' if sys.flags.optimize else 'normal'
    save('verification.'+('extra.' if only else '')+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),'inherited_files':len(sources['dependencies']),'receipts':out,
        'timings':'not hash checked; latest census_timings.json is an observed run'})
    print('PASS',mode,len(scripts),'research receipts',flush=True)
if __name__=='__main__':main()
