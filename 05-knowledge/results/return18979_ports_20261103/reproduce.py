"""LF-stable pinned replay of source-aware repair and composition certificates."""
from pathlib import Path
import hashlib,json,sys,subprocess,platform
import sympy
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PREV=HERE.parent/'return6325_observer_20261103'
SCRIPTS=('grandparent_repair','apex_carry','square_splices','hostiles','moment_search','funnel_search')
EXTRA=('return6325_observer_20261103/return6325.json','return6325_observer_20261103/completion_hostile.json',
 'return6325_observer_20261103/PROOFS.md','return6325_observer_20261103/HANDOFF.md',
 'return6325_observer_20261103/ready_filter.cpp','return6325_completion_20261103/CLONE_LAW.md',
 'return6325_completion_20261103/RANK_DEPTH.md','session173_modular_clocks_20261103/squares/NOTE.md',
 'digit_memory_all_bases_20261103/TOURNAMENT.md','feedback_word_hamilton_spectrum_20261020/proofs.md',
 'feedback_word_hamilton_spectrum_20261020/run_density.md','session172_memory_dynamics_20261103/squares/NOTE.md')

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def save(name,d):(HERE/name).write_text(json.dumps(d,indent=2)+'\n',encoding='utf8',newline='\n')

def main():
    if '--pin' in sys.argv:
        old=json.loads((PREV/'sources.json').read_text());paths={d['path'] for d in old['dependencies']}
        paths.update(p.relative_to(ROOT).as_posix() for p in PREV.glob('*.py'))
        paths.update('05-knowledge/results/'+p for p in EXTRA)
        # Include the direct formula/grammar imports even if the earlier packet's
        # broader pins happen not to include a newly selected file.
        for dirname in ('polarized_path_calculus_20261103','two_entry_funnels_20261103','ready_return_pentagon_20261103'):
            paths.update(p.relative_to(ROOT).as_posix() for p in (HERE.parent/dirname).glob('*.py'))
        inherited_changes=[{'path':d['path'],'startup_sha256':d['sha256'],'current_sha256':digest(ROOT/d['path'])} for d in old['dependencies'] if digest(ROOT/d['path'])!=d['sha256']]
        need(all(x['path'].endswith('.md') for x in inherited_changes),'incoming computational input changed: review before pinning')
        save('sources.json',{'status':'PINNED_INHERITANCE_AND_SELECTED_FRESH_REPLAY','startup_commit':'26fe3955ae64','incoming_commits':['eb88b1beb3','976872ac5c','bced438f00','d53992645e'],
          'hash_scheme':'SHA256 after LF newline normalization','dependencies':[{'path':p,'sha256':digest(ROOT/p)} for p in sorted(paths)],
          'incoming_document_changes':inherited_changes,'hostile81_27':{'status':'FRESH_REPLAY','receipt':'hostiles.json'},'one_birth55_orientations':old['one_birth55_orientations']})
        save('manifest.json',{'code':[{'file':p.name,'sha256':digest(p)} for p in sorted(HERE.glob('*.py'))],
          'receipts':[{'file':k+'.json','sha256':digest(HERE/(k+'.json'))} for k in SCRIPTS]+[{'file':'search_n14_b1.json','sha256':digest(HERE/'search_n14_b1.json')}]})
        print('Pinned',len(paths),'source files; incoming documentation changes',len(inherited_changes));return
    sources=json.loads((HERE/'sources.json').read_text());manifest=json.loads((HERE/'manifest.json').read_text())
    for row in sources['dependencies']:need(digest(ROOT/row['path'])==row['sha256'],'changed inherited input: '+row['path'])
    for row in manifest['code']:need(digest(HERE/row['file'])==row['sha256'],'changed code: '+row['file'])
    expected={r['file']:r['sha256'] for r in manifest['receipts']};out=[]
    scripts=list(SCRIPTS)+(['ready18979'] if '--full' in sys.argv else [])
    for k in scripts:
        args=[sys.executable,'-X','utf8','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/(k+'.py'))]
        r=subprocess.run(args,cwd=ROOT,capture_output=True,text=True,encoding='utf8')
        need(r.returncode==0,k+' failed:\n'+r.stdout+r.stderr)
        receipt='search_n14_b1.json' if k=='ready18979' else k+'.json'
        need(digest(HERE/receipt)==expected[receipt],'receipt drift: '+k)
        out.append({'script':k,'stdout':r.stdout.strip(),'sha256':expected[receipt]});print('PASS',k,flush=True)
    # Independent C++ bank and selective recurrence address the same full n14 bank.
    cpp=json.loads((HERE/'search_n14_b1.json').read_text());mom=json.loads((HERE/'moment_search.json').read_text())
    need(cpp['hit'] is None and mom['hit'] is None and cpp['bank_size']==mom['banks'][0]['examined']==46368,'cross-algorithm order14 search control')
    mode='optimized' if sys.flags.optimize else 'normal'
    save('verification.'+mode+'.json',{'status':'PASS','mode':mode,'python':platform.python_version(),'sympy':sympy.__version__,'inherited_files':len(sources['dependencies']),'full_cpp_rerun':'--full' in sys.argv,'receipts':out})
    print('PASS',mode,'all source-preserving certificates',flush=True)
if __name__=='__main__':main()
