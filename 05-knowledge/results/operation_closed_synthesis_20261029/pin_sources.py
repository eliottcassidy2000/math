"""Pin inspected source blobs; validate packet links and deterministic replay."""
from pathlib import Path
import hashlib,json,re,subprocess,tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
START='4f53ff2859d8c79ec417d441b74777d9e6f3b8bc'
INCOMING='75bd8b27b6'
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def digest(b):return hashlib.sha256(b).hexdigest()
refs=set()
for f in ['README.md','THEORY.md','RECOVERY.md']:
    text=(HERE/f).read_text(encoding='utf-8')
    for link in re.findall(r'\]\(([^)]+)\)',text):
        if link.startswith(('http:','https:')):continue
        p=(HERE/link.split('#')[0]).resolve()
        if not p.is_file():raise RuntimeError('missing link '+str(p))
        if not p.is_relative_to(HERE):refs.add(p.relative_to(ROOT).as_posix())
refs.update([
 '04-computation/experiments/growth_gf_runs_fermat_20261029/audit2/alib.py',
 '04-computation/experiments/growth_gf_runs_fermat_20261029/runs.py',
 '05-knowledge/results/four_vertex_glue_fermat_20261029/verify.py',
 '05-knowledge/results/overnight_rank_return_20261029/two_run.json'])
sources=[]
for p in sorted(refs):
    b=git('show',START+':'+p)
    sources.append({'revision':START,'path':p,'sha256':digest(b),'bytes':len(b),
                    'scope':'targeted passages inspected; not a full independent audit'})
for p in [
 '05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md',
 '05-knowledge/hypotheses/HYP-9446-circulant-collatz-dictionary-and-multiplier-descent-open-questions.md',
 '05-knowledge/results/session102_tournament_classification_20261029/switching/SIGNED_CUTS_AND_SIMPLE_SYMMETRY.md']:
    b=git('show',INCOMING+':'+p)
    sources.append({'revision':git('rev-parse',INCOMING).decode().strip(),'path':p,'sha256':digest(b),'bytes':len(b),
                    'scope':'incoming status/change inspected; not a proved dependency'})
replay={}
for p,q in [('tournament_results.json','operation_closed_tournaments_optimized.json'),('audit_results.json','operation_closed_audit_optimized.json')]:
    a=(HERE/p).read_bytes();b=(Path(tempfile.gettempdir())/q).read_bytes()
    if a!=b:raise RuntimeError('normal/optimized mismatch '+p+'; ensure both jobs finished')
    replay[p]={'sha256':digest(a),'bytes':len(a),'normal_optimized_identical':True}
out={'source_revision':START,'sources':sources,'replay':replay,'primary_literature':[
 {'author':'Michael R. Shelor','title':'On Tournaments Having A Unique Hamiltonian Circuit','date':'May 1977',
  'url':'https://paperzz.com/doc/7419910/shelor--m.r.%3B--1977-.on-tournaments-having-a-unique-hamil...',
  'inspected':'primary manuscript text mirror, sections2--5 and Theorem4; PDF access attempts failed'},
 {'authors':'Gutin, Rafiey, Yeo','title':'On n-partite tournaments with unique n-cycle',
  'url':'https://www.cs.rhul.ac.uk/~gutin/paperstsp/ucycle22.pdf',
  'inspected':'introduction and references; Garey attribution confirmation, not a new imported theorem'}]}
(HERE/'source_manifest.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
files=[p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256']
(HERE/'MANIFEST.sha256').write_text(''.join(digest(p.read_bytes())+'  '+p.name+'\n' for p in sorted(files)),encoding='utf-8',newline='\n')
print(json.dumps({'pinned_documents':len(sources),'replay':replay},indent=2))
