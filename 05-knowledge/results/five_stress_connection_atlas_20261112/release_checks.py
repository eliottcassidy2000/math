"""Validate local links/scopes and pin source and artifact hashes at release."""
from pathlib import Path
import re,json,gzip,hashlib,subprocess,ast,os,importlib.metadata as md
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(x,msg):
    if not x:raise RuntimeError(msg)
def native(p):
    p=p.resolve()
    return Path('\\\\?\\'+str(p)) if os.name=='nt' and not str(p).startswith('\\\\?\\') else p
def digest(p):return hashlib.sha256(native(p).read_bytes()).hexdigest()
def run():
    sources=set();links=0
    for p in HERE.glob('*.md'):
        # This index quotes historical Markdown verbatim; embedded links
        # belong to the quoted source, not to this packet's navigation.
        if p.name=='PROMPT-INDEX.md':continue
        s=p.read_text(encoding='utf8')
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',s):
            target=target.strip('<>').split('#',1)[0]
            if not target or '://' in target or target.startswith('mailto:'):continue
            q=(p.parent/target).resolve();need(native(q).exists(),p.name+' missing link '+target);links+=1
            if not q.is_relative_to(HERE):sources.add(q)
    for p in HERE.glob('*.py'):ast.parse(p.read_text(encoding='utf8'))
    for p in HERE.glob('*.json'):json.loads(p.read_text(encoding='utf8'))
    for p in HERE.glob('*.json.gz'):json.loads(gzip.decompress(p.read_bytes()))
    ext=json.loads((HERE/'round_extension.json').read_text());cl=json.loads((HERE/'seed_classification.json').read_text())
    need(ext['levels'][-1]['n']==18,'round extent')
    need(sum(q['unit_marked_isomorphism_types'] for q in ext['levels'])==590,'round marked total')
    dag=json.loads((HERE/'contraction_dag.json').read_text())
    need(dag['states']==161 and dag['edges_between_marked_types']==63,'contraction graph control')
    need(all(q['nonunique_terminal_births']==0 for q in dag['levels']),'finite confluence')
    desc=json.loads(gzip.decompress((HERE/'descendants_nine.json.gz').read_bytes()))
    need(desc['marked_descendant_types']==122 and desc['edges_to_order8_types']==127,'order9 descendant scope')
    bc=json.loads((HERE/'birth_cones.json').read_text())
    need(bc['births']==98 and not bc['overlaps'] and len(bc['disjoint_certificates'])==10,'all-length birth-cone certificate')
    hs=json.loads((HERE/'history_stats.json').read_text())
    need(hs['candidate_text_objects']==78446 and hs['counts']['bytes_scanned']==56587125074,'history scope')
    ps=json.loads((HERE/'history_prompt_harvest_summary.json').read_text())
    need(ps['unique_normalized_contexts']==1138 and ps['contexts_with_five_in_excerpt']==193,'prompt harvest scope')
    docs=json.loads((HERE/'documentation_check_receipt.json').read_text())
    need(docs['status']=='PASS' and docs['exit_code']==0,'saved documentation checker result')
    replay=json.loads((HERE/'replay_receipt.json').read_text())
    need(len(replay['results'])==16,'optimized replay count')
    for row in replay['results']:
        for check in row['checks']:
            if 'normal_lf_sha256' in check:actual=hashlib.sha256((HERE/check['file']).read_bytes().replace(b'\r\n',b'\n')).hexdigest();expected=check['normal_lf_sha256']
            else:actual=digest(HERE/check['file']);expected=check['normal_sha256']
            need(actual==expected,'replay receipt changed: '+check['file'])
    dependencies=[
      'modular_centers_collatz_pell_20261111/tournament_checks.py',
      'modular_centers_collatz_pell_20261111/content_checks.py',
      'endpoint_coverage_compiler_20261111/inflation_bank.py',
      'endpoint_coverage_compiler_20261111/inflation_bank_7.json',
      'seed_recursion_births_20261111/cones.py',
      'seed_recursion_births_20261111/port_transport.py',
      'metagraph_aspect_sequences_20261111/records.jsonl.gz',
      'session251_creative_completion_20261111/ternary/compiler.py']
    sources.update(HERE.parent/q for q in dependencies)
    pins=[]
    for p in sorted(sources):
        rel=p.relative_to(ROOT).as_posix();need(native(p).is_file(),'source file '+rel)
        result=subprocess.run(['git','rev-parse','HEAD:'+rel],cwd=ROOT,capture_output=True,text=True)
        pins.append({'path':rel,'worktree_sha256':digest(p),'HEAD_blob':result.stdout.strip() if result.returncode==0 else None,'reading_scope':'Linked source or computation dependency; HISTORY and AUDIT specify the inspected portions. This entry does not assert full reading or revalidation.'})
    pdfs=[]
    for name,scope in [
      ('The_Multiset_Dimension_of_Graphs_Extremal_Values_a.pdf','Opening/title and relevant sections14-16, selected later material; finite local automaton reconstructed, no billion-graph or Lean audit.'),
      ('2607.12873v1.pdf','Opening/title and record-composition formula; existing THM-5137 descendant and correction lineage consulted.')]:
        p=Path(r'C:\Users\Eliott\Downloads')/name
        pdfs.append({'filename':name,'sha256':digest(p),'reading_scope':scope})
    source={'baseline':'7510e38f9c77','checkout_parent_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),'sources':pins,'supplied_pdfs':pdfs,'external_primary_routes':[
      'https://doi.org/10.1002/zamm.19330130506',
      'https://www.sciencedirect.com/science/article/pii/S0065215608703004',
      'https://link.springer.com/article/10.1007/s00707-008-0118-3',
      'https://doi.org/10.1007/s00707-010-0370-1',
      'https://ir.cwi.nl/pub/6805',
      'https://oeis.org/A000016','https://oeis.org/A000568',
      'https://oeis.org/A135244','https://oeis.org/A002496']}
    (HERE/'source_manifest.json').write_text(json.dumps(source,indent=2)+'\n',encoding='utf8')
    receipt={'status':'PASS','local_markdown_links':links,'source_pins':len(pins),'python':os.sys.version,'packages':{q:md.version(q) for q in ['sympy','numpy','numba']},'checks':['Packet Python syntax and JSON/gzip parsing','All current packet navigation links resolve; original links inside verbatim PROMPT-INDEX excerpts excluded','Declared census and prompt-harvest extents agree with receipts','Source/PDF hashes recorded'],'matplotlib_note':'Figure generated with Matplotlib in a dedicated temporary package directory; numerical checker environment unchanged.'}
    (HERE/'release_checks_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    files=[]
    for p in sorted(HERE.iterdir()):
        if not p.is_file() or p.name in ['release_manifest.json','RELEASE.md']:continue
        record={'file':p.name,'sha256':digest(p),'bytes':p.stat().st_size}
        if p.suffix in ['.md','.py','.json','.txt','.svg']:record['lf_sha256']=hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        files.append(record)
    (HERE/'release_manifest.json').write_text(json.dumps({'excludes':['release_manifest.json','RELEASE.md'],'byte_convention':'sha256 and bytes record release worktree bytes. lf_sha256 normalizes CRLF to LF for text so Git/platform newline conversion is not a mathematical change.','files':files},indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
if __name__=='__main__':run()
