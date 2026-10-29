"""Four bounded, reproducible waves over tracked mathematical Markdown."""
from pathlib import Path
from random import Random
import subprocess,re,json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
WAVES={
 '1_literal_and_local':{
  'quoted_coefficients':r'(?<![\w-])(2508951|199148647|4372|96256)(?!\d)',
  'Ramanujan_Weber':r'Weber|class invariant|Ramanujan.s constant',
  'quoted_radicals':r'(?:sqrt|square root|surds?).{0,35}(?:22|37|65|130)|(?:22|37|65|130).{0,35}(?:sqrt|surds?)'},
 '2_mechanisms':{
  'Fricke_moonshine':r'Fricke|McKay.Thompson|Hauptmodul',
  'class_groups_trace':r'class number|class group|quadratic.unit|Pell|Dickson|Chebyshev',
  'quotient_transport':r'component.order|common transporter|free frame|Frobenius.*LCA|endpoint.*orbit'},
 '3_sidecars_and_counterexamples':{
  'topology_weight':r'homotopy|independence complex|return.Gram|signless Laplacian',
  'old_phase_losses':r'branch selection|conjugate.*trace|trace.*conjugate|hidden.*phase|phase.*hidden|relative.*coset',
  'arithmetic_interfaces':r'continued fraction|Kronecker|norm.*shell|norm.*form|tetrahedral.*boundary'},
 '4_orthogonal_old_corners':{
  'trace_operators':r'Fekete|Legendre walk|Wronskian|Macaulay|Pluecker|Plücker',
  'finite_field_and_codes':r'Fano|Golay|Gleason|eta.*24|24.*eta',
  'oddness_variants':r'Lucas.*parity|parity.*Lucas|odd.*class group|signed.*copies|copies.*signed'}
}
def main():
 paths=subprocess.check_output(['git','ls-files','*.md'],cwd=ROOT,text=True).splitlines()
 paths=[x for x in paths if x.startswith(('01-canon/theorems/','05-knowledge/results/','07-reflections/'))
        and not any(y in x.lower() for y in ('/inputs/','/copies/','/producer/','/vendor/'))
        and not x.startswith(HERE.relative_to(ROOT).as_posix()+'/')]
 docs=[]
 for name in paths:
  f=ROOT/name
  if f.exists():docs.append((name,f.read_text(encoding='utf-8',errors='replace')))
 result={'snapshot':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
         'corpus':'tracked Markdown in theorem/results/reflection roots; nested inputs/copies/producer/vendor excluded',
         'documents':len(docs),'waves':{}}
 for wi,(wave,terms) in enumerate(WAVES.items()):
  record={}
  for label,pattern in terms.items():
   regex=re.compile(pattern,re.I);hits=[]
   for name,text in docs:
    locations=[(i,line.strip()[:260]) for i,line in enumerate(text.splitlines(),1)
               if regex.search(line) and not (label=='quoted_coefficients' and 'THM-' in line)]
    if locations:hits.append({'path':name,'matches':len(locations),'first':locations[:2]})
   rng=Random(20261028+101*wi+len(label));sample=rng.sample(hits,min(8,len(hits)))
   record[label]={'regex':pattern,'matched_documents':len(hits),'random_tail_sample':sample,'hits':hits}
  result['waves'][wave]=record
  print(json.dumps({'wave':wave,'documents':{k:v['matched_documents'] for k,v in record.items()}},sort_keys=True),flush=True)
 (HERE/'occurrences.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
if __name__=='__main__':main()
