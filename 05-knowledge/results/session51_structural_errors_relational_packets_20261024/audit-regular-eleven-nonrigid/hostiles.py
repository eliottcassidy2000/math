#!/usr/bin/env python3
import argparse,importlib.util,json,pathlib,random,subprocess,tempfile
HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent',HERE/'independent.py');ind=importlib.util.module_from_spec(spec);spec.loader.exec_module(ind)
def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--output',required=True);a=p.parse_args();binary=str(pathlib.Path(a.binary).resolve());data=json.loads((HERE/'normal.json').read_text());graphs=[tuple(x['rows']) for x in data['representatives']];results=[]
 def relabel(r,p):
  out=[0]*len(r)
  for i in range(len(r)):
   for j in range(len(r)):out[p[i]]|=((r[i]>>j)&1)<<p[j]
  return tuple(out)
 rng=random.Random(51120);p=list(range(11));rng.shuffle(p);renamed=relabel(graphs[0],p)
 ind.need(ind.isomorphisms(graphs[0],renamed,True),'positive relabel isomorphism')
 nonself=next(r for r in graphs if not ind.isomorphisms(r,tuple(x^(((1<<11)-1)^(1<<i)) for i,x in enumerate(r)),True))
 results.append(dict(name='ordinary-isomorphism-does-not-fold-converse',accepted=True))
 with tempfile.TemporaryDirectory(prefix='nonrigid-audit-hostile-') as d:
  d=pathlib.Path(d)
  def run(name,records,expected):
   path=d/'representatives.txt';path.write_text(str(len(records))+'\n'+'\n'.join(' '.join(map(str,r)) for r in records)+'\n')
   c=subprocess.run([binary,str(path),str(HERE/'sources/parents-ten.txt'),str(HERE/'sources/states.tsv'),str(d/'result.json')],text=True,capture_output=True)
   ind.need((c.returncode==0)==expected,'incorrect hostile acceptance '+name)
   results.append(dict(name=name,expected_acceptance=expected,returncode=c.returncode,diagnostic=c.stderr.strip() or c.stdout.strip()))
  run('pristine-independent-nonrigid-list',graphs,True)
  run('arbitrarily-relabelled-class', [renamed]+graphs[1:],True)
  run('missing-nonrigid-class',graphs[1:],False)
  run('nonliteral-isomorphic-duplicate',graphs+[renamed],False)
  bad=list(graphs[0]);bad[0]^=2;run('one-sided-arc-toggle',[tuple(bad)]+graphs[1:],False)
  bad=list(graphs[0]);bad[0]^=2;bad[1]^=1;run('valid-tournament-nonregular-arc-flip',[tuple(bad)]+graphs[1:],False)
 out=dict(status='PASS',controls=results);pathlib.Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
