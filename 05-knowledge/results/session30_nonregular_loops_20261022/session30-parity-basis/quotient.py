import itertools,json,collections,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
out=[];summary=[]
for core in json.loads((ROOT/'parity-kernels.json').read_text()):
 edges=[tuple(e) for e in core['edges']];target=sorted(edges)
 autos=[p for p in itertools.permutations(range(core['V'])) if sorted(tuple(sorted((p[a],p[b]))) for a,b in edges)==target]
 terms=collections.defaultdict(int);killed=0
 for t in core['terms']:
  reps=[]
  for p in autos:
   triples=[];sign=1
   for (a,b),i in zip(edges,t['indices']):
    a,b=p[a],p[b]
    if a>b:
     a,b=b,a
     if i>=2:sign=-sign
    triples.append((a,b,i))
   triples.sort();reps.append((tuple(x[2] for x in triples),sign))
  key=min(r[0] for r in reps);signs={s for k,s in reps if k==key}
  if len(signs)==2:killed+=1;continue
  terms[key]+=t['numerator']*next(iter(signs))
 terms={k:v for k,v in terms.items() if v}
 new={k:v for k,v in core.items() if k not in ('terms','loop_pruned_terms')}
 new.update(terms=[{'indices':list(k),'numerator':v} for k,v in sorted(terms.items())],vertex_automorphisms=len(autos),antisymmetric_monomials_killed=killed)
 out.append(new)
 summary.append({'weight':core['weight'],'edges':edges,'parity_terms':len(core['terms']),'quotient_terms':len(terms),'antisymmetric_monomials_killed':killed})
(ROOT/'quotient-kernels.json').write_text(json.dumps(out,indent=2)+'\n')
report={'status':'constructed; contractions require verification','parity_terms':sum(s['parity_terms'] for s in summary),'quotient_terms':sum(s['quotient_terms'] for s in summary),'antisymmetric_monomials_killed':sum(s['antisymmetric_monomials_killed'] for s in summary),'cores':summary}
(ROOT/'quotient-results.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k!='cores'})
