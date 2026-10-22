import json,itertools,collections,pathlib
ROOT=pathlib.Path(__file__).resolve().parent
source=ROOT.parent/'session30-decorated-engine/symbolic-kernels.json'
change=[{1:1,2:1},{1:1,2:-1},{1:-1,2:-1,3:1},{0:-1}]
out=[];summary=[]
for core in json.loads(source.read_text()):
 terms=collections.defaultdict(int)
 for term in core['terms']:
  for choices in itertools.product(*(change[i].items() for i in term['indices'])):
   idx=tuple(x[0] for x in choices);v=term['numerator']
   for _,c in choices:v*=c
   terms[idx]+=v
 terms={k:v for k,v in terms.items() if v}
 odd={k:v for k,v in terms.items() if sum(i>=2 for i in k)%2}
 if odd:raise RuntimeError(('odd parity',core['edges'],odd))
 nonloop={k:v for k,v in terms.items() if all(a!=b or k[e]<2 for e,(a,b) in enumerate(core['edges']))}
 new={k:v for k,v in core.items() if k not in ('basis','terms')}
 new.update(basis=['G','H','O','P'],terms=[{'indices':list(k),'numerator':v} for k,v in sorted(terms.items())],loop_pruned_terms=[{'indices':list(k),'numerator':v} for k,v in sorted(nonloop.items())])
 out.append(new)
 summary.append(dict(weight=core['weight'],edges=core['edges'],old_terms=len(core['terms']),parity_terms=len(terms),loop_pruned_terms=len(nonloop),dense_size=4**len(core['edges'])))
(ROOT/'parity-kernels.json').write_text(json.dumps(out,indent=2)+'\n')
report={'status':'PASS','cores':len(out),'old_total_terms':sum(x['old_terms'] for x in summary),'new_total_terms':sum(x['parity_terms'] for x in summary),'loop_pruned_total_terms':sum(x['loop_pruned_terms'] for x in summary),'cores_summary':summary}
(ROOT/'transformation-results.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:v for k,v in report.items() if k!='cores_summary'})
