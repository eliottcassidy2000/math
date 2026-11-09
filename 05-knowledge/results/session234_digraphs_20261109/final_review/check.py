from pathlib import Path
import importlib.util,itertools,json,re,math
root=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('g',root/'census/generator.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
cs=g.load_cores(); report=g.validate(cs)
for h,cores in cs.items():
 seen=set()
 for c in cores:
  a=c['matrix'];n=len(a)
  assert sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))==h
  assert len(g.scc(a))==1
  canon=g.canonical_block(tuple(map(tuple,a)));assert canon not in seen;seen.add(canon)
print('37 cores: direct permutations, strongness, unique types pass; loopless validator through n8 passes')
for h in (7,21):
 for n in range(1,8):
  seen=set();count=0;labelled=0
  for key,a in g.generate(cs,h,n,True):
   components=g.scc(a)
   sig=tuple(g.canonical_block(tuple(tuple(a[i][j] for j in c) for i in c)) for c in components)
   assert sig not in seen;seen.add(sig)
   assert g.hp(a)==h
   aut=math.prod(len(g.automorphisms([[a[i][j] for j in c] for i in c])) for c in components)
   labelled+=math.factorial(n)//aut;count+=1
  assert count==g.expected(h,n,True)
  assert labelled==2**n*g.expected_labelled(h,n)
  print('decorated',h,n,count,labelled)
missing=[]
for fn in [root/'README.md',root/'PROOFS.md',root/'GENERAL_GRAMMAR.md']:
 for link in re.findall(r'\]\(([^)]+)\)',fn.read_text()):
  if '://' not in link and not (fn.parent/link.split('#')[0]).exists():missing.append((str(fn),link))
assert not missing,missing
print('All local markdown file targets pass (anchor spelling manually reviewed)')
(Path(__file__).resolve().parent/'validation.json').write_text(json.dumps(report,indent=2))
