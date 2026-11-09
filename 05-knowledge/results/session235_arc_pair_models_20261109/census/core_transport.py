#!/usr/bin/env python3
"""Exact transport of Session234 complete semicomplete core grammar."""
import collections,csv,importlib.util,itertools,json,math,pathlib
from fractions import Fraction
ROOT=pathlib.Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]/'session234_digraphs_20261109/census/generator.py'
spec=importlib.util.spec_from_file_location('s234_generator',SOURCE);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
def paths(a):
 return [list(p) for p in itertools.permutations(range(len(a))) if all(a[i][j] for i,j in zip(p,p[1:]))]
def decode(a):return [[int(a[i][j] and not a[j][i]) for j in range(len(a))] for i in range(len(a))]
def aut(a):
 n=len(a)
 return [list(p) for p in itertools.permutations(range(n)) if all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n))]
cores=[];groups=collections.Counter()
for h,encoded in g.CORE_ROWS.items():
 for k,raw in enumerate(encoded):
  a=[[int(v) for v in row] for row in raw.split('/')];b=decode(a)
  ps=paths(a);pt=paths(b);aa=aut(a);ab=aut(b)
  assert len(ps)==h and aa==ab and len(g.scc(a))==1
  assert g.hp(a)==h and g.hp(b)==len(pt)
  groups[h,len(a),len(pt),len(aa)]+=1
  cores.append({'source_id':f'{h}:{len(a)}:{k}','source_H':h,'n':len(a),'decoded_H':len(pt),'automorphism_order':len(aa),'source_matrix':raw.split('/'),'decoded_matrix':[''.join(map(str,r)) for r in b],'source_paths':ps,'decoded_paths':pt,'automorphisms':aa})
summary=[{'source_H':h,'n':n,'decoded_H':r,'automorphism_order':a,'classes':c,'sum_inverse_automorphisms':str(Fraction(c,a))} for (h,n,r,a),c in sorted(groups.items())]
# Terms x^order/(1-x)^gaps, for isomorphism OGF and labelled EGF.
terms=[]
for h in (7,21):
 bins=collections.defaultdict(lambda:[0,Fraction(0)])
 for c in cores:
  if c['source_H']==h:
   z=bins[c['decoded_H'],c['n'],2];z[0]+=1;z[1]+=Fraction(1,c['automorphism_order'])
 if h==21:
  for a in cores:
   if a['source_H']!=7:continue
   for b in cores:
    if b['source_H']!=3:continue
    z=bins[a['decoded_H']*b['decoded_H'],a['n']+b['n'],3];z[0]+=2;z[1]+=Fraction(2,a['automorphism_order']*b['automorphism_order'])
 for (r,n,k),(u,l) in sorted(bins.items()):terms.append({'source_H':h,'decoded_H':r,'numerator_order':n,'denominator_power':k,'isomorphism_coefficient':u,'labelled_EGF_coefficient':str(l)})
checks=[]
for n in range(1,11):
 for h in (7,21):
  actual=collections.Counter();labelled=collections.Counter()
  for key,a in g.generate(g.load_cores(),h,n):
   r=g.hp(decode(a));actual[r]+=1
   aid={c['source_id']:c['automorphism_order'] for c in cores}
   order=math.prod(1 if v[0]=='singleton' else aid[v[0]] for v in key)
   labelled[r]+=math.factorial(n)//order
  pred=collections.Counter();predlab=collections.Counter()
  for t in terms:
   if t['source_H']!=h or n<t['numerator_order']:continue
   ways=math.comb(n-t['numerator_order']+t['denominator_power']-1,t['denominator_power']-1)
   pred[t['decoded_H']]+=t['isomorphism_coefficient']*ways
   predlab[t['decoded_H']]+=Fraction(t['labelled_EGF_coefficient'])*math.factorial(n)*ways
  assert actual==pred and labelled==predlab,(n,h)
  if n<=6:
   joint=csv.DictReader((ROOT/f'joint{n}.csv').open());finite={int(row['H_oriented']):int(row['count']) for row in joint if int(row['H_semicomplete'])==h}
   assert dict(labelled)==finite,(n,h,labelled,finite)
  checks.append({'n':n,'source_H':h,'class_counts':dict(actual),'labelled_counts':dict(labelled)})
out={'source':'session234_digraphs_20261109/census/generator.py:CORE_ROWS','core_count':len(cores),'grouped_cores':summary,'generating_function_terms':terms,'cores':cores,'generated_grammar_checks':checks}
(ROOT/'core_transport.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'grouped_cores':summary,'generating_function_terms':terms},indent=2))
