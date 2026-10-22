import json,pathlib,itertools
from fractions import Fraction
from experiment import frombits
from find_isospectral import charpoly
ROOT=pathlib.Path(__file__).resolve().parent
rows=[json.loads(x) for x in open(ROOT/'census8_predictions.jsonl')];rho={}
for r in rows:
 cp=charpoly(frombits(8,r['bits']));num=sum(c*8**(8-i) for i,c in enumerate(cp));den=abs(sum(c*8**(8-i)*(-1)**((8-i)//2) for i,c in enumerate(cp) if i%2==0));rho[r['index']]=Fraction(num,den)
pairs=[(a,b) for a,b in itertools.combinations(rows,2) if a['H']!=b['H']];out={}
for key in ['rho','kappa1','kappa2']:
 wrong=[];tied=[];correct=[]
 for a,b in pairs:
  if key=='rho':d=rho[a['index']]-rho[b['index']];tie=d==0
  else:d=a['prediction'][key]-b['prediction'][key];tie=abs(d)<1e-10
  pair=[a['index'],b['index']]
  if tie:tied.append(pair)
  elif d*(a['H']-b['H'])<0:wrong.append(pair)
  else:correct.append(pair)
 out[key]=dict(total=len(pairs),reversed=len(wrong),tied_or_numerically_unresolved=len(tied),correct=len(correct),reversed_pairs=wrong,tied_pairs=tied)
json.dump(out,open(ROOT/'census8_ranking.json','w'),indent=2)
print({k:{a:b for a,b in v.items() if not a.endswith('pairs')} for k,v in out.items()})
