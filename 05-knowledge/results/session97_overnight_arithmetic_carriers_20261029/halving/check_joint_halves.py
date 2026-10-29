"""Independent native-model rational group law; no SymPy dependency."""
import copy,json
from fractions import Fraction as Q
from pathlib import Path
D=Path(__file__).resolve().parent;checks=0

def need(b):
 global checks
 checks+=1
 if not b:raise ValueError('joint half rejected')

def verify(c,k):
 a1,a2,a3,a4,a6=map(Q,c['curve']['ainvs'])
 def on(P):
  if P is None:return True
  x,y=P;return y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
 def add(P,R):
  need(on(P)and on(R))
  if P is None:return R
  if R is None:return P
  x,y=P;u,v=R
  if x==u and y+v+a1*x+a3==0:return None
  if P==R:m=(3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
  else:need(x!=u);m=(v-y)/(u-x)
  nu=y-m*x;xx=m*m+a1*m-a2-x-u;yy=-(m+a1)*xx-nu-a3
  T=(xx,yy);need(on(T));return T
 need(k['curve_id']==c['curve']['id'])
 inds=k['indices_zero_based'];need(len(set(inds))==len(inds));need(k['signs']==[1]*len(inds))
 P=None
 for i in inds:need(0<=i<30);P=add(P,tuple(map(Q,c['curve']['points'][i])))
 x,y=P;need([str(x),str(y+(a1*x+a3)/2)]==k['completed_square_sum'])
 H=tuple(map(Q,k['native_half']));need(on(H));need(add(H,H)==P)
 hx,hy=H;need([str(hx),str(hy+(a1*hx+a3)/2)]==k['completed_square_half'])
 # The opposite half doubles to -P, not P; the retained sign matters.
 neg=(hx,-hy-a1*hx-a3);need(add(neg,neg)!=P)
 mask=sum(1<<i for i in inds)
 for row in c['character_rows']:need((mask&row['bits']).bit_count()%2==0)
 return {'curve':k['curve_id'],'support':inds,'exact_rational_half':True,'native_half':k['native_half']}
cs=json.loads((D/'halving_certificate.json').read_text())['curves'];results=[];rejected=0
for c in cs:
 k=json.loads((D/f'joint_candidate_{c["curve"]["id"]}.json').read_text());results.append(verify(c,k))
 for mode in ('half','index','sign'):
  bad=copy.deepcopy(k)
  if mode=='half':bad['native_half'][1]=str(Q(bad['native_half'][1])+1)
  if mode=='index':bad['indices_zero_based'][0]=0
  if mode=='sign':bad['signs'][0]=-1
  try:verify(c,bad)
  except ValueError:rejected+=1
  else:raise ValueError('hostile accepted')
need(rejected==6)
print(json.dumps({'checks':checks,'hostiles_rejected':rejected,'results':results,'scope':'One certified joint half per curve; not full 2-saturation'},sort_keys=True))
