from inherited import *
import json,collections

def seeds():
 out=[]
 for X,js,r,a,b in PAIRS:out.append((X,remark(a,r),remark(b,r)))
 for X,a,r,b,t,w in [(4225,(2,20,41,67,77,27,167,63),6,(6,12,24,49,195,23,47,111),4,'JO'),(4357,(10,20,65,22,37,15,59),3,(2,48,11,67,45,77,23),2,'ZJO')]:
  a,b=remark(a,r),remark(b,t)
  for c in w:a,b=graph_move(a,c),graph_move(b,c)
  out.append((X,a,b))
 return out

def X(p):return p[0]+p[3]-p[1]+2*((p[1]-p[1]%3)//3)
def app(p,w):
 for c in w:p=move(p,c)
 return p
D=8;found=[];collisions=0;longest=[]
for base,aa,bb in seeds():
 a,b=profile(aa),profile(bb)
 def bank(p):
  out=collections.defaultdict(list)
  def walk(p,w):
   out[X(p)].append((w,p))
   if len(w)==D:return
   for c in 'OZJ':walk(move(p,c),w+c)
  walk(p,'');return out
 A,B=bank(a),bank(b)
 for x in A.keys()&B.keys():
  for wa,pa in A[x]:
   for wb,pb in B[x]:
    if 'J' not in wa+wb:continue
    if not wa or not wb:continue
    collisions+=1
    k=1
    while k<=20:
     pa,pb=app(pa,wa),app(pb,wb);k+=1
     if X(pa)!=X(pb):break
    if k>2:longest.append((k,base,wa,wb,x))
    if k>20:
     found.append((base,wa,wb,x));print('PERSISTENT',found[-1],flush=True)
 print('base',base,'collisions',collisions,'persistent',len(found),flush=True)
print(json.dumps({'depth':D,'collisions':collisions,'persistent':found,'longest':sorted(longest,reverse=True)[:20]},indent=2),file=open('MACRO_SEARCH.json','w'))
