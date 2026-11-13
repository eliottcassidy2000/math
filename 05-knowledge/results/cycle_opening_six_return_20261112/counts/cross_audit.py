"""Independent physical checks of sibling L/T laws and the two ready-return targets."""
from pathlib import Path
from hashlib import sha256
import json
import inflation_probe as ip
from audit_counts import birth,substitute,two_covers,require,valid
HERE=Path(__file__).resolve().parent

def state(rows,r):
 h,s,f=ip.literal(rows)
 old=[i for i in range(len(rows))if i!=r]
 deleted=tuple(sum((rows[a]>>b&1)<<j for j,b in enumerate(old))for a in old)
 d=ip.literal(deleted)[0]if old else 1
 require(sum(v>>r&1 for v in rows)==1,'state requires indegree-one root')
 require(f[r]-h==d-s[r],'indegree-one identity')
 return h-d,f[r]-h,d

def apply_word(rows,r,word):
 for op in word:
  oldr=r;rows=birth(rows,oldr);r=len(rows)-1
  if op=='T':rows=birth(rows,oldr);r=len(rows)-1
 return rows,r

def linear(word,v):
 x,b,d=v
 for op in word:
  if op=='L':x,b,d=x+2*b,d-b,x+d
  else:x,b,d=2*x+2*b,d,2*x+2*b+d
 return x,b,d

def funnel_target(module):
 # C3, module, three zero-tail vertices, one final one-tail vertex, funnel root, apex.
 seed=(2,4,1);m=len(module);n=m+8;root=n-2;last=n-3
 owners=[0,1,2]+[3]*m+[4,5,6,7,8]
 # Build k1 parent in seed,tail,new-root order, then substitute its tail start.
 base=[]
 for a in range(9):
  bits=0
  for b in range(9):
   if a==b:continue
   if a<3 and b<3:edge=seed[a]>>b&1
   elif a==8:edge=(b==0 or b in (4,5,6))
   elif b==8:edge=not(a==0 or a in (4,5,6))
   else:edge=a<b
   bits|=int(edge)<<b
  base.append(bits)
 p,idx=substitute(tuple(base),3,module)
 return birth(p,idx[8])

def main():
 report={'status':'VERIFIED independent physical checks','macro_sources':[],'universal_ready_targets':[]}
 C4=(2,4,9,3);P=(78,12,24,16,35,79,30)
 sources=[((2,4,1),0),(birth(C4,2),4),(birth(C4,3),4),(birth(P,4),7)]
 for rows,r in sources:
  v=state(rows,r);q=2*v[0]+2*v[1]+v[2]
  out={}
  for word in ('L','T','TTTLL','TTLTL'):
   result,rr=apply_word(rows,r,word);got=state(result,rr);require(got==linear(word,v),('L/T physical state',rows,r,word))
   out[word]=list(got)
  require(out['TTTLL']==[24*q,9*q,27*q]and out['TTLTL']==[22*q,12*q,29*q],'macro rank-one output')
  A,ra=apply_word(rows,r,'TTTLL');B,rb=apply_word(rows,r,'TTLTL');n=len(rows)
  changes=[(i,j)for i in range(len(A))for j in range(i+1,len(A))if (A[i]>>j&1)!=(B[i]>>j&1)]
  require(changes==[(n+3,n+5),(n+4,n+5),(n+4,n+6),(n+5,n+6)],'exact four arc switch')
  require(ra==rb==n+7 and all((A[i]>>ra&1)==(B[i]>>rb&1)for i in range(len(A))),'same final physical root and entire star')
  report['macro_sources'].append({'rows':rows,'root':r,'state':v,'q':q,'outputs':out,'changed_arcs':changes})
 Q=(14,12,8,64,175,79,535,367,127,447)
 for name,M in [('TT1',(0,)),('TT2',(2,0)),('C3',(2,4,1)),('C4',C4),('five_core',(22,4,8,3,14))]:
  hm=ip.literal(M)[0];pm=two_covers(M);q,idx=substitute(Q,5,M);f=funnel_target(M)
  require(len(q)==len(f),'same target size')
  both=[]
  for target in (q,f):
   valid(target);h,s,e=ip.literal(target)
   require(h==23*hm+135*pm and e[-1]-h==hm,('ready target',name))
   both.append({'H':h,'defect':e[-1]-h})
  old,idxold=substitute(P,1,M)
  for entry in (idxold[4],len(old),len(old)+1):old=birth(old,entry)
  ho,so,fo=ip.literal(old)
  require(ho==23*hm+135*pm and fo[-1]-ho==5*hm+24*pm,('raw excursion',name))
  base_map={0:0,1:3,2:5,3:6,4:7,5:2,6:4,7:1,8:8,9:9};m=len(M)
  lift={}
  for a,b in base_map.items():
   if a==1:
    for j in range(m):lift[1+j]=3+j
   else:lift[a+(m-1 if a>1 else 0)]=b+(m-1 if b>3 else 0)
  pulled=tuple(sum((f[lift[i]]>>lift[j]&1)<<j for j in range(len(old)))for i in range(len(old)))
  actual={(i,j)for i in range(len(old))for j in range(i+1,len(old))if (old[i]>>j&1)!=(pulled[i]>>j&1)}
  def oi(a):return a+(m-1 if a>1 else 0)
  expect={tuple(sorted((oi(a),oi(b))))for a,b in ((0,4),(0,7),(4,5),(4,7),(4,8),(5,8))}
  expect|={tuple(sorted((1+j,oi(b))))for j in range(m)for b in (4,6,8)}
  require(actual==expect and len(actual)==6+3*m,('named funnel return map',name))
  require(lift[len(old)-1]==len(old)-1 and lift[len(old)-2]==len(old)-2,'root and sole predecessor retained')
  report['universal_ready_targets'].append({'module':name,'H_module':hm,'P_module_2':pm,'targets':both,'funnel_changed_arcs':len(actual),'funnel_root_and_predecessor_fixed':True})
 report['script_sha256']=sha256(Path(__file__).read_bytes().replace(b"\r\n",b"\n")).hexdigest()
 (HERE/'cross_audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: L/T and A/B physical states on four sources, exact four-arc same-root-star switch, two ready-target constructions on five modules including nontransitive controls.')
if __name__=='__main__':main()
