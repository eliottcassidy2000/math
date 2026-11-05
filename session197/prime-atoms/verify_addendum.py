from verify import *
from functools import lru_cache

@lru_cache(None)
def local_data(rows,out):
 rows=list(rows);n=len(rows); ps,H=profiles(rows);full=(1<<n)-1
 rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
 qs=profiles(rev)[0]
 h=H[full];Q=sum(H[S]*H[full^S] for S in range(full+1))
 U=sum(ps[x][2] for x in range(n) if out>>x&1)
 V=sum(qs[x][2] for x in range(n) if out>>x&1)
 E=sum(ps[x][1] for x in range(n) if out>>x&1)
 B=sum(qs[x][1] for x in range(n) if out>>x&1)
 T,r=add_apex(rows,out);hh,ee,RR,dd=profiles(T)[0][r]
 f=RR+ee-hh-dd
 return h,Q,U,V,E,B,f

def fast_profile(parts):
 # Data tuples (H,P2,U,V,E,S,f_local).
 n=len(parts);prefix=[1]*(n+1);suffix=[1]*(n+1);suffixQ=[1]*(n+1)
 for i in range(n):prefix[i+1]=prefix[i]*parts[i][0]
 for i in range(n-1,-1,-1):
  suffix[i]=suffix[i+1]*parts[i][0]
  suffixQ[i]=suffixQ[i+1]*parts[i][1]
 d=prefix[n]
 U=sum(prefix[i]*parts[i][2]*suffixQ[i+1] for i in range(n))
 f=sum(prefix[i]*parts[i][6]*suffix[i+1] for i in range(n))
 for i in range(n):
  between=1
  for j in range(i+1,n):
   f+=prefix[i]*parts[i][2]*between*parts[j][3]*suffix[j+1]
   if j==i+1:f+=prefix[i]*parts[i][4]*parts[j][5]*suffix[j+1]
   between*=parts[j][1]
 e=d-(prefix[n-1]*parts[-1][4] if n else 0)
 return U+e-f,e,d+U,d

def joined_apex(blocks,outs):
 rows=[];rootout=0;offset=0;N=sum(map(len,blocks))
 for block,out in zip(blocks,outs):
  after=((1<<N)-1)^((1<<(offset+len(block)))-1)
  rows.extend((r<<offset)|after for r in block)
  rootout|=out<<offset;offset+=len(block)
 return add_apex(rows,rootout)

def fast_classify(cores,d,f):
 k=(1+math.isqrt(1+8*f))//2;ell=max(1,f).bit_length()-1
 found={};examined=0;matching=0
 b=len(cores); corelen=sum(map(len,cores))
 # Ordered positions of nontrivial components among all components.
 for singles in range(k+ell+1):
  for cpos in itertools.combinations(range(singles+b),b):
   for masks in itertools.product(*[range(1<<len(C)) for C in cores]):
    core_out=sum(m.bit_count() for m in masks)
    if core_out>k:continue
    for w in range(1<<singles):
     if w.bit_count()+core_out>k:continue
     parts=[];blocks=[];outs=[];j=0;l=0
     for pos in range(singles+b):
      if j<b and pos==cpos[j]:C=cores[j];out=masks[j];j+=1
      else:C=[0];out=w>>l&1;l+=1
      blocks.append(C);outs.append(out);parts.append(local_data(tuple(C),out))
     nonempty=[i for i,o in enumerate(outs) if o]
     if not nonempty:continue
     a,z=nonempty[0],nonempty[-1]
     if any(len(C)==1 and not o and not (a<i<z) for i,(C,o) in enumerate(zip(blocks,outs))):continue
     examined+=1
     zz=fast_profile(parts);h,e,R,dd=zz
     require(dd==d,'d mismatch')
     if R+e-h-dd==f and R>d:
      matching+=1;U=R-d
      if U not in found:
       T,r=joined_apex(blocks,outs)
       require(profiles(T)[0][r]==zz,'fast witness mismatch')
       found[U]={'rows':T,'root':r,'profile':list(zz)}
 return examined,matching,found

def verify_formula(rows):
 checks=0
 for r,z in enumerate(profiles(rows)[0]):
  vs=[x for x in range(len(rows)) if x!=r]
  S=induced(rows,vs);parts=[]
  for C in components(S):
   block=induced(S,C)
   out=sum(1<<j for j,x in enumerate(C) if rows[r]>>vs[x]&1)
   parts.append(local_data(tuple(block),out))
  require(fast_profile(parts)==z,('component formula',rows,r,z,fast_profile(parts)))
  checks+=1
 return checks


def main():
 types={x:set() for x in [3,5,7,9,11,13]}
 checks=0; strong6=0; minimum6=10**100
 for n in range(1,7):
  for rows in labelled(n):
   if n<=5:checks+=verify_formula(rows)
   if n>1 and len(components(rows))==1:
    h=profiles(rows)[0][0][0]
    if n==6:
     strong6+=1;minimum6=min(minimum6,h)
     require(h>=15,'order-six strong floor')
    if h in types:types[h].add(canonical(rows))
 require({h:len(v) for h,v in types.items()}=={3:1,5:1,7:0,9:2,11:1,13:1},'small strong cores')
 require(minimum6==15,'sharp finite strong floor')
 require(profiles([36,1,2,7,15,30])[0][0][0]==15 and len(components([36,1,2,7,15,30]))==1,'sharp-floor witness')
 require(profiles([52,1,2,7,14,30])[0][5]==(19,1,73,9),'power-of-two hostile')
 rng=random.Random(197019)
 for n in [6,7,8,9]:
  for _ in range(30):
   rows=[0]*n
   for i,j in itertools.combinations(range(n),2):
    a,b=(i,j) if rng.randrange(2) else (j,i)
    rows[a]|=1<<b
   checks+=verify_formula(rows)
 expected={(9,1):[20,24,32,44],(9,-1):[26,27,34,35,41,45,99],
           (11,1):[29,37,47],(11,-1):[],(13,1):[],(13,-1):[]}
 data={};extensions=0
 for d in [9,11,13]:
  banks=[[list(C)] for C in sorted(types[d])]
  if d==9:banks.append([[4,1,2],[4,1,2]])
  for beta in [1,-1]:
   allfound={};tests=0;matches=0
   for bank in banks:
    n,m,f=fast_classify(bank,d,d-beta);tests+=n;matches+=m
    for u,w in f.items():allfound.setdefault(u,w)
   require(sorted(allfound)==expected[d,beta],('extended table',d,beta,allfound))
   for U,w in allfound.items():
    z=tuple(w['profile'])
    for m in [1,2]:
     T=k_graph(w['rows'],w['root'],m)
     actual=profiles(T)[0][w['root']]
     require(actual==k_profile(z,m),'actual K extension')
     require(actual[0]==U*(1<<m)+beta,'actual unit amplitude')
     extensions+=1
   data[f'{d}:{beta}']={'examined_reduced_cores':tests,'matching_reduced_cores':matches,
                      'U_base':sorted(allfound),'witnesses':allfound}
 bases={(1,1):[1],(1,-1):[3,5],(3,1):[6],(3,-1):[9,15],
        (5,1):[11,12,16],(5,-1):[15,39],(7,1):[],(7,-1):[]}
 bases.update(expected)
 allowed=[];budgets=[]
 for p in [3,5,7,11,13,17,19]:
  for beta in [-1,1]:
   D=(2*p+2*beta-1)//3
   ds=[d for d in range(1,D+1,2) if p in bases[d,beta]]
   budgets.append({'prime_amplitude':p,'offset':beta,'deletion_upper_bound':D,'witness_deletions':ds})
   if ds:allowed.append((p,beta))
 require(allowed==[(3,-1),(5,-1),(11,1)],'complete prime amplitude table')
 require(all(23 not in bases[d,-1] for d in range(1,14,2)),'amplitude23 negative offset')
 print(json.dumps({'status':'PASS','component_formula_root_checks':checks,
                   'labelled_strong_order6_count':strong6,'minimum_H_strong_order6':minimum6,
                   'small_strong_cores':{str(h):sorted(v) for h,v in types.items()},
                   'unit_source_table':data,'actual_K_extension_checks':extensions,
                   'prime_amplitude_budgets':budgets,'allowed_odd_prime_amplitudes_through19':allowed,
                   'scope':'All-order bounded-deletion classifications; full spectrum remains open.'},
                  indent=2,sort_keys=True))


if __name__=='__main__':
 main()
