from pathlib import Path
import json,subprocess,random,hashlib
D=Path(__file__).resolve().parent

def need(a,m):
 if not a: raise RuntimeError(m)
def augment(q,a,c,flips=0,perm=None):
 k=len(q);n=3*k;rows=[0]*(n+2)
 for i in range(n):
  bi,pi=divmod(i,3)
  for j in range(n):
   if i==j:continue
   bj,pj=divmod(j,3)
   bit=((pj-pi)%3==(2 if flips>>bi&1 else 1)) if bi==bj else (q[bi]>>bj&1)
   if bit: rows[i]|=1<<j
  if a>>bi&1:rows[i]|=1<<n
  else:rows[n]|=1<<i
  if c>>bi&1:rows[i]|=1<<(n+1)
  else:rows[n+1]|=1<<i
 rows[n]|=1<<(n+1)
 if perm is not None:
  z=[0]*len(rows)
  for i,r in enumerate(rows):
   for j in range(len(rows)):
    if r>>j&1:z[perm[i]]|=1<<perm[j]
  rows=z
 return rows

def main():
 rng=random.Random(18927);cases=[]
 for k in (1,2,3):
  for code in range(1<<(k*(k-1)//2)):
   q=[0]*k;bit=0
   for i in range(k):
    for j in range(i+1,k):
     if code>>bit&1:q[i]|=1<<j
     else:q[j]|=1<<i
     bit+=1
   for a in range(1<<k):
    for c in range(1<<k):cases.append(dict(k=k,quotient=q.copy(),a=a,c=c))
 for _ in range(32):
  k=4;q=[0]*k
  for i in range(k):
   for j in range(i+1,k):
    if rng.randrange(2):q[i]|=1<<j
    else:q[j]|=1<<i
  cases.append(dict(k=k,quotient=q,a=rng.randrange(16),c=rng.randrange(16)))
 orbit_in=str(len(cases))+'\n'+''.join(' '.join(map(str,[i,x['k'],x['a'],x['c']]+x['quotient']))+'\n' for i,x in enumerate(cases))
 (D/'small_orbit.in').write_text(orbit_in)
 orbit=subprocess.check_output([str(D/'orbit_physical')],input=orbit_in.encode(),stderr=subprocess.DEVNULL).decode();(D/'small_orbit.out').write_text(orbit)
 for i,x in enumerate(cases):
  k=x['k'];perm=list(range(3*k+2));rng.shuffle(perm);x['physical_permutation']=perm;x['sides']=[]
  for side in (0,1):
   q=x['quotient'] if side==0 else [((1<<k)-1)^r^(1<<j) for j,r in enumerate(x['quotient'])]
   flips=rng.randrange(1<<k);rows=augment(q,x['a'],x['c'],flips,perm)
   need(all(not(r>>j&1) and all(((r>>l)&1)+((rows[l]>>j)&1)==1 for l in range(len(rows)) if l!=j) for j,r in enumerate(rows)),'physical tournament')
   x['sides'].append(dict(child_orientation_flips=flips,physical_rows=rows))
 li=str(2*len(cases))+'\n'+''.join(str(len(s['physical_rows']))+' '+' '.join(map(str,s['physical_rows']))+'\n' for x in cases for s in x['sides'])
 (D/'small_literal.in').write_text(li)
 lo=subprocess.check_output([str(D/'literal_physical')],input=li.encode()).decode();(D/'small_literal.out').write_text(lo)
 nums=list(map(int,lo.split()));orows=[list(map(int,l.split())) for l in orbit.splitlines()]
 need(len(nums)==2*len(cases) and len(orows)==len(cases),'complete literal/orbit output')
 for i,x in enumerate(cases):
  need(orows[i][0]==i,'orbit record alignment')
  for side in (0,1):
   x['sides'][side]['H']=nums[2*i+side];x['sides'][side]['mod8']=nums[2*i+side]%8
   need(nums[2*i+side]%8==orows[i][side+1],'physical orbit formula vs literal vertices')
 (D/'small_physical_controls.json').write_text(json.dumps(cases,separators=(',',':'))+'\n')
 receipt=dict(status='PASS',query_pairs=len(cases),complete_quotient_query_pairs_through3=548,random_order4_pairs=32,literal_physical_sides=2*len(cases),maximum_physical_order=14,independent_random_child_converses=True,common_arbitrary_physical_permutations=True,sha256={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['orbit_physical.cpp','literal_physical.cpp','small_physical_controls.json','small_orbit.in','small_orbit.out','small_literal.in','small_literal.out']})
 print(json.dumps(receipt,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
