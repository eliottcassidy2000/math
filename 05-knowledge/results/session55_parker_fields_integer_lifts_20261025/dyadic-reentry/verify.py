#!/usr/bin/env python3
"""Independent small-table verification of the dyadic re-entry theorem."""
import argparse, hashlib, itertools, json, pathlib
HERE=pathlib.Path(__file__).resolve().parent
BASE=HERE.parent/'crt-wiring'
ODD=[1,3,5,7,9,11,13,15,17,19,21,23,25,31,33,39,43,47,49,51,57,67,69,75,93,129,141,147,201]
EXCEPT={25,49,75,147}
PAIRS=list(itertools.combinations(range(9),2)); FULL=(1<<36)-1
LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
def require(ok,msg):
 if not ok: raise RuntimeError(msg)
def mask(a):return sum(1<<k for k,(i,j) in enumerate(PAIRS) if a[i]==a[j])
def squares(n):return sorted({r*r%n for r in range(n)})
def table(n):
 # Independent top-row reconstruction; no imports from the frozen checker.
 Q=squares(n); qs=set(Q); out={}; arrays=[]
 for a,b,c,e in itertools.product(Q,repeat=4):
  S=(a+b+c)%n
  if 3*e%n!=S:continue
  g=(S-c-e)%n;i=(S-a-e)%n;h=(S-b-e)%n
  d=(S-a-g)%n;f=(S-c-i)%n
  z=(a,b,c,d,e,f,g,h,i)
  if all(x in qs for x in z) and all(sum(z[j] for j in line)%n==S for line in LINES):
   arrays.append(z);out.setdefault(mask(z),z)
 return out,arrays
def minimum(masks):return sorted(x for x in masks if not any(y!=x and (x&y)==y for y in masks))
def factor(n):
 ans=[];p=2
 while p*p<=n:
  if n%p==0:
   q=1
   while n%p==0:n//=p;q*=p
   ans.append(q)
  p+=1
 if n>1:ans.append(n)
 return ans
def product_masks(moduli,tables):
 masks={FULL}
 for q in moduli:masks=set(minimum({x&y for x in masks for y in minimum(tables[q])}))
 return sorted(masks)
def build():
 for pin in json.loads((HERE/'SOURCE_PINS.json').read_text())['dependencies']:
  require(hashlib.sha256((HERE.parent/pin['path'].replace('session55-crt-collision-wiring/', 'crt-wiring/').replace('session55-finite-field-geometry/', 'finite-fields/')).read_bytes()).hexdigest()==pin['sha256'],f"source pin mismatch: {pin['path']}")
 source=json.loads((BASE/'classification_candidates.json').read_text())
 qs=sorted(set(q for m in ODD for q in factor(m))|{2,4,8,16,32,64})
 tables={};record={};arrays={}
 for q in qs:
  tables[q],arrays[q]=table(q)
  require(set(minimum(tables[q]))=={int(x['mask'],16) for x in source['local_antichains'][str(q)]},f'frozen local mismatch {q}')
  record[str(q)]={'square_residues':squares(q),'magic_array_count':len(arrays[q]),'all_masks':[hex(x) for x in sorted(tables[q])], 'minimal_masks':[{'mask':hex(x),'array':tables[q][x]} for x in minimum(tables[q])]}
 require(set(tables[4])=={FULL},'mod 4 constant')
 require(set(tables[2])==set(tables[8])==set(tables[16]),'binary mask equivalence')
 require(set(arrays[8])=={tuple(4*x for x in z) for z in arrays[2]}|{(1,)*9},'mod 8 exact arrays')
 require(set(arrays[16])=={tuple(1+8*x for x in z) for z in arrays[2]}|{(0,)*9,(4,)*9},'mod 16 exact arrays')
 odd=[];rows=[]
 for m in ODD:
  om=product_masks(factor(m),tables)
  intersections=[[hex(x&y) for y in om] for x in minimum(tables[2])]
  compatible=any(z=='0x0' for row in intersections for z in row)
  require(compatible==(m in EXCEPT),f'odd compatibility {m}')
  odd.append({'odd_part':m,'minimal_masks':[hex(x) for x in om], 'binary_intersections':intersections,'binary_compatible':compatible})
  for a in range(7):
   n=2**a*m;actual=0 not in product_masks(factor(n),tables)
   predicted=(m in ODD if a in [0,2] else m not in EXCEPT if a in [1,3,4] else m in [1,3])
   require(actual==predicted,f'profile mismatch {n}')
   old=next(x for x in source['rows'] if x['n']==n)
   require(actual==old['parker'],f'frozen global mismatch {n}')
   rows.append({'n':n,'a':a,'m':m,'parker':actual})
 P=sorted(r['n'] for r in rows if r['parker']);ps=set(P)
 counterpairs=sorted((n,d) for n in P for d in range(2,n) if n%d==0 and d not in ps)
 require(counterpairs==[(100,50),(196,98),(300,50),(300,150),(588,98),(588,294)],'divisor counterpairs')
 w=next(x for x in source['rows'] if x['n']==50)['witness']
 z=w['entries'];require(len(set(z))==9 and all(x in squares(50) for x in z),'50 squares')
 require(all(sum(z[j] for j in L)%50==28 for L in LINES),'50 line sums')
 require(all(r*r%50==x for r,x in zip(w['square_roots'],z)),'50 roots')
 require(len({x%2 for x in z})==2,'failed lift parity')
 # Audit the AP masks at the first dyadic threshold and at every exceptional odd candidate >=5.
 A=lambda x,y,z:(y,z,x,x,y,z,z,x,y)
 B=lambda x,y,z:tuple(A(x,y,z)[3*r+2-c] for r in range(3) for c in range(3))
 require(mask(A(1,9,17))&mask(B(1,25,49))==0,'orthogonal AP')
 for a in [5,6]:
  q=2**a;require(all(x in squares(q) for x in [1,9,17]),'dyadic AP')
  require(0 not in tables[q],'dyadic Parker 32/64')
 return {'schema':'dyadic-reentry-v1','source_dependency':'Full classification reductions and prime-field theorem remain external proved dependencies. This checker independently verifies every finite local table used here.','odd_parts':ODD,'exceptions':sorted(EXCEPT),'local_tables':record,'odd_binary_compatibility':odd,'candidate_rows':sorted(rows,key=lambda x:x['n']),'profile_counts':[sum(r['parker'] for r in rows if r['a']==a) for a in range(7)],'parker_count_including_1':len(P),'parker_moduli':P,'proper_divisor_counterpairs':counterpairs,'smallest_loss':{'from':50,'to':100,'witness':w,'line_sum':28,'obstruction':'Any magic square array modulo 100 reduces modulo 4 to a constant, so its parity is constant; this witness has mixed parity and cannot lift.'}}
def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args(); data=build();target=HERE/'certificate.json'
 data=json.loads(json.dumps(data))
 if a.verify:require(data==json.loads(target.read_text()),'certificate mismatch')
 else:target.write_text(json.dumps(data,indent=2)+'\n')
 print('PASS: independent local tables; 203 candidates; profile',data['profile_counts'],'; six divisor counterpairs; explicit failed lift')
if __name__=='__main__':main()
