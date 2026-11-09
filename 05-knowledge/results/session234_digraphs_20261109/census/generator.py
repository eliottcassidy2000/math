#!/usr/bin/env python3
"""Generate every semicomplete H7/H21 isomorphism type from classified SCC cores.
Examples: python generator.py --validate; python generator.py --target 21 --n 8 --limit 3
Loop-decorated generation is opt-in (--with-loops); JSON-lines matrices go to stdout.
Classification completeness uses the strong-tournament lower bound documented in REPORT.md.
"""
import argparse,itertools,json,math,collections,functools
from pathlib import Path
ROOT=Path(__file__).resolve().parent
# Exhaustively classified loopless strong cores; rows separated by /.
CORE_ROWS={3: ['010/001/100', '010/101/100'], 7: ['0110/1011/1001/1000', '0111/0011/1001/1000', '0011/1011/0101/1000'], 21: ['01110/10101/10011/11001/10000', '01111/10101/10011/01001/10000', '01111/10101/01011/01001/10000', '01111/00101/01011/11001/10000', '01100/00111/01011/11001/10000', '00111/10111/01011/01001/10000', '01100/10110/10011/10001/11000', '01111/00110/00011/10001/11000', '01011/00110/10011/10001/11000', '01110/10010/11011/10001/11000', '01111/10010/01011/10001/11000', '01110/00110/01011/10001/11000', '00110/10110/01011/10001/11000', '01011/10100/11011/01001/11000', '00011/10100/11011/11001/11000', '00101/10010/11011/11001/11000', '011110/001111/000111/100011/000001/100000', '011111/001111/000111/100011/000001/100000', '011010/101111/100111/100011/000001/100000', '011011/101111/100111/100011/000001/100000', '001010/101011/110111/110011/000001/100000', '000110/101011/110111/110011/000001/100000', '001011/101011/110111/110011/000001/100000', '000111/101011/110111/110011/000001/100000', '010111/001111/100111/000011/100001/100000', '010111/101111/100111/000011/100001/100000', '010100/001111/100111/100011/100001/100000', '010100/101111/100111/100011/100001/100000', '010010/101011/100111/110011/100001/100000', '010010/101011/110111/110011/100001/100000', '000011/101011/100111/110011/010001/100000', '000011/101011/110111/110011/010001/100000']}

def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for j in range(n):
   if s>>j&1:
    p=s^(1<<j)
    if p:d[s][j]=sum(d[p][i] for i in range(n) if p>>i&1 and a[i][j])
 return sum(d[-1])

def automorphisms(a):
 n=len(a)
 return [p for p in itertools.permutations(range(n)) if all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n))]

def cycle_type(p):
 remaining=set(range(len(p)));lengths=[]
 while remaining:
  i=min(remaining);j=i;length=0
  while j in remaining:remaining.remove(j);length+=1;j=p[j]
  lengths.append(length)
 return tuple(sorted(lengths))

def mask_transform(mask,p):return sum(((mask>>p[i])&1)<<i for i in range(len(p)))

def load_cores():
 raw={h:[[[int(v) for v in row] for row in encoded.split('/')] for encoded in values] for h,values in CORE_ROWS.items()}
 out={}
 for h,ms in raw.items():
  out[h]=[]
  for k,a in enumerate(ms):
   aut=automorphisms(a);masks=sorted({min(mask_transform(m,p) for p in aut) for m in range(1<<len(a))})
   assert hp(a)==h
   out[h].append({'id':f'{h}:{len(a)}:{k}','matrix':a,'aut':aut,'loop_masks':masks})
 return out

def compositions(total,parts):
 if parts==1:
  yield (total,);return
 for x in range(total+1):
  for rest in compositions(total-x,parts-1):yield (x,)+rest

def decorate(c,mask):
 a=[r[:] for r in c['matrix']]
 for i in range(len(a)):a[i][i]=mask>>i&1
 return a

def join(blocks):
 n=sum(len(a) for a in blocks);a=[[0]*n for _ in range(n)];k=0
 for b in blocks:
  for i in range(len(b)):
   for j in range(len(b)):a[k+i][k+j]=b[i][j]
   for j in range(k+len(b),n):a[k+i][j]=1
  k+=len(b)
 return a

def generate(cores,target,n,loops=False):
 families=[(c,) for c in cores[target]]
 if target==21:
  for a in cores[3]:
   for b in cores[7]:families.extend(((a,b),(b,a)))
 for family in families:
  left=n-sum(len(c['matrix']) for c in family)
  if left<0:continue
  mask_lists=[c['loop_masks'] if loops else [0] for c in family]
  for gaps in compositions(left,len(family)+1):
   for masks in itertools.product(*mask_lists):
    for singletons in itertools.product(range(2) if loops else range(1),repeat=left):
     blocks=[];key=[];j=0
     for i,gap in enumerate(gaps):
      for _ in range(gap):blocks.append([[singletons[j]]]);key.append(('singleton',singletons[j]));j+=1
      if i<len(family):blocks.append(decorate(family[i],masks[i]));key.append((family[i]['id'],masks[i]))
     yield tuple(key),join(blocks)

def expected(target,n,loops=False):
 def term(k,c,p):
  return 0 if n<k else c*math.comb(n-k+p-1,p-1)*(2 if loops else 1)**(n-k)
 if target==7:return term(4,48 if loops else 3,2)
 return term(5,480 if loops else 16,2)+term(6,1024 if loops else 16,2)+term(7,1152 if loops else 12,3)

def scc(a):
 n=len(a);r=[row[:] for row in a]
 for i in range(n):r[i][i]=1
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k]&r[k][j]
 left=set(range(n));cs=[]
 while left:
  i=min(left);c=sorted(j for j in left if r[i][j] and r[j][i]);left-=set(c);cs.append(c)
 # Generated forward ordinal sums have SCCs in increasing vertex order.
 return cs

@functools.lru_cache(None)
def canonical_block(a):
 n=len(a)
 return min(''.join(str(a[p[i]][p[j]]) for i in range(n) for j in range(n)) for p in itertools.permutations(range(n)))

def expected_labelled(h,n):
 f=math.factorial(n)
 if h==7:return 0 if n<4 else 3*(n-3)*f
 return (0 if n<5 else 44*(n-4)*f//3)+(0 if n<6 else 16*(n-5)*f)+(0 if n<7 else 4*(n-5)*(n-6)*f)

def validate(cores):
 report={'decorated_core_census':[],'loopless_family_checks':[]}
 for h,cs in cores.items():
  for n in sorted({len(c['matrix']) for c in cs}):
   subset=[c for c in cs if len(c['matrix'])==n];decorated=0;details=[]
   for c in subset:
    aut=c['aut'];orbits={};fixed_total=0
    for mask in range(1<<n):
     orbit=frozenset(mask_transform(mask,p) for p in aut)
     orbits.setdefault(min(orbit),orbit)
    assert len(orbits)==len(c['loop_masks'])
    assert sum(len(o) for o in orbits.values())==2**n
    for p in aut:fixed_total+=sum(mask_transform(m,p)==m for m in range(1<<n))
    assert fixed_total//len(aut)==len(orbits) and fixed_total%len(aut)==0
    for mask in c['loop_masks']:assert hp(decorate(c,mask))==h
    decorated+=len(orbits);details.append({'id':c['id'],'automorphism_order':len(aut),'automorphisms':aut,'automorphism_cycle_types':dict(collections.Counter(str(cycle_type(p)) for p in aut)),'loop_orbit_count':len(orbits),'orbit_sizes':dict(collections.Counter(map(len,orbits.values())))})
   report['decorated_core_census'].append({'target':h,'n':n,'unlooped_cores':len(subset),'decorated_cores':decorated,'details':details})
 for h in (7,21):
  for n in range(1,9):
   seen=set();actual_seen=set();labelled=0
   for key,a in generate(cores,h,n):
    assert key not in seen;seen.add(key);assert hp(a)==h
    assert all(a[i][j] or a[j][i] for i in range(n) for j in range(i+1,n))
    components=scc(a);assert len(components)==len(key)
    actual=tuple(canonical_block(tuple(tuple(a[i][j] for j in c) for i in c)) for c in components)
    assert actual not in actual_seen;actual_seen.add(actual)
    order_by_id={c['id']:len(c['aut']) for cs in cores.values() for c in cs}
    aut=math.prod(1 if x[0]=='singleton' else order_by_id[x[0]] for x in key)
    labelled+=math.factorial(n)//aut
   assert len(seen)==expected(h,n)
   assert labelled==expected_labelled(h,n)
   report['loopless_family_checks'].append({'target':h,'n':n,'isomorphism_classes':len(seen),'labelled_count':labelled,'every_matrix_hp_verified':True,'distinct_ordered_scc_signatures':True})
 return report

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--target',type=int,choices=(7,21));p.add_argument('--n',type=int);p.add_argument('--with-loops',action='store_true');p.add_argument('--limit',type=int);p.add_argument('--validate',action='store_true');args=p.parse_args();cores=load_cores()
 if args.validate:
  report=validate(cores);(ROOT/'generator_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));return
 if args.target is None or args.n is None:p.error('supply --target and --n, or --validate')
 if args.n<1:p.error('--n must be positive')
 records=generate(cores,args.target,args.n,args.with_loops)
 if args.limit is not None:records=itertools.islice(records,args.limit)
 for key,a in records:print(json.dumps({'target':args.target,'n':args.n,'with_loops':args.with_loops,'scc_signature':key,'matrix':a}))
if __name__=='__main__':main()
