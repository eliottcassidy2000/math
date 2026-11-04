from pathlib import Path
import csv,itertools,json,random
D=Path(__file__).resolve().parent;P=D/'input'
TYPES={'cycle_source':([2,4,1,7],3),'cycle_sink':([10,12,9,0],3),'strong_four':([6,4,8,3],5)}
def need(x,m):
 if not x:raise RuntimeError(m)
def build(q,j,typ,a,c):
 k=len(q);cells=[list(range(3*i,3*i+3)) for i in range(k)];cells[j].append(3*k);n=3*k+3;rows=[0]*n
 for i,cell in enumerate(cells):
  template=TYPES[typ][0] if i==j else [2,4,1]
  for t,v in enumerate(cell):
   for u,w in enumerate(cell):
    if template[t]>>u&1:rows[v]|=1<<w
   for l,other in enumerate(cells):
    if q[i]>>l&1:rows[v]|=sum(1<<w for w in other)
   for port,mask in ((n-2,a),(n-1,c)):
    if mask>>i&1:rows[v]|=1<<port
    else:rows[port]|=1<<v
 rows[-2]|=1<<(n-1)
 return rows,cells

def compress(rows,cells,j,typ):
 if typ=='strong_four':blocks=[c for i,c in enumerate(cells) if i!=j]+[[v] for v in cells[j]]
 else:blocks=[c[:3] for c in cells]+[[cells[j][-1]]]
 blocks += [[len(rows)-2],[len(rows)-1]];cycles=sum(len(c)==3 for c in blocks);out=[0]*len(blocks)
 for i,ci in enumerate(blocks):
  for l,cl in enumerate(blocks):
   if l==i:continue
   bits={rows[v]>>w&1 for v in ci for w in cl};need(len(bits)==1,'uniform intergroup relation')
   if 1 in bits:out[i]|=1<<l
  if len(ci)==3:need(all(sum(rows[v]>>w&1 for w in ci)==1 for v in ci),'actual directed triangle')
 return cycles,out

def main():
 W=[{k:int(v) for k,v in x.items()} for x in csv.DictReader((P/'witnesses.tsv').open(),delimiter='\t')]
 by={x['node']:x['rows'] for x in json.load((P/'candidates.json').open())['candidates']}
 miss={(x['node'],x['position']) for x in json.load((P/'initial_misses.json').open())};selected={}
 for w in W:
  if (w['node'],w['position']) in miss:selected[w['node'],w['position']]=w
 for key in ['position','query_index']:
  seen=set()
  for w in W:
   if w[key] not in seen:selected[w['node'],w['position']]=w;seen.add(w[key])
 selected[W[-1]['node'],W[-1]['position']]=W[-1]
 records=[];lines=[]
 for w in selected.values():
  q=by[w['node']]
  for typ,(_,factor) in TYPES.items():
   for opposite in range(2):
    qq=[511^r^(1<<i) for i,r in enumerate(q)] if opposite else q
    rows,cells=build(qq,w['position'],typ,w['a'],w['c']);cycles,group=compress(rows,cells,w['position'],typ)
    idx=len(records);expected=w['normalized_target' if opposite else 'normalized_source']*factor%8
    records.append(dict(index=idx,witness=w,type=typ,opposite=bool(opposite),rows=rows,cycles=cycles,group_rows=group,expected=expected))
    lines.append(' '.join(map(str,[idx,cycles,len(group),0]+group)))
 (D/'large_physical.in').write_text(str(len(lines))+'\n'+'\n'.join(lines)+'\n');(D/'large_physical_records.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
 small=[];lines=[]
 for k in range(1,4):
  for code in range(1<<(k*(k-1)//2)):
   q=[0]*k;b=0
   for i in range(k):
    for j in range(i+1,k):
     if code>>b&1:q[i]|=1<<j
     else:q[j]|=1<<i
     b+=1
   for typ in TYPES:
    for a in range(1<<k):
     for c in range(1<<k):
      j=(a+c)%k;rows,cells=build(q,j,typ,a,c);cycles,groups=compress(rows,cells,j,typ);idx=len(small)
      # Every seventh case is also relabeled physically before literal count.
      if idx%7==0:
       perm=list(range(len(rows)));random.Random(idx).shuffle(perm)
       rel=[0]*len(rows)
       for i,r in enumerate(rows):
        for t in range(len(rows)):
         if r>>t&1:rel[perm[i]]|=1<<perm[t]
       rows=rel
      small.append(dict(index=idx,k=k,code=code,type=typ,position=j,a=a,c=c,physical_rows=rows,cycles=cycles,group_rows=groups))
      lines.append(' '.join(map(str,[2*idx,cycles,len(groups),0]+groups)));lines.append(' '.join(map(str,[2*idx+1,0,len(rows),1]+rows)))
 (D/'small_physical.in').write_text(str(len(lines))+'\n'+'\n'.join(lines)+'\n');(D/'small_physical_records.json').write_text(json.dumps(small,separators=(',',':'))+'\n')
 print(json.dumps(dict(large_placements=len(selected),large_graphs=len(records),small_graphs=len(small))))
if __name__=='__main__':main()
