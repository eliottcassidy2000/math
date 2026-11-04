from pathlib import Path
import csv,itertools,json,random
from model import D,P,FIVE_A,FIVE_B,FOURS,CYCLE,inflate,ambient,compress,profile,need

def readw():
 out=[]
 for r in csv.DictReader((P/'root_witnesses.tsv').open(),delimiter='\t'):
  out.append({k:(v if k=='type' else json.loads(v) if k=='positions' else int(v)) for k,v in r.items()})
 return out

def main():
 W=readw();by={(w['node'],w['configuration']):w for w in W};quotients={r['node']:r['rows'] for r in json.load((P/'candidates.json').open())['candidates']};selected={}
 for w in W:selected.setdefault(('config',w['configuration']),w)
 for miss in json.load((P/'root_initial_misses.json').open()):selected['miss',miss['node'],miss['configuration']]=by[miss['node'],miss['configuration']]
 for w in W:selected.setdefault(('query',w['query_index']),w)
 specs=[]
 for w in selected.values():
  if w['type']=='two_four':local={j:FOURS[0] for j in w['positions']}
  else:local={w['positions'][0]:FIVE_A if w['type']=='five_6' else FIVE_B}
  specs.append((w,local))
 # All nontransitive five types and all nine ordered four-type combinations.
 five=json.load((P/'five_child_records.json').open())
 for r in five:
  if not r['nontransitive']:continue
  j=r['code']%9;config=1+j if r['normalized_mod8'][1]==6 else 46+j
  specs.append((next(w for w in W if w['configuration']==config),{j:r['rows']}))
 pairw=next(w for w in W if w['type']=='two_four')
 for left in FOURS:
  for right in FOURS:specs.append((pairw,{pairw['positions'][0]:left,pairw['positions'][1]:right}))
 records=[];lines=[]
 for w,local in specs:
  factor=1
  for j in range(9):factor=factor*profile(local.get(j,CYCLE))[0]%8
  for opposite in range(2):
   q=quotients[w['node']];qq=[511^r^(1<<i) for i,r in enumerate(q)] if opposite else q
   rows,cells=inflate(qq,local,w['a'],w['c']);cycles,groups=compress(rows,cells);idx=len(records);expected=factor*w['normalized_target' if opposite else 'normalized_source']%8
   records.append(dict(index=idx,kind='root',witness=w,locals=local,opposite=bool(opposite),physical_rows=rows,cycles=cycles,group_rows=groups,expected=expected))
   lines.append(' '.join(map(str,[idx,cycles,len(groups),0]+groups)))
 amb=[{k:int(v) for k,v in r.items()} for r in csv.DictReader((P/'witnesses.tsv').open(),delimiter='\t')];ch={};
 for w in amb:ch.setdefault((w['a'],w['c']),w)
 ch[('extra',)]=next(w for w in amb if w['node']==56186)
 for w in ch.values():
  for opposite in range(2):
   q=quotients[w['node']];qq=[511^r^(1<<i) for i,r in enumerate(q)] if opposite else q
   rows,cells=ambient(qq,w['a'],w['c']);cycles,groups=compress(rows,cells);idx=len(records);expected=w['ambient_target' if opposite else 'ambient_source']
   records.append(dict(index=idx,kind='ambient',witness=w,opposite=bool(opposite),physical_rows=rows,cycles=cycles,group_rows=groups,expected=expected));lines.append(' '.join(map(str,[idx,cycles,len(groups),0]+groups)))
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
   for a in range(1<<k):
    for c in range(1<<k):
     specs=[{(a+c)%k:r['rows']} for r in five]
     if k>=2:specs += [{0:left,1:right} for left in FOURS for right in FOURS]
     for local in specs:
      rows,cells=inflate(q,local,a,c);cycles,groups=compress(rows,cells);idx=len(small)
      if idx%17==0:
       perm=list(range(len(rows)));random.Random(idx).shuffle(perm);rel=[0]*len(rows)
       for i,r in enumerate(rows):
        for j in range(len(rows)):
         if r>>j&1:rel[perm[i]]|=1<<perm[j]
       rows=rel
      small.append(dict(index=idx,k=k,code=code,a=a,c=c,locals=local,physical_rows=rows,cycles=cycles,group_rows=groups))
      lines.append(' '.join(map(str,[2*idx,cycles,len(groups),0]+groups)));lines.append(' '.join(map(str,[2*idx+1,0,len(rows),1]+rows)))
 (D/'small_physical.in').write_text(str(len(lines))+'\n'+'\n'.join(lines)+'\n');(D/'small_physical_records.json').write_text(json.dumps(small,separators=(',',':'))+'\n')
 print(json.dumps(dict(root_pairs=sum(r['kind']=='root' for r in records)//2,ambient_pairs=sum(r['kind']=='ambient' for r in records)//2,large_graphs=len(records),small_graphs=len(small))))
if __name__=='__main__':main()
