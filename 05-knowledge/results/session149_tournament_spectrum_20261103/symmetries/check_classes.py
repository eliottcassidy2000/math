#!/usr/bin/env python3
"""Literal exhaustive class and converse census through order five."""
import itertools,json,math
from pathlib import Path

def pairs(n): return list(itertools.combinations(range(n),2))
def beats(mask,n,a,b):
    if a==b: return False
    i,j=sorted((a,b)); k=pairs(n).index((i,j)); return bool((mask>>k)&1) if a==i else not bool((mask>>k)&1)
def relabel(mask,n,p):
    return sum(int(beats(mask,n,p[i],p[j]))<<k for k,(i,j) in enumerate(pairs(n)))
def main():
 result={}
 for n in (3,4,5):
    ps=list(itertools.permutations(range(n))); limit=1<<(n*(n-1)//2); seen=set(); rows=[]
    for m in range(limit):
     if m in seen:continue
     orbit={relabel(m,n,p) for p in ps}; seen.update(orbit); c=min(orbit); conv=min(relabel(c^(limit-1),n,p) for p in ps)
     H=sum(all(beats(c,n,p[i],p[i+1]) for i in range(n-1)) for p in ps)
     aut=sum(relabel(c,n,p)==c for p in ps)
     score=sorted(sum(beats(c,n,i,j) for j in range(n)) for i in range(n))
     rows.append({'mask':c,'converse':conv,'paths':H,'aut':aut,'path_orbits':H//aut,'outdegrees':score})
    rows.sort(key=lambda a:a['mask'])
    if len(seen)!=limit: raise RuntimeError('Incomplete labelled census')
    idx={r['mask']:i for i,r in enumerate(rows)}
    for row in rows:
      mate=rows[idx[row['converse']]]
      if mate['converse']!=row['mask'] or mate['paths']!=row['paths'] or row['paths']%row['aut']!=0:
       raise RuntimeError('Converse or free-action check failed')
    expected={3:(2,2),4:(4,2),5:(12,8)}[n]
    if (len(rows),sum(r['mask']==r['converse'] for r in rows))!=expected: raise RuntimeError('Known count mismatch')
    result[n]={'rows':rows,'class_count':len(rows),'fixed':sum(r['mask']==r['converse'] for r in rows),'artificial_pairs':[[rows[i]['mask'],rows[i+1]['mask']] for i in range(0,len(rows),2)]}
    print('n=',n,'classes=',len(rows),'fixed=',result[n]['fixed'])
    for row in rows:print(row)
 Path(__file__).with_name('class_census.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
