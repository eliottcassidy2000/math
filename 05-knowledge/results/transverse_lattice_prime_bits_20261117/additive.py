"""Optimal finite additive discriminator of order-eight same-P pixels."""
from collections import defaultdict
from itertools import combinations
import csv, io, subprocess, argparse
from summary import N8, jets, need

def main():
 p=argparse.ArgumentParser();p.add_argument('repo'); p.add_argument('binary')
 p.add_argument('--ref',default='HEAD');args=p.parse_args()
 table=subprocess.check_output(['git','-C',args.repo,'show',args.ref+':'+N8+'/n8_classes.csv'])
 tab=list(csv.DictReader(io.StringIO(table.decode())))
 inp=''.join(' '.join(str(int(r['code'])>>(8*i)&255) for i in range(8))+'\n' for r in tab)
 raw=subprocess.check_output([args.binary,'8','rows'],input=inp.encode())
 bins=defaultdict(set)
 for l in raw.splitlines():
  _,_,*G=map(int,l.split()); need(len(G)==5)
  bins[jets(8,G)[0]].add((G[2],G[4]))
 dif={(a-c,b-d) for pts in bins.values() for (a,b),(c,d) in combinations(sorted(pts),2)}
 need(all((4*a+b)%29 for a,b in dif))
 count=0
 for d in range(1,29):
  for a in range(1,d+1):
   if d%a:continue
   c=d//a
   for b in range(a):
    count+=1
    need(any(y%c==0 and (x-b*(y//c))%a==0 for x,y in dif))
 print('all order-eight fibers: minimal additive quotient index 29;',
       count, 'smaller HNF lattices ruled out,',len(dif),'directed difference vectors checked')
if __name__=='__main__':main()
