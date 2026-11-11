#!/usr/bin/env python3
"""Exact fixed-core transfer: bit j of a color says core j beats spine vertex."""
from itertools import permutations,product
from functools import lru_cache
from pathlib import Path
import json,argparse

@lru_cache(None)
def machine(core):
 k=len(core); states=[]
 for p in permutations(range(k)):
  states.extend((p,s) for s in product(*([range(3)]*k+[range(2)])))
 index={x:i for i,x in enumerate(states)}; init=[0]*len(states); accept=[]
 for i,(p,s) in enumerate(states):
  if not any(s):init[i]=1
  accept.append(all(s[j]!=2 and (s[j] or (core[p[j-1]]>>p[j]&1)) for j in range(1,k)) and s[0]!=2)
 transitions=[]
 for color in range(1<<k):
  edges=[]
  for i,(p,s) in enumerate(states):
   for slot in range(k+1):
    if not s[slot] and slot and not(color>>p[slot-1]&1):continue
    t=list(s);t[slot]=1 if slot==k or not(color>>p[slot]&1) else 2
    edges.append((i,index[p,tuple(t)]))
  transitions.append(edges)
 return init,list(map(int,accept)),transitions

def step(v,edges):
 out=[0]*len(v)
 for i,j in edges:out[j]+=v[i]
 return tuple(out)

def value(word,core=(2,0)):
 v,a,M=machine(core)
 for c in word:v=step(v,M[c])
 return sum(x*y for x,y in zip(v,a))

def graph(word,core=(2,0)):
 k=len(core);m=len(word);r=list(core)+[0]*m
 for i,c in enumerate(word):
  for j in range(k):
   if c>>j&1:r[j]|=1<<(k+i)
   else:r[k+i]|=1<<j
  for j in range(i+1,m):r[k+i]|=1<<(k+j)
 return r

def count(rows):
 n=len(rows);dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,x in dp[mask].items():
   nxt=rows[i]&~mask
   while nxt:
    b=nxt&-nxt;j=b.bit_length()-1;nxt-=b
    d=dp[mask|b];d[j]=d.get(j,0)+x
 return sum(dp[-1].values())

def controls():
 checked=0
 for k,core,bound in [(1,(0,),7),(2,(2,0),5),(3,(2,4,1),3)]:
  for m in range(bound+1):
   for word in product(range(1<<k),repeat=m):
    a=value(word,core);b=count(graph(word,core))
    if a!=b:raise RuntimeError((word,core,a,b))
    checked+=1
 return {'independent_subset_dp_controls':checked,'core_orders':[1,2,3]}

def cpp():
 v,a,M=machine((2,0));D=len(v)
 s='#include <array>\n#include <vector>\n#include <iostream>\n#include <cstdint>\n#include <algorithm>\nusing namespace std;\n'
 s+=f'const int D={D}; using V=array<uint64_t,D>;\n'
 s+='const V initial={'+','.join(map(str,v))+'};\n'
 s+='const V accept={'+','.join(map(str,a))+'};\n'
 s+='const vector<pair<int,int>> E[4]={'+','.join('{'+','.join('{%s,%s}'%e for e in es)+'}' for es in M)+'};\n'
 s+=r'''
int depth,limit;vector<int> cost;vector<uint64_t> witness;vector<uint64_t> distinct;uint64_t visited=0;
void search(const V& v,int m,uint64_t word){
 uint64_t h=0;for(int i=0;i<D;i++)h+=v[i]*accept[i];visited++;
 if(h<=uint64_t(limit)&&cost[h]<0){cost[h]=m+2;witness[h]=word;distinct[m]++;}
 if(m==depth)return;
 for(int c=0;c<4;c++){V w={};for(auto e:E[c])w[e.second]+=v[e.first];search(w,m+1,(word<<2)|c);}
}
int main(int argc,char**argv){depth=argc>1?atoi(argv[1]):10;limit=argc>2?atoi(argv[2]):100000;
 if(depth<0||depth>14||limit<1){cerr<<"supported census depth: 0..14; limit must be positive\n";return 2;}
 cost.assign(limit+1,-1);witness.resize(limit+1);distinct.resize(depth+1);
 // Iterative deepening means recorded costs are minimum within this word language.
 for(int d=0;d<=depth;d++){int saved=depth;depth=d;search(initial,0,0);depth=saved;
 int total=0,first=0;for(int h=1;h<=limit;h+=2){if(cost[h]>=0)total++;else if(!first&&h!=7&&h!=21)first=h;}
 cerr<<d<<" "<<total<<" "<<first<<"\n";}
 cout<<"{\"max_spine\":"<<depth<<",\"limit\":"<<limit<<",\"visited_prefixes_with_repetition\":"<<visited<<",\"witnesses\":[";
 bool first=true;for(int h=1;h<=limit;h+=2)if(cost[h]>=0){if(!first)cout<<",";first=false;cout<<"["<<h<<","<<cost[h]<<","<<witness[h]<<"]";}cout<<"]}\n";
}
'''
 return s.rstrip()+'\n'
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cpp',action='store_true');a=p.parse_args()
 print(cpp() if a.cpp else json.dumps(controls(),indent=2),end='' if a.cpp else '\n')
