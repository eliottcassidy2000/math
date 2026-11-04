#!/usr/bin/env python3
from pathlib import Path
from itertools import product,permutations
import json,subprocess,tempfile
from automaton import kernel
from cover_automaton import cover_kernel
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def small_profile(T):
 n=len(T);q=[0]*6
 for word in permutations(range(n)):
  for cuts in range(1<<(n-1)):
   if all(cuts>>j&1 or T[word[j]]>>word[j+1]&1 for j in range(n-1)):q[1+cuts.bit_count()]+=1
 return q
def physical(blocks,a,c,rev):
 sizes=list(map(len,blocks));offs=[0,sizes[0],sizes[0]+sizes[1]];n=sum(sizes);R=[0]*(n+2);R[n]=1<<(n+1)
 for i,B in enumerate(blocks):
  nxt=(i+(2 if rev else 1))%3
  for j,r in enumerate(B):
   v=offs[i]+j;R[v]=(r<<offs[i])|(((1<<sizes[nxt])-1)<<offs[nxt])
   for mask,p in[(a,n),(c,n+1)]:
    if mask>>i&1:R[v]|=1<<p
    else:R[p]|=1<<v
 return R

def main():
 blocks=[[0],[2,0],[6,4,0],[2,4,1]];profiles=list(map(small_profile,blocks));K=[cover_kernel(1,3,r) for r in(False,True)]
 for rev in(False,True):need(K[rev][1]==kernel(1,3,rev),'one-path coefficientwise automata agree')
 graphs=[];records=[]
 for ids in product(range(4),repeat=3):
  B=[blocks[i] for i in ids];P=[profiles[i] for i in ids];expected=[]
  for rev in(False,True):
   graphs.append(physical(B,1,3,rev));expected.append([sum(co*P[0][e[0]]*P[1][e[1]]*P[2][e[2]] for e,co in K[rev][h].items()) for h in range(1,6)])
  records.append(dict(child_indices=ids,expected_exact_augmented_profiles=expected))
 with tempfile.TemporaryDirectory(prefix='cover-physical-') as td:
  exe=Path(td)/'cover';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'cover_count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';actual=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(actual)==len(graphs)*5,'complete physical cover coordinates')
 for k,z in enumerate(records):
  vals=[actual[10*k:10*k+5],actual[10*k+5:10*k+10]];need(vals==z['expected_exact_augmented_profiles'],'literal cover convolution equals run-cut automaton');z['literal_exact_augmented_profiles']=vals
 return dict(status='PASS',actual_child_rows=blocks,complete_small_profiles=profiles,actual_augmented_graphs=len(graphs),integer_cover_coordinates=len(actual),records=records)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
