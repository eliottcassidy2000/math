#!/usr/bin/env python3
import json,itertools
from independent_guard_leakage import rank_mod
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def det(a):
 return sum((-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))*__import__('math').prod(a[i][p[i]] for i in range(len(p))) for p in itertools.permutations(range(len(a))))
def main():
 B=[[2*(i==j)+int((i-j)%5 in(1,4)) for j in range(5)] for i in range(5)]
 P=[[1,0,0,1],[-1,0,1,0],[0,-1,-1,0],[0,0,0,-1],[0,1,0,0]]
 G=[[1,1,0,0],[1,2,0,0],[0,0,1,1],[0,0,1,2]]
 ck(mm(B,P)==mm(P,G),'integral intertwining')
 ck(all(sum(c)==0 for c in zip(*P)),'augmentation')
 ck(det(P[:4])==-1,'unimodular augmentation')
 T=[[1]+r for r in P];ck(det(T)==-5,'ambient index')
 ck([r[0]%5 for r in mm(P,[[2],[1],[3],[-1]])]==[1]*5,'constant in augmentation')
 D=[[B[i][j]+int(i==j) for j in range(5)] for i in range(5)];D2=mm(D,D);D3=mm(D2,D)
 ck(all(x%5==1 for r in D2 for x in r),'square allones')
 ck(all(x%5==0 for r in D3 for x in r),'cube zero')
 ck([rank_mod(x,5) for x in(D,D2,D3)]==[3,1,0],'ambient ranks')
 ck([rank_mod(mm(x,P),5) for x in(D,D2)]==[2,0],'augmentation ranks')
 print(json.dumps({'status':'PASS','checks':N},indent=2))
if __name__=='__main__':main()
