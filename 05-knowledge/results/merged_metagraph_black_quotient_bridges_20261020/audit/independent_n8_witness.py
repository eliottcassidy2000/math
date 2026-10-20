#!/usr/bin/env python3
"""Brute-force all 8! permutations for two order-eight blue-cap hostiles."""
from itertools import permutations
from math import comb
import json
from pathlib import Path
OUT=Path(__file__).resolve().parent
rows=[]
for mask in [40077,40109]:
 n=8
 tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
 A=[[0]*n for _ in range(n)]
 for i in range(n-1):A[i][i+1]=1
 for k,(x,y) in enumerate(tiles):
  i,j=n-x,n-y
  if mask&(1<<k):A[j][i]=1
  else:A[i][j]=1
 H=Aut=blue_paths=0;blue_words=[]
 for p in permutations(range(n)):
  if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n)):Aut+=1
  if not all(A[p[i]][p[i+1]] for i in range(n-1)):continue
  H+=1
  if all(A[p[i]][p[j]]==A[p[n-1-j]][p[n-1-i]] for i in range(n) for j in range(n)):
   blue_paths+=1
   blue_words.append(p)
 scores=list(map(sum,A));c3=comb(n,3)-sum(comb(s,2) for s in scores)
 assert H%Aut==0 and blue_paths%Aut==0
 row={'mask':mask,'n':n,'outneighbors':[[j for j in range(n) if A[i][j]] for i in range(n)],'scores':scores,'c3':c3,'HP':H,'Aut':Aut,'fibre':H//Aut,'blue_HP':blue_paths,'blue_fibre':blue_paths//Aut,'claimed_cap':3**((n-2)//2)}
 assert row['blue_fibre']==39>row['claimed_cap']==27
 assert (H,Aut,c3)==((637,1,20) if mask==40077 else (621,3,20))
 rows.append(row);print(json.dumps(row),flush=True)
(OUT/'independent_n8_witness.json').write_text(json.dumps(rows,indent=2)+'\n')
