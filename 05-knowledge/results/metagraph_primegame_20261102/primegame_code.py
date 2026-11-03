#!/usr/bin/env python3
"""Exact arithmetic-return code; formal circulation and executable words differ."""
from fractions import Fraction
from itertools import combinations,product
from pathlib import Path
from math import prod,isqrt
import json
HERE=Path(__file__).resolve().parent
PRIMES=[2,3,5,7,11,13,17,19,23,29]
FRACTIONS=[(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
LETTERS='ABCDEFGHIJKLMN'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def exponents(n):
 out=[]
 for p in PRIMES:
  a=0
  while n%p==0:n//=p;a+=1
  out.append(a)
 need(n==1,'supported factors');return out
def rank(A,p=None):
 a=[[Fraction(x) if p is None else x%p for x in row] for row in A];r=0
 for j in range(len(a[0])):
  hit=next((i for i in range(r,len(a)) if a[i][j]),None)
  if hit is None:continue
  a[r],a[hit]=a[hit],a[r];v=a[r][j]
  a[r]=[x/v for x in a[r]] if p is None else [x*pow(v,-1,p)%p for x in a[r]]
  for i in range(len(a)):
   if i!=r:
    t=a[i][j];a[i]=[x-t*y for x,y in zip(a[i],a[r])]
    if p:a[i]=[x%p for x in a[i]]
  r+=1
  if r==len(a):break
 return r
def vector(g,j,t,m):return [g+2*t+m,g+t+m,t,g,g+t+2*m,g+t+2*m,g,t,0,j,t+j,t+m,m,0]
def step(n):
 for i,(a,b) in enumerate(FRACTIONS):
  if n%b==0:return LETTERS[i],n//b*a
 raise RuntimeError('no enabled fraction')
def det3(c):
 a,b,d=c
 return a[0]*(b[1]*d[2]-b[2]*d[1])-a[1]*(b[0]*d[2]-b[2]*d[0])+a[2]*(b[0]*d[1]-b[1]*d[0])
def main():
 cols=[[u-v for u,v in zip(exponents(a),exponents(b))] for a,b in FRACTIONS]
 A=list(map(list,zip(*cols)))
 B=[vector(1,0,0,0),vector(0,1,0,0),vector(0,0,1,0),vector(0,0,0,1)]
 need(rank(A)==10 and rank(B)==4,'rational ranks')
 for p in [2,3,5,7,11,13]:need(rank(A,p)==10 and rank(B,p)==4,'finite field ranks')
 need(all(sum(x*y for x,y in zip(row,b))==0 for row in A for b in B),'kernel basis')
 entries=sorted(set(a for a,b in FRACTIONS)|set(b for a,b in FRACTIONS))
 incidence=[[int(a==v)-int(b==v) for a,b in FRACTIONS] for v in entries]
 valuation=list(map(list,zip(*(exponents(v) for v in entries))))
 need(len(entries)==19 and rank(incidence)==13,'entry graph rank')
 need([[sum(valuation[i][k]*incidence[k][j] for k in range(19)) for j in range(14)] for i in range(10)]==A,'valuation factors through entry incidence')
 kernel_images=[[sum(row[i]*b[i] for i in range(14)) for row in incidence] for b in B]
 need(all(v==0 for v in kernel_images[1]) and rank([kernel_images[i] for i in [0,2,3]])==3,'reciprocal quotient has three arithmetic directions')
 probe_groups={}
 for index,letter in enumerate(LETTERS):
  if letter in 'JK':continue
  form=tuple(B[i][index] for i in [0,2,3])
  if any(form):probe_groups.setdefault(form,[]).append(letter)
 need(len(probe_groups)==7,'exactly seven distinct nonzero quotient coordinate observations')
 selected=list('DCEMALB');atoms=[7,2,5,19,89,59,29]
 G=[[b[LETTERS.index(k)] for k in selected] for b in [B[0],B[2],B[3]]]
 columns=list(map(list,zip(*G)))
 labels=[sum((x%2)<<i for i,x in enumerate(c)) for c in columns]
 need(labels==list(range(1,8)),'all Fano points')
 def word(u):return tuple(sum(u[i]*G[i][j] for i in range(3))%2 for j in range(7))
 simplex={word(u) for u in product([0,1],repeat=3)}
 need(len(simplex)==8 and sorted(sum(w) for w in simplex)==[0]+[4]*7,'simplex parameters')
 accepted=[];code_weights={};equivalences=0
 for x in product([0,1],repeat=7):
  even_all=all(sum(a*b for a,b in zip(x,w))%2==0 for w in simplex)
  sigma=prod((p+1) for p,b in zip(atoms,x) if b);square=isqrt(sigma)**2==sigma
  need(even_all==square,'square-sigma iff universal formal-return parity');equivalences+=1
  if even_all:accepted.append(x);code_weights[sum(x)]=code_weights.get(sum(x),0)+1
 need(code_weights=={0:1,3:7,4:7,7:1},'Hamming parameters')
 circuits=[]
 for w in accepted:
  if any(w) and not any(any(v) and v!=w and all(a<=b for a,b in zip(v,w)) for v in accepted):circuits.append([selected[i] for i,b in enumerate(w) if b])
 need(len(circuits)==14,'all circuits, not just seven lines')
 lines=[]
 for ids in combinations(range(7),3):
  if labels[ids[0]]^labels[ids[1]]^labels[ids[2]]==0:
   lines.append(dict(instructions=[selected[i] for i in ids],integer_columns=[columns[i] for i in ids],determinant=det3([columns[i] for i in ids])))
 cycle=[];n=4301
 for _ in range(6):k,n=step(n);cycle.append(k)
 need(n==4301 and ''.join(cycle)=='GBDEFA','actual six-cycle control')
 need(step(11)==('K',13) and step(13)==('J',11),'reciprocal control')
 # The actual periodic words inherited from Session129 use only g and j.
 actual_accept=sum(sum(x[i]*(B[0][LETTERS.index(k)]%2) for i,k in enumerate(selected))%2==0 for x in product([0,1],repeat=7))
 need(actual_accept==64,'priority loses two parity tests')
 output=dict(status='PASS',entry_values=entries,entry_incidence_rank=13,arithmetic_quotient_rank=3,quotient_probe_groups=[dict(form=list(k),letters=v) for k,v in sorted(probe_groups.items())],prime_order=PRIMES,fractions=FRACTIONS,increment_matrix=A,kernel_basis_gjtm=B,
             selected_letters=selected,arithmetic_primes=atoms,integer_probe_matrix=G,parity_column_labels=labels,
             simplex_words=sorted(simplex),hamming_weight_enumerator=code_weights,circuits=circuits,
             fano_lines_with_integer_lifts=lines,exact_square_sigma_equivalences=equivalences,
             actual_cycle_word=''.join(cycle),actual_periodic_parity_acceptances=actual_accept,
             hostile=dict(selection='C',arithmetic_prime=2,sigma=3,formal_t_count=1,actual_periodic_count=0),
             scope='Formal arithmetic return multisets; actual periodic words rely on the inherited complete priority classification')
 (HERE/'primegame_code.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n');print('PASS rank10/kernel4; Fano7; Hamming16; actual-periodic64')
if __name__=='__main__':main()
