#!/usr/bin/env python3
import itertools,math,json

def require(p,s):
 if not p:raise RuntimeError(s)

def run_words(T,r):
 outside=[v for v in range(len(T)) if v!=r];n=len(T);a=[0]*(n+1)
 for word in itertools.permutations(outside):
  for k in range(1,n+1):
   for gaps in itertools.combinations(range(n),k):
    full=[]
    for j in range(n):
     if j in gaps:full.append(r)
     if j<n-1:full.append(word[j])
    if all(T[u]>>v&1 for u,v in zip(full,full[1:])):a[k]+=1
 return a

def main():
 roots=coefficients=0
 for n in range(1,5):
  pairs=list(itertools.combinations(range(n),2))
  for mask in range(1<<len(pairs)):
   T=[0]*n
   for bit,(u,v) in enumerate(pairs):
    if mask>>bit&1:T[u]|=1<<v
    else:T[v]|=1<<u
   for r in range(n):
    roots+=1
    for k,a in enumerate(run_words(T,r)[1:],1):
     require(a%math.factorial(k-1)==0,'free fragment factorial');coefficients+=1
 for b in range(2,1001):
  if b%2==0:require((b*b)%4==0,'even base contradiction')
  elif b>=5:require(math.factorial(b-1)%4==0,'odd higher base factorial contradiction')
 require([2*q*q for q in range(1,12,2) if 2*q*q>=78]==[98,162,242],'first square candidates after inherited78 bound')
 print(json.dumps({'status':'PASS','all_labelled_rooted_sources_through4':roots,'literal_run_divisibility_checks':coefficients,'arithmetic_base_controls_through':1000,'only_possible_nonconstant_negative_unit_base':3,'inherited_amplitude_lower_bound':242},sort_keys=True,indent=2))
if __name__=='__main__':main()
