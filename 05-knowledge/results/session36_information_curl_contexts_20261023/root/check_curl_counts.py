#!/usr/bin/env python3
import json, math, itertools
from collections import defaultdict
from check_curl_decoder import matrix,key,hp,check

def logarithmic_counts(a):
    c=[0]*len(a)
    for n in range(1,len(a)):
        c[n]=a[n]-sum(math.comb(n-1,j-1)*c[j]*a[n-j] for j in range(1,n))
    return c

def main():
    a=[1]+[2**math.comb(n,2) for n in range(1,13)]
    b=[1]+[math.factorial(n)*2**((n-1)*(n-2)//2) for n in range(1,13)]
    c=logarithmic_counts(a); d=logarithmic_counts(b); rows=[]
    for n in range(1,7):
        classes={}
        for mask in range(1<<math.comb(n,2)):
            M=matrix(n,mask); K=key(M); H=hp(M)
            check(K not in classes or classes[K]==H,'H not constant')
            classes[K]=H
        check(len(classes)==c[n],'unweighted log mismatch')
        check(sum(classes.values())==d[n],'weighted log mismatch')
        rows.append({'n':n,'classes':c[n],'sum_H':d[n],'H_histogram':dict(sorted(__import__('collections').Counter(classes.values()).items()))})
    print(json.dumps({'status':'PASS','unweighted_to_12':c[1:],'weighted_to_12':d[1:],'exact_labelled_checks':rows},sort_keys=True,indent=2))
if __name__=='__main__': main()
