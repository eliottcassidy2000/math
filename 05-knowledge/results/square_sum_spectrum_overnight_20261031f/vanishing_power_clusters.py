#!/usr/bin/env python3
"""Independent finite audit of the bounded-span lemma for minimal signed powers."""
import itertools,json,math

def subset_sums(exponents):
    out={0}
    for e in exponents:out|={v+(1<<e) for v in out}
    return out

def main():
    maximum_terms=8;maximum_exponent=6;groups={}
    for size in range(1,maximum_terms):
        for es in itertools.combinations_with_replacement(range(maximum_exponent+1),size):
            total=sum(1<<e for e in es)
            groups.setdefault(total,[]).append((es,subset_sums(es)))
    equal=minimal=0;maxspan={}
    for total,rows in groups.items():
        for positive,psum in rows:
            for negative,nsum in rows:
                m=len(positive)+len(negative)
                if m>maximum_terms:continue
                equal+=1
                if psum&nsum!={0,total}:continue
                minimal+=1;span=max(positive+negative)-min(positive+negative)
                if span>m-2:raise RuntimeError('minimal vanishing span')
                maxspan[m]=max(maxspan.get(m,0),span)
    sharp=[]
    for m in range(2,17):
        positive=(0,) if m==2 else (0,0)+tuple(range(1,m-2));negative=(m-2,)
        total=sum(1<<e for e in positive)
        if total!=1<<negative[0] or subset_sums(positive)&subset_sums(negative)!={0,total}:
            raise RuntimeError('sharp carry-chain witness')
        sharp.append([m,max(positive+negative)-min(positive+negative)])
    c=9;bounds=[]
    for r in range(1,9):
        M=r+r*(r-1)//2+1+(c+1).bit_count();B=M+(c+1).bit_length()-1+4
        bounds.append({'max_apex_losses':r,'term_bound':M,'rhs_bound':B,
                       'exponent_upper_bound':math.factorial(r+1)*(1<<(r+1))*B})
    print(json.dumps({'status':'passed','universe':{'maximum_terms':maximum_terms,'exponents':[0,maximum_exponent]},
                      'equal_signed_multisets':equal,'minimal_vanishing_multisets':minimal,'maximum_span_by_terms':maxspan,
                      'sharp_carry_chains':sharp,'deficit9_coarse_effective_bounds':bounds,
                      'scope':'Finite audit of the span lemma; the effective fan bound uses the all-size affine classification and the explicit polyhedral argument.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
