#!/usr/bin/env python3
"""Session31 exact/numerical controls for the tree-series mirror.
Python 3.11+ standard library only. Numerical checks do not prove transcendence.
Run: python3 verify_mirror_tree.py --output mirror_verification.json
"""
from __future__ import annotations
import argparse, itertools, json, math
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path


def coefficients(mmax):
    a=[F(1)]
    for m in range(mmax):
        a.append(F(m+2, 2*(m+1))*sum(a[i]*a[m-i] for i in range(m+1)))
    for m in range(1,mmax+1):
        assert a[m]==F((m+1)**(m-1),math.factorial(m))
    return a


def forest_count(m):
    # Each vertex chooses a parent in {0,..,m} other than itself;
    # 0 is an added distinguished root. Reject every directed cycle.
    count=0
    for parents in itertools.product(*[[v for v in range(m+1) if v!=i] for i in range(1,m+1)]):
        good=True
        for start in range(1,m+1):
            seen=set(); v=start
            while v:
                if v in seen:
                    good=False; break
                seen.add(v); v=parents[v-1]
            if not good: break
        count += good
    return count


def mirror(n):
    if n<2: raise ValueError('No distinct positive finite mirror for integer n<2.')
    if n==2: return D(4)
    target=D(n).ln()/n
    lo,hi=D(1),D(1).exp()
    for _ in range(230):
        mid=(lo+hi)/2
        if mid.ln()/mid < target: lo=mid
        else: hi=mid
    return (lo+hi)/2


def series(z,M):
    # a_(m+1)/a_m=(1+1/(m+1))^m, independent of the convolution check.
    term=D(1); total=term
    for m in range(M):
        term *= z*(D(m+2)/D(m+1))**m
        total += term
    nextterm=term*z*(D(M+2)/D(M+1))**M
    bound=nextterm/(1-D(1).exp()*z)
    return total,bound


def hp(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            value=dp[mask][v]
            if not value: continue
            avail=out[v]&~mask
            while avail:
                bit=avail&-avail; avail-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])


def out_from_mask(n,mask):
    out=[0]*n
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v=(i,j) if (mask>>k)&1 else (j,i)
        out[u] |= 1<<v
    return out


def strong(out):
    n=len(out)
    for s in range(n):
        seen=1<<s; frontier=seen
        while frontier:
            bit=frontier&-frontier; frontier-=bit
            nxt=out[bit.bit_length()-1]&~seen
            seen|=nxt; frontier|=nxt
        if seen!=(1<<n)-1: return False
    return True


def directed_forest_count(out):
    n=len(out); count=0
    choices=[[-1]+[j for j in range(n) if out[i]>>j&1] for i in range(n)]
    for parents in itertools.product(*choices):
        good=True
        for start in range(n):
            seen=set(); v=start
            while v!=-1:
                if v in seen: good=False; break
                seen.add(v); v=parents[v]
            if not good: break
        count+=good
    return count


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='mirror_verification.json')
    args=parser.parse_args(); report={}
    a=coefficients(40)
    assert a[1:7]==[F(1),F(3,2),F(8,3),F(125,24),F(54,5),F(16807,720)]
    report['coefficient_checks']={'independent_convolution_through':40,'prefix':[str(x) for x in a[:9]]}
    report['rooted_forests']={}
    for m in range(1,7):
        count=forest_count(m); assert count==(m+1)**(m-1)
        report['rooted_forests'][m]={'exhaustive_parent_functions':m**m,'accepted':count}
    report['mirror_checks']=[];report['branching_checks']=[]
    with localcontext() as ctx:
        ctx.prec=90
        for n,M in [(2,800),(3,10000),(4,800),(5,400),(7,200),(21,80),(63,50),(100,40),(1000,25)]:
            nD=D(n); lam=nD.ln(); z=lam/nD; r=mirror(n)
            assert abs(r.ln()/r-z)<D('1e-65')
            principal,bound=series(z,M)
            expected=D(2) if n==2 else r
            error=expected-principal
            assert error>=-D('1e-65') and error<=bound+D('1e-65')
            q=D(1) if n==2 else r/nD
            mu=lam*q
            assert abs(q-(lam*(q-1)).exp())<D('1e-65')
            assert 0<mu<1
            if n>=3:
                assert abs(mu-r.ln())<D('1e-65')
                assert abs(mu*(-mu).exp()-z)<D('1e-65')
                assert abs(r/nD-q)<D('1e-65')
            # Borel point masses and offspring tilt, k=1..80 and j=0..25.
            max_borel_err=D(0);max_offspring_err=D(0)
            for k in range(1,81):
                p=(-lam*k).exp()*(lam*k)**(k-1)/D(math.factorial(k))
                dual=(-mu*k).exp()*(mu*k)**(k-1)/D(math.factorial(k))
                max_borel_err=max(max_borel_err,abs(p/q-dual))
            for j in range(26):
                p=(-lam).exp()*lam**j/D(math.factorial(j))
                tilted=p*q**j/q
                dual=(-mu).exp()*mu**j/D(math.factorial(j))
                max_offspring_err=max(max_offspring_err,abs(tilted-dual))
            assert max_borel_err<D('1e-65') and max_offspring_err<D('1e-65')
            assert abs(principal/nD-q)<=bound/nD+D('1e-65')
            report['mirror_checks'].append({'n':n,'mirror':str(r),'terms_M':M,'series_value':str(principal),'series_error':str(error),'proved_real_tail_bound_numerically_evaluated':str(bound)})
            report['branching_checks'].append({'n':n,'lambda':str(lam),'extinction_q':str(q),'conditional_mean_offspring_mu':str(mu),'conditional_mean_total_progeny':str(1/(1-mu)),'max_borel_identity_residual':str(max_borel_err),'max_offspring_identity_residual':str(max_offspring_err)})
        assert mirror(2)/2>D(1) # wrong-branch "extinction probability" hostile.
        report['hostile_n2']={'mirror':4,'principal_series':2,'actual_extinction_probability':1,'wrong_mirror_over_n':2}
        # Pointwise asymptotic and absolutely-summable correction diagnostics.
        report['asymptotic_checks']=[]
        for n in [100,1000,10000,1000000]:
            z=D(n).ln()/n;r=mirror(n)
            ratio=(r-z.exp())/(z*z)
            report['asymptotic_checks'].append({'n':n,'(r_n-exp(z))/z^2':str(ratio)})
            assert ratio>1
        # Rational-root finite positive control, not a transcendence proof.
    pairs=[(n,m) for n in range(2,31) for m in range(2,31) if n!=m and n**m==m**n]
    assert pairs==[(2,4),(4,2)]
    report['integer_mirror_pairs_through_30']=pairs
    try: mirror(1)
    except ValueError: report['n1_distinct_mirror_rejected']=True
    else: raise AssertionError('n=1 must be rejected')
    cycle=[2,4,1]
    assert hp(cycle)==3 and directed_forest_count(cycle)==7
    report['forest_not_path_hostile']={'tournament':'0->1->2->0','Hamiltonian_paths':3,'rooted_directed_forests':7}
    report['labelled_tournament_spectra_through_6']={}
    for n in range(1,7):
        spectrum=set(); strongs=set(); total=1<<(n*(n-1)//2)
        for mask in range(total):
            out=out_from_mask(n,mask); h=hp(out);spectrum.add(h)
            if strong(out): strongs.add(h)
        assert 7 not in spectrum and 21 not in spectrum
        report['labelled_tournament_spectra_through_6'][n]={'tournaments':total,'H_spectrum':sorted(spectrum),'strong_H_spectrum':sorted(strongs)}
    # Exhibit the failure of any "all multiples of 7 are forbidden" inference.
    witness=None
    for word in range(1<<7):
        out=[0]*8
        for i in range(1,8):
            for j in range(i+1,8): out[i] |= 1<<j
            if word>>(i-1)&1: out[0] |= 1<<i
            else: out[i] |= 1
        if hp(out)==63:
            witness={'apex_word_low_bit_is_vertex_1':word,'outmasks':out,'Hamiltonian_paths':63};break
    assert witness is not None
    report['attained_multiple_7_hostile']=witness
    report['status']='PASS: exact finite checks plus high-precision numerical diagnostics; proofs are in README.md.'
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(report['status']);print('Wrote',args.output)
    print('Cayley forest counts:',[x['accepted'] for x in report['rooted_forests'].values()])
    print('Integer mirror pairs:',pairs)
    print('H=63 witness:',witness)
    print('Strong spectra:',{n:r['strong_H_spectrum'] for n,r in report['labelled_tournament_spectra_through_6'].items()})
if __name__=='__main__': main()
