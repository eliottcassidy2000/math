#!/usr/bin/env python3
"""Median wall-clock comparison; no performance assertion. No omitted failures."""
import argparse,json,platform,statistics,time
from baseline87 import degree_bucket_dp
from baseline98 import fixed_degree
from integer_response import count,count_packed

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeat',type=int,default=3)
    parser.add_argument('--orders',nargs='+',type=int,default=[40,60,100,150,200])
    args=parser.parse_args();rows=[]
    for n in args.orders:
        answer=None;row={'n':n,'methods':{}}
        methods={
          'session87_all_coefficients':lambda:degree_bucket_dp(n)[0][n],
          'session98_fraction_projective':lambda:fixed_degree(n)[0],
          'integer_tuple_uncached':lambda:count(n,1,False)[0],
          'integer_tuple_cached':lambda:count(n,1,True)[0],
          'integer_tuple_three_tail':lambda:count(n,3,True)[0],
          'packed_full_factorial':lambda:count_packed(n,False)[0],
          'packed_odd_local':lambda:count_packed(n,True)[0]}
        for name,fn in methods.items():
            times=[]
            for _ in range(args.repeat):
                start=time.perf_counter();value=fn();times.append(time.perf_counter()-start)
                if answer is not None and value!=answer:raise RuntimeError('benchmark disagreement')
                answer=value
            row['methods'][name]={'seconds':times,'median':statistics.median(times)}
        rows.append(row)
    print(json.dumps({'python':platform.python_version(),'platform':platform.platform(),'repeat':args.repeat,'rows':rows},indent=2))
if __name__=='__main__':main()
