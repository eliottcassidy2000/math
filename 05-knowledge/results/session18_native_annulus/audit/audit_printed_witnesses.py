#!/usr/bin/env python3
"""Parse the producer's printed numerical witnesses and authenticate receipts."""
import argparse,hashlib,json,re
from pathlib import Path

def need(ok,why):
    if not ok:raise RuntimeError(why)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--audit',type=Path,required=True);args=p.parse_args();root=args.package
    text=(root/'README.md').read_text();ann=json.loads((root/'annulus_receipt.json').read_text());ray=json.loads((root/'ray_receipt.json').read_text());audit=json.loads(args.audit.read_text())
    table=[list(map(int,m)) for m in re.findall(r'^\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|$',text,re.M)]
    expected=[[r['N'],r['max_increment'],r['first_peak_depth'],r['budget']] for r in ray['results']]
    need(table==expected,'printed six-ray peak table')
    small=list(map(int,re.search(r'old energy=(\d+), Q\(H\)=(\d+), mixed=(-?\d+),\narrived-new norm=(\d+), moving-new norm=(\d+),\nnew/new cross=(\d+), total cross=(\d+),\nDelta=(\d+) >153\*7=(\d+)',text).groups()))
    row=audit['N7_first_failure'];need(small==row[1:8]+[row[2]+row[3],1071],'printed N7 obstruction')
    old=list(map(int,re.search(r'old arrived sources are(\d+) and(\d+); its new arrived sources are(\d+) and(\d+)',text).groups()));need(old==[9,6371,203889,407779],'printed terminal-nine sources')
    nonmono=list(map(int,re.search(r'G_4096=(-?\d+) but G_4097=(-?\d+), an increase of(\d+)',text).groups()));need(nonmono==[-201674,-201616,58],'printed nonmonotonicity')
    ratio=list(map(int,re.search(r'maximum ratio of full-budget use at N=(\d+), with peak increment(\d+) at depth(\d+), against budget(\d+)',text).groups()));need(ratio==ann['worst_increment_ratio']+[153*ann['worst_increment_ratio'][0]],'printed finite sweep maximum')
    # Validate its consequence object in every archived row, without claiming
    # an independent re-execution of the producer's full contiguous sweep.
    records=ann['records'];need([r[0] for r in records]==list(range(4096,16385)),'contiguous saved record coverage')
    need(all(r[3]==153*r[0]-r[1] and r[4]<=0 and r[3]>=0 for r in records),'saved finite consequences')
    worst=records[0]
    for r in records[1:]:
        if r[1]*worst[0]>worst[1]*r[0]:worst=r
    need(worst[:3]==ann['worst_increment_ratio'],'exact-rational saved maximum')
    parsed_logs=[]
    for log,receipt in [('normal.log','annulus_receipt.json'),('optimized.log','annulus_receipt.json'),('rays-normal.log','ray_receipt.json'),('rays-optimized.log','ray_receipt.json'),('recurrence-normal.log','recurrence_receipt.json'),('recurrence-optimized.log','recurrence_receipt.json')]:
        j=json.loads((root/log).read_text());r=json.loads((root/receipt).read_text());need(j['status']=='PASS' and j['receipt_sha256']==digest(root/receipt) and j['checks']==r['checks'],'printed log '+log);parsed_logs.append(log)
    print(json.dumps({'status':'PASS','printed_ray_rows':table,'printed_N7_witness':small,'printed_terminal_nine_sources':old,'printed_nonmonotonicity':nonmono,'printed_sweep_maximum':ratio,'saved_records_checked':len(records),'logs_parsed':parsed_logs},sort_keys=True))
if __name__=='__main__':main()
