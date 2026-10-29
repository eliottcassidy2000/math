"""Independent grouped-core coloured kernel versus the parked response DP."""
import json
import baseline87 as b
import parked_counts as p

def verify():
    values,stats=b.palette_prefixes(30,(1,2,3))
    for q in (1,2,3):
        if p.prefix(30,q)[0]!=values[q]:
            raise ValueError(('independent coloured kernel',q))
    return dict(status='PASS',colors=[1,2,3],full_prefix_through=30,
                method='grouped-core kernel with analytic 3-tail; parked response DP',stats=stats)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
