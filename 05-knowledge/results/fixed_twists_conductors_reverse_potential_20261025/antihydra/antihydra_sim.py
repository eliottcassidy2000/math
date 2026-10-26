# Antihydra: a0=8, b0=0; even a: (3a/2, b+2); odd a: halt if b==0 else ((3a-1)/2, b-1).  Track b, odd-to-odd valuations.
import sys, time, json
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
a, b = 8, 0
minb, argmin = 0, 0
vals = {}          # histogram of v = v_2(3a-1) at odd a (odd-to-odd valuation)
odd_count = even_count = 0
last_odd_step = None
walk_check = True  # check b == 5 + sum(2v_i - 3) at odd terms after the first odd term
B = None
t0 = time.time()
step = 0
first_odd_b = None
while step < N:
    if a & 1:
        if b == 0:
            print(json.dumps({'halted_at_step': step, 'a_bits': a.bit_length()})); sys.exit()
        v = ((3*a - 1) & -(3*a - 1)).bit_length() - 1   # v_2(3a-1)
        vals[v] = vals.get(v, 0) + 1
        if B is None:
            B = b - 1; first_odd_b = b - 1
        else:
            B = B + 2*v_prev - 3
            if B != b - 1: walk_check = False
        v_prev = v
        a = (3*a - 1) >> 1; b -= 1; odd_count += 1
    else:
        a = (3*a) >> 1; b += 2; even_count += 1
    step += 1
    if b < minb: minb, argmin = b, step
print(json.dumps({'steps': N, 'odd': odd_count, 'even': even_count, 'final_b': b, 'min_b': minb, 'argmin_b': argmin,
                  'first_odd_b': first_odd_b, 'walk_identity_holds': walk_check, 'a_bits': a.bit_length(),
                  'valuation_histogram': dict(sorted(vals.items())), 'seconds': round(time.time()-t0,1)}))
