"""Exhaustive short-word tests of the exact non-3 macro inverse normal form."""
from itertools import product
from pathlib import Path
import json
import digital
BASE=Path(__file__).resolve().parent
roots=json.loads((BASE/'actual_seed_roots.json').read_text())
states=[(r['state'][0],r['state'][1],r['state'][2],r['state'][0]-r['state'][3]) for r in roots]
words=[w for k in range(3) for w in product((0,1),repeat=k)]
count=0
for seed in states:
    for pre in words:
        for post in words:
            for r in (0,1,2):
                for middle in product(('C1','C2'),repeat=r):
                    blocks=pre+middle+post
                    while blocks and blocks[-1]==1:blocks=blocks[:-1]
                    final=digital.run_blocks(seed,blocks)
                    h,F,q,g=final
                    if F%3==0:continue
                    if blocks and q>=4*F:raise RuntimeError(('final bound',seed,blocks,final))
                    if q%3**r:raise RuntimeError(('shadow q not integral',seed,blocks,final))
                    decodings=digital.macro_decodings(seed,q//3**r,r)
                    if not any(digital.run_blocks(seed,d)==final for d in decodings):
                        raise RuntimeError(('inverse misses actual state',seed,blocks,final,decodings))
                    count+=1
print('all_273_actual_seeds_pre_post_lengths_0_to_2_and_0_1_2_triangle_macros',count,'PASS')
# Soundness and a finite completeness sanity check, not a spectrum census.
for state in states[:20]:
    for target in range(1,32,2):
        if target%3==0:continue
        answer=digital.inverse_non3(target,state)
        if answer is not None and digital.run_blocks(state,answer)[1]!=target:
            raise RuntimeError(('inverse unsound',state,target,answer))
print('small_target_exact_inverse_soundness_PASS')
print('PASS_FINAL_BOUND_NOT_ARBITRARY_BLOCK_PRUNING')
