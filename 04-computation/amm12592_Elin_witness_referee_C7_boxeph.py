"""LANE D2 (E-lin) -- C7: independent-referee verification of a
linear-slack rule-A closure.

Blocks are reconstructed from the certified fast engine (full_blocks) at
(R, D0) = (512, 32) [eps = 1/16] and (256, 16); they are then verified by
the INDEPENDENT hostile-referee implementation
(amm12592_independent_witness_referee_boxeph.py): its own admissible(),
its own block_poly/padd/peq polynomial arithmetic, its own exact
Lucas/Fibonacci floor engine for the profile, its own qpow target.
Nothing from the engine's toolbox is used on the verification side.
"""
import os, json, importlib.util, contextlib, io

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_Elin_witnessC7_boxeph.json")

spec = importlib.util.spec_from_file_location(
    "fast", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fast = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fast)

spec2 = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec2)
with contextlib.redirect_stdout(io.StringIO()):
    spec2.loader.exec_module(iref)

report = {}
for R, D0 in ((256, 16), (512, 32)):
    blocks, res = fast.full_blocks(R, D0)
    assert blocks is not None and res["outcome"]["status"] == "CLOSED"
    profile = [iref.floor_gamma_star(R + i) + D0 for i in range(R)]
    bad = []
    for i, (d, row) in enumerate(zip(profile, blocks)):
        ok, why = iref.admissible(row, d)
        if not ok:
            bad.append((i, why))
    S = []
    for i, (d, row) in enumerate(zip(profile, blocks)):
        S = iref.padd(S, iref.pshift(iref.block_poly(row, d), i))
    ident = iref.peq(S, iref.qpow(R - 1))
    rec = {"R": R, "D0": D0, "n_blocks": len(blocks),
           "admissible_all": len(bad) == 0, "adm_bad": bad[:5],
           "identity": bool(ident),
           "unit_ok": all(row[0] in (1, -1) for row in blocks)}
    report[f"R{R}_D0{D0}"] = rec
    print(f"C7 R={R} D0={D0}: admissible={rec['admissible_all']} "
          f"identity={rec['identity']} unit={rec['unit_ok']}", flush=True)
    assert rec["admissible_all"] and rec["identity"]

report["ALL_PASS"] = True
json.dump(report, open(OUT, "w"), indent=1)
print("C7 ALL_PASS ->", OUT, flush=True)
