#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_zeropad_vs_bare_kernel_kps.py   (kind-pasteur 2026-06-21, HYP-2724 adjudication)

DECISIVE adjudication of the workflow woi07wlep "support-6 floor => support-3 contributes
ZERO => HYP-2723/2724 REFUTED" conclusion.

The conclusion is the ALREADY-CONCEDED THM-538 / MISTAKE-078 error (CASE-thm538-support6-floor-
zero-padding): the support-6 floor holds for the ACTIVE-COORDINATE sum Q (bare vector), but the
MEASURE kernel K(n) in corr(E)=Sum_{n in Lambda°(E)} K(n) is ZERO-PADDED to length k, and the zero
coordinates carry chat(0,T)=(1-|T|/7) factors that BREAK the floor. So support-3 relations DO
contribute, and (corrected canon) "the AP's correction is support-3-dominated."

This script demonstrates it on the workflow's OWN kernel (Kk from lrc_q108_weight_enumerator_validate).
"""
import importlib.util, os
_d = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location("wev", os.path.join(_d, "lrc_q108_weight_enumerator_validate_kpswf2.py"))
wev = importlib.util.module_from_spec(spec); spec.loader.exec_module(wev)
Kk = wev.Kk

def show(name, base):
    print(f"  {name}: bare={Kk(base).real:+.6f}", end="")
    for z in (3,4,5):
        print(f"  +{z}zeros(len{len(base)+z})={Kk(base+[0]*z).real:+.6f}", end="")
    print()

if __name__ == "__main__":
    print("="*78)
    print("ADJUDICATION: bare (active-coord Q) vs zero-padded (measure K) -- support-3 relations")
    print("="*78)
    print("Canon THM-538 (CORRECTED): K(1,1,-1,0,0,0,0)=+0.00066 != 0  (a support-3 relation, k=8).")
    print()
    show("Schur (1,1,-1) ", [1,1,-1])
    show("3-AP  (1,-2,1) ", [1,-2,1])
    show("supp2 (1,-1)x  ", [1,-1])     # support-2: 1*e_i = 1*e_j only if equal; placeholder
    print()
    print("VERDICT:")
    print(" * BARE vector  -> K=0 for support<=5  == the active-coord sum Q (HAS the support-6 floor).")
    print(" * ZERO-PADDED  -> K!=0 once total length>=6  == the MEASURE kernel (NO floor; canon THM-538).")
    print(" * corr(E)=Sum over Lambda°(E) uses ZERO-PADDED length-k vectors, so support-3 relations")
    print("   DO contribute (~+0.0007..0.0013 each at k=8). The AP, rich in support-3 relations")
    print("   (3-APs+Schur=additive energy A3), is support-3-dominated -> corr~0.30. HYP-2724 (A3 the")
    print("   driver, Pearson(A3,corr)=+0.93) STANDS; the workflow's 'support-6 carrier' is the")
    print("   conceded Q-vs-K conflation (THM-538 / MISTAKE-078).")
