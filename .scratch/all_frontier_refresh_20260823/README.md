# All-frontier refresh artifacts

Scratch-only packet for the 2026-08-23 bounded frontier refresh.  `REPORT.md`
contains the synthesis and truth/status firewall.  The four scripts have
frozen matching outputs and make no canon claim.

## Files

- `rooted_switching_probe.py/.out` — exhaustive strong-tournament switching
  hostile through order six, with independent Pfaffian reconstruction.
- `lrc_two_cube_address_probe.py/.out` — complete THM-3743 support-two address
  census and nonprimitive two-cube control through pair sum 1000.
- `planar_surface_sidecar_probe.py/.out` — exact Hensel-resonance,
  r-independent `z^2`, and cubic branch-fibre controls.
- `mixed_anchored_step3_probe.py/.out` — exact THM-3795-to-THM-3798
  map/loss audit with the mandatory nodal anchors restored, for both the
  THM-3803 affine control and the next adjacent quadratic correction.
- `REPORT.md` — broad inheritance matrix and three typed connection audits.

## Reproduction

```powershell
python3 -B .scratch/all_frontier_refresh_20260823/rooted_switching_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/rooted_switching_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/lrc_two_cube_address_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/lrc_two_cube_address_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/planar_surface_sidecar_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/planar_surface_sidecar_probe.py
python3 -B .scratch/all_frontier_refresh_20260823/mixed_anchored_step3_probe.py
python3 -B -O .scratch/all_frontier_refresh_20260823/mixed_anchored_step3_probe.py
```

At packet creation, the eight script/output SHA-256 values were:

```text
3ba234f1a4da9f0691238720e90e4b15ad404bb9c940c5e99822772421f0338c  rooted_switching_probe.py
99804febaab6f6c691f7a29f13412146c4a3c42af2f54dbe2bf1e052b4d37501  rooted_switching_probe.out
58ea4bf57caef4ce2f7fc26cf6cadeac1793cbe431d02bd7c9c5d0502531c524  lrc_two_cube_address_probe.py
8b2a69da5d5cffa451239919fafa34f38c1d5c3708dbb1e622fca802b9bbf10d  lrc_two_cube_address_probe.out
8fe865771f6c564c92f5827fbe747ca117b538dee38a95e18e2ca827ecdc59ac  planar_surface_sidecar_probe.py
d0e08977bfca288d7bcb45e4f0482359de363072e44661392955c6e053bb3e93  planar_surface_sidecar_probe.out
c9ddfa2c1d627102e8dabc7d95e6cb642711577331603b4b1a530b56a5905d67  mixed_anchored_step3_probe.py
342b9d2e71b016c4c0ef8efa61bcb2669047c1b239e18f1dc624cd372cd57c63  mixed_anchored_step3_probe.out
```
