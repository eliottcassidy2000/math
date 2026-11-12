# Restoration provenance

The runtime restored an older filesystem before the final Session281 publication. Checkpoints 9b364259d2c4be4d5b39ad228aa012b6de7e06d9 and 6301e35a651b89d6a48ab516f043f932f2bd820e survived remotely. Final files were reconstructed from the workers' retained derivations and coordinator record. This release distinguishes newly replayed checks from earlier unavailable receipts; it does not fabricate old source hashes or imply the lost files were byte-for-byte recovered.

Chan–Pak completeness remains open. Reconstructed files are reviewed and replayed before release.
