# Checkpoint2 verification additions

Use the commands in [VERIFICATION.md](VERIFICATION.md); the newest versioned
manifest is selected automatically. All foundational files and manifests remain
unchanged and are included in the new inventory.

There are22 declared positive replay jobs, each run normally and with `-O`, plus
five negative-control families in both modes. The complete bounded two-step
search is replayed by `replay_depth_search.py`: it copies the unchanged producer
script into a fresh temporary directory, checks its generated JSON against both
stdout and the frozen receipt, and removes only that temporary directory. This
prevents a producer output file from modifying the frozen packet. The third-step
search and independent saved-receipt verifiers are also replayed. These are exact
finite domains, not exhaustive native reachability classifications.

The producer explicitly excluded three convenience depth-search duplicates.
Their frozen auditor copies remain audit evidence; no unlisted native copies
are added. All source snapshots remain inert. Existing core manifests, the new
native manifests/hash lists, reviewed-copy equality, Git origins and every
outer-packet byte are actively checked. Normal/optimized equivalence is required,
not inferred from an assertion-free source inspection.
