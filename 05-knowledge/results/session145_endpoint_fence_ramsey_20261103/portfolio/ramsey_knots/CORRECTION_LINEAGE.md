# Tangle scope correction

The first sealed edition used the broad phrase “actual tangles” in Section C after establishing its matching formula for disjoint paths. That phrase did not explicitly exclude closed components already present in the tangle.

This edition restricts that statement to exactly r arcs with no pre-existing closed components. For a tangle with m additional closed components, the correct total is m+c(P,Q), so m must also be retained. A disjoint extra circle leaves the endpoint matching unchanged and increases the result by one. The two braid examples have no pre-existing closed components and remain unchanged. The independent checker now constructs matching-closure graphs with zero through three additional disjoint circles and verifies the corrected formula and the extra-circle hostile.

The original sealed edition is preserved without alteration. Its manifest SHA256 is `817bcc5b94af4857857e8b355a32b3457ddaec8ed482135995ea68f9f0bc68ce`.

Inherited Markdown snapshots now use .md.txt filenames. Their contents, source revisions, Git blob identities and SHA256 hashes are unchanged; only their archive paths and the corresponding source-pin entries differ. The new manifest records the revised package.
