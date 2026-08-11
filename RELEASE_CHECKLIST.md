# Release checklist

## Scientific payload

- [x] Proof report separates rigorous, numerical, and catalog-conditional claims.
- [x] Machine-readable initial box and period interval are included.
- [x] Full source, mesh, center enclosures, outer C1 tubes, and certificate are included.
- [x] Original proof payload has a SHA-256 manifest.
- [x] Archived-log certificate check reproduces `PROOF_OK 1`.
- [x] Full regeneration command is documented.
- [ ] Obtain at least one external Level 2 rerun.
- [ ] Preferably obtain one Level 3 independent implementation.

## Release engineering

- [x] Ubuntu 24.04 CI workflow added.
- [x] Container recipe added.
- [x] Release-level structural and checksum script added.
- [x] Independent-verification issue template added.
- [x] Deterministic archive builder added.
- [ ] Replace the release-candidate tag after external verification.

## Publication metadata

- [x] Record OpenAI Codex / GPT-5.6 Sol (`gpt-5.6-sol`) as an AI computational contributor in `AUTHORSHIP.md` and `DISCLOSURE.md`.
- [x] Record `Rexyysilent` as the approved public creator identity.
- [x] Apply Apache-2.0 to code and CC-BY-4.0 to documentation/data.
- [x] Add the public repository URL.
- [x] Complete `.zenodo.json` for the future archive.
- [x] Complete `CITATION.cff` for GitHub and Zenodo.
- [ ] Add ORCID, affiliation, and contact details if the creator chooses to publish them.
- [x] Prepare the `v1.0.0-rc1` GitHub pre-release; publishing the containing tag completes distribution.
- [ ] Create a signed/tagged `v1.0.0` GitHub release.
- [ ] Archive that exact tag with Zenodo and record the version DOI.
- [ ] Add the DOI to the manuscript, README, citation file, and release notes.

## Manuscript

- [ ] Convert `proof/proof_report.md` into reviewed LaTeX/PDF.
- [x] Add a release-level AI/tool-use contribution record and disclosure draft.
- [ ] Adapt the disclosure and authorship statement to the target venue's current policy.
- [ ] Have a domain expert check the reversibility argument and novelty wording.
- [ ] Submit the preprint and cite the immutable version DOI and commit/tag.
