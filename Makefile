.PHONY: check checksums quick full archive

check:
	bash scripts/check_release.sh

checksums:
	cd proof && sha256sum -c checksums.sha256
	sha256sum -c RELEASE_SHA256SUMS

quick:
	cd proof && bash reproduction/scripts/verify_existing.sh

full:
	cd proof && JOBS=$${JOBS:-4} bash reproduction/scripts/rerun_all.sh

archive:
	python3 scripts/build_archive.py
