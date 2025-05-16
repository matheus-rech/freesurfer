#!/usr/bin/env bash

source "$(dirname $0)/../test.sh"

test_command infant_recon_all --s infantfs_test_subj -i "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/infantfs_test_subj.mgz" --age 9 --outdir "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/recon"

compare_vol "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/recon/infantfs_test_subj/mri/norm.mgz" "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/infantfs_test_subj_norm.mgz"
compare_vol "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/recon/infantfs_test_subj/mri/aseg.mgz" "${FSTEST_TESTDATA_DIR}/infantfs_test_subj/infantfs_test_subj_aseg.mgz"
