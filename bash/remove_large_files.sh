#!/bin/bash
git ls-tree -r --name-only --full-tree HEAD | \
grep -E "data/infernal/g1_intron_IMGVR_phagescope.aout|data/infernal/g1_intron_IMGVR.tblout|data/infernal/g1_intron_IMGVR_phagescope.out|data/infernal/g1_intron_IMGVR.aout|data/infernal/g1_intron_IMGVR.out" | \
xargs git rm --cached

