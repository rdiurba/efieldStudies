#!/bin/bash
INDIR=/pnfs/dune/persistent/users/rdiurba/laserFiles/golden/
FCL=laseranalysispdhd.fcl

for run in $(ls $INDIR | grep "\.root$" | grep -v "_hist_" | sed -E 's/^np04hd_raw_run([0-9]+)_.*/\1/' | sort -u) ; do
  src=()
  for f in $(ls $INDIR/np04hd_raw_run${run}_*.root | grep -v "_hist_") ; do src+=(-s "$f") ; done
  echo "run $run"
  lar -c $FCL "${src[@]}" -T laseranalysis_run${run}.root > /dev/null 2>&1
  root -l -b -q "add_laser_truth.C(\"laseranalysis_run${run}.root\")"
done
