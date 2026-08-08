#!/bin/bash
# Threshold-25 tower for F(15).  Establishes that no order-15 signing has
# M <= 25; with the parity law and the witness of maximum 27 this gives F(15) = 27.
#
# Thresholds by parity: order 12 -> M <= 24, order 13 -> M <= 24, order 14 -> M <= 25,
# order 15 -> M <= 25.  Level 12 comes from a complete geng pass over the
# 1,018,997,864 graphs on 11 vertices; every later level is a one-vertex expansion.
#
# Level 12->13 is split into $PASSES passes over the same seed.  Pass p keeps only the
# extensions whose canonical form hashes to p, so the raw file that has to be sorted is
# 1/$PASSES of the whole.  That is the only concession to disk: the union of the passes
# is exactly the deduped level.
#
# Needs nauty (geng, labelg) on PATH, numpy, and roughly 15 GB of scratch with
# PASSES=8.  Expect 15 to 20 hours on 16 cores.  Run it from this directory:
#   PASSES=8 JOBS=16 bash tower25.sh
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
PASSES=${PASSES:-8}
JOBS=${JOBS:-16}

if [ ! -s T12_thresh24.g6 ]; then
  mkdir -p T12run24
  seq 0 63 | xargs -P $JOBS -I{} sh -c 'KEEP_PREFIX=T12run24/h python3 sample_T.py 11 24 64 $1 > T12run24/log_$1.txt 2>&1' _ {}
  cat T12run24/h_*.g6 > T12_thresh24.g6
fi
echo "|T12(M<=24)| = $(wc -l < T12_thresh24.g6)"

mkdir -p t25a
: > T13_thresh24.g6
for p in $(seq 0 $((PASSES-1))); do
  find t25a -type f -delete
  seq 0 $((JOBS-1)) | xargs -P $JOBS -I{} sh -c \
    'python3 expand_prune.py T12_thresh24.g6 11 24 /dev/stdout --res $1 --mod '"$JOBS"' --K 128 2> t25a/log$1.txt | labelg -q | python3 hashfilter.py '"$PASSES $p"' | sort -u -T /tmp > t25a/r$1.g6' _ {}
  sort -m -u t25a/r*.g6 >> T13_thresh24.g6
  echo "pass $p done, running |T13| = $(wc -l < T13_thresh24.g6)"
done
sort -u -o T13_thresh24.g6 T13_thresh24.g6
echo "|T13(M<=24)| = $(wc -l < T13_thresh24.g6)"

mkdir -p t25b; find t25b -type f -delete
seq 0 $((JOBS-1)) | xargs -P $JOBS -I{} sh -c \
  'python3 expand_prune.py T13_thresh24.g6 12 25 t25b/r$1.g6 --res $1 --mod '"$JOBS"' --K 128 > t25b/log$1.txt 2>&1' _ {}
cat t25b/r*.g6 | labelg -q | sort -u > T14_thresh25.g6
echo "|T14(M<=25)| = $(wc -l < T14_thresh25.g6)"

mkdir -p t25c; find t25c -type f -delete
seq 0 $((JOBS-1)) | xargs -P $JOBS -I{} sh -c \
  'python3 expand_prune.py T14_thresh25.g6 13 25 t25c/r$1.g6 --res $1 --mod '"$JOBS"' --K 128 > t25c/log$1.txt 2>&1' _ {}
cat t25c/r*.g6 > T15_raw25.g6
if [ -s T15_raw25.g6 ]; then
  labelg -q T15_raw25.g6 | sort -u > T15_thresh25.g6
  echo "|T15(M<=25)| = $(wc -l < T15_thresh25.g6)   =>  F(15) = 25"
else
  echo "T15(M<=25) EMPTY  =>  F(15) = 27 (witness already verified)"
fi
echo TOWER25_DONE
