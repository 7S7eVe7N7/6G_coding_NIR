#!/bin/bash
OUT_FILE="polar_mmtc_results.txt"
rm -f $OUT_FILE
for ebno in 0 0.5 1 1.5 2 2.5 3 3.5 4; do
    esno=$(echo "$ebno + 10*l(256/512)/l(10) + 10*l(2)/l(10)" | bc -l)
    echo "Eb/N0 = $ebno dB" >> $OUT_FILE
    ./bin/aff3ct -C "POLAR" -K 256 -N 512 -m $esno -M $esno -e 200 --n-frames 200000 --dec-type "SCL" --list-size 8 --sim-stats 1 >> $OUT_FILE 2>&1
    echo "----------------------------------------" >> $OUT_FILE
done
