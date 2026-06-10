#!/bin/bash
OUT_FILE="polar_urllc_results.txt"
rm -f $OUT_FILE
for ebno in $(seq 0 0.5 4); do
    esno=$(echo "$ebno + 10*l(512/1024)/l(10) + 10*l(2)/l(10)" | bc -l)
    echo "Eb/N0 = $ebno dB" >> $OUT_FILE
    ./bin/aff3ct -C "POLAR" -K 512 -N 1024 -m $esno -M $esno -e 200 --n-frames 200000 --dec-type "SCL" --list-size 8 --sim-stats 1 >> $OUT_FILE 2>&1
    echo "----------------------------------------" >> $OUT_FILE
done
