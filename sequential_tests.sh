#!/bin/bash

#make parallel
make sequential

list_b=(50 80 100 150 180)
list_sizes=(50 100 150 200 250 300)

#echo "B;T;Ag;Size;Time" > results/parallel_tests.csv
echo "B;T;Ag;Size;Time" > results/seq_tests.csv

for b in ${list_b[@]}; do
    for size in ${list_sizes[@]}; do
	    ag=$((10 * $b))
        echo -n "$b;$b;$ag;$size;"
        #./exe/parallel $b $b $ag $size | rev | cut -d' ' -f1 | rev
	    ./exe/sequential $b $b $ag $size | rev | cut -d' ' -f1 | rev
    done
done >> results/seq_tests.csv #results/parallel_tests.csv
