make sequential

list_b=(150 300 450 600 750 1000)
list_sizes=(200 400 600 800 1000 1200)

echo "B;T;Ag;Size;Time" > results/seq_tests.csv

for b in ${list_b[@]}; do
    for size in ${list_sizes[@]}; do
        ag=$((10 * $b))
        echo -n "$b;$b;$ag;$size;"
        ./exe/sequential $b $b $ag $size | rev | cut -d' ' -f1 | rev
    done
done >> results/seq_tests.csv