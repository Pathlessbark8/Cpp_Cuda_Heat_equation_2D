#!/bin/sh

SECONDS=0
loop_cnt=0
rm -rf ./a.out
# rm -rf ./cpp.perf
rm -rf ./cuda.perf

ntimes=$1
while [ $loop_cnt -lt $ntimes ]
do
    echo "Running Perf test... ($loop_cnt/$ntimes)"
    # ./cpp.out >> cpp.perf 
    # P1=$!
    ./cuda.out >> cuda.perf 
    # P2=$!
    # wait $P1 $P2
    loop_cnt=`expr $loop_cnt + 1`
done
duration=$SECONDS
echo "Script took $(($duration / 60)) minutes and $(($duration % 60)) seconds."