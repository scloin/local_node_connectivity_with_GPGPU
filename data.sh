#!/bin/bash

home=/home/deepl/sooho/local_node_connectivity_with_GPGPU/code

#default value
N=20
P=0.4
S=20

while getopts n:s:p: opt # if first argument str length is not 0
do
  #echo $1
  case ${opt} in
    n) N="$OPTARG";;
    s) S="$OPTARG";;
    p) P="$OPTARG";;
    *) echo "Option.."
    esac
done

for i in {0..4}
do

NN=$((20000 * i + N))

#echo "=====${NN} ${P} ${S}\n"

python3 ${home}/python/gengraph.py ${NN} ${P} ${S}

echo "=====${NN} ${P} ${S}\n"
${home}/exec/frontier_thd_test
time0=$({ time ${home}/exec/frontier_thd_test; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
echo "fr_t  : ${time0}"

time2=$({ time python3 ${home}/python/pyver_t.py; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
echo "py    : ${time2}"

time3=$({ time ${home}/exec/no_frontier_test; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
echo "nfr_t : ${time3}"


time3=$({ time ${home}/exec/justgpu_test; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
echo "gpu_t : ${time3}"

# time3=$({ time ${home}/python/pyver_t.bin; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
# echo "py_c  : ${time3}"

# time4=$({ time python3 ${home}/python/pyver_t.py >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
# echo "${time4}"

# DD1=$(diff -y --suppress-common-lines result/no_frontier.txt result/pyver.txt | wc -l)
# echo "diff no_frontier: ${DD1}"

# DD1=$(diff -y --suppress-common-lines result/frontier.txt result/pyver.txt | wc -l)
# echo "diff between pyver: ${DD1}"

done
echo ""
