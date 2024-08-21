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
    esac
done

for i in {1..5}
do

NN=$((100 * i + N))
# NN=$N

python3 ${home}/python/gengraph.py ${NN} ${P} ${S}

echo "=====${NN} ${P} ${S}"

# time0=$({ time ${home}/exec/frontier_thd >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
# echo "VTG ${time0}"

# time2=$({ time python3 ${home}/python/pyver.py >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
# echo "CPU ${time2}"

# time3=$({ time ${home}/exec/no_frontier >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
# echo "FTG ${time3}"

time3=$({ time ${home}/exec/justgpu >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
echo "GPU ${time3}"

echo " "

done
# python3 ${home}/python/sort.py ${N}

# DD1=$(diff -y --suppress-common-lines result/no_frontier.txt result/pyver.txt | wc -l)
# echo "diff no_frontier: ${DD1}"

# DD1=$(diff -y --suppress-common-lines result/pyver.txt result/frontier.txt | wc -l)
# echo "diff between pyver: ${DD1}"
