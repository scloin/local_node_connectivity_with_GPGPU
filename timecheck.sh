home=/home/deepl/sooho/test/K-com/test/code

#default value
N=20
P=0.4
S=20

while getopts n:s:p:c: opt # if first argument str length is not 0
do
  #echo $1
  case ${opt} in
    n) N="$OPTARG";;
    s) S="$OPTARG";;
    p) P="$OPTARG";;
    esac
done

python ${home}/python/gengraph.py ${N} ${P} ${S}

time0=$({ time ${home}/exec/frontier_thd >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time1=$({ time ${home}/exec/frontier >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time2=$({ time python ${home}/python/pyver.py >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time3=$({ time ${home}/exec/no_frontier >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')

echo "no_frontier   (${time3})"
echo "frontier_thd  (${time0})"
echo "frontier      (${time1})"
echo "python        (${time2})"

python ${home}/python/sort.py ${N}

DD1=$(diff -y --suppress-common-lines result/no_frontier.txt result/original.txt | wc -l)
echo "diff no_frontier: ${DD1}"

DD1=$(diff -y --suppress-common-lines result/pyver.txt result/original.txt | wc -l)
echo "diff between pyver: ${DD1}"
