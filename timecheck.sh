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

time0=$({ time ${home}/exec/revise4t >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time1=$({ time ${home}/exec/revise4 >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time2=$({ time python ${home}/python/pyver.py >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')

echo "ver addthread  (${time0})"
echo "ver origin     (${time1})"
echo "ver python     (${time2})"

python ${home}/python/sort.py ${N}

DD1=$(diff -y --suppress-common-lines result/addthread.txt result/original.txt | wc -l)
echo "diff: ${DD1}"

DD1=$(diff -y --suppress-common-lines result/pyver.txt result/original.txt | wc -l)
echo "diff between pyver: ${DD1}"
