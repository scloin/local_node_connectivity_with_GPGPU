home=/home/deepl/sooho/test/K-com/test

while getopts n: opt # if first argument str length is not 0
do
  #echo $1
  case ${opt} in
    n) N="$OPTARG";;
    esac
done

time0=$({ time ${home}/exec/revise4t >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time1=$({ time ${home}/exec/revise4 >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')
time2=$({ time python ${home}/pyver.py >/dev/null; } 2>&1 | grep real | grep -o '[[:digit:]].*$')

echo "revise4t   (${time0})"
echo "revise4    (${time1})"
echo "python     (${time2})"

python sort.py ${N}

DD1=$(diff -y --suppress-common-lines result/thread.txt result/test2.txt | wc -l)
echo "diff: ${DD1}"

DD1=$(diff -y --suppress-common-lines result/test1.txt result/test2.txt | wc -l)
echo "diff between pyver: ${DD1}"
