import re
import sys

if len(sys.argv) < 1:
    print("120, 27(default)")
    num=20
    #sys.exit()
else:
    num =int(sys.argv[1])
    #print(N,S)

af=[]
af.append("./result/frontier.txt")
af.append("./result/no_frontier.txt")

for fs in af:
    f=open(fs, 'r')
    array = [[1 for col in range(row,num,1)] for row in range(1,num,1)]

    line =f.read()
    f.close
    import re

    regex = r"^\W(\d*)\W*(\d*)\W*(\d*)$"

    matches = re.finditer(regex, line, re.MULTILINE)

    for match in (matches):
        group0 = int(match.group(1))
        group1 = int(match.group(2))
        group2 = int(match.group(3))
        array[group0][group1-group0-1]=group2
    
    f=open(fs, 'w')
    for i,arr in enumerate(array,0):
        for j,val in enumerate(arr,0):
            f.write("[%d, %d] %d\n" %(i,j+i+1,val))
    f.close