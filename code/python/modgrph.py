import numpy as np

def reader(A: str):
    return list(map(int,A.split(',')))

fp = open("./graph/fb_forum/fb-forum.edges", "r"); fp.readline()
# _, N, E = reader(fp.readline())
N=899; E=33719
ori = [[] for _ in range(N)]
mod = []
for i in range(E):
    l = reader(fp.readline())
    ori[l[1]-1]+=[l[0]-1]
    ori[l[0]-1]+=[l[1]-1]
fp.close

fp1 = open("./graph/forum.txt", 'w')
dest=''
edges='0'
num=0
for i in (ori):
    dest+=' '.join(list(map(str,i)))
    dest+=' '
    num+=len(i)
    edges+=' %d' %(num)
E*=2

fp1.writelines([str(E)+"\n", str(N)+"\n", dest+"\n", edges])
fp1.close()

#list2 = np.concatenate(list1).tolist()