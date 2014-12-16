import os
import matplotlib.pyplot as plt

x=[]
y=[]
t=[]

# Tree:
#   A           C
#    \         /
#     \       /
#      \i___j/
#      /     \ 
#     B       D

for bl in range(1,30):
    Ai=bl/10.0
    Bi=0.03
    Cj=bl/10.0
    Dj=0.03
    ij=0.03
    L=10000
    nexfi="sim.nex"

    pairwise_true={}
    pairwise_true[('A','B')]=Ai+Bi
    pairwise_true[('A','C')]=Ai+ij+Cj
    pairwise_true[('A','D')]=Ai+ij+Dj
    pairwise_true[('B','C')]=Bi+ij+Cj
    pairwise_true[('C','D')]=Dj+Cj

    for i in range(5):
#            os.system('python gap_sim.py -ta {} -tb {} -tc {} -td {} -ti {} -l {} > {}'.format(Ai,Bi,Cj,Dj,ij,L,nexfi))
            newfi=open('sim.tre','w')
            newfi.write('((A:{},B:{}):{},(C:{},D:{}):0.1);'.format(Ai,Bi,ij,Cj,Dj))
            newfi.close()
            os.system('seq-gen -mGTR -on < {} > {}'.format('sim.tre',nexfi))
            os.system('paup -n paupblock0 > out.p')
            fi=open('out.p').readlines()
            pinv=[]
            for lin in fi:
                if lin.startswith("  Estimated value of proportion of invariable sites = "):
                    pinvar=lin.split()[-1]
            pinv.append(pinvar)
            pb1=open("paupblock1","w")
            pb1.writelines(["execute {};".format(nexfi),"set criterion=distance;\n","dset dist=jc pinv={};".format(pinvar),"[!pairwise corrected dists]","showdist;","[!end pairwise corrected dists]"])
            pb1.close()
            os.system('paup -n paupblock1 > out1.p')
            out1=open('out1.p').readlines()
            taxa=['A','B','C','D']
            pairwise_p={}
            for i, lin in enumerate(out1):
                if "1         2         3         4" in lin or lin.startswith("             1        2        3        4\n"):
                    sp=i
            for lin in out1[sp+1:sp+5]:
                lii=lin.split()
                tax=lii[1]
                for i, item in enumerate(lii[2:]):
                    if item != '-':
        #            pairwise_p[(tax, taxa[i])]=float(item)
                        try:
                            pairwise_p[(taxa[i], tax)]=float(item)
                        except:
                            pairwise_p[(taxa[i], tax)]=float(item[1:])
            trf=open('tmp.nex').readlines()
            import re
            alt=[]
            for lin in trf:
                lii=re.split(r'(\s+)', lin)
                for i,item in enumerate(lii):
                  if 'e-' in item:
                    if ';' in item:
                        lii[i]='0;'
                    else:
                        lii[i]='0'
                alt.append(''.join(lii))
            open('tmp_no_e.nex','w').writelines(alt)
            os.system('paup -n paupblock2 > out2.p')
            fi2=open('out2.p').readlines()
            taxa=['A','B','C','D']
            pairwise_ML={}
            for i, lin in enumerate(fi2):
                if lin.startswith("Patristic distance matrix"):
                    startpoint=i
            for lin in fi2[startpoint+6:startpoint+10]:
                lii=lin.split()
                tax=lii[1]
                for i, item in enumerate(lii[2:]):
                    if item != '-':
                        try:
                            pairwise_ML[(tax, taxa[i])]=float(item)
                        except:
                            pairwise_ML[(tax, taxa[i])]=float(item[1:])
            for item in pairwise_p.keys():
                if item not in pairwise_ML.keys():
                    print(item)
            for key in pairwise_true:
                x.append(pairwise_p[key])
                y.append(pairwise_ML[key])
                t.append(pairwise_true[key])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    p = ax.plot(t, x, 'ro')
    ax.set_xlabel('pairwise true')
    ax.set_ylabel('pairwise p')
    plt.savefig('true_p.pdf')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p = ax.plot(x, y, 'ro')
    ax.set_xlabel('pairwise p')
    ax.set_ylabel('pairwise ML')
    plt.savefig('p_ML.pdf')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    p = ax.plot(t, y, 'ro')
    ax.set_xlabel('true t')
    ax.set_ylabel('pairwise ML')
    plt.savefig('true_ML.pdf')