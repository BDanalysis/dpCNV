"dc=0.5,prob=0.01"

'跑真实数据，对预处理部分的代码作很多修改，不一定适合仿真数据'
import numpy as np
import pysam
import math
import sys
import matplotlib.pyplot as plt
from numba import njit
from sklearn.metrics import euclidean_distances
from scipy.stats import multivariate_normal
from sklearn import preprocessing
import rpy2.robjects as robjects
import os
import pandas as pd
import datetime
import subprocess

chro_num = 168
print('12314444444444444444444444444')
def get_chrlist(filename):
    samfile = pysam.AlignmentFile(filename, "rb")
    List = samfile.references	#List=[1,2,3,,,,21,x,y,mt,nc....]  length=113
    chrList = np.full(len(List), 0)
    for i in range(len(List)):
        chr = str(List[i]).strip('chr')
        #print('chr:',chr)
        if chr.isdigit():
            chrList[i] = int(chr)
    chrList = np.delete(chrList,np.where(chrList==0))   # 刪除掉不属于22条染色体的元素
    #print('chrList:',chrList)
    chrList=np.array([chro_num])
    print('chrList:',chrList)   #[21]
    return chrList

def get_RC(filename, chrList, ReadCount):
    samfile = pysam.AlignmentFile(filename, "rb")
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            #print('chr:',type(chr))
            #print('chr:',chr)
            if chr.isdigit():
                num = np.argwhere(chrList == int(chr))[0][0]#####################################33
                #print('num',num)
                posList = line.positions
                #print('position:',posList)
                #print('-----------')
                #print('posList',posList)	
                ReadCount[num][posList] += 1
    print('len(num)',num.size)	#num就是染色体的个数
    print('lenposlist',len(posList))	#每个read的长度
    return ReadCount

def read_ref_file(filename, ref, num):
    # read reference file
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            for line in f:
                linestr = line.strip()
                ref[num] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref


def ReadDepth(ReadCount, binNum, ref):
    RD = np.full(binNum, 0.0)
    GC = np.full(binNum, 0)
    pos = np.arange(1, binNum+1)
    for i in range(binNum):
        RD[i] = np.mean(ReadCount[i * binSize:(i + 1) * binSize])
        cur_ref = ref[i * binSize:(i + 1) * binSize]
        N_count = cur_ref.count('N') + cur_ref.count('n')
        if N_count == 0:
            gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
        else:
            RD[i] = -10000
            gc_count = 0
        GC[i] = int(round(gc_count / binSize, 3) * 1000)

    index = RD > 0
    RD = RD[index]
    GC = GC[index]
    pos = pos[index]
    RD = gc_correct(RD, GC)
    
    #print('maxRD:',max(RD),' minRD:',min(RD))
    #plot(pos,RD)
    return pos, RD


def gc_correct(RD, GC):
    # correcting gc bias
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean
    return RD


def RDsimilar(RD):
    simRD = np.full(len(RD), 0.0)
    for k in range(1, len(RD)-1):
        dot = 1 + abs(seg_rd[k+1] - seg_rd[k]) * abs(seg_rd[k] - seg_rd[k-1])
        data = (1 + (pow(seg_rd[k+1] - seg_rd[k], 2))) * (1 + (pow(seg_rd[k] - seg_rd[k-1], 2)))
        simRD[k] = dot/data
    simRD[0] = simRD[1]
    simRD[-1] = simRD[-2]
    return simRD


def read_seg_file(filename):
    segRD = []
    with open(filename, 'r') as f:
        line = f.readline()
        for line in f:
            linestr = line.strip()
            linestrlist = linestr.split('\t')
            segRD.append(float(linestrlist[1]))
    segRD = np.array(segRD)
    return segRD

#没用到这个函数？？
def write(data, data1, data2):
    output = open("segrd.txt", "w")
    for i in range(len(data)):
        output.write(str(data[i]) + '\t' + str(data1[i]) + '\t' + str(data2[i]) + '\n')


def dis_matrix(RD):
    # calculating euclidean_distances matrix

    RD = RD.astype(np.float)
    pos = np.arange(1, len(RD)+1)
    #pos = simRD.astype(np.float)
    nr_min = np.min(RD)
    nr_max = np.max(RD)
    newpos = (pos - min(pos)) / (max(pos) - min(pos)) * (nr_max - nr_min) + nr_min
    newpos = newpos.astype(np.float)

    RD = RD.astype(np.float)
    newpos = newpos.astype(np.float)
    rd = np.c_[RD, newpos]
    print("dis_matrix calculating")
    dis_matrix = euclidean_distances(rd, rd)
    return dis_matrix


def density(dis_matrix):
    num = len(dis_matrix[0])
    ds = np.full(num, 0)
    for i in range(num):
        ds[i] = np.sum(dis_matrix[i] < dc)
    return ds


def distance(dis_matrix, ds):
    center = []
    num = len(dis_matrix[0])
    dt = np.full(num, 0.0)
    for i in range(num):
        index = ds > ds[i]
        if np.sum(index) == 0:
            center.append(i)
            dt[i] = np.max(dis_matrix[i])
        else:
            dt[i] = np.min(dis_matrix[i][index])
    #print('dis_matrix',dis_matrix)
    print('center:',center)
    print('len(center):',len(center))
    return dt, center


def get_dc_matrix(dis_matrix):
    num = len(dis_matrix[0])
    dc_matrix = []
    pos = np.arange(num)
    for i in range(num):
        index = dis_matrix[i] < dc
        dc_matrix.append(pos[index])
    return dc_matrix


def read_RC_file(filename, data):
    with open(filename, 'r') as f:
        for line in f:
            linestr = line.strip()
            linestrlist = linestr.split('\t')
            pos = int(linestrlist[1])
            data[pos] = int(linestrlist[3])
    return data


def segment(pos, segrd):
    start = []
    end = []
    seg_rd = []
    i = 0
    j = 1
    while j < len(segrd):
        if j == len(segrd) - 1:
            start.append(int(pos[i]))
            end.append(int(pos[j-1]))
            seg_rd.append(float(segrd[i]))
            j += 1
        else:
            if segrd[i] == segrd[j]:
                j += 1
            else:
                start.append(int(pos[i]))
                end.append(int(pos[j-1]))
                seg_rd.append(float(segrd[i]))
                i = j
    return start, end, seg_rd


def write_CNV(filename, chr, start, end, type):
    count = 0
    output = open(filename, "w")
    for i in range(len(start)):
        if type[i] == 2:
            count += 1
            output.write('chr' + str(chr) + '\t' + str(start[i] + 1) + '\t' + str(end[i]) + '\t' + 'gain' + '\n')
        elif type[i] == 1:
            count += 1
            output.write('chr' + str(chr) + '\t' + str(start[i] + 1) + '\t' + str(end[i]) + '\t' + 'loss' + '\n')
    print('detected count:',count)


def plot(pos, data):
    plt.scatter(pos, data, s=3, c="black")
    #plt.scatter(pos1, data1, s=3, c="red")
    plt.xlabel("ds")
    plt.ylabel("rd")
    plt.show()


def caculating_CNV(dc_m, center, ds, dt, start, end, rd):
    normRD = np.mean(rd[center])
    num = len(rd)
    flag = np.full(num, 0)
    CNV_start = []
    CNV_end = []
    CNV_type = []
    D_max_value = 0

    ds_score = preprocessing.scale(ds)
    dt_score = preprocessing.scale(dt)

    for i in range(len(center)):
        pos = dc_m[center[i]]
        flag[pos] = 1
        D_value = abs(rd[pos] - normRD)
        if max(D_value) > D_max_value:
            D_max_value = max(D_value)

    rd_value = abs(rd - normRD)
    mean_ds = np.mean(ds_score)
    mean_dt = np.mean(dt_score)
    mu = np.array([mean_ds, mean_dt])
    sigma = np.cov(ds_score, dt_score)

    print('mu',mu)
    print('sigma',sigma)
    #print('ds_score',ds_score)
    #print('dt_score',dt_score)
    print('mean_ds',mean_ds)

    var = multivariate_normal(mean=mu, cov=sigma)

    for i in range(num):
        prob = var.pdf([ds_score[i], dt_score[i]])
        if prob < 0.01 and ds_score[i] < mean_ds:
            if rd[i] < normRD:
                type = 1
                CNV_type.append(int(1))
            else:
                type = 2
                CNV_type.append(int(2))
            CNV_start.append(start[i] * binSize)
            CNV_end.append(end[i] * binSize + binSize)
        #print(start[i] * binSize+ 1, end[i] * binSize + binSize, prob, ds_score[i], rd[i])

    CNVstart = np.array(CNV_start)
    CNVend = np.array(CNV_end)
    CNVtype = np.array(CNV_type)

    for i in range(len(CNVtype)-1):
        if CNVstart[i+1] <= CNVend[i] and CNVtype[i] == CNVtype[i+1]:
            CNVstart[i+1] = CNVstart[i]
            CNVtype[i] = 0
    index = CNVtype > 0
    CNVstart = CNVstart[index]
    CNVend = CNVend[index]
    CNVtype = CNVtype[index]

    return CNVstart, CNVend, CNVtype

starttime = datetime.datetime.now()
# get params
#path = "/Volumes/TOSHIBA/NGS_data/SInC/10x/"
#outpath = "/Volumes/TOSHIBA/NGS_data/SInC/10x/dp_result/"
#bamName = sys.argv[1]
#bam = path + bamName
bam = sys.argv[1]
num_ = ''
for each_str in range(len(bam)):
    if bam[each_str].isdigit():
        num_ += bam[each_str]

chro_num = int(num_)
#bam = "sim1_6x_150_sort.bam"
binSize = 1000
dc = 0.5		#参数，在0.5-1.5之间调整
CNVfile = bam + '_result.txt'
chrList = get_chrlist(bam)	
chrNum = len(chrList)

refList = [[] for i in range(chrNum)]
#print('chrNum:',chrNum)	
for i in range(chrNum):
    #reference = ["chr1.fa","chr2.fa","chr3.fa","chr4.fa","chr5.fa","chr6.fa","chr7.fa","chr8.fa","chr9.fa","chr10.fa","chr11.fa","chr12.fa","chr13.fa","chr14.fa","chr15.fa","chr16.fa","chr17.fa","chr18.fa","chr19.fa","chr20.fa","chr21.fa","chr22.fa",]

    reference = ["chr"+str(chro_num)+".fa"]
    print("ref:",reference)
    refList = read_ref_file("/home/ty/桌面/真实数据/hg37/"+reference[i], refList, i)
    #refList = read_ref_file(reference[i], refList, i)

chrLen = np.full(chrNum, 0)
for i in range(chrNum):
    chrLen[i] = len(refList[i])
    print('chrLen:',chrLen[i])		#48129895
print('maxChrLen',np.max(chrLen))  	#48129895
print("Read bam file:", bam)
ReadCount = np.full((chrNum, np.max(chrLen)), 0)
ReadCount = get_RC(bam, chrList, ReadCount)
print('ReadCount:',ReadCount[0].size)	#48129895
#plot(np.arange(ReadCount[0].size),ReadCount)	

print('ReadCount,',ReadCount)	
for i in range(chrNum):
    binNum = int(chrLen[i]/binSize)+1
    pos, RD = ReadDepth(ReadCount[0], binNum, refList[i])

    #plot(pos,RD)
    v = robjects.FloatVector(RD)
    m = robjects.r['matrix'](v, ncol=1)
    robjects.r.source("segment.R")
    robjects.r.segment(m)

    #subprocess.call('Rscript segment.R',shell=True)
    segFile = "seg.txt"
    segRD = read_seg_file(segFile)
    start, end, seg_rd = segment(pos, segRD)

    #plot(np.arange(len(seg_rd)),seg_rd)
    start = np.array(start)
    end = np.array(end)
    seg_rd = np.array(seg_rd)
    #simRD = RDsimilar(seg_rd)
    dis_m = dis_matrix(seg_rd)
    print("calculate density")
    ds = density(dis_m)
    print("calculate distance")
    dt, center = distance(dis_m, ds)
    dc_matrix = get_dc_matrix(dis_m)
    #plot(ds,dt)
    CNVstart, CNVend, CNVtype = caculating_CNV(dc_matrix, center, ds, dt, start, end, seg_rd)
    write_CNV(CNVfile, chrList[i], CNVstart, CNVend, CNVtype)
    #subprocess.call('rm RD seg.txt',shell=True)
    endtime = datetime.datetime.now()
    print("running time: " + str((endtime - starttime).seconds) + " seconds")
    print('----------------------------------------------------')
