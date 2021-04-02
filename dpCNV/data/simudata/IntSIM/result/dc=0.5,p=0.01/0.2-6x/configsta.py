import subprocess

#params list
binLen = 1000

bam = ['sim%d_6_6100_read.sort.bam_result.txt' % (i) for i in range(1, 51)]

for i in range(len(bam)):
    subprocess.call('python sta_score.py ' + bam[i], shell = True)
subprocess.call('python cal_mean.py', shell = True)
