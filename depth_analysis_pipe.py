import os 
import pandas as pd
import datetime
start_time = datetime.datetime.now()
depth = 4000000
os.system('mkdir ./tmp')
os.system('ls *star_gene_exon_tagged_TagIntronic_clean.bam > ./tmp/file_list_1.txt')
file_list = pd.read_csv('./tmp/file_list_1.txt', sep='\t', names=['file'])
bam_lis = list(file_list['file'])

for i in bam_lis:
    print(i, 'analysis start!')
    os.system('python depth_analysis.py -b {} -d {} -p {}'.format(i, depth, i.split('_')[0]))

print('all done')
os.system('python /home/songjia/bigdisk/xx/dingding.py')