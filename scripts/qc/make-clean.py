#!/usr/bin/env python

import os
import glob
from tqdm import tqdm

idx = {}
path_in = "/mnt/neuro-genomic-1-ro/glass/RNAseq/fastq"

for _ in tqdm(sorted(glob.glob(path_in + "/*.fastq.gz"))):
    sid = _.split('_')[1]
    res = _.split('_')[-1].replace('.fastq.gz','')

    if sid not in idx:
        idx[sid] = {'R1':[] , 'R2': []}

    idx[sid][res].append(_)


for key in idx:
    print("")
    
    if(len(idx[key]['R1']) != len(idx[key]['R2']) ):
        print ("ERR: " + key )
        import sys
        sys.exit(1)
    
    for i in range(len(idx[key]['R1'])):
        if idx[key]['R1'][i].replace("_R1","_R2") != idx[key]['R2'][i]:
            print ("ERR: " + str(idx[key]) )
            import sys
            sys.exit(1)
        else:
            r1_in = idx[key]['R1'][i]
            r2_in = idx[key]['R2'][i]
            
            basename = os.path.basename(r1_in).replace('_R1.fastq.gz','')
            #print(basename)

            log_json = '/mnt/neuro-genomic-1-rw/glass/RNAseq/fastq-clean/' + basename + '_fastp_log.json'
            log_html = '/mnt/neuro-genomic-1-rw/glass/RNAseq/fastq-clean/' + basename + '_fastp_log.html'

            r1_out = '/mnt/neuro-genomic-1-rw/glass/RNAseq/fastq-clean/' + basename + '_trimmed_R1.fastq.gz'
            r2_out = '/mnt/neuro-genomic-1-rw/glass/RNAseq/fastq-clean/' + basename + '_trimmed_R2.fastq.gz'

            #print(r1_out)
            #fastp -y -x
            # -i "/data/users/youri/mnt/neurogen-ro/gsam/RNA/Liege_part_2/"$sid"R1_001.fastq.gz"
            # -I "/data/users/youri/mnt/neurogen-ro/gsam/RNA/Liege_part_2/"$sid"R2_001.fastq.gz"
            # -o "fastq-clean/"$sid"001.clean_yx_R1.fastq.gz"
            # -O "fastq-clean/"$sid"001.clean_yx_R2.fastq.gz"
            # -j "log/"$sid"001.xy.json" -h "log/"$sid"001.xy.html"
            # --trim_poly_g -l 35
            out = ['nice', '/home/youri/projects/glass/bin/fastp/fastp',
                '-y',
                '-x',
                '-p',
                '--thread', '4',
                '-i', '"'+r1_in+'"',
                '-I', '"'+r2_in+'"',
                '-o', '"'+r1_out+'"',
                '-O', '"'+r2_out+'"',
                '-j', '"'+log_json+'"',
                '-h', '"'+log_html+'"',
                '--trim_poly_g',
                '-l', '35'
            ]
            
            print(' '.join(out))
            

