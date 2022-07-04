#!/bin/bash

### adapt the bed files with minor annotation fields
# python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed
# python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed
# python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed
# python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed



cat /tmp/vcf-test/219-R3-I1_intersect_fun.vcf | \
    bcftools annotate \
    -a tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed \
    -h ~/projects/glass/assets/PON_LoFreq_2_0_025.annots.hdr \
    -c CHROM,FROM,TO,PON_LoFreq_2_0_025 \
    - | \
 \
    bcftools annotate \
    -a tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed \
    -h ~/projects/glass/assets/PON_LoFreq_2_0_05.annots.hdr \
    -c CHROM,FROM,TO,PON_LoFreq_2_0_05 \
    - | \
 \
    bcftools annotate \
    -a tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed \
    -h ~/projects/glass/assets/PON_LoFreq_3_0_025.annots.hdr \
    -c CHROM,FROM,TO,PON_LoFreq_3_0_025 \
    - | \
 \
    bcftools annotate \
    -a tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed \
    -h ~/projects/glass/assets/PON_LoFreq_3_0_05.annots.hdr \
    -c CHROM,FROM,TO,PON_LoFreq_3_0_05 \
    > /tmp/vcf-test/test4.vcf


