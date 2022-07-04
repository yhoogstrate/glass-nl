#!/usr/bin/env python

"""
# adapt the bed files with minor annotation fields:


python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.025minAF.bed

python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min2_pos_0.05minAF.bed

python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.025minAF.bed

python3 ./scripts/py/add_column.py "PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed" data/glass/WES/Assets/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed > tmp/PON_n207_LoFreq_Hyperexome_WES_min3_pos_0.05minAF.bed

"""


import click
from tqdm import tqdm

@click.command()
@click.argument('Column_value_to_append')
@click.argument('fh_in', type=click.File('r'))
def main(column_value_to_append, fh_in):
  for line in tqdm(fh_in):
    print(line.strip() + "\t" + column_value_to_append)

if __name__ == '__main__':
    main()


