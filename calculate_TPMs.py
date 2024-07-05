import pandas as pd
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE', description='Calculate TPMs from featureCounts file')
parser.add_argument('-i', '--input', metavar='', action='store', type=str, required=True, help='Input featureCounts text file (TSV)')

args = parser.parse_args()
fc_file = Path(args.input)
fc_df = pd.read_csv(f"{fc_file}", delim_whitespace=True, skiprows=1)

rm_cols = ['Chr', 'Start', 'End', 'Strand']
cols = [col for col in fc_df.columns if col not in rm_cols]
fc_df = fc_df.filter(cols)
fc_df = fc_df.set_index('Geneid')
rpk_df = fc_df
rpk_df.iloc[:,1:] = rpk_df.iloc[:,1:].div(rpk_df.Length, axis=0)
tpm_df = rpk_df
scale_dict = {}
for col in [col for col in rpk_df.columns if col != "Length"]:
    norm_tm = rpk_df[col].sum()
    tpm_df[col] = rpk_df[col].mul(1000000).div(norm_tm)

print(f'TPM file saved to: {fc_file.stem}_tpms.tsv')
tpm_df.to_csv(f'{fc_file.stem}_tpms.tsv', sep='\t', header=True, index=True)
