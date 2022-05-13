import click
import pandas as pd


fam_col_names = ['fid', 'iid', 'pid', 'mid', 'sex', 'batch']

sqc_col_names = [x.lower() for x in list(
    ['affymetrix.1', 'affymetrix.2', 'genotyping.array', 'Batch',
     'Plate.Name', 'Well', 'Cluster.CR', 'dQC', 'Internal.Pico..ng.uL.',
     'Submitted.Gender', 'Inferred.Gender', 'X.intensity', 'Y.intensity',
     'Submitted.Plate.Name', 'Submitted.Well', 'sample.qc.missing.rate',
     'heterozygosity', 'heterozygosity.pc.corrected',
     'het.missing.outliers', 'putative.sex.chromosome.aneuploidy',
     'in.kinship.table', 'excluded.from.kinship.inference',
     'excess.relatives', 'in.white.British.ancestry.subset',
     'used.in.pca.calculation'] +
    [f'{a}{b}' for a, b in zip(['pc'] * 40, list(range(1, 40 + 1)))] +
    ['in.Phasing.Input.chr1_22', 'in.Phasing.Input.chrX',
     'in.Phasing.Input.chrXY'])]


@click.command()
@click.option('--fam', help='Path to fam file')
@click.option('--sqc', help='Path to ukb_sqc_v2.txt')
@click.option('--out', help='Path to output file')
def write_sqc_exclusions(fam, sqc, out):
    sqc = pd.read_csv(sqc, names=sqc_col_names, sep=r'\s+')
    # sqc.insert(0, "sample_index", range(1, 1 + len(sqc)))
    fam = pd.read_csv(fam, names=fam_col_names, sep=r'\s+')

    sqc['eid'] = fam['iid']

    sqc = sqc.loc[
        (sqc['submitted.gender'].ne(sqc['inferred.gender'])) |
        (sqc['het.missing.outliers'] == 1) |
        (sqc['putative.sex.chromosome.aneuploidy'] == 1) |
        (sqc['excess.relatives'] == 1),
        ['eid', 'eid']
    ]

    sqc.to_csv(out, sep=" ", header=False, index=False)


if __name__ == '__main__':
    write_sqc_exclusions()
