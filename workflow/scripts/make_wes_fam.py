import click
import pandas as pd


@click.command()
@click.option('-f', '--fam', help='Path to WES reference fam file')
@click.option('-b', '--bridge', help='Path to bridging file')
@click.option('-o', '--out', help='Path to output file')
def make_wes_fam(bridge, fam, out):
    """Writes WES fam with application specific IDs"""
    fam_col_names = ['fid', 'iid', 'pid', 'mid', 'sex', 'batch']
    fam = pd.read_csv(fam, names=fam_col_names, sep=r'\s+')
    bridge = pd.read_csv(bridge, header=0, names=['eid_new', 'eid_18177'])

    lookup = bridge.set_index('eid_18177').to_dict()['eid_new']

    fam['fid'] = (fam['fid']
                  .map(lookup)
                  .fillna(-fam.index.to_series().astype(int))
                  .astype(int))

    fam['iid'] = fam['fid']
    fam.to_csv(out, sep=' ', header=None, index=None)


if __name__ == '__main__':
    make_wes_fam()
