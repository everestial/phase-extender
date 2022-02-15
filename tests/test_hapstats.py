import os
import shutil
import filecmp
import tempfile
import pandas as pd
from hapstats import compute_haplotype_stats
from tutils import replace_mkdir, is_same



df1 = pd.read_csv('tests/inputs/eg1/input_haplotype_file.txt', sep='\t')
df2 = pd.read_csv('tests/inputs/eg2/haplotype_file_test02.txt', sep='\t', na_values= '.', keep_default_na= False)

def test_small_file(tmp_path):
    ## TODO: tmp_path and tmpdir is not working in this scenario
    with tempfile.TemporaryDirectory() as outputdir1:
        compute_haplotype_stats(df1, 'ms02g', 'initial',str(outputdir1))
        assert is_same(outputdir1,'tests/outdir/plotseg1')
        
def test_hapstats_large_file(tmp_path):
    outputdir2 = tmp_path.mkdir(exist_ok = True)
    with tempfile.TemporaryDirectory() as outputdir2:
        compute_haplotype_stats(df2, 'ms02g', 'initial',outputdir2)
        assert is_same(outputdir2, 'tests/outdir/plotseg2')




