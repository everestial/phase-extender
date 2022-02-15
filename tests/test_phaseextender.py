import filecmp
import os
import shutil
import tempfile


from phaser import phase_converter
from tutils import replace_mkdir, is_same

# initialize the default parameters
soi ='ms02g'
addmissingsites='no'
bed_file = None
maxed_as = '*'
hapstats='no'
input_file='tests/inputs/eg1/input_haplotype_file.txt'
lods=10
nt=2
numHets=40
snp_threshold=3
use_sample='all'
writelod='no'
refhap = None


def test_first_example():
    with tempfile.TemporaryDirectory() as tmpdir1:
        phase_converter(soi, replace_mkdir(tmpdir1), nt, input_file, lods, snp_threshold, numHets,  maxed_as, bed_file, refhap, use_sample ,hapstats, writelod, addmissingsites )
        assert is_same('tests/outdir/eg1out/', tmpdir1) is True
    # assert comparison.diff_files is None and  comparison.right_only is None and comparison.left_only is None, **kwargs


def test_second_example():
    writelod= 'yes'
    with tempfile.TemporaryDirectory() as temp_dir:
        phase_converter(soi, temp_dir, nt, input_file, lods, snp_threshold, numHets,  maxed_as, bed_file, refhap, use_sample ,hapstats, writelod, addmissingsites )
        assert is_same('tests/outdir/eg2out', temp_dir ) is True



def test_third_example():
    lods= 10
    nt= 1
    numHets= 25
    input_file='tests/inputs/eg2/haplotype_file_test02.txt'
    with tempfile.TemporaryDirectory() as temp_dir1:
        phase_converter(soi, temp_dir1, nt, input_file, lods, snp_threshold, numHets,  maxed_as, bed_file, refhap, use_sample ,'yes', writelod, addmissingsites )
        assert is_same('tests/outdir/eg3out', temp_dir1 ) is True

def test_fourth_example():
    lods= 5
    nt= 4
    numHets= 40
    use_sample = 'ms01e,ms02g,MA605,Sp76'
    input_file='tests/inputs/eg2/haplotype_file_test02.txt'
    with tempfile.TemporaryDirectory() as temp_dir2:
        phase_converter(soi, temp_dir2, nt, input_file, lods, snp_threshold, numHets,  maxed_as, bed_file, refhap, use_sample ,hapstats, writelod, addmissingsites )
        assert is_same('tests/outdir/eg4out', temp_dir2 ) is True

def test_fifth_example():
    lods= 5
    nt= 4
    numHets= 25
    use_sample = 'all'
    input_file ='tests/inputs/eg2/haplotype_file_test02.txt'
    refhap = 'tests/inputs/eg2/refPanel_lyrata_test02.txt'
    bed_file = 'tests/inputs/eg2/bed_boundries.bed'
    hapstats = 'yes'
    with tempfile.TemporaryDirectory() as temp_dir3:
        phase_converter(soi, temp_dir3, nt, input_file, lods, snp_threshold, numHets,  maxed_as, bed_file, refhap, use_sample ,hapstats, writelod, addmissingsites )
        assert is_same('tests/outdir/eg5out', temp_dir3 ) is True



def replace_mkdir(dir_name):
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name, ignore_errors=False, onerror=None)
    os.makedirs(dir_name, exist_ok=True)
    return dir_name


