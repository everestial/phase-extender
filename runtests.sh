

python3 phase_extender.py --input eg1/input_haplotype_file.txt --output outdir/eg1out --SOI ms02g --lods 10


# Namespace(SOI='ms02g', addMissingSites='no', bed='', culLH='maxPd', hapStats='no',
# input='eg1/input_haplotype_file.txt', lods='10', nt=1, numHets=40,
#  output='outdir/eg1out', refHap='', snpTh=3, useSample='all', writeLOD='no')


python3 phase_extender.py --input eg1/input_haplotype_file.txt --SOI ms02g --writeLOD yes
# Namespace(SOI='ms02g', addMissingSites='no', bed='', culLH='maxPd', hapStats='no',
# input='eg1/input_haplotype_file.txt', lods=5, nt=1, numHets=40, 
# output='outdir/eg2out', refHap='', snpTh=3, useSample='all', writeLOD='yes')


python3 phase_extender.py  --nt 2 --input eg2/haplotype_file_test02.txt  --output outdir/eg3out --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --lods 10
# Namespace(SOI='ms02g', addMissingSites='no', bed='', culLH='maxPd', hapStats='yes',
# input='eg2/haplotype_file_test02.txt', lods='10', nt='2', numHets='25',
# output='outdir/eg3out', refHap='', snpTh=3, useSample='all', writeLOD='no')



python3 phase_extender.py  --input eg2/haplotype_file_test02.txt  --output outdir/eg4out --SOI ms02g --useSample ms01e,ms02g,MA605,Sp76

# Namespace(SOI='ms02g', addMissingSites='no', bed='', culLH='maxPd', hapStats='no', 
# input='eg2/haplotype_file_test02.txt', lods=5, nt=1, numHets=40, 
# output='outdir/eg4out', refHap='', snpTh=3, useSample='ms01e,ms02g,MA605,Sp76', writeLOD='no')

python3 phase_extender.py --nt 1 --input eg2/haplotype_file_test02.txt --output outdir/eg5out --SOI ms02g --numHets 25 --culLH maxPd --hapStats yes --refHap eg2/refPanel_lyrata_test02.txt --bed eg2/bed_boundries.bed

# Namespace(SOI='ms02g', addMissingSites='no', bed='eg2/bed_boundries.bed', culLH='maxPd',
#  hapStats='yes', input='eg2/haplotype_file_test02.txt', lods=5, nt='1', numHets='25',
#   output='outdir/eg5out', 
# refHap='eg2/refPanel_lyrata_test02.txt', snpTh=3, useSample='all', writeLOD='no')