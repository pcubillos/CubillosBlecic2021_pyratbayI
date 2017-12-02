# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# Clone (download) the necessary code:
cd $topdir
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout cea5ca0  # FINDME: Update when final
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 4ba3633  # FINDME: Update when final
make

# Generate filter files:
cd $topdir/code
$topdir/code/make_filters.py > $topdir/code/filter_info.txt

# Download Exomol data:
cd $topdir/inputs/opacity
wget -i wget_exomol_NH3.txt
wget -i wget_exomol_CH4.txt
wget -i wget_exomol_HCN.txt

# Download HITRAN/HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_H2O_CO2.txt
unzip '*.zip'
rm -f *.zip

# Download CO data:
cd $topdir/inputs/opacity
wget http://iopscience.iop.org/0067-0049/216/1/15/suppdata/apjs504015_data.tar.gz
tar -xvzf apjs504015_data.tar.gz
rm -f apjs504015_data.tar.gz ReadMe Table_S1.txt Table_S2.txt \
      Table_S3.txt Table_S6.par

# Download TiO data:
cd $topdir/inputs/opacity
wget http://kurucz.harvard.edu/molecules/tio/tioschwenke.bin
wget http://kurucz.harvard.edu/molecules/tio/tiopart.dat


# Generate partition-function files for H2O and NH3:
cd $topdir/run01
$topdir/code/pf_tips_H2O-NH3.py

# Generate partition-function files for HCN, CH4, and TiO:
cd $topdir/run01
$topdir/pyratbay/scripts/PFformat_Exomol.py        \
      $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
      $topdir/inputs/opacity/1H-13C-14N__Larner.pf
$topdir/pyratbay/scripts/PFformat_Exomol.py \
      $topdir/inputs/opacity/12C-1H4__YT10to10.pf
$topdir/pyratbay/scripts/PFformat_Schwenke_TiO.py \
      $topdir/inputs/opacity/tiopart.dat


# Compress LBL databases:
cd $topdir/run01
$topdir/repack/repack.py repack_H2O.cfg
$topdir/repack/repack.py repack_HCN.cfg
$topdir/repack/repack.py repack_NH3.cfg
$topdir/repack/repack.py repack_CH4.cfg
$topdir/repack/repack.py repack_TiO.cfg


# Make TLI files:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c tli_hitemp_CO2.cfg
$topdir/pyratbay/pbay.py -c tli_hitemp_CO.cfg
$topdir/pyratbay/pbay.py -c tli_hitemp_H2O.cfg
$topdir/pyratbay/pbay.py -c tli_exomol_HCN.cfg
$topdir/pyratbay/pbay.py -c tli_exomol_NH3.cfg
$topdir/pyratbay/pbay.py -c tli_exomol_CH4.cfg


# Make atmospheric files:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c atm_1xsolar_1000K.cfg
$topdir/pyratbay/pbay.py -c atm_1xsolar_1500K.cfg
$topdir/pyratbay/pbay.py -c atm_1xsolar_2000K.cfg
$topdir/pyratbay/pbay.py -c atm_1xsolar_2500K.cfg
$topdir/pyratbay/pbay.py -c atm_uniform.cfg

# Make nominal opacity file (H2O CO CO2 CH4 HCN NH3):
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c opacity_nominal_0.3-33um.cfg
$topdir/pyratbay/pbay.py -c opacity_nominal_0.3-5.5um.cfg
$topdir/pyratbay/pbay.py -c opacity_nominal_1.0-5.5um.cfg  # TBD

# Run MCMC retrievals:
cd $topdir/run02_HATP01b/
$topdir/pyratbay/pbay.py -c mcmc_hatp01b.cfg

# :::  OK  :::

# Run transmission retrievals:
cd $topdir/run02_HATP01b
$topdir/pyratbay/pbay.py -c mcmc_hatp01b.cfg

cd $topdir/run03_HATP11b
$topdir/pyratbay/pbay.py -c mcmc_hatp11b.cfg

cd $topdir/run04_HATP12b
$topdir/pyratbay/pbay.py -c mcmc_hatp12b.cfg

cd $topdir/run05_HATP26b
$topdir/pyratbay/pbay.py -c mcmc_hatp26b.cfg

cd $topdir/run06_WASP006b
$topdir/pyratbay/pbay.py -c mcmc_wasp006b.cfg

cd $topdir/run07_WASP012b
$topdir/pyratbay/pbay.py -c mcmc_wasp012b.cfg

cd $topdir/run08_WASP017b
$topdir/pyratbay/pbay.py -c mcmc_wasp017b.cfg

cd $topdir/run09_WASP019b
$topdir/pyratbay/pbay.py -c mcmc_wasp019b.cfg

cd $topdir/run10_WASP031b
$topdir/pyratbay/pbay.py -c mcmc_wasp031b.cfg

cd $topdir/run11_WASP039b
$topdir/pyratbay/pbay.py -c mcmc_wasp039b.cfg

cd $topdir/run15_WASP121b
$topdir/pyratbay/pbay.py -c mcmc_wasp121b.cfg


# Figure 3:
cd $topdir/run01/
$topdir/fig3.py

# Figures 4 and 5:
cd $topdir/run02/
$topdir/fig4_5.py

