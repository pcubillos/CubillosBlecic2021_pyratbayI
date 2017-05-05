# Clone code
# From the directory where this file is located, execute:
topdir=`pwd`
cd $topdir
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 33f4c70  # FINDME: Update when final
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 7da48c3  # FINDME: Update when final
make


# Download Exomol data:
cd $topdir/inputs/opacity
wget -i wget_exomol_NH3.txt
wget -i wget_exomol_HCN.txt
bzip2 -d *.bz2

# Download HITRAN/HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_H2O_CO2_CO_CH4.txt
unzip '*.zip'
rm -f *.zip


# Generate partition-function files for H2O and NH3:
cd $topdir/run01
$topdir/code/pf_tips_H2O-NH3.py

# Compress LBL databases:
cd $topdir/run01
$topdir/repack/repack.py repack_H2O.cfg
$topdir/repack/repack.py repack_HCN.cfg
$topdir/repack/repack.py repack_NH3.cfg

# :::  OK  :::

# Generate HCN partition function:  
cd $topdir/run01
$topdir/pyratbay/scripts/PFformat_Exomol.py \
      $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
      $topdir/inputs/opacity/1H-13C-14N__Larner.pf


# Make TLI files:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c tli_byte_NH3.cfg
$topdir/pyratbay/pbay.py -c tli_exomol_HCN.cfg
$topdir/pyratbay/pbay.py -c tli_hitemp_CO2.cfg
$topdir/pyratbay/pbay.py -c tli_hitemp_CO.cfg
$topdir/pyratbay/pbay.py -c tli_hitemp_H2O.cfg
$topdir/pyratbay/pbay.py -c tli_hitran_CH4.cfg


# Make atmospheric files:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c atm_000.10xsolar.cfg
$topdir/pyratbay/pbay.py -c atm_001.00xsolar.cfg
$topdir/pyratbay/pbay.py -c atm_100.00xsolar.cfg


# Make H2O opacity file:
cd $topdir/run01/
$topdir/pyratbay/pbay.py -c opacity_H2O.cfg


# Run MCMC for solar abundance model:
cd $topdir/run02/
$topdir/pyratbay/pbay.py -c mcmc_000.1xsolar.cfg
$topdir/pyratbay/pbay.py -c mcmc_001.0xsolar.cfg
$topdir/pyratbay/pbay.py -c mcmc_100.0xsolar.cfg


# Figure 3:
cd $topdir/run01/
$topdir/fig3.py

# Figures 4 and 5:
cd $topdir/run02/
$topdir/fig4_5.py
