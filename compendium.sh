# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# Clone (download) the necessary code:
cd $topdir
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 8198c1e
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 4ba3633
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

# Generate partition-function files for H2O and NH3:
cd $topdir/run01
python $topdir/code/pf_tips_H2O-NH3.py

# Generate partition-function files for HCN, CH4, and TiO:
cd $topdir/run01
python $topdir/pyratbay/scripts/PFformat_Exomol.py  \
       $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
       $topdir/inputs/opacity/1H-13C-14N__Larner.pf
python $topdir/pyratbay/scripts/PFformat_Exomol.py \
       $topdir/inputs/opacity/12C-1H4__YT10to10.pf
python $topdir/pyratbay/scripts/PFformat_Schwenke_TiO.py \
       $topdir/inputs/opacity/tiopart.dat


# Compress LBL databases:
cd $topdir/run01
python $topdir/repack/repack.py repack_H2O.cfg
python $topdir/repack/repack.py repack_HCN.cfg
python $topdir/repack/repack.py repack_NH3.cfg
python $topdir/repack/repack.py repack_CH4.cfg


# Make TLI files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c tli_hitemp_CO2.cfg
python $topdir/pyratbay/pbay.py -c tli_Li_CO.cfg
python $topdir/pyratbay/pbay.py -c tli_hitemp_H2O.cfg
python $topdir/pyratbay/pbay.py -c tli_exomol_HCN.cfg
python $topdir/pyratbay/pbay.py -c tli_exomol_NH3.cfg
python $topdir/pyratbay/pbay.py -c tli_exomol_CH4.cfg


# Make atmospheric files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c atm_1xsolar_1000K.cfg
python $topdir/pyratbay/pbay.py -c atm_1xsolar_1500K.cfg
python $topdir/pyratbay/pbay.py -c atm_1xsolar_2000K.cfg
python $topdir/pyratbay/pbay.py -c atm_1xsolar_2500K.cfg
python $topdir/pyratbay/pbay.py -c atm_uniform.cfg

# Make nominal opacity file (H2O CO CO2 CH4 HCN NH3):
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c opacity_nominal_1.0-5.5um.cfg


# Run MCMC transmission retrievals:
cd $topdir/run02_HATP11b
python $topdir/pyratbay/pbay.py -c mcmc_hatp11b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp11b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp11b_wm-00000-c.cfg

cd $topdir/run03_HATP32b
python $topdir/pyratbay/pbay.py -c mcmc_hatp32b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp32b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp32b_wm-00000-c.cfg

cd $topdir/run04_HATP38b
python $topdir/pyratbay/pbay.py -c mcmc_hatp38b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp38b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_hatp38b_wm-00000-c.cfg

cd $topdir/run05_WASP043b
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_w0-000h0-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-000h0-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-000h0-c.cfg

cd $topdir/run06_WASP063b
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_w0-000h0-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-000h0-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-000h0-c.cfg

cd $topdir/run07_WASP067b
python $topdir/pyratbay/pbay.py -c mcmc_wasp067b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp067b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp067b_wm-00000-c.cfg

cd $topdir/run08_WASP101b
python $topdir/pyratbay/pbay.py -c mcmc_wasp101b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp101b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp101b_wm-00000-c.cfg

cd $topdir/run09_WASP107b
python $topdir/pyratbay/pbay.py -c mcmc_wasp107b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp107b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp107b_wm-00000-c.cfg


# Figure 3:
cd $topdir/run01/
python $topdir/fig3.py

# Figures 4 and 5:
cd $topdir/run02/
python $topdir/fig4_5.py

