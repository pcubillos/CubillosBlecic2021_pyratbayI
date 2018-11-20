# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# Clone (download) the necessary code:
cd $topdir
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 7bf813d
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 4ba3633
make

# Generate filter files:
cd $topdir/code
$topdir/code/make_filters.py > $topdir/code/filter_info.txt

# Download HCN Exomol data:
cd $topdir/inputs/opacity
wget -i wget_exomol_HCN.txt

# Download H2O HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp_H2O.txt
unzip '*.zip'
rm -f *.zip

# Generate partition-function file for H2O:
cd $topdir/run01
python $topdir/code/pf_tips_H2O.py

# Generate partition-function file for HCN:
cd $topdir/run01
python $topdir/pyratbay/scripts/PFformat_Exomol.py  \
       $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
       $topdir/inputs/opacity/1H-13C-14N__Larner.pf


# Compress LBL databases:
cd $topdir/run01
python $topdir/repack/repack.py repack_H2O.cfg
python $topdir/repack/repack.py repack_HCN.cfg


# Make TLI files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c tli_hitemp_H2O.cfg
python $topdir/pyratbay/pbay.py -c tli_exomol_HCN.cfg


# Make atmospheric files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c atm_uniform.cfg

# Make opacity files:
cd $topdir/run01/
python $topdir/pyratbay/pbay.py -c opacity_H2O_1.0-2.0um.cfg
python $topdir/pyratbay/pbay.py -c opacity_H2O_1.0-5.5um.cfg
python $topdir/pyratbay/pbay.py -c opacity_H2O-HCN_1.0-2.0um.cfg


# Flat-curve fit to the data:
cd $topdir/code/
python $topdir/code/flatfit.py > flatfit.dat


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
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_wm-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp043b_w0-000h0-c.cfg

cd $topdir/run06_WASP063b
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_w0-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-00000-0.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_wm-00000-c.cfg
python $topdir/pyratbay/pbay.py -c mcmc_wasp063b_w0-000h0-c.cfg
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

