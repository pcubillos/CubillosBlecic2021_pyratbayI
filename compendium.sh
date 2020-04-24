# Define topdir (in your top working directory) to make your life easier:
topdir=`pwd`

# Clone (download) the necessary code:
pip install mc3>=3.0.0
pip install pyratbay>=0.9.0a5

# Install petitRADTRANS:
# TBD

# Install taurex:
cd $topdir
git clone https://github.com/ucl-exoplanets/TauREx3_public.git taurex
cd taurex
python setup.py develop

# Install rate:
cd $topdir
git clone https://github.com/pcubillos/rate


# Generate filter files:
cd $topdir
python code/make_filters.py

# Download exomol repack data:
# TBD

# Download HITEMP data:
cd $topdir/inputs/opacity
wget -i wget_hitemp_CO2.txt
unzip '*.zip'
rm -f *.zip

# Download CO data:
cd $topdir/inputs/opacity
wget https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2
bzip2 -d 05_HITEMP2019.par.bz2

# Download and format HITRAN CIA data:
cd $topdir/inputs/opacity
wget https://hitran.org/data/CIA/H2-H2_2011.cia
cd $topdir/run_setup
pbay -cs hitran ../inputs/opacity/H2-H2_2011.cia 2 10

# Generate partition-function files:
cd $topdir/run_setup
pbay -pf exomol ../inputs/opacity/1H-12C-14N__Harris.pf \
                ../inputs/opacity/1H-13C-14N__Larner.pf
pbay -pf exomol ../inputs/opacity/1H2-16O__POKAZATEL.pf
pbay -pf exomol ../inputs/opacity/12C-1H4__YT10to10.pf

# Make TLI files:
cd $topdir/run_setup
pbay -c tli_exomol_H2O.cfg
pbay -c tli_exomol_CH4.cfg
pbay -c tli_exomol_HCN.cfg
pbay -c tli_hitemp_CO.cfg
pbay -c tli_hitemp_CO2.cfg

# Make atmospheric files:
cd $topdir/run_setup
pbay -c atm_uniform.cfg

# Make opacity files:
cd $topdir/run_setup
pbay -c opacity_H2O_0.5-10.0um.cfg
pbay -c opacity_HCN_0.5-10.0um.cfg
pbay -c opacity_CH4_0.5-10.0um.cfg
pbay -c opacity_CO2_0.5-10.0um.cfg
pbay -c opacity_CO_0.5-10.0um.cfg

# Opacity comparison:
cd $topdir/run_validation_opacities
python ../code/fig_validation_opacities.py

# Forward-model comparison:
cd $topdir/run_validation_forward_model
python ../code/fig_validation_forward_model.py


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
cd $topdir/run_setup
python $topdir/fig3.py

# Figures 4 and 5:
cd $topdir/run02/
python $topdir/fig4_5.py

