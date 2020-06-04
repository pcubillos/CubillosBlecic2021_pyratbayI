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
# Patch constants:
cp code/taurex_patch_constants.py taurex/taurex/constants.py
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
wget -i wget_repack_exomol_H2O-CH4-HCN.txt

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
cd $topdir/run_benchmark_opacities
python ../code/fig_benchmark_opacities.py

# Forward-model comparison:
cd $topdir/run_benchmark_forward_model
python ../code/fig_benchmark_forward_model.py

# Retrieval comparison:
cd $topdir
python code/setup_taurex_ariel.py
python code/taurex_ariel_sim.py

cd $topdir/run_benchmark_retrieval
# TBD: Run benchmark retrievals script
python ../code/fig_benchmark_retrieval.py


# Flat-curve fit to the data:
cd $topdir/code/
python $topdir/code/flat_fit.py > stats/flat_fit.dat


# Run MCMC transmission retrievals:
cd $topdir
# TBD: Run WFC3 retrievals script


# Figure 3:
cd $topdir/run_setup
python $topdir/fig3.py

# Figures 4 and 5:
cd $topdir/run02/
python $topdir/fig4_5.py

