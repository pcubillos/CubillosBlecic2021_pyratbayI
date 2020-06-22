# Note about the crazy names:
# The appendix denotes the fitting parameters:
# - zero means no fit, non-zero means fit
# - the order is wmdmha-mc, which stands for:
#   (w)ater, carbon (m)onoxide, carbon (d)ioxide, (m)ethane,
#   (h)ydrogen cyanide, (a)monnia - (m)ean molecular mass, (c)louds

topdir=`pwd`

cd $topdir/run_GJ-0436b_tsiaras
pbay -c mcmc_GJ-0436b_tsiaras_w00000-mc.cfg
pbay -c mcmc_GJ-0436b_tsiaras_w00000-0c.cfg
pbay -c mcmc_GJ-0436b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-03b_tsiaras/
pbay -c mcmc_HATP-03b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-03b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-03b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-17b_tsiaras/
pbay -c mcmc_HATP-17b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-17b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-17b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-18b_tsiaras/
pbay -c mcmc_HATP-18b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-18b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-18b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-32b_damiano/
pbay -c mcmc_HATP-32b_damiano_w00000-mc.cfg
pbay -c mcmc_HATP-32b_damiano_w00000-0c.cfg
pbay -c mcmc_HATP-32b_damiano_w00000-m0.cfg

cd $topdir/run_HATP-32b_tsiaras/
pbay -c mcmc_HATP-32b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-32b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-32b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-38b_bruno/
pbay -c mcmc_HATP-38b_bruno_w00000-mc.cfg
pbay -c mcmc_HATP-38b_bruno_w00000-0c.cfg
pbay -c mcmc_HATP-38b_bruno_w00000-m0.cfg

cd $topdir/run_HATP-38b_tsiaras/
pbay -c mcmc_HATP-38b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-38b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-38b_tsiaras_w00000-m0.cfg

cd $topdir/run_HATP-41b_tsiaras/
pbay -c mcmc_HATP-41b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HATP-41b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HATP-41b_tsiaras_w00000-m0.cfg

cd $topdir/run_HD-149026b_tsiaras/
pbay -c mcmc_HD-149026b_tsiaras_w00000-mc.cfg
pbay -c mcmc_HD-149026b_tsiaras_w00000-0c.cfg
pbay -c mcmc_HD-149026b_tsiaras_w00000-m0.cfg

cd $topdir/run_K2-018b_benneke
pbay -c mcmc_K2-018b_benneke_wmdm00-mc.cfg
pbay -c mcmc_K2-018b_benneke_wmdm00-0c.cfg
pbay -c mcmc_K2-018b_benneke_wmdm00-m0.cfg
pbay -c mcmc_K2-018b_benneke_w00000-mc.cfg
pbay -c mcmc_K2-018b_benneke_w00000-0c.cfg
pbay -c mcmc_K2-018b_benneke_w00000-m0.cfg

cd $topdir/run_K2-018b_tsiaras
pbay -c mcmc_K2-018b_tsiaras_w00000-mc.cfg
pbay -c mcmc_K2-018b_tsiaras_w00000-0c.cfg
pbay -c mcmc_K2-018b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-029b_tsiaras/
pbay -c mcmc_WASP-029b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-029b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-029b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-043b_stevenson/
pbay -c mcmc_WASP-043b_stevenson_wmdm00-mc.cfg
pbay -c mcmc_WASP-043b_stevenson_wmdm00-0c.cfg
pbay -c mcmc_WASP-043b_stevenson_wmdm00-m0.cfg
pbay -c mcmc_WASP-043b_stevenson_w00000-mc.cfg
pbay -c mcmc_WASP-043b_stevenson_w00000-0c.cfg
pbay -c mcmc_WASP-043b_stevenson_w00000-m0.cfg

cd $topdir/run_WASP-043b_tsiaras/
pbay -c mcmc_WASP-043b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-043b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-043b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-063b_kilpatrick/
pbay -c mcmc_WASP-063b_kilpatrick_w00000-mc.cfg
pbay -c mcmc_WASP-063b_kilpatrick_w00000-0c.cfg
pbay -c mcmc_WASP-063b_kilpatrick_w00000-m0.cfg
pbay -c mcmc_WASP-063b_kilpatrick_w000h0-mc.cfg
pbay -c mcmc_WASP-063b_kilpatrick_w000h0-m0.cfg
pbay -c mcmc_WASP-063b_kilpatrick_w000h0-0c.cfg

cd $topdir/run_WASP-063b_tsiaras/
pbay -c mcmc_WASP-063b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-063b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-063b_tsiaras_w00000-m0.cfg
pbay -c mcmc_WASP-063b_tsiaras_w000h0-mc.cfg

cd $topdir/run_WASP-067b_bruno/
pbay -c mcmc_WASP-067b_bruno_w00000-mc.cfg
pbay -c mcmc_WASP-067b_bruno_w00000-0c.cfg
pbay -c mcmc_WASP-067b_bruno_w00000-m0.cfg

cd $topdir/run_WASP-067b_tsiaras/
pbay -c mcmc_WASP-067b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-067b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-067b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-069b_tsiaras/
pbay -c mcmc_WASP-069b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-069b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-069b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-074b_tsiaras/
pbay -c mcmc_WASP-074b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-074b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-074b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-080b_tsiaras/
pbay -c mcmc_WASP-080b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-080b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-080b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-101b_tsiaras/
pbay -c mcmc_WASP-101b_tsiaras_w00000-mc.cfg
pbay -c mcmc_WASP-101b_tsiaras_w00000-0c.cfg
pbay -c mcmc_WASP-101b_tsiaras_w00000-m0.cfg

cd $topdir/run_WASP-101b_wakeford/
pbay -c mcmc_WASP-101b_wakeford_w00000-mc.cfg
pbay -c mcmc_WASP-101b_wakeford_w00000-0c.cfg
pbay -c mcmc_WASP-101b_wakeford_w00000-m0.cfg

cd $topdir/run_WASP-103b_kreidberg/
pbay -c mcmc_WASP-103b_kreidberg_wmdm00-mc.cfg
pbay -c mcmc_WASP-103b_kreidberg_wmdm00-0c.cfg
pbay -c mcmc_WASP-103b_kreidberg_wmdm00-m0.cfg
pbay -c mcmc_WASP-103b_kreidberg_w00000-mc.cfg
pbay -c mcmc_WASP-103b_kreidberg_w00000-0c.cfg
pbay -c mcmc_WASP-103b_kreidberg_w00000-m0.cfg

cd $topdir/run_XO-1b_tsiaras/
pbay -c mcmc_XO-1b_tsiaras_w00000-mc.cfg
pbay -c mcmc_XO-1b_tsiaras_w00000-0c.cfg
pbay -c mcmc_XO-1b_tsiaras_w00000-m0.cfg
