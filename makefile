RSCRIPT = Rscript

# if you wildcard the all-target, then nothing will happen if the target doesn't
# exist (no target). hard code the target.
WRITEUP = presanis-conflict.pdf

TEX_FILES = $(wildcard tex-input/*.tex) \
	$(wildcard tex-input/*/*.tex) \
	$(wildcard tex-input/*/*/*.tex)

# base files
RDS = rds
SCRIPTS = scripts
PLOTS = plots
PLOT_SETTINGS = scripts/common/plot-settings.R

# hiv-example
HIV_EXAMPLE = hiv-example
RDS_HIV = $(RDS)/$(HIV_EXAMPLE)
SCRIPTS_HIV = $(SCRIPTS)/$(HIV_EXAMPLE)
PLOTS_HIV = $(PLOTS)/$(HIV_EXAMPLE)
STAN_FILES = $(SCRIPTS)/stan-files

DATA_HIV = $(RDS_HIV)/hiv-data.rds
PARS_HIV = $(RDS_HIV)/pars-of-interest.rds

ALL_PLOTS = $(PLOTS_HIV)/prior-post-compare.pdf \
	$(PLOTS_HIV)/posterior-qq-plot.pdf \
	$(PLOTS_HIV)/p12-prior-post-compare.pdf \
  $(PLOTS_HIV)/p12-only-melding-dists.pdf 

all : $(WRITEUP) 

# knitr is becoming more picky about encoding, specify UTF-8 input
$(WRITEUP) : $(wildcard *.rmd) $(TEX_FILES) $(ALL_PLOTS)
	$(RSCRIPT) -e "rmarkdown::render(input = Sys.glob('*.rmd'), encoding = 'UTF-8')"

$(DATA_HIV) : $(SCRIPTS_HIV)/01-data.R
	$(RSCRIPT) $<

$(PARS_HIV) : $(DATA_HIV)

$(RDS_HIV)/full-model-fit.rds : $(SCRIPTS_HIV)/02-full-model.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn.stan $(PARS_HIV)
	$(RSCRIPT) $<

$(RDS_HIV)/prior-samples.rds : $(SCRIPTS_HIV)/03-sample-prior.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn-prior.stan $(PARS_HIV)
	$(RSCRIPT) $<

$(RDS_HIV)/big-submodel-samples.rds : $(SCRIPTS_HIV)/06-sample-big-submodel.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn-big-submodel.stan $(PARS_HIV)
	$(RSCRIPT) $<

$(RDS_HIV)/small-submodel-samples.rds : $(SCRIPTS_HIV)/07-sample-small-submodel.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn-small-submodel.stan $(PARS_HIV)
	$(RSCRIPT) $<

$(RDS_HIV)/stage-one-samples.rds : $(SCRIPTS_HIV)/08-stage-one-targeter.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn-stage-one-target-big.stan $(RDS_HIV)/prior-samples.rds $(DATA_HIV) $(PARS_HIV)
	$(RSCRIPT) $<

$(RDS_HIV)/big-sub-prior-wsre-est.rds : $(SCRIPTS_HIV)/09-wsre-prior-marginal.R $(STAN_FILES)/hiv-ev-sythn-prior-wsre.stan
	$(RSCRIPT) $<

$(RDS_HIV)/stage-one-wsre-samples.rds : $(SCRIPTS_HIV)/10-stage-one-wsre-est.R $(DATA_HIV) $(STAN_FILES)/hiv-ev-sythn-big-submodel.stan $(PARS_HIV) $(RDS_HIV)/big-sub-prior-wsre-est.rds
	$(RSCRIPT) $<

$(PLOTS_HIV)/prior-post-compare.pdf : $(SCRIPTS_HIV)/04-output-plots.R $(PLOT_SETTINGS) $(RDS_HIV)/prior-samples.rds $(RDS_HIV)/full-model-fit.rds $(PARS_HIV) $(RDS_HIV)/big-submodel-samples.rds $(RDS_HIV)/small-submodel-samples.rds $(RDS_HIV)/stage-one-samples.rds $(RDS_HIV)/stage-one-wsre-samples.rds $(RDS_HIV)/stage-two-samples.rds $(RDS_HIV)/stage-two-wsre-samples.RDS
	$(RSCRIPT) $<

$(PLOTS_HIV)/p12-prior-post-compare.pdf : $(PLOTS_HIV)/prior-post-compare.pdf

$(PLOTS_HIV)/p12-only-melding-dists.pdf : $(PLOTS_HIV)/prior-post-compare.pdf

# this is a bit of a kludge because I `source` 11- as part of 12-
$(RDS_HIV)/stage-two-samples.rds : $(SCRIPTS_HIV)/12-stage-two-melded-posterior.R $(SCRIPTS_HIV)/11-stage-two-priors.R $(PARS_HIV) $(DATA_HIV) $(RDS_HIV)/stage-one-samples.rds $(RDS_HIV)/prior-samples.rds
	$(RSCRIPT) $<

$(RDS_HIV)/stage-two-wsre-samples.rds : $(SCRIPTS_HIV)/13-stage-two-melded-posterior-wsre.R $(SCRIPTS_HIV)/11-stage-two-priors.R $(PARS_HIV) $(DATA_HIV) $(RDS_HIV)/stage-one-wsre-samples.rds $(RDS_HIV)/prior-samples.rds $(RDS_HIV)/big-sub-prior-wsre-est.rds
	$(RSCRIPT) $<

$(PLOTS_HIV)/posterior-qq-plot.pdf : $(SCRIPTS_HIV)/14-posterior-qq-compare.R $(PLOT_SETTINGS) $(RDS_HIV)/stage-two-samples.rds $(RDS_HIV)/stage-two-wsre-samples.rds $(RDS_HIV)/full-model-fit.rds
	$(RSCRIPT) $<
