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

ALL_PLOTS = $(PLOTS_HIV)/prior-post-compare.pdf

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

$(PLOTS_HIV)/prior-post-compare.pdf : $(SCRIPTS_HIV)/04-output-plots.R $(PLOT_SETTINGS) $(RDS_HIV)/prior-samples.rds $(RDS_HIV)/full-model-fit.rds $(PARS_HIV)
	$(RSCRIPT) $<