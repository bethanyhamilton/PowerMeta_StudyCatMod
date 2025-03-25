FROM rocker/r-base:4.4.2

# Make a directory in the container
RUN mkdir -p /home/r-environment/output

# Install R packages
RUN install2.r --error \
    tidyr \
    dplyr \
    purrr \
    metafor \
    clubSandwich \
    mvtnorm 

COPY functionsdatagen.R /home/r-environment/functionsdatagen.R

RUN chmod +rx /home/r-environment/functionsdatagen.R

COPY run_sim_study.R /home/r-environment/run_sim_study.R

RUN chmod +rx /home/r-environment/run_sim_study.R

COPY dat_kjN_mathdat.rds /home/r-environment/dat_kjN_mathdat.rds

# Run the R script
CMD R -e "source('/home/r-environment/run_sim_study.R')" --args PYL_ID 1 