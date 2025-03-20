FROM rocker/r-base:4.4.2

MAINTAINER "Bethany Hamilton" bethanybhat@gmail.com

# Install R packages
RUN install2.r --error \
    tidyr \
    dplyr \
    purrr \
    metafor \
    clubSandwich \
    mvtnorm 

COPY functionsdatagen.R /SimFunctions/functionsdatagen.R

RUN chmod +rx /SimFunctions/functionsdatagen.R

COPY run_sim_study.R /SimFunctions/run_sim_study.R

RUN chmod +rx /SimFunctions/run_sim_study.R

COPY dat_kjN_mathdat.rds /SimFunctions/dat_kjN_mathdat.rds

RUN chmod +rx /SimFunctions/dat_kjN_mathdat.rds

COPY pyl_id_values_test.csv /SimFunctions/pyl_id_values_test.csv

RUN chmod +rx /SimFunctions/pyl_id_values_test.csv

RUN export PATH=/SimFunctions:$PATH

ENV PATH="/SimFunctions:$PATH"