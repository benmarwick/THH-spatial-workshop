FROM rocker/binder:4.4.1

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libglpk-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

COPY readme.qmd /home/rstudio/readme.qmd
COPY install.R /tmp/install.R

WORKDIR /home/rstudio

RUN Rscript /tmp/install.R

USER rstudio
