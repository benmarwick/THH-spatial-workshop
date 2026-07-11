# pull base image 
FROM rocker/binder:4.4

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

RUN adduser "$NB_USER" sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
USER ${NB_USER}

# --- Install R Packages ---
# Copy the installation script into the image and run it as root
COPY install.R /tmp/install.R

# --- Copy RStudio preferences ---
# Ensure the config directory exists and copy the preferences file
RUN mkdir -p /home/${NB_USER}/.config/rstudio/
COPY rstudio-prefs.json /home/${NB_USER}/.config/rstudio/rstudio-prefs.json

RUN Rscript /tmp/install.R

# ---  Copy GitHub files into the container ---
# Copy all files from your repo into the home directory
COPY  readme.qmd /home/${NB_USER}/

#  permissions so the binder user owns everything
USER root
RUN chown -R ${NB_USER}:${NB_USER} /home/${NB_USER}
USER ${NB_USER}
