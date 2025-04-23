FROM r-base:4.4.3
RUN mkdir /app && wget -P /app https://raw.githubusercontent.com/morgannprice/BobaseqFitness/refs/heads/main/bobaseq.R
RUN R -e "install.packages('optparse', lib='/usr/lib/R/site-library', repos='http://cran.us.r-project.org')"
COPY ./build_r_image.R ./plot_fitness.R /app/
RUN chmod +x /app/*.R
