FROM rhub/r-minimal:4.5.0
RUN mkdir /app && wget -P /app https://raw.githubusercontent.com/morgannprice/BobaseqFitness/refs/heads/main/bobaseq.R
RUN R -e 'install.packages("optparse", repos="http://cloud.r-project.org")'
COPY ./build_r_image.R /app/
RUN chmod +x /app/*.R
