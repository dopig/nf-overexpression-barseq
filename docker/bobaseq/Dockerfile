FROM python:3.12.9-bookworm

# Install system dependencies
RUN apt-get update && apt-get install -y minimap2=2.24* vsearch=2.22.1*
    # minimap2=2.24-r1122 \
    # vsearch=2.22.1-1

# Install usearch (adjust version as needed)
RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
    gunzip usearch11.0.667_i86linux32.gz && \
    chmod +x usearch11.0.667_i86linux32 && \
    mv usearch11.0.667_i86linux32 /usr/local/bin/usearch

# Install python dependencies.
RUN pip install pandas==2.2.3 matplotlib==3.10.1

# # Clone Boba-seq repository
RUN git clone https://github.com/OGalOz/Boba-seq.git

# Pandas append is gone after 1.4 and installing older pandas didn't work well
RUN sed -i 's/pos_bc = pos_bc.append(df_bc_keep\[\["query", "qhi_x", "qlo_y"\]\])/pos_bc = pd.concat(\[pos_bc, df_bc_keep\[\["query", "qhi_x", "qlo_y"\]\]\])/g' /Boba-seq/src/step2.py
RUN sed -i 's/pos_ins = pos_ins.append(df_ins_keep\[\["query", "qhi_x", "qlo_y"\]\])/pos_ins = pd.concat(\[pos_ins, df_ins_keep\[\["query", "qhi_x", "qlo_y"\]\]\])/g' /Boba-seq/src/step2.py
# Update for current Matplotlib
RUN sed -i 's/"seaborn-whitegrid"/"seaborn-v0_8-whitegrid"/g' /Boba-seq/src/step6.py


# RUN sed -i 's/shutil.rmtree(op_dir, ignore_errors=True)/shutil.rmtree(op_dir, ignore_errors=False)/g' /Boba-seq/src/util_mod.py

# # Set working directory
WORKDIR /work

ENTRYPOINT ["python3", "/Boba-seq/src/run_steps.py"]

## Below seems smoother, but can't get conda usearch to work

# FROM  mambaorg/micromamba:latest

# COPY ./environment.yml /environment.yml

# RUN micromamba create -y -f /environment.yml && \
#     micromamba clean --all --yes

# WORKDIR /Boba-seq

# RUN micromamba run -n bobaseq-env git clone https://github.com/OGalOz/Boba-seq.git /Boba-seq

# ENTRYPOINT ["micromamba", "run", "-n", "bobaseq-env"]
