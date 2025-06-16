# Lightweight, reproducible container for the whole pipeline
FROM mambaorg/micromamba:1.5.6

# 1. install conda environment
COPY environment.yml /env.yml
RUN micromamba env create -f /env.yml && \
    micromamba clean --all --yes

# 2. activate env + set PATH
ENV PATH=/opt/conda/envs/preclinical-radiomics/bin:$PATH
WORKDIR /workspace

# 3. copy the source (so `pip install -e .` works in CI too, optional)
# COPY . /workspace

CMD ["bash"]
