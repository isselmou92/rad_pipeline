# ---------------------------------------------------------------------------
# Lightweight, reproducible container for the whole pipeline (pip-only)
# ---------------------------------------------------------------------------
FROM python:3.9 AS base

# 1. system deps required by SimpleITK, Pyradiomics, etc.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        && \
    rm -rf /var/lib/apt/lists/*

# 2. copy dependency lists first (layer caching)
WORKDIR /workspace
COPY requirements.txt pyproject.toml ./

# 3. install Python deps (wheel, build back-end, runtime libs)
RUN python -m pip install --upgrade pip wheel build && \
    python -m pip install --no-cache-dir -r requirements.txt && \
    python -m pip install --no-cache-dir .  \
    && python -m pip cache purge          # keep image slim

# 4. copy the rest of the source tree (for editable workflows / notebooks)
COPY src ./src

# 5. default shell
CMD [ "bash" ]
