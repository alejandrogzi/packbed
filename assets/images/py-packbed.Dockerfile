# syntax=docker/dockerfile:1

FROM ghcr.io/pyo3/maturin:latest AS build
WORKDIR /io
COPY packbed/ packbed/
COPY py-packbed/ py-packbed/
RUN maturin build --release -m py-packbed/Cargo.toml --out /dist

FROM python:3.12-slim
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm /tmp/*.whl
ENTRYPOINT ["python3"]
