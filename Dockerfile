FROM python:3.11-slim-bookworm

# ViennaRNA 2.5.1 is supplied by Debian bookworm's non-free component.
RUN sed -i 's/Components: main/Components: main non-free/' /etc/apt/sources.list.d/debian.sources \
    && apt-get update \
    && apt-get install --yes --no-install-recommends vienna-rna=2.5.1+dfsg-1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements-dev.lock pyproject.toml README.md LICENSE MANIFEST.in ./
COPY src ./src
COPY tests ./tests
COPY examples ./examples
RUN python -m pip install --no-cache-dir -r requirements-dev.lock \
    && python -m pip install --no-cache-dir --no-deps .

ENV PAYLOAD_COMPILER_REQUIRE_REAL_RNAFOLD=1
CMD ["pytest", "-q"]
