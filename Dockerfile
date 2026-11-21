# reproducible environment: docker build -t loyscope . && docker run --rm -v $PWD/results:/loyscope/results loyscope make figures
FROM ubuntu:24.04
RUN apt-get update && apt-get install -y --no-install-recommends build-essential curl ca-certificates libhdf5-dev pkg-config \
    python3 python3-pip git gh librsvg2-bin fonts-liberation && rm -rf /var/lib/apt/lists/*
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.99.0 --profile minimal
ENV PATH=/root/.cargo/bin:$PATH
WORKDIR /loyscope
COPY requirements.txt .
RUN pip3 install --break-system-packages -r requirements.txt
COPY . .
RUN cargo build --release && cargo test --release
CMD ["make", "figures"]
