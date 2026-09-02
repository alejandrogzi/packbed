# syntax=docker/dockerfile:1

FROM rust:1-alpine AS build
RUN apk add --no-cache musl-dev
WORKDIR /build
COPY packbed/ .
RUN cargo build --release

FROM alpine:3.20
RUN apk add --no-cache procps
COPY --from=build /build/target/release/packbed /usr/local/bin/packbed
ENTRYPOINT ["/usr/local/bin/packbed"]
