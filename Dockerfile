FROM rust:slim-bullseye
RUN apt-get update && apt-get install -y libgmp-dev bc
RUN mkdir /QV-net
WORKDIR /QV-net
COPY ./ .

ENTRYPOINT [ "/bin/bash" ]