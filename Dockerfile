FROM python:3.9

ARG MULTIQC_VERSION=v1.14

LABEL author="David Brawand" \
      description="MultiQC ${MULTIQC_VERSION} with SEGLH plugin" \
      maintainer="dbrawand@nhs.net"

# Install MultiQC
RUN git clone https://github.com/ewels/MultiQC.git --branch ${MULTIQC_VERSION} && \
    cd MultiQC && \
    git checkout  && \
    python setup.py install

# Add plugin and install
COPY . multiqc_plugins
RUN cd multiqc_plugins && \
    python setup.py install

# uncomment this is the dockerfile is distributed outside of this repository 
#RUN git clone https://github.com/moka-guys/multiqc_plugins.git --branch main && \
#    cd multiqc_plugins && \
#    python setup.py install

RUN mkdir -p /data
WORKDIR /data

ENTRYPOINT [ "multiqc" ]
CMD [ "." ]
