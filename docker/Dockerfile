FROM ferreol/mira:latest

MAINTAINER Ferréol Soulez <ferreol.soulez@univ-lyon1.fr>

RUN echo http://nl.alpinelinux.org/alpine/edge/testing >> /etc/apk/repositories && \
    apk add --update --no-cache --virtual .build-deps wget  &&\
    cd /opt/yorick/i && wget https://raw.githubusercontent.com/kluskaj/mira-sparco-multi/master/mira2_plugin_sparcomulti.i && \
    apk del .build-deps  && \
    rm -rf /var/cache/apk/* && rm -rf /usr/local/share/doc

ENTRYPOINT ["/usr/local/bin/ymira","-plugin=sparcomulti"]
