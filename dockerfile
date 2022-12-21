FROM archlinux

RUN yes | pacman -Syu gcc make arb parallel

COPY . /MAXDICUT

WORKDIR /MAXDICUT

RUN make
