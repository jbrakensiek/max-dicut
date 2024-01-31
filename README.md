This repository is a companion to the paper

Separating MAX 2-AND, MAX DI-CUT and MAX CUT by Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

-----

**Running in Docker**

To run this in docker, first install docker:

https://docs.docker.com/get-docker/

Then, try the following in a terminal (cross-platform):

- `docker build -t "maxdicut" .`
- `docker run -it maxdicut`

This opens up an interactive shell in the docker. Then, run

- `pacman -Syu gcc make flint parallel`

to install dependencies. Now, run

- `make`
- `mkdir out`
- `time ./par.sh 0.5 data/candidate_dicut.txt`

Where `0.5` should be replaced with the ratio to be verified and `candidate_dicut.txt` can also be `candidate_2and.txt` (or `improved_dicut.txt` and `improved_2and.txt`)

**Distribution File Format**

Each distribution file has the following format:

- The first line has two space-separated integers: (1) the number of break points in the function and (2) the number of functions.

- The second line has a list of the breakpoints.

- Each subsequent line has each function in the distribution. The first number is the probability weight for that function in the distribution, each subsequent number is the value of the function at the breakpoints (to then be linearly interpolated).

Note that the total sum of the weights need not be 1.
