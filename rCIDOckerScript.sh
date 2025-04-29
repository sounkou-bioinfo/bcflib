# launch docker using 'r-bspm:20.04' from rocker

docker run --rm -ti -v ${PWD}:/work -w /work rocker/r-bspm:20.04

# fetch the script
docker run --rm -ti -v ${PWD}:/work -w /work wget https://eddelbuettel.github.io/r-ci/run.sh && chmod 0755 run.sh

# bootstrap
docker run --rm -ti -v ${PWD}:/work -w /work ./run.sh bootstrap

# install just depends
docker run --rm -ti -v ${PWD}:/work -w /work ./run.sh install_deps

# test (with just depends)
export _R_CHECK_FORCE_SUGGESTS_="false"; ./run.sh run_tests

# alternate: install all
./run.sh install_all

# test (with all)
export _R_CHECK_FORCE_SUGGESTS_="true"; ./run.sh run_tests
