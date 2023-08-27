#! /usr/bin/bash

set -eu -o pipefail

[[ -x $1 ]]
[[ $1 = /* ]]
[[ $2 =~ ^[[:digit:]]+$ ]]

N=$2
NPROC="$( nproc )"
WORKSPACE="$( mktemp -d --tmpdir 'freq.XXXXXX' )"

function on_exit {
    set +e
    set -v
    popd >/dev/null

    if [[ $( cat /sys/devices/system/cpu/smt/control ) == "off" ]]
    then
        echo on | sudo tee /sys/devices/system/cpu/smt/control >/dev/null
    fi
    if command -v tuna >/dev/null
    then
        sudo tuna include --cpus=0-$NPROC
        >&2 echo "cpus 0-$NPROC are included"
    fi
    #echo powersave | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor >/dev/null

    rm -r "$WORKSPACE"
}
trap on_exit EXIT

pushd "$WORKSPACE"

#echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor >/dev/null
if [[ NPROC > 1 ]]
then
    if [[ $( cat /sys/devices/system/cpu/smt/control ) == "on" ]]
    then
        echo off | sudo tee /sys/devices/system/cpu/smt/control >/dev/null
    fi
    if command -v tuna >/dev/null
    then
        sudo tuna isolate --cpus=1-$NPROC
        >&2 echo "cpus 1-$NPROC are isolated"
    fi
fi

for (( i = 0 ; i < N ; ++i ))
do
    time LC_ALL=C taskset --cpu-list 1-$NPROC "$1"
done
