#/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/../../../manageUtils.sh

mirroredProject chemistry

BASE=$HGROOT/programs/data/chemistry

case "$1" in
mirror)
  syncHg  
;;

esac

