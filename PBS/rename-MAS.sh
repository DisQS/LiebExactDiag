#!/bin/bash
find . -name "*-????.raw" | awk -F'[.-]' -v mvCmd='mv -vu "%s" "%s"\n'     '{ num=sprintf("%05d", $(NF-1));
       old=$0;
       sub(/[0-9]+'.raw'/,num".raw");
       printf (mvCmd,old,$0);
     }' > rename-SLA.sh
#cat rename-SLA.sh
source ./rename-SLA.sh

