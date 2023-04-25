#!/bin/bash

while getopts ":o:s:a:g:x:c:b:t:f:d:m:" arg;
    do
        case $arg in
             o)
                OUTPUT_DIR=$OPTARG
                ;;
             s)
                SAMPLE_LIST=$OPTARG
                ;;
             a)
                SIF=$OPTARG
                ;;
             g)
                GTF=$OPTARG
                ;;
             x)
                HISAT2_INDEX=$OPTARG
                ;;
             c)
                CORES=$OPTARG
                ;;
             b)
                EXBIND=$OPTARG
        ;;
             t)
                TYPE=$OPTARG
        ;;   
             f)
                FASTP=$OPTARG
                ;;
             d)
                DRY=$OPTARG
                ;;
         m)
            MODEL=$OPTARG
        ;;
             ?)
                echo "Invalid option: -$OPTARG"
                exit 1
                ;;
        esac
    done
