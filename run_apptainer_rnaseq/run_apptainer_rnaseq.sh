#!/bin/bash
# script_version=0.3.0
# image_version=0.3.0

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

# 用于检测更新
FileName=$(basename $0)
basedir=`cd $(dirname $0); pwd -P`

mkdir -p "/home/$(id -un)/.config/rnaseq_apptainer/"

export NEW_STARTSCRIPT="/home/$(id -un)/.config/rnaseq_apptainer/run_apptainer_rnaseq.sh"
export NEW_STARTSCRIPT_PATH="https://raw.githubusercontent.com/luotao-hust/Bioinfomatics-scripts/main/run_apptainer_rnaseq/run_apptainer_rnaseq.sh"

wget -q -4 --timeout=6 --tries=2 ${NEW_STARTSCRIPT_PATH} -O ${NEW_STARTSCRIPT}
NEW_VERSION_SCRIPT=`sed -n '2p' ${NEW_STARTSCRIPT} | awk -F "=" '{print $2}'`
CURRENT_VERSION_SCRIPT=`sed -n '2p' "${basedir}/${FileName}" | awk -F "=" '{print $2}'`
NEW_VERSION_IMAGE=`sed -n '3p' ${NEW_STARTSCRIPT} | awk -F "=" '{print $2}'`
CURRENT_VERSION_IMAGE=`sed -n '3p' "${basedir}/${FileName}" | awk -F "=" '{print $2}'`

### update startscript

if [[ "$(printf '%s\n' "${CURRENT_VERSION_IMAGE}" "${NEW_VERSION_IMAGE}" | sort -rV | head -n1)" != "${CURRENT_VERSION_IMAGE}" ]]; then
    echo -e "\e[31mWarning: \e[0m New image version is available!"
fi

if [[ "$(printf '%s\n' "${CURRENT_VERSION_SCRIPT}" "${NEW_VERSION_SCRIPT}" | sort -rV | head -n1)" != "${CURRENT_VERSION_SCRIPT}" ]]; then
    while true
        do
            echo -e "\e[31mWarning: \e[0m New version is available!"
            read -r -p "Current version = ${CURRENT_VERSION_SCRIPT} , new version = ${NEW_VERSION_SCRIPT}; Do you want to update it (y/n)? " input
            case $input in
                [yY][eE][sS]|[yY])
                    mv ${basedir}/${FileName} "/home/$(id -un)/.config/rnaseq_apptainer/run_apptainer_rnaseq_bk.sh"
                    cp ${NEW_VERSION_SCRIPT} ${basedir}/${FileName}
                    chmod 755 ${basedir}/${FileName}
                    echo -e "\033[34mINFO:\033[0m  Finishing script update. Please restart!"
                    ;;
                [nN][oO]|[nN])
                    break
                    ;;

                *)
                    echo -e "\e[31m Warning: \e[0m Invalid input..."
                    sleep 3
                    ;;
            esac
        done
fi

# 用于检测更新

INPUT_FASTQ_FILE01=`head -n 1 ${SAMPLE_LIST} | awk -F ";" '{print $2}' `
INPUT_DIR=$(python -c 'print("/"+"/".join("'$INPUT_FASTQ_FILE01'".split("/")[1:3]))')
mkdir -p ${OUTPUT_DIR}
SAMPLE_LIST_DIR=${SAMPLE_LIST%/*}
USER_NAME=$(id -un)

if test -z $HISAT2_INDEX
then
    export APPTAINER_BIND="${INPUT_DIR}:${INPUT_DIR},${OUTPUT_DIR}:${OUTPUT_DIR},${SAMPLE_LIST_DIR}:${SAMPLE_LIST_DIR},/tmp:${OUTPUT_DIR},${EXBIND}"
else
    TOTAL_WORK_DIR=$(python -c 'print("/"+"/".join("'$HISAT2_INDEX'".split("/")[1:3]))')
    export APPTAINER_BIND="${INPUT_DIR}:${INPUT_DIR},${OUTPUT_DIR}:${OUTPUT_DIR},${SAMPLE_LIST_DIR}:${SAMPLE_LIST_DIR},${TOTAL_WORK_DIR}:${TOTAL_WORK_DIR},/tmp:${OUTPUT_DIR},${EXBIND}"
fi

if test -z $DRY
then
    DRY_RUN="F"
else
    DRY_RUN="T"
fi

apptainer exec ${SIF}  \
    bash -c "cd /biosoftware/script/ && python3 parse_submit.py \
            -s ${SAMPLE_LIST} \
            -o ${OUTPUT_DIR} \
            -g ${GTF} \
            -x ${HISAT2_INDEX} \
            -c ${CORES} \
            -t ${TYPE} \
            -u ${USER_NAME} \
            -f ${FASTP} \
            -d ${DRY_RUN} \
            -m ${MODEL}
            "
