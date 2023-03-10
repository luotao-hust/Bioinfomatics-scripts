#!/bin/bash
# VERSION=1.0.0

export TOP_PID=$$
trap 'exit 1' TERM

while getopts ":b:p:" arg;
    do
        case $arg in
             b)
                EXBIND=$(python -c 'print(",".join([i+":"+i for i in "'$OPTARG'".split(",") ]))')
                ;;
             p)
                PORT=$OPTARG
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

mkdir -p "/home/$(id -un)/.config/jupyterlab_apptainer/"

export NEW_STARTSCRIPT="/home/$(id -un)/.config/jupyterlab_apptainer/run_apptainer_jupyterlab_latest.sh"
export NEW_STARTSCRIPT_PATH="https://raw.githubusercontent.com/luotao-hust/Bioinfomatics-scripts/main/jupyterlab_apptainer/run_apptainer_jupyterlab_beta.sh"

wget -q --timeout=6 --tries=2 ${NEW_STARTSCRIPT_PATH} -O ${NEW_STARTSCRIPT}
NEW_VERSION=`sed -n '2p' ${NEW_STARTSCRIPT} | awk -F "=" '{print $2}'`
CURRENT_VERSION=`sed -n '2p' "${basedir}/${FileName}" | awk -F "=" '{print $2}'`

### update startscript
if [[ "${CURRENT_VERSION}" == "${NEW_VERSION}" ]]; then
    :
fi
if [[ "$(printf '%s\n' "${CURRENT_VERSION}" "${NEW_VERSION}" | sort -rV | head -n1)" != "${CURRENT_VERSION}" ]]; then
    echo -e "\e[31mWarning: \e[0m New version is available!"
    mv ${basedir}/${FileName} "/home/$(id -un)/.config/jupyterlab_apptainer/run_apptainer_jupyterlab_bck.sh"
    cp ${NEW_STARTSCRIPT} ${basedir}/${FileName}
    chmod 755 ${basedir}/${FileName}
    echo -e "\033[34mINFO:\033[0m  Finishing script update. Please restart!"
    kill -s TERM $TOP_PID
fi

# 用于检测更新

JUPYTERLAB_TMP=${basedir}/tmp_${PORT}
mkdir -p -m 700 ${JUPYTERLAB_TMP}/run ${JUPYTERLAB_TMP}/tmp

# 默认状态下镜像存放在下面的地址(郭老师的服务器上都有);所以在服务器之外的地方使用，需要自己拷贝镜像文件。
RSTUDIO_SIF="/home/luot/software/apptainer_rstudio/apptainer_rstudio_latest.sif"
# RSTUDIO_SIF="/***/rstudio_latest.sif" # 镜像的地址,可自己修改

export APPTAINER_BIND="${JUPYTERLAB_TMP}/run:/run,${JUPYTERLAB_TMP}/tmp:/tmp,$EXBIND"

apptainer exec ${RSTUDIO_SIF} python3  -m jupyterlab --notebook-dir / --port ${PORT} --ip=0.0.0.0 --no-browser
