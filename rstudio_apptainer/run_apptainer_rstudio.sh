#!/bin/bash
# VERSION=1.0.1

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

# 默认状态下镜像存放在下面的地址(郭老师的服务器上都有);所以在服务器之外的地方使用，需要自己拷贝镜像文件。
RSTUDIO_SIF="/home/luot/software/apptainer_rstudio/apptainer_rstudio_latest.sif"
# RSTUDIO_SIF="/***/rstudio_latest.sif" # 镜像的地址,可自己修改

RSTUDIO_TMP=`pwd`/tmp_${PORT}
SERVER_IP=$(python -c "import socket;print([(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1])")
echo -e "\033[34mINFO:\033[0m    Rstudio server is running at \033[32m http://${SERVER_IP}:${PORT}\033[0m  Using your web browser with your server username to log in !"

export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_PORT=${PORT}
export APPTAINERENV_USER=$(id -un)

mkdir -p -m 700 \
        ${RSTUDIO_TMP}/run \
        ${RSTUDIO_TMP}/tmp/rstudio-server \
        ${RSTUDIO_TMP}/var/lib/rstudio-server \
        ${RSTUDIO_TMP}/.config \
        ${RSTUDIO_TMP}/.local \
        ${RSTUDIO_TMP}/.cache

UUID_PATH=${RSTUDIO_TMP}/tmp/rstudio-server/secure-cookie-key

if [ ! -f "${UUID_PATH}" ]; then
    uuidgen | sed 's/-//g' > ${UUID_PATH}
    chmod 0600 ${UUID_PATH}
fi

cat > ${RSTUDIO_TMP}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

cat > ${RSTUDIO_TMP}/rsession.conf <<END
# R Session Configuration File
END

export APPTAINER_BIND="${RSTUDIO_TMP}/run:/run,${RSTUDIO_TMP}/tmp:/tmp,${RSTUDIO_TMP}/database.conf:/etc/rstudio/database.conf,${RSTUDIO_TMP}/rsession.conf:/etc/rstudio/rsession.conf,${RSTUDIO_TMP}/var/lib/rstudio-server:/var/lib/rstudio-server,${RSTUDIO_TMP}/.cache:/home/${APPTAINERENV_USER}/.cache,${RSTUDIO_TMP}/.config:/home/${APPTAINERENV_USER}/.config,${RSTUDIO_TMP}/.local/:/home/${APPTAINERENV_USER}/.local,$EXBIND"

# PASSWD CHECK
PASSWORD_FILE="/home/${APPTAINERENV_USER}/.config/rstudio_apptainer/passwd"
if [ -f "$PASSWORD_FILE" ]; then
    export APPTAINERENV_PASSWORD=`cat ${PASSWORD_FILE}`
else
    mkdir -p "/home/${APPTAINERENV_USER}/.config/rstudio_apptainer"
	read  -s  -p "Please input an new password for rstudio server:" passwd_input
	password=`echo -n ${passwd_input} | openssl dgst -sha256`
    password=${password: -64}
	echo ${password} > ${PASSWORD_FILE}
	echo
	export APPTAINERENV_PASSWORD=${password}
fi


apptainer instance start ${RSTUDIO_SIF} rstudio_${PORT}