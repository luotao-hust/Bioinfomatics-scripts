#!/bin/bash
# VERSION=1.0.3
export TOP_PID=$$
trap 'exit 1' TERM

while getopts ":b:p:" arg;
    do
        case $arg in
             b)
                EXBIND=$(python -c 'print(",".join([i+":"+i for i in "'$OPTARG'".split(",") ]))')
		RAW_BIND=$OPTARG
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

mkdir -p "/home/$(id -un)/.config/rstudio_apptainer"
# PASSWD CHECK
PASSWORD_FILE="/home/$(id -un)/.config/rstudio_apptainer/passwd"
if [ -f "$PASSWORD_FILE" ]; then
    export APPTAINERENV_PASSWORD=`cat ${PASSWORD_FILE}`
else
   
	while true
		do
			read  -s  -p "Please input an new password for rstudio server: " passwd_input
			echo
			read  -s  -p "Verify password: " passwd_va
			echo
			if [[ ${passwd_va} == ${passwd_input} ]]; 
			then
				break
			else
				echo -e "\e[31mWarning: \e[0m The passwords do not match!"
				sleep 3
			fi
		done
	password=`echo -n ${passwd_input} | openssl dgst -sha256`
	password=${password: -64}
	touch ${PASSWORD_FILE}
	echo ${password} > ${PASSWORD_FILE}
	echo
	export APPTAINERENV_PASSWORD=${password}
fi
# 用于检测更新
FileName=$(basename $0)
basedir=`cd $(dirname $0); pwd -P`


export NEW_STARTSCRIPT="/home/$(id -un)/.config/rstudio_apptainer/run_apptainer_rstudio_latest.sh"
export NEW_STARTSCRIPT_PATH="https://raw.githubusercontent.com/luotao-hust/Bioinfomatics-scripts/main/rstudio_apptainer/run_apptainer_rstudio.sh"

wget -q --timeout=6 --tries=2 ${NEW_STARTSCRIPT_PATH} -O ${NEW_STARTSCRIPT}
NEW_VERSION=`sed -n '2p' ${NEW_STARTSCRIPT} | awk -F "=" '{print $2}'`
CURRENT_VERSION=`sed -n '2p' "${basedir}/${FileName}" | awk -F "=" '{print $2}'`

### update startscript
if [[ "${CURRENT_VERSION}" == "${NEW_VERSION}" ]]; then
    :
fi
if [[ "$(printf '%s\n' "${CURRENT_VERSION}" "${NEW_VERSION}" | sort -rV | head -n1)" != "${CURRENT_VERSION}" ]]; then
    while true
        do
            echo -e "\e[31mWarning: \e[0m New version is available!"
            read -r -p "Current version = ${CURRENT_VERSION} , new version = ${NEW_VERSION}; Do you want to update it (y/n)? " input
            case $input in
                [yY][eE][sS]|[yY])
                    mv ${basedir}/${FileName} "/home/$(id -un)/.config/rstudio_apptainer/run_apptainer_rstudio_bk.sh"
                    cp ${NEW_STARTSCRIPT} ${basedir}/${FileName}
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





# 用于检测镜像更新 

# 默认状态下镜像存放在下面的地址(郭老师的服务器上都有);所以在服务器之外的地方使用，需要自己拷贝镜像文件。
OR_RSTUDIO_SIF="/home/luot/software/apptainer_rstudio/apptainer_jupyterlab_rstudio-server.sif"
RSTUDIO_SIF=${basedir}/apptainer_jupyterlab_rstudio-server/
# RSTUDIO_SIF="/***/rstudio_latest.sif" # 镜像的地址,可自己修改

if [ -d "${RSTUDIO_SIF}" ];then
    :
    else        
        echo -e "\033[34mINFO:\033[0m Extracting files.........."
        apptainer build --sandbox ${RSTUDIO_SIF} ${OR_RSTUDIO_SIF}
fi

###################################################################################################


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


apptainer instance start ${RSTUDIO_SIF} rstudio_${PORT}
